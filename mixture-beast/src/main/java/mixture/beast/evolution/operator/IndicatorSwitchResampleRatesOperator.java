package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import mixture.beast.evolution.mixture.RelaxedRatesPriorSVS;

import java.util.Arrays;

@Description("Bridge move for SVS: flip indicator (0<->1) AND resample the entire rawRates vector "
        + "from the target model's prior (UC or AC) given current hyperparameters. "
        + "This helps the chain move between the UC-mode and AC-mode basins.")
public class IndicatorSwitchResampleRatesOperator extends Operator {

    public final Input<RelaxedRatesPriorSVS> priorInput = new Input<>(
            "prior",
            "RelaxedRatesPriorSVS instance that defines UC/AC densities and holds references to tree/rates/indicator/hypers.",
            Input.Validate.REQUIRED
    );

    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator", "SVS indicator", Input.Validate.REQUIRED);

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates", "raw branch rates", Input.Validate.REQUIRED);


    private RelaxedRatesPriorSVS prior;

    private Tree tree;
    private RealParameter rates;
    private IntegerParameter indicator;

    private RealParameter ucldStdev;
    private RealParameter sigma2;
    private RealParameter rootLogRate; // may be null
    private double minDt;

    private int nNodes;
    private int rootNr;
    private int[] idxMap;

    @Override
    public void initAndValidate() {
        prior = priorInput.get();
        indicator = indicatorInput.get();
        rates = ratesInput.get();
        // Grab the SAME objects used by the prior
        tree = prior.treeInput.get();
        rates = prior.ratesInput.get();
        indicator = prior.indicatorInput.get();
        ucldStdev = prior.ucldStdevInput.get();
        sigma2 = prior.sigma2Input.get();
        rootLogRate = prior.rootLogRateInput.get(); // optional
        minDt = prior.minBranchLengthInput.get();

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();

        if (indicator.getDimension() != 1) {
            throw new IllegalArgumentException("IndicatorSwitchResampleRatesOperator: indicator must have dimension=1.");
        }
        if (rates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("IndicatorSwitchResampleRatesOperator: rates must have dimension nodeCount-1.");
        }

        buildIndexMap();
    }

    private void buildIndexMap() {
        idxMap = new int[nNodes];
        Arrays.fill(idxMap, -2);

        int k = 0;
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            int nr = node.getNr();
            if (nr == rootNr) {
                idxMap[nr] = -1;
            } else {
                idxMap[nr] = k++;
            }
        }

        if (k != nNodes - 1) {
            throw new IllegalStateException("Expected to map " + (nNodes - 1) + " non-root nodes but mapped " + k);
        }
    }

    @Override
    public double proposal() {
        final int oldK = indicator.getValue(0);
        if (!(oldK == 0 || oldK == 1)) {
            return Double.NEGATIVE_INFINITY;
        }
        final int newK = 1 - oldK;

        // log prior density for the CURRENT state under its CURRENT model
        final double logOld = (oldK == 0) ? prior.logPriorUCOnly() : prior.logPriorACOnly();
        if (Double.isInfinite(logOld) || Double.isNaN(logOld)) {
            return Double.NEGATIVE_INFINITY;
        }

        // Propose new rawRates from the TARGET model prior
        final boolean ok = (newK == 0) ? sampleUC() : sampleAC();
        if (!ok) {
            return Double.NEGATIVE_INFINITY;
        }

        // Flip indicator
        indicator.setValue(0, newK);

        // log prior density for the PROPOSED state under the TARGET model
        final double logNew = (newK == 0) ? prior.logPriorUCOnly() : prior.logPriorACOnly();
        if (Double.isInfinite(logNew) || Double.isNaN(logNew)) {
            return Double.NEGATIVE_INFINITY;
        }

        // Hastings ratio:
        // q(new|old) uses the TARGET prior density; q(old|new) uses the OLD prior density.
        // So logHR = log q(old|new) - log q(new|old) = logOld - logNew
        return logOld - logNew;
    }

    /** Sample i.i.d LogNormal on rates with E[r]=1 (matches prior implementation). */
    private boolean sampleUC() {
        final double s = ucldStdev.getValue();
        if (!(s > 0.0)) return false;

        final double s2 = s * s;
        final double mu = -0.5 * s2;

        for (int i = 0; i < rates.getDimension(); i++) {
            final double z = Randomizer.nextGaussian();
            final double logR = mu + s * z;
            final double r = Math.exp(logR);
            if (!(r > 0.0) || Double.isInfinite(r) || Double.isNaN(r)) return false;
            rates.setValue(i, r);
        }
        return true;
    }

    /** Sample mean-corrected AC increments along the tree (matches prior implementation). */
    private boolean sampleAC() {
        final double s2 = sigma2.getValue();
        if (!(s2 > 0.0)) return false;
        if (!(minDt > 0.0)) return false;

        final double rootLog = (rootLogRate == null ? 0.0 : rootLogRate.getValue());
        Node root = tree.getRoot();
        return sampleACFrom(root, rootLog, s2);
    }

    private boolean sampleACFrom(Node parent, double logPar, double s2) {
        for (int c = 0; c < parent.getChildCount(); c++) {
            Node child = parent.getChild(c);

            final double dt = child.getLength();
            if (!(dt > minDt)) {
                // If your trees can have extremely short branches, set minBranchLength very small in the prior.
                return false;
            }

            final double var = s2 * dt;
            if (!(var > 0.0)) return false;

            final double mean = logPar - 0.5 * var;
            final double z = Randomizer.nextGaussian();
            final double logChi = mean + Math.sqrt(var) * z;

            final double rChi = Math.exp(logChi);
            if (!(rChi > 0.0) || Double.isInfinite(rChi) || Double.isNaN(rChi)) return false;

            final int idx = idxMap[child.getNr()];
            if (idx < 0) return false;
            rates.setValue(idx, rChi);

            if (!sampleACFrom(child, logChi, s2)) return false;
        }
        return true;
    }
}
