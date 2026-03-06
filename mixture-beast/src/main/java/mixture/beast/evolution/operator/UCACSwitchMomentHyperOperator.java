package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.util.Arrays;

@Description("UC<->AC switch operator that also proposes sigma2 and ucldStdev using fast moment/MLE-like "
        + "estimates from the current shared rates. Helps big trees avoid 'prior mismatch' when switching.")
public class UCACSwitchMomentHyperOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);
    public final Input<RealParameter> ratesInput = new Input<>("rates", "shared positive branch rates (non-root nodes)", Input.Validate.REQUIRED);
    public final Input<IntegerParameter> indicatorInput = new Input<>("indicator", "0=UC, 1=AC", Input.Validate.REQUIRED);

    public final Input<RealParameter> ucldStdevInput = new Input<>("ucldStdev", "UC lognormal stdev (sigma on log scale)", Input.Validate.REQUIRED);
    public final Input<RealParameter> sigma2Input = new Input<>("sigma2", "AC Brownian variance per unit time", Input.Validate.REQUIRED);
    public final Input<RealParameter> rootLogRateInput = new Input<>("rootLogRate", "optional root log-rate anchor (default 0)", Input.Validate.OPTIONAL);

    public final Input<Double> minBranchLengthInput = new Input<>("minBranchLength", "min dt used for AC moment estimate", 1e-12);

    public final Input<Double> ucldLogSdInput = new Input<>("ucldLogSd", "proposal SD on log(ucldStdev)", 0.6);
    public final Input<Double> sigma2LogSdInput = new Input<>("sigma2LogSd", "proposal SD on log(sigma2)", 0.6);

    public final Input<Double> minPositiveInput = new Input<>("minPositive", "floor for estimated centers", 1e-12);

    private static final double LOG_SQRT_2PI = 0.5 * Math.log(2.0 * Math.PI);

    private Tree tree;
    private RealParameter rates;
    private IntegerParameter indicator;
    private RealParameter ucldStdev;
    private RealParameter sigma2;
    private RealParameter rootLogRate;

    private int nNodes;
    private int rootNr;
    private int[] idxMap; // nodeNr -> ratesIndex, root -> -1

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();
        indicator = indicatorInput.get();
        ucldStdev = ucldStdevInput.get();
        sigma2 = sigma2Input.get();
        rootLogRate = rootLogRateInput.get();

        if (indicator.getDimension() != 1) throw new IllegalArgumentException("indicator dimension must be 1");
        if (ucldStdev.getDimension() != 1) throw new IllegalArgumentException("ucldStdev dimension must be 1");
        if (sigma2.getDimension() != 1) throw new IllegalArgumentException("sigma2 dimension must be 1");
        if (rootLogRate != null && rootLogRate.getDimension() != 1) throw new IllegalArgumentException("rootLogRate dimension must be 1");

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();
        if (rates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("rates must have dimension nodeCount-1");
        }
        rebuildIdxMap_likeYourPrior();
    }

    private void rebuildIdxMap_likeYourPrior() {
        idxMap = new int[nNodes];
        Arrays.fill(idxMap, -2);

        int k = 0;
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            int nr = node.getNr();
            if (nr < 0 || nr >= nNodes) throw new IllegalArgumentException("Node nr out of range: " + nr);
            if (nr == rootNr) idxMap[nr] = -1;
            else idxMap[nr] = k++;
        }
        if (k != nNodes - 1) throw new IllegalStateException("Expected " + (nNodes - 1) + " non-root nodes, got " + k);
    }

    private void ensureMapsUpToDate() {
        int newN = tree.getNodeCount();
        int newRoot = tree.getRoot().getNr();
        if (newN != nNodes || newRoot != rootNr) {
            nNodes = newN;
            rootNr = newRoot;
            if (rates.getDimension() != nNodes - 1) throw new IllegalStateException("rates dimension inconsistent with tree");
            rebuildIdxMap_likeYourPrior();
        }
    }

    private double rootLog() {
        return (rootLogRate == null ? 0.0 : rootLogRate.getValue(0));
    }


    private static double logLogNormalPdf(final double value, final double meanLog, final double sdLog) {
        if (!(value > 0.0) || !(sdLog > 0.0)) return Double.NEGATIVE_INFINITY;
        final double x = Math.log(value);
        final double z = (x - meanLog) / sdLog;
        return -x - Math.log(sdLog) - LOG_SQRT_2PI - 0.5 * z * z;
    }

    /** Fast estimate of UC sigma (ucldStdev) from current x=log(r): MLE under mu=-0.5*s^2. */
    private double estimateUCLDStdevFromRates() {
        final int n = rates.getDimension();
        double sumX2 = 0.0;
        for (int i = 0; i < n; i++) {
            final double r = rates.getValue(i);
            if (!(r > 0.0)) return Double.NaN;
            final double x = Math.log(r);
            sumX2 += x * x;
        }
        final double meanX2 = sumX2 / n;
        // v = s^2 = 2( sqrt(1 + mean(x^2)) - 1 )
        final double vHat = 2.0 * (Math.sqrt(1.0 + meanX2) - 1.0);
        return Math.sqrt(Math.max(vHat, minPositiveInput.get()));
    }

    /** Fast estimate of AC sigma2 from increments d = logChi - logPar (MLE-like with mean correction). */
    private double estimateSigma2FromRatesAndTree() {
        final double minDt = minBranchLengthInput.get();
        if (!(minDt > 0.0)) return Double.NaN;

        double Sdt = 0.0;
        double S = 0.0;
        int m = 0;

        final double rootLog = rootLog();

        for (int i = 0; i < nNodes; i++) {
            final Node node = tree.getNode(i);
            if (node.isRoot()) continue;

            final double dt = node.getLength();
            if (!(dt > minDt)) return Double.NaN;

            final int idxChi = idxMap[node.getNr()];
            if (idxChi < 0) return Double.NaN;

            final double rChi = rates.getValue(idxChi);
            if (!(rChi > 0.0)) return Double.NaN;
            final double logChi = Math.log(rChi);

            final Node parent = node.getParent();
            final double logPar;
            if (parent.isRoot()) {
                logPar = rootLog;
            } else {
                final int idxPar = idxMap[parent.getNr()];
                if (idxPar < 0) return Double.NaN;

                final double rPar = rates.getValue(idxPar);
                if (!(rPar > 0.0)) return Double.NaN;
                logPar = Math.log(rPar);
            }

            final double d = logChi - logPar;

            Sdt += dt;
            S += (d * d) / dt;
            m++;
        }

        if (!(Sdt > 0.0) || m <= 0) return Double.NaN;

        // Solve: 0.25*Sdt*v^2 + m*v - S = 0, v>0
        final double disc = (double) m * (double) m + Sdt * S;
        if (!(disc >= 0.0)) return Double.NaN;

        final double vHat = 2.0 * (-m + Math.sqrt(disc)) / Sdt;
        return Math.max(vHat, minPositiveInput.get());
    }

    @Override
    public double proposal() {
        ensureMapsUpToDate();

        final double ucldCenter = estimateUCLDStdevFromRates();
        final double sigma2Center = estimateSigma2FromRatesAndTree();
        if (!Double.isFinite(ucldCenter) || !Double.isFinite(sigma2Center)) return Double.NEGATIVE_INFINITY;

        final double uSd = ucldLogSdInput.get();
        final double s2Sd = sigma2LogSdInput.get();
        if (!(uSd > 0.0) || !(s2Sd > 0.0)) return Double.NEGATIVE_INFINITY;

        final double oldU = ucldStdev.getValue(0);
        final double oldS2 = sigma2.getValue(0);
        if (!(oldU > 0.0) || !(oldS2 > 0.0)) return Double.NEGATIVE_INFINITY;

        // Propose new hyperparameters around centers (independence sampler)
        final double newU = ucldCenter * Math.exp(uSd * Randomizer.nextGaussian());
        final double newS2 = sigma2Center * Math.exp(s2Sd * Randomizer.nextGaussian());

        final double meanLogU = Math.log(ucldCenter);
        final double meanLogS2 = Math.log(sigma2Center);

        final double logQForward =
                logLogNormalPdf(newU, meanLogU, uSd) +
                        logLogNormalPdf(newS2, meanLogS2, s2Sd);

        final double logQReverse =
                logLogNormalPdf(oldU, meanLogU, uSd) +
                        logLogNormalPdf(oldS2, meanLogS2, s2Sd);

        final double logH = logQReverse - logQForward;

        ucldStdev.startEditing(this);
        ucldStdev.setValue(0, newU);

        sigma2.startEditing(this);
        sigma2.setValue(0, newS2);

        indicator.startEditing(this);
        final int k = indicator.getValue(0);
        indicator.setValue(0, (k == 0 ? 1 : 0));

        return logH;
    }
}
