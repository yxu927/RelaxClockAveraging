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

/**
 * Bridge-style cross-model switch for RelaxedRatesPriorSVS.
 *
 * Key idea:
 *  - DO NOT resample the full rawRates vector (keeps likelihood unchanged).
 *  - Flip indicator (0<->1).
 *  - Propose BOTH hyperparameters (ucldStdev and sigma2) so the move is reversible
 *    in the full state space.
 *    * hyper becoming ACTIVE in the proposed model: propose from a pseudo-prior centered
 *      at a moment estimate computed from current rates.
 *    * hyper becoming INACTIVE: propose with a small log-random-walk around its current value.
 */
@Description("Bridge-style switch for SVS indicator inside RelaxedRatesPriorSVS that keeps rawRates fixed "
        + "and proposes hyperparameters (sigma2/ucldStdev) to improve UC<->AC mixing.")
public class IndicatorSwitchBridgeHyperOperator extends Operator {

    public final Input<RelaxedRatesPriorSVS> priorInput = new Input<>(
            "prior",
            "RelaxedRatesPriorSVS instance that defines UC/AC densities and holds references to tree/rates/indicator/hypers.",
            Input.Validate.REQUIRED
    );

    // keep same XML interface as old operator for easy swapping
    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator",
            "SVS indicator (dimension=1). Should match prior.indicatorInput.",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "raw branch rates (dimension=nodeCount-1). Should match prior.ratesInput.",
            Input.Validate.REQUIRED
    );

    // --- tuning knobs ---
    public final Input<Double> ucldPseudoLogSdInput = new Input<>(
            "ucldPseudoLogSd",
            "Log-space SD for pseudo-prior proposal of ucldStdev when UC becomes active.",
            0.35
    );

    public final Input<Double> sigma2PseudoLogSdInput = new Input<>(
            "sigma2PseudoLogSd",
            "Log-space SD for pseudo-prior proposal of sigma2 when AC becomes active.",
            0.75
    );

    public final Input<Double> ucldInactiveLogSdInput = new Input<>(
            "ucldInactiveLogSd",
            "Log-space SD for random-walk update of ucldStdev when UC is inactive (i_relax=1).",
            0.10
    );

    public final Input<Double> sigma2InactiveLogSdInput = new Input<>(
            "sigma2InactiveLogSd",
            "Log-space SD for random-walk update of sigma2 when AC is inactive (i_relax=0).",
            0.10
    );

    public final Input<Double> minHatInput = new Input<>(
            "minHat",
            "Lower bound to avoid log(0) when estimating hyperparameters from current rates.",
            1e-8
    );

    private static final double LOG_2PI = Math.log(2.0 * Math.PI);

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

        tree = prior.treeInput.get();
        ucldStdev = prior.ucldStdevInput.get();
        sigma2 = prior.sigma2Input.get();
        rootLogRate = prior.rootLogRateInput.get(); // optional
        minDt = prior.minBranchLengthInput.get();

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();

        if (indicator.getDimension() != 1) {
            throw new IllegalArgumentException("IndicatorSwitchBridgeHyperOperator: indicator must have dimension=1.");
        }
        if (rates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("IndicatorSwitchBridgeHyperOperator: rates must have dimension nodeCount-1.");
        }

        final double uP = ucldPseudoLogSdInput.get();
        final double sP = sigma2PseudoLogSdInput.get();
        final double uI = ucldInactiveLogSdInput.get();
        final double sI = sigma2InactiveLogSdInput.get();
        if (!(uP > 0 && sP > 0 && uI > 0 && sI > 0)) {
            throw new IllegalArgumentException("All log-SD tuning inputs must be > 0.");
        }

        buildIndexMap();

        // sanity checks: avoid mismatched objects
        if (prior.indicatorInput.get() != null && prior.indicatorInput.get() != indicator) {
            throw new IllegalArgumentException("IndicatorSwitchBridgeHyperOperator: 'indicator' does not match prior.indicatorInput.");
        }
        if (prior.ratesInput.get() != null && prior.ratesInput.get() != rates) {
            throw new IllegalArgumentException("IndicatorSwitchBridgeHyperOperator: 'rates' does not match prior.ratesInput.");
        }
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

        // moment estimates from CURRENT rates (rates are NOT changed here)
        final double minHat = minHatInput.get();
        final double uHat = Math.max(minHat, estimateUcldStdevFromRates());
        final double s2Hat = Math.max(minHat, estimateSigma2FromRatesAndTree());

        final double logUHat = Math.log(uHat);
        final double logS2Hat = Math.log(s2Hat);

        // current hyper values
        final double uOld = ucldStdev.getValue();
        final double s2Old = sigma2.getValue();
        if (!(uOld > 0.0) || !(s2Old > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        // draw proposals:
        // newK=0 (UC active): ucldStdev from pseudo-prior; sigma2 inactive RW around s2Old
        // newK=1 (AC active): sigma2 from pseudo-prior; ucldStdev inactive RW around uOld
        final double uNew;
        final double s2New;

        if (newK == 0) {
            uNew = drawLogNormalAround(logUHat, ucldPseudoLogSdInput.get());
            s2New = drawLogNormalAround(Math.log(s2Old), sigma2InactiveLogSdInput.get());
        } else {
            s2New = drawLogNormalAround(logS2Hat, sigma2PseudoLogSdInput.get());
            uNew = drawLogNormalAround(Math.log(uOld), ucldInactiveLogSdInput.get());
        }

        if (!(uNew > 0.0) || !(s2New > 0.0) || Double.isInfinite(uNew) || Double.isInfinite(s2New) || Double.isNaN(uNew) || Double.isNaN(s2New)) {
            return Double.NEGATIVE_INFINITY;
        }

        // Forward density q(new|old)
        final double logqForwardU;
        final double logqForwardS2;
        if (newK == 0) {
            logqForwardU  = logLogNormalDensity(uNew,  logUHat,  ucldPseudoLogSdInput.get());
            logqForwardS2 = logLogNormalDensity(s2New, Math.log(s2Old), sigma2InactiveLogSdInput.get());
        } else {
            logqForwardS2 = logLogNormalDensity(s2New, logS2Hat, sigma2PseudoLogSdInput.get());
            logqForwardU  = logLogNormalDensity(uNew,  Math.log(uOld), ucldInactiveLogSdInput.get());
        }
        if (Double.isInfinite(logqForwardU) || Double.isInfinite(logqForwardS2)) {
            return Double.NEGATIVE_INFINITY;
        }

        // Reverse density q(old|new) (reverse target is oldK)
        final double logqReverseU;
        final double logqReverseS2;
        if (oldK == 0) {
            // reverse target UC active: ucldStdev pseudo-prior; sigma2 inactive RW around s2New
            logqReverseU  = logLogNormalDensity(uOld,  logUHat, ucldPseudoLogSdInput.get());
            logqReverseS2 = logLogNormalDensity(s2Old, Math.log(s2New), sigma2InactiveLogSdInput.get());
        } else {
            // reverse target AC active: sigma2 pseudo-prior; ucldStdev inactive RW around uNew
            logqReverseS2 = logLogNormalDensity(s2Old, logS2Hat, sigma2PseudoLogSdInput.get());
            logqReverseU  = logLogNormalDensity(uOld,  Math.log(uNew), ucldInactiveLogSdInput.get());
        }
        if (Double.isInfinite(logqReverseU) || Double.isInfinite(logqReverseS2)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double logHR = (logqReverseU + logqReverseS2) - (logqForwardU + logqForwardS2);

        // apply changes to the state
        indicator.startEditing(this);
        sigma2.startEditing(this);
        ucldStdev.startEditing(this);

        indicator.setValue(0, newK);
        sigma2.setValue(0, s2New);
        ucldStdev.setValue(0, uNew);

        return logHR;
    }

    private double drawLogNormalAround(final double muLog, final double sdLog) {
        final double z = Randomizer.nextGaussian();
        return Math.exp(muLog + sdLog * z);
    }

    private static double logLogNormalDensity(final double x, final double muLog, final double sdLog) {
        if (!(x > 0.0) || !(sdLog > 0.0) || Double.isNaN(x) || Double.isInfinite(x)) {
            return Double.NEGATIVE_INFINITY;
        }
        final double lx = Math.log(x);
        final double z = (lx - muLog) / sdLog;
        return -Math.log(x) - Math.log(sdLog) - 0.5 * LOG_2PI - 0.5 * z * z;
    }

    // UC: ucldStdev ~ SD(log rates)
    private double estimateUcldStdevFromRates() {
        final int dim = rates.getDimension();
        double sum = 0.0;
        double sumSq = 0.0;
        for (int i = 0; i < dim; i++) {
            final double r = rates.getValue(i);
            if (!(r > 0.0) || Double.isInfinite(r) || Double.isNaN(r)) {
                continue;
            }
            final double lr = Math.log(r);
            sum += lr;
            sumSq += lr * lr;
        }
        final double n = (double) dim;
        final double mean = sum / n;
        final double var = Math.max(0.0, (sumSq / n) - mean * mean);
        return Math.sqrt(var);
    }

    // AC: sigma2 ~ sum(delta^2)/sum(dt), delta=log(r_child)-log(r_parent)
    private double estimateSigma2FromRatesAndTree() {
        final double rootLog = (rootLogRate == null ? 0.0 : rootLogRate.getValue());

        double sumDelta2 = 0.0;
        double sumDt = 0.0;

        for (int i = 0; i < nNodes; i++) {
            final Node child = tree.getNode(i);
            if (child.getNr() == rootNr) {
                continue;
            }
            final Node parent = child.getParent();
            if (parent == null) {
                continue;
            }

            final double dt = child.getLength();
            if (!(dt > 0.0) || Double.isInfinite(dt) || Double.isNaN(dt)) {
                continue;
            }
            if (minDt > 0.0 && dt <= minDt) {
                continue;
            }

            final double logPar;
            if (parent.getNr() == rootNr) {
                logPar = rootLog;
            } else {
                final int pIdx = idxMap[parent.getNr()];
                if (pIdx < 0) {
                    continue;
                }
                final double rPar = rates.getValue(pIdx);
                if (!(rPar > 0.0) || Double.isInfinite(rPar) || Double.isNaN(rPar)) {
                    continue;
                }
                logPar = Math.log(rPar);
            }

            final int cIdx = idxMap[child.getNr()];
            if (cIdx < 0) {
                continue;
            }
            final double rChi = rates.getValue(cIdx);
            if (!(rChi > 0.0) || Double.isInfinite(rChi) || Double.isNaN(rChi)) {
                continue;
            }
            final double logChi = Math.log(rChi);

            final double delta = logChi - logPar;
            sumDelta2 += delta * delta;
            sumDt += dt;
        }

        if (!(sumDt > 0.0)) {
            return 0.0;
        }
        return sumDelta2 / sumDt;
    }
}
