package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import mixture.beast.evolution.util.BranchRateIndexHelper;

@Description("AC-only non-centered hyper move: changes sigma2 while keeping the latent AC increments fixed, "
        + "and reconstructs the entire shared rate vector accordingly. "
        + "Useful when the bottleneck is within-relax mixing under the autocorrelated branch of the mixture.")
public class ACSigma2NonCenteredOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "shared positive branch rates (non-root nodes)",
            Input.Validate.REQUIRED
    );

    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator",
            "0=UC, 1=AC",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "AC Brownian variance per unit time",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "optional root log-rate anchor (default 0)",
            Input.Validate.OPTIONAL
    );

    public final Input<Double> minBranchLengthInput = new Input<>(
            "minBranchLength",
            "minimum branch length allowed for AC mapping",
            1e-12
    );

    public final Input<Double> windowInput = new Input<>(
            "window",
            "log-scale proposal window: eps ~ Uniform(-window, +window)",
            0.15
    );

    /** If true (default), reject move when indicator!=1. For mixture runs, prefer autoOptimize=false in XML. */
    public final Input<Boolean> rejectIfNotACInput = new Input<>(
            "rejectIfNotAC",
            "reject move when indicator!=1",
            true
    );

    private Tree tree;
    private RealParameter rates;
    private IntegerParameter indicator;
    private RealParameter sigma2;
    private RealParameter rootLogRate;

    private BranchRateIndexHelper.Mapping mapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();
        indicator = indicatorInput.get();
        sigma2 = sigma2Input.get();
        rootLogRate = rootLogRateInput.get();

        if (indicator.getDimension() != 1) {
            throw new IllegalArgumentException("indicator dimension must be 1");
        }
        if (sigma2.getDimension() != 1) {
            throw new IllegalArgumentException("sigma2 dimension must be 1");
        }
        if (rootLogRate != null && rootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate dimension must be 1");
        }

        BranchRateIndexHelper.validateRatesDimension(tree, rates, "ACSigma2NonCenteredOperator");
        mapping = BranchRateIndexHelper.buildDeterministic(tree);
    }

    private void ensureMappingUpToDate() {
        mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, rates, "ACSigma2NonCenteredOperator");
    }

    private double rootLog() {
        return (rootLogRate == null ? 0.0 : rootLogRate.getValue(0));
    }

    @Override
    public double proposal() {
        ensureMappingUpToDate();

        final int k = indicator.getValue(0);
        if (k != 1) {
            return rejectIfNotACInput.get() ? Double.NEGATIVE_INFINITY : 0.0;
        }

        final double minDt = minBranchLengthInput.get();
        final double window = windowInput.get();
        if (!(minDt > 0.0) || !(window > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double oldSigma2 = sigma2.getValue(0);
        if (!(oldSigma2 > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double eps = (Randomizer.nextDouble() * 2.0 - 1.0) * window;
        final double newSigma2 = oldSigma2 * Math.exp(eps);
        if (!(newSigma2 > 0.0) || Double.isInfinite(newSigma2) || Double.isNaN(newSigma2)) {
            return Double.NEGATIVE_INFINITY;
        }

        final int nEdges = rates.getDimension();
        final double[] xOld = new double[nEdges];
        final double[] u = new double[nEdges];
        final double[] xNew = new double[nEdges];

        for (int i = 0; i < nEdges; i++) {
            final double r = rates.getValue(i);
            if (!(r > 0.0)) {
                return Double.NEGATIVE_INFINITY;
            }
            xOld[i] = Math.log(r);
        }

        try {
            fillUFromLogRatesAC(tree.getRoot(), rootLog(), xOld, u, oldSigma2, minDt);
            fillLogRatesFromU(tree.getRoot(), rootLog(), u, xNew, newSigma2, minDt);
        } catch (ArithmeticException bad) {
            return Double.NEGATIVE_INFINITY;
        }

        double sumDelta = 0.0;
        for (int i = 0; i < nEdges; i++) {
            if (!Double.isFinite(xNew[i])) {
                return Double.NEGATIVE_INFINITY;
            }
            sumDelta += (xNew[i] - xOld[i]);
        }

        rates.startEditing(this);
        for (int i = 0; i < nEdges; i++) {
            rates.setValue(i, Math.exp(xNew[i]));
        }

        sigma2.startEditing(this);
        sigma2.setValue(0, newSigma2);

        return sumDelta + (0.5 * nEdges + 1.0) * eps;
    }

    private void fillUFromLogRatesAC(final Node parent,
                                     final double logPar,
                                     final double[] x,
                                     final double[] uOut,
                                     final double sigma2Value,
                                     final double minDt) {
        final int cc = parent.getChildCount();
        for (int c = 0; c < cc; c++) {
            final Node child = parent.getChild(c);
            final int idx = mapping.idxForNode(child);
            if (idx < 0) {
                throw new ArithmeticException("Bad mapping for child");
            }

            final double dt = child.getLength();
            if (!(dt > minDt)) {
                throw new ArithmeticException("dt too small");
            }

            final double var = sigma2Value * dt;
            if (!(var > 0.0)) {
                throw new ArithmeticException("var<=0");
            }

            final double mean = logPar - 0.5 * var;
            uOut[idx] = (x[idx] - mean) / Math.sqrt(var);

            fillUFromLogRatesAC(child, x[idx], x, uOut, sigma2Value, minDt);
        }
    }

    private void fillLogRatesFromU(final Node parent,
                                   final double logPar,
                                   final double[] u,
                                   final double[] xOut,
                                   final double sigma2Value,
                                   final double minDt) {
        final int cc = parent.getChildCount();
        for (int c = 0; c < cc; c++) {
            final Node child = parent.getChild(c);
            final int idx = mapping.idxForNode(child);
            if (idx < 0) {
                throw new ArithmeticException("Bad mapping for child");
            }

            final double dt = child.getLength();
            if (!(dt > minDt)) {
                throw new ArithmeticException("dt too small");
            }

            final double var = sigma2Value * dt;
            if (!(var > 0.0)) {
                throw new ArithmeticException("var<=0");
            }

            final double logChild = logPar - 0.5 * var + Math.sqrt(var) * u[idx];
            xOut[idx] = logChild;

            fillLogRatesFromU(child, logChild, u, xOut, sigma2Value, minDt);
        }
    }
}
