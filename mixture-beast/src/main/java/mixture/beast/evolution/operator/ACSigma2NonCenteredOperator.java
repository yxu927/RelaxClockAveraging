package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import mixture.beast.evolution.util.BranchRateIndexHelper;

@Description("AC-only non-centered hyper move: changes sigma2 while keeping the latent AC increments fixed, "
        + "and reconstructs the entire shared rate vector accordingly. "
        + "Useful when the bottleneck is within-relax mixing under the autocorrelated branch of the mixture.")
public class ACSigma2NonCenteredOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "Legacy shared positive branch rates (non-root nodes).",
            Input.Validate.OPTIONAL
    );

    public final Input<RealVectorParam<?>> ratesVectorInput = new Input<>(
            "ratesVector",
            "BEAST3 typed mutable shared positive branch rates (non-root nodes).",
            Input.Validate.OPTIONAL
    );

    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator",
            "Legacy scalar indicator; 0=UC, 1=AC.",
            Input.Validate.OPTIONAL
    );

    public final Input<IntScalarParam<?>> indicatorScalarInput = new Input<>(
            "indicatorScalar",
            "BEAST3 typed mutable scalar indicator; 0=UC, 1=AC.",
            Input.Validate.OPTIONAL
    );

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "Legacy AC Brownian variance per unit time.",
            Input.Validate.OPTIONAL
    );

    public final Input<RealScalarParam<?>> sigma2ScalarInput = new Input<>(
            "sigma2Scalar",
            "BEAST3 typed mutable AC Brownian variance per unit time.",
            Input.Validate.OPTIONAL
    );

    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "Legacy optional root log-rate anchor (default 0).",
            Input.Validate.OPTIONAL
    );

    public final Input<RealScalarParam<?>> rootLogRateScalarInput = new Input<>(
            "rootLogRateScalar",
            "BEAST3 typed optional root log-rate anchor (default 0).",
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
    private RealParameter legacyRates;
    private RealVectorParam<?> typedRates;
    private IntegerParameter legacyIndicator;
    private IntScalarParam<?> typedIndicator;
    private RealParameter legacySigma2;
    private RealScalarParam<?> typedSigma2;
    private RealParameter legacyRootLogRate;
    private RealScalarParam<?> typedRootLogRate;

    private BranchRateIndexHelper.Mapping mapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        legacyRates = ratesInput.get();
        typedRates = ratesVectorInput.get();
        legacyIndicator = indicatorInput.get();
        typedIndicator = indicatorScalarInput.get();
        legacySigma2 = sigma2Input.get();
        typedSigma2 = sigma2ScalarInput.get();
        legacyRootLogRate = rootLogRateInput.get();
        typedRootLogRate = rootLogRateScalarInput.get();

        requireExactlyOne(legacyRates, typedRates, "rates", "ratesVector");
        requireExactlyOne(legacyIndicator, typedIndicator, "indicator", "indicatorScalar");
        requireExactlyOne(legacySigma2, typedSigma2, "sigma2", "sigma2Scalar");
        requireAtMostOne(legacyRootLogRate, typedRootLogRate, "rootLogRate", "rootLogRateScalar");

        if (legacyIndicator != null && legacyIndicator.getDimension() != 1) {
            throw new IllegalArgumentException("ACSigma2NonCenteredOperator: indicator dimension must be 1");
        }
        if (legacySigma2 != null && legacySigma2.getDimension() != 1) {
            throw new IllegalArgumentException("ACSigma2NonCenteredOperator: sigma2 dimension must be 1");
        }
        if (legacyRootLogRate != null && legacyRootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("ACSigma2NonCenteredOperator: rootLogRate dimension must be 1");
        }

        validateOrExpandRatesDimension();
        mapping = BranchRateIndexHelper.buildDeterministic(tree);
    }

    private static void requireExactlyOne(final Object legacy,
                                          final Object typed,
                                          final String legacyName,
                                          final String typedName) {
        if (legacy == null && typed == null) {
            throw new IllegalArgumentException("ACSigma2NonCenteredOperator: either "
                    + legacyName + " or " + typedName + " must be specified.");
        }
        requireAtMostOne(legacy, typed, legacyName, typedName);
    }

    private static void requireAtMostOne(final Object legacy,
                                         final Object typed,
                                         final String legacyName,
                                         final String typedName) {
        if (legacy != null && typed != null) {
            throw new IllegalArgumentException("ACSigma2NonCenteredOperator: specify only one of "
                    + legacyName + " or " + typedName + ".");
        }
    }

    private void ensureMappingUpToDate() {
        if (legacyRates != null) {
            mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, legacyRates, "ACSigma2NonCenteredOperator");
        } else {
            mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, typedRates, "ACSigma2NonCenteredOperator");
        }
    }

    private void validateOrExpandRatesDimension() {
        if (legacyRates != null) {
            BranchRateIndexHelper.validateRatesDimension(tree, legacyRates, "ACSigma2NonCenteredOperator");
        } else {
            BranchRateIndexHelper.validateRatesDimension(tree, typedRates, "ACSigma2NonCenteredOperator");
        }
    }

    private int rateDimension() {
        return legacyRates != null ? legacyRates.getDimension() : typedRates.size();
    }

    private double rateValue(final int i) {
        return legacyRates != null ? legacyRates.getValue(i) : typedRates.get(i);
    }

    private void setRateValue(final int i, final double value) {
        if (legacyRates != null) {
            legacyRates.setValue(i, value);
        } else {
            typedRates.set(i, value);
        }
    }

    private int indicatorValue() {
        return legacyIndicator != null ? legacyIndicator.getValue(0) : typedIndicator.get();
    }

    private double sigma2Value() {
        return legacySigma2 != null ? legacySigma2.getValue(0) : typedSigma2.get();
    }

    private void setSigma2Value(final double value) {
        if (legacySigma2 != null) {
            legacySigma2.setValue(0, value);
        } else {
            typedSigma2.set(value);
        }
    }

    private double rootLog() {
        if (legacyRootLogRate != null) {
            return legacyRootLogRate.getValue(0);
        }
        if (typedRootLogRate != null) {
            return typedRootLogRate.get();
        }
        return 0.0;
    }

    @Override
    public double proposal() {
        ensureMappingUpToDate();

        final int k = indicatorValue();
        if (k != 1) {
            return rejectIfNotACInput.get() ? Double.NEGATIVE_INFINITY : 0.0;
        }

        final double minDt = minBranchLengthInput.get();
        final double window = windowInput.get();
        if (!(minDt > 0.0) || !(window > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double oldSigma2 = sigma2Value();
        if (!(oldSigma2 > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double eps = (Randomizer.nextDouble() * 2.0 - 1.0) * window;
        final double newSigma2 = oldSigma2 * Math.exp(eps);
        if (!(newSigma2 > 0.0) || Double.isInfinite(newSigma2) || Double.isNaN(newSigma2)) {
            return Double.NEGATIVE_INFINITY;
        }

        final int nEdges = rateDimension();
        final double[] xOld = new double[nEdges];
        final double[] u = new double[nEdges];
        final double[] xNew = new double[nEdges];

        for (int i = 0; i < nEdges; i++) {
            final double r = rateValue(i);
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

        if (legacyRates != null) {
            legacyRates.startEditing(this);
        }
        for (int i = 0; i < nEdges; i++) {
            setRateValue(i, Math.exp(xNew[i]));
        }

        if (legacySigma2 != null) {
            legacySigma2.startEditing(this);
        }
        setSigma2Value(newSigma2);

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
