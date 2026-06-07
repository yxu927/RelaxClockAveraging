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
import mixture.beast.evolution.util.BranchRateIndexHelper;

@Description("UC<->AC switch operator that deterministically maps the shared rate vector through "
        + "a latent standard-normal vector u (non-centered bridge). "
        + "Greatly improves mixing of the indicator on large trees.")
public class UCACSwitchBridgeOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);
    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "Legacy shared positive branch rates (non-root nodes).",
            Input.Validate.OPTIONAL);
    public final Input<RealVectorParam<?>> ratesVectorInput = new Input<>(
            "ratesVector",
            "BEAST3 typed mutable shared positive branch rates (non-root nodes).",
            Input.Validate.OPTIONAL);
    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator",
            "Legacy scalar indicator; 0=UC, 1=AC.",
            Input.Validate.OPTIONAL);
    public final Input<IntScalarParam<?>> indicatorScalarInput = new Input<>(
            "indicatorScalar",
            "BEAST3 typed mutable scalar indicator; 0=UC, 1=AC.",
            Input.Validate.OPTIONAL);

    public final Input<RealParameter> ucldStdevInput = new Input<>(
            "ucldStdev",
            "Legacy UC lognormal stdev (sigma on log scale).",
            Input.Validate.OPTIONAL);
    public final Input<RealScalarParam<?>> ucldStdevScalarInput = new Input<>(
            "ucldStdevScalar",
            "BEAST3 typed mutable UC lognormal stdev (sigma on log scale).",
            Input.Validate.OPTIONAL);
    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "Legacy AC Brownian variance per unit time.",
            Input.Validate.OPTIONAL);
    public final Input<RealScalarParam<?>> sigma2ScalarInput = new Input<>(
            "sigma2Scalar",
            "BEAST3 typed mutable AC Brownian variance per unit time.",
            Input.Validate.OPTIONAL);
    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "Legacy optional root log-rate anchor (default 0).",
            Input.Validate.OPTIONAL);
    public final Input<RealScalarParam<?>> rootLogRateScalarInput = new Input<>(
            "rootLogRateScalar",
            "BEAST3 typed optional root log-rate anchor (default 0).",
            Input.Validate.OPTIONAL);

    public final Input<Double> minBranchLengthInput = new Input<>("minBranchLength", "min dt allowed for AC mapping", 1e-12);

    private Tree tree;
    private RealParameter legacyRates;
    private RealVectorParam<?> typedRates;
    private IntegerParameter legacyIndicator;
    private IntScalarParam<?> typedIndicator;
    private RealParameter legacyUcldStdev;
    private RealScalarParam<?> typedUcldStdev;
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
        legacyUcldStdev = ucldStdevInput.get();
        typedUcldStdev = ucldStdevScalarInput.get();
        legacySigma2 = sigma2Input.get();
        typedSigma2 = sigma2ScalarInput.get();
        legacyRootLogRate = rootLogRateInput.get();
        typedRootLogRate = rootLogRateScalarInput.get();

        requireExactlyOne(legacyRates, typedRates, "rates", "ratesVector");
        requireExactlyOne(legacyIndicator, typedIndicator, "indicator", "indicatorScalar");
        requireExactlyOne(legacyUcldStdev, typedUcldStdev, "ucldStdev", "ucldStdevScalar");
        requireExactlyOne(legacySigma2, typedSigma2, "sigma2", "sigma2Scalar");
        requireAtMostOne(legacyRootLogRate, typedRootLogRate, "rootLogRate", "rootLogRateScalar");

        if (legacyIndicator != null && legacyIndicator.getDimension() != 1) {
            throw new IllegalArgumentException("UCACSwitchBridgeOperator: indicator dimension must be 1");
        }
        if (legacyUcldStdev != null && legacyUcldStdev.getDimension() != 1) {
            throw new IllegalArgumentException("UCACSwitchBridgeOperator: ucldStdev dimension must be 1");
        }
        if (legacySigma2 != null && legacySigma2.getDimension() != 1) {
            throw new IllegalArgumentException("UCACSwitchBridgeOperator: sigma2 dimension must be 1");
        }
        if (legacyRootLogRate != null && legacyRootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("UCACSwitchBridgeOperator: rootLogRate dimension must be 1");
        }

        validateOrExpandRatesDimension();
        mapping = BranchRateIndexHelper.buildDeterministic(tree);
    }

    private static void requireExactlyOne(final Object legacy,
                                          final Object typed,
                                          final String legacyName,
                                          final String typedName) {
        if (legacy == null && typed == null) {
            throw new IllegalArgumentException("UCACSwitchBridgeOperator: either "
                    + legacyName + " or " + typedName + " must be specified.");
        }
        requireAtMostOne(legacy, typed, legacyName, typedName);
    }

    private static void requireAtMostOne(final Object legacy,
                                         final Object typed,
                                         final String legacyName,
                                         final String typedName) {
        if (legacy != null && typed != null) {
            throw new IllegalArgumentException("UCACSwitchBridgeOperator: specify only one of "
                    + legacyName + " or " + typedName + ".");
        }
    }

    private void ensureMappingUpToDate() {
        if (legacyRates != null) {
            mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, legacyRates, "UCACSwitchBridgeOperator");
        } else {
            mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, typedRates, "UCACSwitchBridgeOperator");
        }
    }

    private void validateOrExpandRatesDimension() {
        if (legacyRates != null) {
            BranchRateIndexHelper.validateRatesDimension(tree, legacyRates, "UCACSwitchBridgeOperator");
        } else {
            BranchRateIndexHelper.validateRatesDimension(tree, typedRates, "UCACSwitchBridgeOperator");
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

    private void setIndicatorValue(final int value) {
        if (legacyIndicator != null) {
            legacyIndicator.setValue(0, value);
        } else {
            typedIndicator.set(value);
        }
    }

    private double ucldStdevValue() {
        return legacyUcldStdev != null ? legacyUcldStdev.getValue(0) : typedUcldStdev.get();
    }

    private double sigma2Value() {
        return legacySigma2 != null ? legacySigma2.getValue(0) : typedSigma2.get();
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

    private static final class Sum {
        double sumLogVar = 0.0;
    }

    @Override
    public double proposal() {
        ensureMappingUpToDate();

        final int k = indicatorValue();
        if (k == 0) {
            return proposeUCtoAC();
        } else if (k == 1) {
            return proposeACtoUC();
        }
        return Double.NEGATIVE_INFINITY;
    }

    private double proposeUCtoAC() {
        final double s = ucldStdevValue();
        final double sig2 = sigma2Value();
        final double minDt = minBranchLengthInput.get();
        if (!(s > 0.0) || !(sig2 > 0.0) || !(minDt > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final int nEdges = rateDimension();
        final double muUC = -0.5 * s * s;

        final double[] xOld = new double[nEdges];
        final double[] u = new double[nEdges];

        for (int i = 0; i < nEdges; i++) {
            final double r = rateValue(i);
            if (!(r > 0.0)) {
                return Double.NEGATIVE_INFINITY;
            }
            final double x = Math.log(r);
            xOld[i] = x;
            u[i] = (x - muUC) / s;
        }

        final double[] xNew = new double[nEdges];
        final Sum sum = new Sum();

        try {
            fillLogRatesFromU(tree.getRoot(), rootLog(), u, xNew, sig2, minDt, sum);
        } catch (ArithmeticException bad) {
            return Double.NEGATIVE_INFINITY;
        }

        double sumDelta = 0.0;
        for (int i = 0; i < nEdges; i++) {
            sumDelta += (xNew[i] - xOld[i]);
        }
        final double logH = sumDelta + 0.5 * sum.sumLogVar - nEdges * Math.log(s);

        if (legacyRates != null) {
            legacyRates.startEditing(this);
        }
        for (int i = 0; i < nEdges; i++) {
            setRateValue(i, Math.exp(xNew[i]));
        }

        if (legacyIndicator != null) {
            legacyIndicator.startEditing(this);
        }
        setIndicatorValue(1);

        return logH;
    }

    private double proposeACtoUC() {
        final double s = ucldStdevValue();
        final double sig2 = sigma2Value();
        final double minDt = minBranchLengthInput.get();
        if (!(s > 0.0) || !(sig2 > 0.0) || !(minDt > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final int nEdges = rateDimension();
        final double muUC = -0.5 * s * s;

        final double[] xOld = new double[nEdges];
        for (int i = 0; i < nEdges; i++) {
            final double r = rateValue(i);
            if (!(r > 0.0)) {
                return Double.NEGATIVE_INFINITY;
            }
            xOld[i] = Math.log(r);
        }

        final double[] u = new double[nEdges];
        final Sum sum = new Sum();

        try {
            fillUFromLogRatesAC(tree.getRoot(), rootLog(), xOld, u, sig2, minDt, sum);
        } catch (ArithmeticException bad) {
            return Double.NEGATIVE_INFINITY;
        }

        final double[] xNew = new double[nEdges];
        for (int i = 0; i < nEdges; i++) {
            xNew[i] = muUC + s * u[i];
        }

        double sumDelta = 0.0;
        for (int i = 0; i < nEdges; i++) {
            sumDelta += (xNew[i] - xOld[i]);
        }
        final double logH = sumDelta + nEdges * Math.log(s) - 0.5 * sum.sumLogVar;

        if (legacyRates != null) {
            legacyRates.startEditing(this);
        }
        for (int i = 0; i < nEdges; i++) {
            setRateValue(i, Math.exp(xNew[i]));
        }

        if (legacyIndicator != null) {
            legacyIndicator.startEditing(this);
        }
        setIndicatorValue(0);

        return logH;
    }

    private void fillLogRatesFromU(final Node parent,
                                   final double logPar,
                                   final double[] u,
                                   final double[] xOut,
                                   final double sigma2Value,
                                   final double minDt,
                                   final Sum sum) {
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

            sum.sumLogVar += Math.log(var);

            final double logChild = logPar - 0.5 * var + Math.sqrt(var) * u[idx];
            xOut[idx] = logChild;

            fillLogRatesFromU(child, logChild, u, xOut, sigma2Value, minDt, sum);
        }
    }

    private void fillUFromLogRatesAC(final Node parent,
                                     final double logPar,
                                     final double[] x,
                                     final double[] uOut,
                                     final double sigma2Value,
                                     final double minDt,
                                     final Sum sum) {
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

            sum.sumLogVar += Math.log(var);

            final double mean = logPar - 0.5 * var;
            uOut[idx] = (x[idx] - mean) / Math.sqrt(var);

            fillUFromLogRatesAC(child, x[idx], x, uOut, sigma2Value, minDt, sum);
        }
    }
}
