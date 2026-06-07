package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;

@Description("UC-only non-centered hyper move: changes ucldStdev while keeping latent z_i fixed, "
        + "and reconstructs the shared rate vector accordingly. "
        + "Useful when the bottleneck is within-relax mixing under the UCLD branch of the mixture.")
public class UCLDStdevNonCenteredOperator extends Operator {

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

    public final Input<RealParameter> ucldStdevInput = new Input<>(
            "ucldStdev",
            "Legacy UC lognormal stdev (sigma on log scale).",
            Input.Validate.OPTIONAL
    );

    public final Input<RealScalarParam<?>> ucldStdevScalarInput = new Input<>(
            "ucldStdevScalar",
            "BEAST3 typed mutable UC lognormal stdev (sigma on log scale).",
            Input.Validate.OPTIONAL
    );

    public final Input<Double> windowInput = new Input<>(
            "window",
            "log-scale proposal window: eps ~ Uniform(-window, +window)",
            0.2
    );

    /** If true (default), reject move when indicator!=0. For mixture runs, prefer autoOptimize=false in XML. */
    public final Input<Boolean> rejectIfNotUCInput = new Input<>(
            "rejectIfNotUC",
            "reject move when indicator!=0",
            true
    );

    private RealParameter legacyRates;
    private RealVectorParam<?> typedRates;
    private IntegerParameter legacyIndicator;
    private IntScalarParam<?> typedIndicator;
    private RealParameter legacyUcldStdev;
    private RealScalarParam<?> typedUcldStdev;

    @Override
    public void initAndValidate() {
        legacyRates = ratesInput.get();
        typedRates = ratesVectorInput.get();
        legacyIndicator = indicatorInput.get();
        typedIndicator = indicatorScalarInput.get();
        legacyUcldStdev = ucldStdevInput.get();
        typedUcldStdev = ucldStdevScalarInput.get();

        requireExactlyOne(legacyRates, typedRates, "rates", "ratesVector");
        requireExactlyOne(legacyIndicator, typedIndicator, "indicator", "indicatorScalar");
        requireExactlyOne(legacyUcldStdev, typedUcldStdev, "ucldStdev", "ucldStdevScalar");

        if (rateDimension() < 1) {
            throw new IllegalArgumentException("UCLDStdevNonCenteredOperator: rates dimension must be >= 1");
        }
        if (legacyIndicator != null && legacyIndicator.getDimension() != 1) {
            throw new IllegalArgumentException("UCLDStdevNonCenteredOperator: indicator dimension must be 1");
        }
        if (legacyUcldStdev != null && legacyUcldStdev.getDimension() != 1) {
            throw new IllegalArgumentException("UCLDStdevNonCenteredOperator: ucldStdev dimension must be 1");
        }
    }

    private static void requireExactlyOne(final Object legacy,
                                          final Object typed,
                                          final String legacyName,
                                          final String typedName) {
        if (legacy == null && typed == null) {
            throw new IllegalArgumentException("UCLDStdevNonCenteredOperator: either "
                    + legacyName + " or " + typedName + " must be specified.");
        }
        if (legacy != null && typed != null) {
            throw new IllegalArgumentException("UCLDStdevNonCenteredOperator: specify only one of "
                    + legacyName + " or " + typedName + ".");
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

    private double ucldStdevValue() {
        return legacyUcldStdev != null ? legacyUcldStdev.getValue(0) : typedUcldStdev.get();
    }

    private void setUcldStdevValue(final double value) {
        if (legacyUcldStdev != null) {
            legacyUcldStdev.setValue(0, value);
        } else {
            typedUcldStdev.set(value);
        }
    }

    @Override
    public double proposal() {
        final int k = indicatorValue();
        if (k != 0) {
            return rejectIfNotUCInput.get() ? Double.NEGATIVE_INFINITY : 0.0;
        }

        final double window = windowInput.get();
        if (!(window > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double oldS = ucldStdevValue();
        if (!(oldS > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double eps = (Randomizer.nextDouble() * 2.0 - 1.0) * window;
        final double newS = oldS * Math.exp(eps);
        if (!(newS > 0.0) || Double.isInfinite(newS) || Double.isNaN(newS)) {
            return Double.NEGATIVE_INFINITY;
        }

        final int nEdges = rateDimension();
        final double oldVar = oldS * oldS;
        final double newVar = newS * newS;

        final double[] rNew = new double[nEdges];
        double sumDelta = 0.0;

        for (int i = 0; i < nEdges; i++) {
            final double r = rateValue(i);
            if (!(r > 0.0)) {
                return Double.NEGATIVE_INFINITY;
            }

            final double xOld = Math.log(r);
            final double z = (xOld + 0.5 * oldVar) / oldS;
            final double xNew = -0.5 * newVar + newS * z;

            final double rCandidate = Math.exp(xNew);
            if (!(rCandidate > 0.0) || Double.isInfinite(rCandidate) || Double.isNaN(rCandidate)) {
                return Double.NEGATIVE_INFINITY;
            }

            rNew[i] = rCandidate;
            sumDelta += (xNew - xOld);
        }

        if (legacyRates != null) {
            legacyRates.startEditing(this);
        }
        for (int i = 0; i < nEdges; i++) {
            setRateValue(i, rNew[i]);
        }

        if (legacyUcldStdev != null) {
            legacyUcldStdev.startEditing(this);
        }
        setUcldStdevValue(newS);

        return sumDelta + (nEdges + 1.0) * eps;
    }
}
