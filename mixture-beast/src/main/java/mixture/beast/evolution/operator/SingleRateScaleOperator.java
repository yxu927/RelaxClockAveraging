package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;

@Description("Scale a single randomly chosen element of a positive rate vector: r <- r * exp(eps). "
        + "This is typically better than additive random-walk for positive parameters.")
public class SingleRateScaleOperator extends Operator {

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "Legacy positive rate vector parameter. Kept for backwards compatibility.",
            Input.Validate.OPTIONAL
    );

    public final Input<RealVectorParam<?>> ratesVectorInput = new Input<>(
            "ratesVector",
            "BEAST3 typed mutable positive rate vector parameter.",
            Input.Validate.OPTIONAL
    );

    public final Input<Double> windowInput = new Input<>(
            "window",
            "eps ~ Uniform(-window, +window) on log scale",
            0.3
    );

    private RealParameter legacyRates;
    private RealVectorParam<?> typedRates;

    @Override
    public void initAndValidate() {
        legacyRates = ratesInput.get();
        typedRates = ratesVectorInput.get();

        if (legacyRates == null && typedRates == null) {
            throw new IllegalArgumentException("SingleRateScaleOperator: either rates or ratesVector must be specified.");
        }
        if (legacyRates != null && typedRates != null) {
            throw new IllegalArgumentException("SingleRateScaleOperator: specify only one of rates or ratesVector.");
        }
        if (rateDimension() < 1) {
            throw new IllegalArgumentException("SingleRateScaleOperator: rates dimension must be >= 1");
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

    @Override
    public double proposal() {
        final double window = windowInput.get();
        if (!(window > 0.0)) return Double.NEGATIVE_INFINITY;

        final int dim = rateDimension();
        final int i = Randomizer.nextInt(dim);

        final double r = rateValue(i);
        if (!(r > 0.0)) return Double.NEGATIVE_INFINITY;

        final double eps = (Randomizer.nextDouble() * 2.0 - 1.0) * window;
        final double m = Math.exp(eps);

        final double rNew = r * m;
        if (!(rNew > 0.0) || Double.isNaN(rNew) || Double.isInfinite(rNew)) {
            return Double.NEGATIVE_INFINITY;
        }

        setRateValue(i, rNew);

        // scale proposal Hastings: log(q(old|new)/q(new|old)) = eps
        return eps;
    }
}
