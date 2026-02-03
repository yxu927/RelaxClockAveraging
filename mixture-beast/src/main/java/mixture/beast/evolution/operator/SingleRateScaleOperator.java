package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Scale a single randomly chosen element of a positive rate vector: r <- r * exp(eps). "
        + "This is typically better than additive random-walk for positive parameters.")
public class SingleRateScaleOperator extends Operator {

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "positive rate vector parameter",
            Input.Validate.REQUIRED
    );

    public final Input<Double> windowInput = new Input<>(
            "window",
            "eps ~ Uniform(-window, +window) on log scale",
            0.3
    );

    private RealParameter rates;

    @Override
    public void initAndValidate() {
        rates = ratesInput.get();
        if (rates.getDimension() < 1) {
            throw new IllegalArgumentException("SingleRateScaleOperator: rates dimension must be >= 1");
        }
    }

    @Override
    public double proposal() {
        final double window = windowInput.get();
        if (!(window > 0.0)) return Double.NEGATIVE_INFINITY;

        final int dim = rates.getDimension();
        final int i = Randomizer.nextInt(dim);

        final double r = rates.getValue(i);
        if (!(r > 0.0)) return Double.NEGATIVE_INFINITY;

        final double eps = (Randomizer.nextDouble() * 2.0 - 1.0) * window;
        final double m = Math.exp(eps);

        final double rNew = r * m;
        if (!(rNew > 0.0) || Double.isNaN(rNew) || Double.isInfinite(rNew)) {
            return Double.NEGATIVE_INFINITY;
        }

        rates.setValue(i, rNew);

        // scale proposal Hastings: log(q(old|new)/q(new|old)) = eps
        return eps;
    }
}
