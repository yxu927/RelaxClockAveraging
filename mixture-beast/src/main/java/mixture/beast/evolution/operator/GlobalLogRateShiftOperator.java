package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Shift rootLogRate and all node log-rates by the same additive delta.")
public class GlobalLogRateShiftOperator extends Operator {

    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "root log-rate (dimension=1)",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> nodeLogRatesInput = new Input<>(
            "nodeRates",
            "log-rates for non-root nodes",
            Input.Validate.REQUIRED
    );

    public final Input<Double> windowSizeInput = new Input<>(
            "windowSize",
            "uniform shift size in log space",
            0.1
    );

    private RealParameter rootLogRate;
    private RealParameter nodeLogRates;

    @Override
    public void initAndValidate() {
        rootLogRate = rootLogRateInput.get();
        nodeLogRates = nodeLogRatesInput.get();

        if (rootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate must have dimension=1.");
        }
    }

    /** Hook for deterministic testing */
    protected double drawDelta(final double w) {
        return (Randomizer.nextDouble() - 0.5) * 2.0 * w;
    }

    @Override
    public double proposal() {
        final double w = windowSizeInput.get();
        final double delta = drawDelta(w);

        rootLogRate.startEditing(this);
        nodeLogRates.startEditing(this);

        rootLogRate.setValue(0, rootLogRate.getValue(0) + delta);

        final int dim = nodeLogRates.getDimension();
        for (int i = 0; i < dim; i++) {
            nodeLogRates.setValue(i, nodeLogRates.getValue(i) + delta);
        }

        return 0.0;
    }
}
