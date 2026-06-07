package mixture.beast.evolution.operator;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.RealScalarParam;


public class AlphaAnnealingOperator extends Operator {

    public final Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Legacy annealing exponent parameter (dimension 1).",
            Validate.OPTIONAL);

    public final Input<RealScalarParam<?>> alphaScalarInput = new Input<>(
            "alphaScalar",
            "BEAST3 typed mutable annealing exponent scalar.",
            Validate.OPTIONAL);

    public final Input<Double> alphaStartInput = new Input<>(
            "alphaStart",
            "Nominal starting value alpha_start (used to compute step size). "
                    + "Should match the initial value of the alpha parameter.",
            1.0);

    public final Input<Double> alphaEndInput = new Input<>(
            "alphaEnd",
            "Final value alpha_end after annealing (typically 0.0).",
            0.0);

    public final Input<Integer> alphaStepsInput = new Input<>(
            "alphaSteps",
            "Approximate number of annealing steps over which alpha moves "
                    + "from alphaStart to alphaEnd. If 0, alpha is never changed.",
            0);

    private RealParameter legacyAlpha;
    private RealScalarParam<?> typedAlpha;
    private double alphaStart;
    private double alphaEnd;
    private int    alphaSteps;
    private double delta; // step size


    @Override
    public void initAndValidate() {

        legacyAlpha = alphaInput.get();
        typedAlpha = alphaScalarInput.get();

        if (legacyAlpha == null && typedAlpha == null) {
            throw new IllegalArgumentException(
                    "AlphaAnnealingOperator: either alpha or alphaScalar must be specified.\n"
                            + "For legacy XML add alpha='@alpha'. For BEAST3 typed XML add alphaScalar='@alpha'.");
        }
        if (legacyAlpha != null && typedAlpha != null) {
            throw new IllegalArgumentException(
                    "AlphaAnnealingOperator: specify only one of alpha or alphaScalar.");
        }

        if (legacyAlpha != null && legacyAlpha.getDimension() != 1) {
            throw new IllegalArgumentException(
                    "AlphaAnnealingOperator: alpha parameter must have dimension 1.");
        }

        alphaStart = alphaStartInput.get();
        alphaEnd   = alphaEndInput.get();
        alphaSteps = alphaStepsInput.get();

        if (alphaSteps < 0) {
            throw new IllegalArgumentException("alphaSteps must be ≥ 0.");
        }
        if (alphaStart < alphaEnd) {
            throw new IllegalArgumentException("require alphaStart ≥ alphaEnd.");
        }

        final double current = alphaValue();

        if (Math.abs(current - alphaStart) > 1e-8) {
            Log.warning(String.format(
                    "AlphaAnnealingOperator: current alpha=%.6f but alphaStart=%.6f; using current.",
                    current, alphaStart));
            alphaStart = current;
        }

        if (alphaSteps == 0 || alphaStart == alphaEnd || current <= alphaEnd) {
            delta = 0.0;
        } else {
            delta = (alphaStart - alphaEnd) / (double) alphaSteps;
        }
    }

    @Override
    public double proposal() {

        if (delta <= 0.0) {
            return 0.0;
        }

        final double a = alphaValue();

        if (a <= alphaEnd) {
            return 0.0;
        }

        double newA = a - delta;
        if (newA < alphaEnd) newA = alphaEnd;

        setAlphaValue(newA);

        return 0.0;
    }

    private double alphaValue() {
        return legacyAlpha != null ? legacyAlpha.getArrayValue(0) : typedAlpha.get();
    }

    private void setAlphaValue(final double value) {
        if (legacyAlpha != null) {
            legacyAlpha.setValue(0, value);
        } else {
            typedAlpha.set(value);
        }
    }
}
