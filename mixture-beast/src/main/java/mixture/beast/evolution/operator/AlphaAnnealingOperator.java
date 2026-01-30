package mixture.beast.evolution.operator;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;


public class AlphaAnnealingOperator extends Operator {

    public final Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Annealing exponent parameter (dimension 1).",
            Validate.REQUIRED);

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

    private RealParameter alphaParam;
    private double alphaStart;
    private double alphaEnd;
    private int    alphaSteps;
    private double delta; // step size


    @Override
    public void initAndValidate() {

        alphaParam = alphaInput.get();

        if (alphaParam == null) {
            throw new IllegalArgumentException(
                    "AlphaAnnealingOperator: input 'alpha' is not set.\n"
                            + "Fix your XML by adding: alpha='@alpha' (or <alpha idref='alpha'/>),\n"
                            + "and make sure a RealParameter with id='alpha' exists INSIDE <state> as a stateNode.");
        }

        if (alphaParam.getDimension() != 1) {
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

        final double current = alphaParam.getArrayValue(0);

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

        final double a = alphaParam.getArrayValue(0);

        if (a <= alphaEnd) {
            return 0.0;
        }

        double newA = a - delta;
        if (newA < alphaEnd) newA = alphaEnd;

        alphaParam.setValue(0, newA);

        return 0.0;
    }
}
