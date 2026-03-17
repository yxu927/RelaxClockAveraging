package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("UC-only non-centered hyper move: changes ucldStdev while keeping latent z_i fixed, "
        + "and reconstructs the shared rate vector accordingly. "
        + "Useful when the bottleneck is within-relax mixing under the UCLD branch of the mixture.")
public class UCLDStdevNonCenteredOperator extends Operator {

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

    public final Input<RealParameter> ucldStdevInput = new Input<>(
            "ucldStdev",
            "UC lognormal stdev (sigma on log scale)",
            Input.Validate.REQUIRED
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

    private RealParameter rates;
    private IntegerParameter indicator;
    private RealParameter ucldStdev;

    @Override
    public void initAndValidate() {
        rates = ratesInput.get();
        indicator = indicatorInput.get();
        ucldStdev = ucldStdevInput.get();

        if (rates.getDimension() < 1) {
            throw new IllegalArgumentException("rates dimension must be >= 1");
        }
        if (indicator.getDimension() != 1) {
            throw new IllegalArgumentException("indicator dimension must be 1");
        }
        if (ucldStdev.getDimension() != 1) {
            throw new IllegalArgumentException("ucldStdev dimension must be 1");
        }
    }

    @Override
    public double proposal() {
        final int k = indicator.getValue(0);
        if (k != 0) {
            return rejectIfNotUCInput.get() ? Double.NEGATIVE_INFINITY : 0.0;
        }

        final double window = windowInput.get();
        if (!(window > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double oldS = ucldStdev.getValue(0);
        if (!(oldS > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double eps = (Randomizer.nextDouble() * 2.0 - 1.0) * window;
        final double newS = oldS * Math.exp(eps);
        if (!(newS > 0.0) || Double.isInfinite(newS) || Double.isNaN(newS)) {
            return Double.NEGATIVE_INFINITY;
        }

        final int nEdges = rates.getDimension();
        final double oldVar = oldS * oldS;
        final double newVar = newS * newS;

        final double[] rNew = new double[nEdges];
        double sumDelta = 0.0;

        for (int i = 0; i < nEdges; i++) {
            final double r = rates.getValue(i);
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

        rates.startEditing(this);
        for (int i = 0; i < nEdges; i++) {
            rates.setValue(i, rNew[i]);
        }

        ucldStdev.startEditing(this);
        ucldStdev.setValue(0, newS);

        return sumDelta + (nEdges + 1.0) * eps;
    }
}
