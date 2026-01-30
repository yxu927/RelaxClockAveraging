package mixture.beast.evolution.clockmodel.auto;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.List;
import java.util.Random;

@Description("Categorical distribution over K categories. "
        + "If 'parameter' (x) is multidimensional, each entry contributes independently to the total logP. "
        + "By default categories are zero-indexed [0..K-1]; set oneIndexed=\"true\" to interpret x in [1..K].")
public class CategoricalDistribution extends Distribution {

    // Inputs
    public final Input<RealParameter> pInput =
            new Input<>("p", "Probability vector over K categories; entries must be non-negative and sum to 1.", Validate.REQUIRED);

    public final Input<IntegerParameter> xInput =
            new Input<>("parameter", "Integer category (or vector of categories).", Validate.REQUIRED);

    public final Input<Boolean> oneIndexedInput =
            new Input<>("oneIndexed", "Interpret categories as 1..K instead of 0..K-1 (default false).", false);

    private static final double SUM_TOL = 1e-8;

    @Override
    public double calculateLogP() {
        this.logP = 0.0;

        final RealParameter p = pInput.get();
        final IntegerParameter x = xInput.get();
        final boolean oneIndexed = oneIndexedInput.get() != null && oneIndexedInput.get();

        final int K = p.getDimension();
        if (K < 1) {
            return Double.NEGATIVE_INFINITY;
        }

        // Check probabilities and sum-to-one
        double sum = 0.0;
        for (int i = 0; i < K; i++) {
            double pi = p.getArrayValue(i);
            if (pi < 0.0 || Double.isNaN(pi) || Double.isInfinite(pi)) {
                return Double.NEGATIVE_INFINITY;
            }
            sum += pi;
        }
        if (Math.abs(sum - 1.0) > SUM_TOL) {
            return Double.NEGATIVE_INFINITY;
        }

        final int D = x.getDimension();
        for (int d = 0; d < D; d++) {
            int xv = (int) Math.round(x.getArrayValue(d));

            int idx = oneIndexed ? (xv - 1) : xv;
            if (idx < 0 || idx >= K) {
                return Double.NEGATIVE_INFINITY;
            }

            double px = p.getArrayValue(idx);
            if (px <= 0.0) {
                return Double.NEGATIVE_INFINITY;
            }
            this.logP += Math.log(px);
        }

        return this.logP;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void initAndValidate() {
        final RealParameter p = pInput.get();
        final IntegerParameter x = xInput.get();

        if (p.getDimension() < 1) {
            throw new RuntimeException("p must have dimension >= 1.");
        }
        if (x.getDimension() < 1) {
            throw new RuntimeException("parameter (x) must have dimension >= 1.");
        }
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }
}
