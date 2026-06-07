package mixture.beast.evolution.mixture;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.type.IntScalar;
import beast.base.spec.type.IntVector;
import beast.base.spec.type.RealVector;

import java.util.List;
import java.util.Random;

@Description("Categorical distribution over K categories. "
        + "If 'parameter' (x) is multidimensional, each entry contributes independently to the total logP. "
        + "By default categories are zero-indexed [0..K-1]; set oneIndexed=\"true\" to interpret x in [1..K].")
public class CategoricalDistribution extends Distribution {

    // Inputs
    public final Input<RealParameter> pInput =
            new Input<>("p",
                    "Legacy probability vector over K categories; entries must be non-negative and sum to 1.",
                    Validate.OPTIONAL);

    public final Input<RealVector> pVectorInput =
            new Input<>("pVector",
                    "BEAST3 typed probability vector over K categories; entries must be non-negative and sum to 1.",
                    Validate.OPTIONAL);

    public final Input<IntegerParameter> xInput =
            new Input<>("parameter", "Legacy integer category or vector of categories.", Validate.OPTIONAL);

    public final Input<IntScalar> xScalarInput =
            new Input<>("parameterScalar", "BEAST3 typed scalar category.", Validate.OPTIONAL);

    public final Input<IntVector> xVectorInput =
            new Input<>("parameterVector", "BEAST3 typed vector of categories.", Validate.OPTIONAL);

    public final Input<Boolean> oneIndexedInput =
            new Input<>("oneIndexed", "Interpret categories as 1..K instead of 0..K-1 (default false).", false);

    private static final double SUM_TOL = 1e-8;

    private RealParameter legacyP;
    private RealVector typedP;
    private IntegerParameter legacyX;
    private IntScalar typedXScalar;
    private IntVector typedXVector;

    @Override
    public double calculateLogP() {
        this.logP = 0.0;

        final boolean oneIndexed = oneIndexedInput.get() != null && oneIndexedInput.get();

        final int K = pDimension();
        if (K < 1) {
            return Double.NEGATIVE_INFINITY;
        }

        // Check probabilities and sum-to-one
        double sum = 0.0;
        for (int i = 0; i < K; i++) {
            double pi = pValue(i);
            if (pi < 0.0 || Double.isNaN(pi) || Double.isInfinite(pi)) {
                return Double.NEGATIVE_INFINITY;
            }
            sum += pi;
        }
        if (Math.abs(sum - 1.0) > SUM_TOL) {
            return Double.NEGATIVE_INFINITY;
        }

        final int D = xDimension();
        for (int d = 0; d < D; d++) {
            int xv = xValue(d);

            int idx = oneIndexed ? (xv - 1) : xv;
            if (idx < 0 || idx >= K) {
                return Double.NEGATIVE_INFINITY;
            }

            double px = pValue(idx);
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
        legacyP = pInput.get();
        typedP = pVectorInput.get();
        legacyX = xInput.get();
        typedXScalar = xScalarInput.get();
        typedXVector = xVectorInput.get();

        requireExactlyOne(legacyP, typedP, "p", "pVector");

        final int categoryInputCount = (legacyX == null ? 0 : 1)
                + (typedXScalar == null ? 0 : 1)
                + (typedXVector == null ? 0 : 1);
        if (categoryInputCount != 1) {
            throw new IllegalArgumentException("CategoricalDistribution: exactly one of parameter, parameterScalar, or parameterVector must be specified.");
        }

        if (pDimension() < 1) {
            throw new RuntimeException("p must have dimension >= 1.");
        }
        if (xDimension() < 1) {
            throw new RuntimeException("parameter (x) must have dimension >= 1.");
        }
    }

    private static void requireExactlyOne(final Object legacy,
                                          final Object typed,
                                          final String legacyName,
                                          final String typedName) {
        if (legacy == null && typed == null) {
            throw new IllegalArgumentException("CategoricalDistribution: either " + legacyName
                    + " or " + typedName + " must be specified.");
        }
        if (legacy != null && typed != null) {
            throw new IllegalArgumentException("CategoricalDistribution: specify only one of " + legacyName
                    + " or " + typedName + ".");
        }
    }

    private int pDimension() {
        return legacyP != null ? legacyP.getDimension() : typedP.size();
    }

    private double pValue(final int i) {
        return legacyP != null ? legacyP.getArrayValue(i) : typedP.get(i);
    }

    private int xDimension() {
        if (legacyX != null) {
            return legacyX.getDimension();
        }
        if (typedXScalar != null) {
            return 1;
        }
        return typedXVector.size();
    }

    private int xValue(final int d) {
        if (legacyX != null) {
            return (int) Math.round(legacyX.getArrayValue(d));
        }
        if (typedXScalar != null) {
            return typedXScalar.get();
        }
        return typedXVector.get(d);
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
