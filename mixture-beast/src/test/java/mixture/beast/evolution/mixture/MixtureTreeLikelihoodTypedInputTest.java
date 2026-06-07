package mixture.beast.evolution.mixture;

import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import org.junit.Test;

import java.util.List;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class MixtureTreeLikelihoodTypedInputTest {

    private static final double EPS = 1.0e-10;

    @Test
    public void typedInitRejectsWeightsDimensionMismatch() {
        final MixtureTreeLikelihood mix = shellTwoComponents();
        mix.weightsVectorInput.setValue(weights(1.0), mix);

        assertThrows(IllegalArgumentException.class, mix::initAndValidate);
    }

    @Test
    public void typedInitRejectsNegativeWeights() {
        final MixtureTreeLikelihood mix = shellTwoComponents();
        mix.weightsVectorInput.setValue(weights(1.1, -0.1), mix);

        assertThrows(IllegalArgumentException.class, mix::initAndValidate);
    }

    @Test
    public void typedWeightsThatDoNotSumToOneAreWarnedButNotRejected() {
        final MixtureTreeLikelihood mix = shellTwoComponents();
        mix.weightsVectorInput.setValue(weights(0.2, 0.2), mix);

        mix.initAndValidate();

        assertTrue(Double.isFinite(mix.calculateLogP()));
    }

    @Test
    public void typedCalculateLogPWithoutAlphaIsWeightedLogSumExp() {
        final MixtureTreeLikelihood mix = mixtureTyped(
                weights(0.25, 0.75),
                null,
                -10.0,
                -12.0
        );

        final double expected = logSumExp(Math.log(0.25) - 10.0, Math.log(0.75) - 12.0);

        assertEquals(expected, mix.calculateLogP(), EPS);
    }

    @Test
    public void typedZeroWeightComponentIsIgnoredWhenAlphaAbsent() {
        final CountingGenericTreeLikelihood skipped = constantLogL(Double.NaN);
        final CountingGenericTreeLikelihood used = constantLogL(-3.0);

        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        addSubLikelihood(mix, skipped);
        addSubLikelihood(mix, used);
        mix.weightsVectorInput.setValue(weights(0.0, 1.0), mix);
        mix.initAndValidate();

        assertEquals(-3.0, mix.calculateLogP(), EPS);
        assertEquals(0, skipped.callCount);
        assertEquals(1, used.callCount);
    }

    @Test
    public void typedNonFinitePositiveWeightComponentIsIgnoredForMixtureTermWhenAlphaAbsent() {
        final MixtureTreeLikelihood mix = mixtureTyped(
                weights(0.5, 0.5),
                null,
                Double.NEGATIVE_INFINITY,
                -2.0
        );

        assertEquals(Math.log(0.5) - 2.0, mix.calculateLogP(), EPS);
    }

    @Test
    public void typedAllMixtureTermsInvalidReturnsNegativeInfinity() {
        final MixtureTreeLikelihood mix = mixtureTyped(
                weights(0.5, 0.5),
                null,
                Double.NEGATIVE_INFINITY,
                Double.NaN
        );

        assertEquals(Double.NEGATIVE_INFINITY, mix.calculateLogP(), 0.0);
    }

    @Test
    public void typedAlphaZeroBehavesLikeNoCoupling() {
        final MixtureTreeLikelihood noAlpha = mixtureTyped(
                weights(0.25, 0.75),
                null,
                -10.0,
                -12.0
        );
        final MixtureTreeLikelihood alphaZero = mixtureTyped(
                weights(0.25, 0.75),
                alpha(0.0),
                -10.0,
                -12.0
        );

        assertEquals(noAlpha.calculateLogP(), alphaZero.calculateLogP(), EPS);
    }

    @Test
    public void typedPositiveAlphaAddsAlphaTimesSumLogLikelihoods() {
        final double alpha = 0.2;
        final MixtureTreeLikelihood mix = mixtureTyped(
                weights(0.25, 0.75),
                alpha(alpha),
                -10.0,
                -12.0
        );

        final double logMix = logSumExp(Math.log(0.25) - 10.0, Math.log(0.75) - 12.0);
        final double expected = logMix + alpha * (-10.0 - 12.0);

        assertEquals(expected, mix.calculateLogP(), EPS);
    }

    @Test
    public void typedPositiveAlphaRequiresAllComponentLogLikelihoodsFinite() {
        final MixtureTreeLikelihood mix = mixtureTyped(
                weights(0.5, 0.5),
                alpha(0.2),
                -1.0,
                Double.NEGATIVE_INFINITY
        );

        assertEquals(Double.NEGATIVE_INFINITY, mix.calculateLogP(), 0.0);
    }

    @Test
    public void legacyAndTypedWeightsCannotBothBeSpecified() {
        final MixtureTreeLikelihood mix = shellTwoComponents();
        mix.weightsInput.setValue(new RealParameter("0.5 0.5"), mix);
        mix.weightsVectorInput.setValue(weights(0.5, 0.5), mix);

        assertThrows(IllegalArgumentException.class, mix::initAndValidate);
    }

    @Test
    public void missingWeightsIsRejected() {
        final MixtureTreeLikelihood mix = shellTwoComponents();

        assertThrows(IllegalArgumentException.class, mix::initAndValidate);
    }

    @Test
    public void legacyAndTypedAlphaCannotBothBeSpecified() {
        final MixtureTreeLikelihood mix = shellTwoComponents();
        mix.weightsVectorInput.setValue(weights(0.5, 0.5), mix);
        mix.alphaInput.setValue(new RealParameter("0.2"), mix);
        mix.alphaScalarInput.setValue(alpha(0.2), mix);

        assertThrows(IllegalArgumentException.class, mix::initAndValidate);
    }

    @Test
    public void typedAndLegacyMixtureEquivalentWithoutAlpha() {
        final MixtureTreeLikelihood legacy = mixtureLegacy(new double[]{0.25, 0.75}, null, -10.0, -12.0);
        final MixtureTreeLikelihood typed = mixtureTyped(weights(0.25, 0.75), null, -10.0, -12.0);

        assertEquals(legacy.calculateLogP(), typed.calculateLogP(), EPS);
    }

    @Test
    public void typedAndLegacyMixtureEquivalentWithAlpha() {
        final MixtureTreeLikelihood legacy = mixtureLegacy(new double[]{0.25, 0.75}, 0.2, -10.0, -12.0);
        final MixtureTreeLikelihood typed = mixtureTyped(weights(0.25, 0.75), alpha(0.2), -10.0, -12.0);

        assertEquals(legacy.calculateLogP(), typed.calculateLogP(), EPS);
    }

    @Test
    public void requiresRecalculationArgumentsConditionsAndSampleRemain() {
        final MixtureTreeLikelihood mix = mixtureTyped(weights(0.5, 0.5), null, -1.0, -2.0);

        assertTrue(mix.requiresRecalculation());
        assertEquals(List.of(), mix.getArguments());
        assertEquals(List.of(), mix.getConditions());
        mix.sample((State) null, new Random(1L));
    }

    private static MixtureTreeLikelihood shellTwoComponents() {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        addSubLikelihood(mix, constantLogL(-1.0));
        addSubLikelihood(mix, constantLogL(-2.0));
        return mix;
    }

    private static MixtureTreeLikelihood mixtureTyped(final RealVectorParam<Real> weights,
                                                      final RealScalarParam<Real> alpha,
                                                      final double... logL) {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        for (final double value : logL) {
            addSubLikelihood(mix, constantLogL(value));
        }
        mix.weightsVectorInput.setValue(weights, mix);
        if (alpha != null) {
            mix.alphaScalarInput.setValue(alpha, mix);
        }
        mix.initAndValidate();
        return mix;
    }

    private static MixtureTreeLikelihood mixtureLegacy(final double[] weights,
                                                       final Double alpha,
                                                       final double... logL) {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        for (final double value : logL) {
            addSubLikelihood(mix, constantLogL(value));
        }
        mix.weightsInput.setValue(realParameter(weights), mix);
        if (alpha != null) {
            mix.alphaInput.setValue(new RealParameter(Double.toString(alpha)), mix);
        }
        mix.initAndValidate();
        return mix;
    }

    private static void addSubLikelihood(final MixtureTreeLikelihood mix,
                                         final GenericTreeLikelihood likelihood) {
        mix.subLikelihoodsInput.get().add(likelihood);
    }

    private static RealVectorParam<Real> weights(final double... values) {
        return new RealVectorParam<>(values, Real.INSTANCE);
    }

    private static RealScalarParam<Real> alpha(final double value) {
        return new RealScalarParam<>(value, Real.INSTANCE);
    }

    private static RealParameter realParameter(final double[] values) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < values.length; i++) {
            if (i > 0) {
                sb.append(' ');
            }
            sb.append(values[i]);
        }
        return new RealParameter(sb.toString());
    }

    private static CountingGenericTreeLikelihood constantLogL(final double value) {
        return new CountingGenericTreeLikelihood(value);
    }

    private static double logSumExp(final double a, final double b) {
        if (a > b) {
            return a + Math.log1p(Math.exp(b - a));
        }
        return b + Math.log1p(Math.exp(a - b));
    }

    public static final class CountingGenericTreeLikelihood extends GenericTreeLikelihood {
        private final double value;
        private int callCount = 0;

        private CountingGenericTreeLikelihood(final double value) {
            this.value = value;
        }

        @Override
        public void initAndValidate() {
            // Do not require data/tree/siteModel for this constant test stub.
        }

        @Override
        public double calculateLogP() {
            callCount++;
            logP = value;
            return logP;
        }
    }
}
