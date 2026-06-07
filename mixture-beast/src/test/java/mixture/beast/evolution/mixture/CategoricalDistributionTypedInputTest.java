package mixture.beast.evolution.mixture;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Int;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.IntVectorParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;

public class CategoricalDistributionTypedInputTest {

    private static final double EPS = 1.0e-12;

    @Test
    public void typedZeroIndexedScalarCategoryReturnsLogProbability() {
        final CategoricalDistribution cat = categoricalWithScalar(p(0.25, 0.75), xScalar(1), false);

        assertEquals(Math.log(0.75), cat.calculateLogP(), EPS);
    }

    @Test
    public void typedOneIndexedScalarCategoryReturnsLogProbability() {
        final CategoricalDistribution cat = categoricalWithScalar(p(0.25, 0.75), xScalar(2), true);

        assertEquals(Math.log(0.75), cat.calculateLogP(), EPS);
    }

    @Test
    public void typedMultidimensionalCategoriesSumIndependentLogProbabilities() {
        final CategoricalDistribution cat = categoricalWithVector(p(0.2, 0.3, 0.5),
                xVector(0, 2, 1, 2), false);

        final double expected = Math.log(0.2) + Math.log(0.5) + Math.log(0.3) + Math.log(0.5);

        assertEquals(expected, cat.calculateLogP(), EPS);
    }

    @Test
    public void typedProbabilityVectorMustSumToOne() {
        final CategoricalDistribution cat = categoricalWithScalar(p(0.2, 0.2), xScalar(0), false);

        assertEquals(Double.NEGATIVE_INFINITY, cat.calculateLogP(), 0.0);
    }

    @Test
    public void typedNegativeProbabilityReturnsNegativeInfinity() {
        final CategoricalDistribution cat = categoricalWithScalar(p(1.1, -0.1), xScalar(0), false);

        assertEquals(Double.NEGATIVE_INFINITY, cat.calculateLogP(), 0.0);
    }

    @Test
    public void typedSelectedZeroProbabilityReturnsNegativeInfinity() {
        final CategoricalDistribution cat = categoricalWithScalar(p(0.0, 1.0), xScalar(0), false);

        assertEquals(Double.NEGATIVE_INFINITY, cat.calculateLogP(), 0.0);
    }

    @Test
    public void typedOutOfRangeZeroIndexedCategoryReturnsNegativeInfinity() {
        assertEquals(Double.NEGATIVE_INFINITY,
                categoricalWithScalar(p(0.25, 0.75), xScalar(-1), false).calculateLogP(), 0.0);
        assertEquals(Double.NEGATIVE_INFINITY,
                categoricalWithScalar(p(0.25, 0.75), xScalar(2), false).calculateLogP(), 0.0);
    }

    @Test
    public void typedOutOfRangeOneIndexedCategoryReturnsNegativeInfinity() {
        assertEquals(Double.NEGATIVE_INFINITY,
                categoricalWithScalar(p(0.25, 0.75), xScalar(0), true).calculateLogP(), 0.0);
        assertEquals(Double.NEGATIVE_INFINITY,
                categoricalWithScalar(p(0.25, 0.75), xScalar(3), true).calculateLogP(), 0.0);
    }

    @Test
    public void legacyAndTypedPInputsCannotBothBeSpecified() {
        final CategoricalDistribution cat = new CategoricalDistribution();
        cat.pInput.setValue(new RealParameter("0.25 0.75"), cat);
        cat.pVectorInput.setValue(p(0.25, 0.75), cat);
        cat.xScalarInput.setValue(xScalar(1), cat);

        assertThrows(IllegalArgumentException.class, cat::initAndValidate);
    }

    @Test
    public void missingPInputIsRejected() {
        final CategoricalDistribution cat = new CategoricalDistribution();
        cat.xScalarInput.setValue(xScalar(1), cat);

        assertThrows(IllegalArgumentException.class, cat::initAndValidate);
    }

    @Test
    public void legacyAndTypedParameterInputsCannotBothBeSpecified() {
        final CategoricalDistribution cat = new CategoricalDistribution();
        cat.pVectorInput.setValue(p(0.25, 0.75), cat);
        cat.xInput.setValue(new IntegerParameter("1"), cat);
        cat.xScalarInput.setValue(xScalar(1), cat);

        assertThrows(IllegalArgumentException.class, cat::initAndValidate);
    }

    @Test
    public void scalarAndVectorTypedParameterInputsCannotBothBeSpecified() {
        final CategoricalDistribution cat = new CategoricalDistribution();
        cat.pVectorInput.setValue(p(0.25, 0.75), cat);
        cat.xScalarInput.setValue(xScalar(1), cat);
        cat.xVectorInput.setValue(xVector(1), cat);

        assertThrows(IllegalArgumentException.class, cat::initAndValidate);
    }

    @Test
    public void missingParameterInputIsRejected() {
        final CategoricalDistribution cat = new CategoricalDistribution();
        cat.pVectorInput.setValue(p(0.25, 0.75), cat);

        assertThrows(IllegalArgumentException.class, cat::initAndValidate);
    }

    @Test
    public void typedAndLegacyCategoricalEquivalentForScalar() {
        final CategoricalDistribution legacy = legacyCategorical("0.25 0.75", "1", false);
        final CategoricalDistribution typed = categoricalWithScalar(p(0.25, 0.75), xScalar(1), false);

        assertEquals(legacy.calculateLogP(), typed.calculateLogP(), EPS);
    }

    @Test
    public void sampleUnsupportedAndMetadataNullRemainOnTypedPath() {
        final CategoricalDistribution cat = categoricalWithScalar(p(0.25, 0.75), xScalar(1), false);

        assertThrows(UnsupportedOperationException.class, () -> cat.sample(null, new Random(1L)));
        assertNull(cat.getArguments());
        assertNull(cat.getConditions());
    }

    private static CategoricalDistribution categoricalWithScalar(final RealVectorParam<Real> p,
                                                                final IntScalarParam<Int> x,
                                                                final boolean oneIndexed) {
        final CategoricalDistribution cat = new CategoricalDistribution();
        cat.initByName(
                "pVector", p,
                "parameterScalar", x,
                "oneIndexed", oneIndexed
        );
        return cat;
    }

    private static CategoricalDistribution categoricalWithVector(final RealVectorParam<Real> p,
                                                                final IntVectorParam<Int> x,
                                                                final boolean oneIndexed) {
        final CategoricalDistribution cat = new CategoricalDistribution();
        cat.initByName(
                "pVector", p,
                "parameterVector", x,
                "oneIndexed", oneIndexed
        );
        return cat;
    }

    private static CategoricalDistribution legacyCategorical(final String p,
                                                            final String x,
                                                            final boolean oneIndexed) {
        final CategoricalDistribution cat = new CategoricalDistribution();
        cat.initByName(
                "p", new RealParameter(p),
                "parameter", new IntegerParameter(x),
                "oneIndexed", oneIndexed
        );
        return cat;
    }

    private static RealVectorParam<Real> p(final double... values) {
        return new RealVectorParam<>(values, Real.INSTANCE);
    }

    private static IntScalarParam<Int> xScalar(final int value) {
        return new IntScalarParam<>(value, Int.INSTANCE);
    }

    private static IntVectorParam<Int> xVector(final int... values) {
        return new IntVectorParam<>(values, Int.INSTANCE);
    }
}
