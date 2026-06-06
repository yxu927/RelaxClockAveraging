package mixture.beast.evolution.mixture;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;

public class CategoricalDistributionCharacterizationTest {

    private static final double EPS = 1.0e-12;

    @Test
    public void zeroIndexedSingleCategoryReturnsLogProbability() {
        final CategoricalDistribution cat = categorical("0.25 0.75", "1", false);

        assertEquals(Math.log(0.75), cat.calculateLogP(), EPS);
    }

    @Test
    public void oneIndexedSingleCategoryReturnsLogProbability() {
        final CategoricalDistribution cat = categorical("0.25 0.75", "2", true);

        assertEquals(Math.log(0.75), cat.calculateLogP(), EPS);
    }

    @Test
    public void multidimensionalCategoriesSumIndependentLogProbabilities() {
        final CategoricalDistribution cat = categorical("0.2 0.3 0.5", "0 2 1 2", false);

        final double expected = Math.log(0.2) + Math.log(0.5) + Math.log(0.3) + Math.log(0.5);

        assertEquals(expected, cat.calculateLogP(), EPS);
    }

    @Test
    public void probabilityVectorMustSumToOne() {
        final CategoricalDistribution cat = categorical("0.2 0.2", "0", false);

        assertEquals(Double.NEGATIVE_INFINITY, cat.calculateLogP(), 0.0);
    }

    @Test
    public void negativeProbabilityReturnsNegativeInfinity() {
        final CategoricalDistribution cat = categorical("1.1 -0.1", "0", false);

        assertEquals(Double.NEGATIVE_INFINITY, cat.calculateLogP(), 0.0);
    }

    @Test
    public void selectedZeroProbabilityReturnsNegativeInfinity() {
        final CategoricalDistribution cat = categorical("0.0 1.0", "0", false);

        assertEquals(Double.NEGATIVE_INFINITY, cat.calculateLogP(), 0.0);
    }

    @Test
    public void outOfRangeZeroIndexedCategoryReturnsNegativeInfinity() {
        assertEquals(Double.NEGATIVE_INFINITY, categorical("0.25 0.75", "-1", false).calculateLogP(), 0.0);
        assertEquals(Double.NEGATIVE_INFINITY, categorical("0.25 0.75", "2", false).calculateLogP(), 0.0);
    }

    @Test
    public void outOfRangeOneIndexedCategoryReturnsNegativeInfinity() {
        assertEquals(Double.NEGATIVE_INFINITY, categorical("0.25 0.75", "0", true).calculateLogP(), 0.0);
        assertEquals(Double.NEGATIVE_INFINITY, categorical("0.25 0.75", "3", true).calculateLogP(), 0.0);
    }

    @Test
    public void sampleThrowsUnsupportedOperationException() {
        final CategoricalDistribution cat = categorical("0.25 0.75", "1", false);

        assertThrows(UnsupportedOperationException.class, () -> cat.sample(null, new Random(1L)));
    }

    @Test
    public void getArgumentsAndConditionsCurrentlyReturnNull() {
        final CategoricalDistribution cat = categorical("0.25 0.75", "1", false);

        assertNull(cat.getArguments());
        assertNull(cat.getConditions());
    }

    private static CategoricalDistribution categorical(final String p,
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
}
