package mixture.lphybeast.tobeast;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TypedParameterUtilsTest {

    @Test
    public void positiveRealVectorFromLegacyScalarBroadcastsToTargetDimension() {
        final RealParameter legacy = new RealParameter("1.25");
        legacy.setID("rawRates");

        final RealVectorParam<?> typed = TypedParameterUtils.positiveRealVectorFrom(
                legacy, "fallbackRates", 4, 1.0, true);

        assertEquals("rawRates", typed.getID());
        assertEquals(4, typed.size());
        for (int i = 0; i < typed.size(); i++) {
            assertEquals(1.25, typed.get(i), 1.0e-12);
        }
    }

    @Test
    public void positiveRealVectorFromNullUsesFallbackValue() {
        final RealVectorParam<?> typed = TypedParameterUtils.positiveRealVectorFrom(
                null, "rawRates", 3, 1.0, true);

        assertEquals("rawRates", typed.getID());
        assertEquals(3, typed.size());
        for (int i = 0; i < typed.size(); i++) {
            assertEquals(1.0, typed.get(i), 1.0e-12);
        }
    }

    @Test
    public void realVectorFromConstantValuesIsUnestimated() {
        final RealParameter legacy = new RealParameter("0.25 0.75");
        legacy.setID("weights");

        final RealVectorParam<?> typed = TypedParameterUtils.realVectorFrom(
                legacy, "fallbackWeights", new double[]{0.25, 0.75}, false);

        assertEquals("weights", typed.getID());
        assertEquals(2, typed.size());
        assertEquals(0.25, typed.get(0), 1.0e-12);
        assertEquals(0.75, typed.get(1), 1.0e-12);
        assertTrue(!typed.isEstimated());
    }

    @Test
    public void intScalarFromLegacyIntegerPreservesValueAndEstimateFlag() {
        final IntegerParameter legacy = new IntegerParameter("1");
        legacy.setID("indicator");

        final IntScalarParam<?> typed = TypedParameterUtils.intScalarFrom(
                legacy, "fallbackIndicator", 0, true);

        assertEquals("indicator", typed.getID());
        assertEquals(1, typed.get());
    }
}
