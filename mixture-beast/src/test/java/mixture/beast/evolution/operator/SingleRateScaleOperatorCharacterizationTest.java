package mixture.beast.evolution.operator;

import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class SingleRateScaleOperatorCharacterizationTest {

    private static final double EPS = 1.0e-10;

    @Test
    public void invalidWindowReturnsNegativeInfinityAndDoesNotMutateRates() {
        final RealParameter rates = realParameter("1.0 2.0 3.0");
        final double[] before = copyValues(rates);
        final SingleRateScaleOperator op = operator(rates, 0.0);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void negativeWindowReturnsNegativeInfinityAndDoesNotMutateRates() {
        final RealParameter rates = realParameter("1.0 2.0 3.0");
        final double[] before = copyValues(rates);
        final SingleRateScaleOperator op = operator(rates, -0.2);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void allNonPositiveRatesReturnNegativeInfinityAndDoNotMutateRates() {
        Randomizer.setSeed(11L);

        final RealParameter rates = realParameter("0.0 0.0 0.0");
        final double[] before = copyValues(rates);
        final SingleRateScaleOperator op = operator(rates, 0.3);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void validProposalChangesExactlyOneRateAndReturnsLogScaleFactor() {
        Randomizer.setSeed(101L);

        final RealParameter rates = realParameter("1.0 2.0 4.0 8.0");
        final double[] before = copyValues(rates);
        final double window = 0.37;
        final SingleRateScaleOperator op = operator(rates, window);

        final double hr = op.proposal();

        assertTrue(Double.isFinite(hr));
        assertTrue(Math.abs(hr) <= window + EPS);

        final double[] after = copyValues(rates);
        assertEquals(1, changedCount(before, after));

        final int idx = firstChangedIndex(before, after);
        assertEquals(Math.log(after[idx] / before[idx]), hr, EPS);

        for (int i = 0; i < before.length; i++) {
            if (i != idx) {
                assertEquals(before[i], after[i], EPS);
            }
        }
    }

    @Test
    public void defaultWindowIsPointThreeByObservedHastingsBound() {
        Randomizer.setSeed(202L);

        final RealParameter rates = realParameter("1.0 2.0 4.0 8.0");
        final SingleRateScaleOperator op = operatorWithDefaultWindow(rates);

        final double hr = op.proposal();

        assertTrue(Double.isFinite(hr));
        assertTrue(Math.abs(hr) <= 0.3 + EPS);
    }

    @Test
    public void repeatedValidProposalsKeepRatesPositiveAndReturnFiniteHastingsTerms() {
        Randomizer.setSeed(303L);

        final RealParameter rates = realParameter("1.0 2.0 4.0 8.0");
        final SingleRateScaleOperator op = operator(rates, 0.2);

        for (int i = 0; i < 100; i++) {
            final double hr = op.proposal();
            assertTrue(Double.isFinite(hr));
            for (int j = 0; j < rates.getDimension(); j++) {
                assertTrue(rates.getValue(j) > 0.0);
            }
        }
    }

    private static SingleRateScaleOperator operator(final RealParameter rates, final double window) {
        final SingleRateScaleOperator op = new SingleRateScaleOperator();
        op.ratesInput.setValue(rates, op);
        op.windowInput.setValue(window, op);
        op.initAndValidate();
        return op;
    }

    private static SingleRateScaleOperator operatorWithDefaultWindow(final RealParameter rates) {
        final SingleRateScaleOperator op = new SingleRateScaleOperator();
        op.ratesInput.setValue(rates, op);
        op.initAndValidate();
        return op;
    }

    private static RealParameter realParameter(final String value) {
        return new RealParameter(value);
    }

    private static double[] copyValues(final RealParameter parameter) {
        final double[] values = new double[parameter.getDimension()];
        for (int i = 0; i < values.length; i++) {
            values[i] = parameter.getValue(i);
        }
        return values;
    }

    private static int changedCount(final double[] before, final double[] after) {
        int count = 0;
        for (int i = 0; i < before.length; i++) {
            if (Math.abs(before[i] - after[i]) > EPS) {
                count++;
            }
        }
        return count;
    }

    private static int firstChangedIndex(final double[] before, final double[] after) {
        for (int i = 0; i < before.length; i++) {
            if (Math.abs(before[i] - after[i]) > EPS) {
                return i;
            }
        }
        throw new IllegalArgumentException("No changed index found");
    }

    private static void assertArrayEquals(final double[] expected, final double[] observed, final double eps) {
        assertEquals(expected.length, observed.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], observed[i], eps);
        }
    }
}
