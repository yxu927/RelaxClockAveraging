package mixture.beast.evolution.operator;

import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class SingleRateScaleOperatorTypedInputTest {

    private static final double EPS = 1.0e-10;

    @Test
    public void typedValidProposalChangesExactlyOneRateAndReturnsLogScaleFactor() {
        Randomizer.setSeed(101L);

        final RealVectorParam<PositiveReal> rates = positiveRates(1.0, 2.0, 4.0, 8.0);
        final double[] before = copyValues(rates);
        final double window = 0.37;
        final SingleRateScaleOperator op = typedOperator(rates, window);

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
    public void typedDefaultWindowIsPointThreeByObservedHastingsBound() {
        Randomizer.setSeed(202L);

        final RealVectorParam<PositiveReal> rates = positiveRates(1.0, 2.0, 4.0, 8.0);
        final SingleRateScaleOperator op = typedOperatorWithDefaultWindow(rates);

        final double hr = op.proposal();

        assertTrue(Double.isFinite(hr));
        assertTrue(Math.abs(hr) <= 0.3 + EPS);
    }

    @Test
    public void typedRepeatedValidProposalsKeepRatesPositiveAndReturnFiniteHastingsTerms() {
        Randomizer.setSeed(303L);

        final RealVectorParam<PositiveReal> rates = positiveRates(1.0, 2.0, 4.0, 8.0);
        final SingleRateScaleOperator op = typedOperator(rates, 0.2);

        for (int i = 0; i < 100; i++) {
            final double hr = op.proposal();
            assertTrue(Double.isFinite(hr));
            for (int j = 0; j < rates.size(); j++) {
                assertTrue(rates.get(j) > 0.0);
            }
        }
    }

    @Test
    public void typedWindowZeroReturnsNegativeInfinityAndDoesNotMutateRates() {
        final RealVectorParam<PositiveReal> rates = positiveRates(1.0, 2.0, 3.0);
        final double[] before = copyValues(rates);
        final SingleRateScaleOperator op = typedOperator(rates, 0.0);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void typedNegativeWindowReturnsNegativeInfinityAndDoesNotMutateRates() {
        final RealVectorParam<PositiveReal> rates = positiveRates(1.0, 2.0, 3.0);
        final double[] before = copyValues(rates);
        final SingleRateScaleOperator op = typedOperator(rates, -0.2);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void typedNonPositiveSelectedRateReturnsNegativeInfinityAndDoesNotMutateRates() {
        Randomizer.setSeed(11L);

        final RealVectorParam<Real> rates = realRates(0.0, 0.0, 0.0);
        final double[] before = copyValues(rates);
        final SingleRateScaleOperator op = typedOperator(rates, 0.3);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void legacyAndTypedRatesCannotBothBeSpecified() {
        final SingleRateScaleOperator op = new SingleRateScaleOperator();
        op.ratesInput.setValue(new RealParameter("1.0 2.0 3.0"), op);
        op.ratesVectorInput.setValue(positiveRates(1.0, 2.0, 3.0), op);

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void missingRatesIsRejected() {
        final SingleRateScaleOperator op = new SingleRateScaleOperator();

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void typedAndLegacyOperatorsEquivalentForSameRandomDraw() {
        final long seed = 707L;
        final double window = 0.61;

        final RealParameter legacyRates = new RealParameter("1.0 2.0 4.0 8.0");
        final SingleRateScaleOperator legacyOperator = legacyOperator(legacyRates, window);
        final double[] legacyBefore = copyValues(legacyRates);

        Randomizer.setSeed(seed);
        final double legacyHr = legacyOperator.proposal();

        final RealVectorParam<PositiveReal> typedRates = positiveRates(1.0, 2.0, 4.0, 8.0);
        final SingleRateScaleOperator typedOperator = typedOperator(typedRates, window);
        final double[] typedBefore = copyValues(typedRates);

        Randomizer.setSeed(seed);
        final double typedHr = typedOperator.proposal();

        assertArrayEquals(legacyBefore, typedBefore, EPS);
        assertEquals(legacyHr, typedHr, EPS);
        assertArrayEquals(copyValues(legacyRates), copyValues(typedRates), EPS);
        assertEquals(firstChangedIndex(legacyBefore, copyValues(legacyRates)),
                firstChangedIndex(typedBefore, copyValues(typedRates)));
    }

    private static SingleRateScaleOperator typedOperator(final RealVectorParam<?> rates,
                                                        final double window) {
        final SingleRateScaleOperator op = new SingleRateScaleOperator();
        op.ratesVectorInput.setValue(rates, op);
        op.windowInput.setValue(window, op);
        op.initAndValidate();
        return op;
    }

    private static SingleRateScaleOperator typedOperatorWithDefaultWindow(final RealVectorParam<?> rates) {
        final SingleRateScaleOperator op = new SingleRateScaleOperator();
        op.ratesVectorInput.setValue(rates, op);
        op.initAndValidate();
        return op;
    }

    private static SingleRateScaleOperator legacyOperator(final RealParameter rates,
                                                         final double window) {
        final SingleRateScaleOperator op = new SingleRateScaleOperator();
        op.ratesInput.setValue(rates, op);
        op.windowInput.setValue(window, op);
        op.initAndValidate();
        return op;
    }

    private static RealVectorParam<PositiveReal> positiveRates(final double... values) {
        return new RealVectorParam<>(values, PositiveReal.INSTANCE);
    }

    private static RealVectorParam<Real> realRates(final double... values) {
        return new RealVectorParam<>(values, Real.INSTANCE);
    }

    private static double[] copyValues(final RealParameter parameter) {
        final double[] values = new double[parameter.getDimension()];
        for (int i = 0; i < values.length; i++) {
            values[i] = parameter.getValue(i);
        }
        return values;
    }

    private static double[] copyValues(final RealVectorParam<?> parameter) {
        final double[] values = new double[parameter.size()];
        for (int i = 0; i < values.length; i++) {
            values[i] = parameter.get(i);
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
