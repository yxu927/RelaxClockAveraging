package mixture.beast.evolution.operator;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import mixture.beast.evolution.mixture.RelaxedRatesPriorSVS;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class IndicatorGibbsOperatorCharacterizationTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";
    private static final String POSITIVE_RATES = "0.8 1.1 1.4 0.9";

    @Test
    public void initRejectsIndicatorDimensionGreaterThanOne() {
        final IntegerParameter indicator = integerParameter("0 1");
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter(POSITIVE_RATES), integerParameter("0"),
                realParameter("0.5"), realParameter("0.0"), realParameter("0.2"));

        final IndicatorGibbsOperator op = new IndicatorGibbsOperator();
        op.indicatorInput.setValue(indicator, op);
        op.priorInput.setValue(prior, op);

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void initRejectsInvalidPOne() {
        assertThrows(IllegalArgumentException.class,
                () -> operator(integerParameter("0"), validPrior(integerParameter("0")), 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> operator(integerParameter("0"), validPrior(integerParameter("0")), 1.0));
        assertThrows(IllegalArgumentException.class,
                () -> operator(integerParameter("0"), validPrior(integerParameter("0")), -0.1));
        assertThrows(IllegalArgumentException.class,
                () -> operator(integerParameter("0"), validPrior(integerParameter("0")), 1.1));
    }

    @Test
    public void defaultPOneIsOneHalfAndProposalUsesIt() {
        Randomizer.setSeed(12345L);

        final IntegerParameter indicator = integerParameter("0");
        final RelaxedRatesPriorSVS prior = validPrior(indicator);
        final IndicatorGibbsOperator op = new IndicatorGibbsOperator();
        op.indicatorInput.setValue(indicator, op);
        op.priorInput.setValue(prior, op);
        op.initAndValidate();

        final int oldK = indicator.getValue(0);
        final Conditional conditional = conditional(prior, 0.5);
        final double hr = op.proposal();
        final int newK = indicator.getValue(0);

        assertBinary(newK);
        assertEquals(conditional.expectedLogHastings(oldK, newK), hr, EPS);
    }

    @Test
    public void invalidOldIndicatorReturnsNegativeInfinityAndLeavesIndicatorUnchanged() {
        final IntegerParameter indicator = integerParameter("2");
        final RelaxedRatesPriorSVS prior = validPrior(indicator);
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertEquals(2, indicator.getValue(0).intValue());
    }

    @Test
    public void bothPriorWeightsNegativeInfinityReturnsNegativeInfinityAndLeavesIndicatorUnchanged() {
        final IntegerParameter indicator = integerParameter("0");
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter("0.0 0.0 0.0 0.0"), indicator,
                realParameter("0.5"), realParameter("0.0"), realParameter("0.2"));
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertEquals(0, indicator.getValue(0).intValue());
    }

    @Test
    public void proposalFromOldZeroReturnsExpectedGibbsHastingsTerm() {
        Randomizer.setSeed(101L);

        final IntegerParameter indicator = integerParameter("0");
        final RelaxedRatesPriorSVS prior = validPrior(indicator);
        final double pOne = 0.5;
        final IndicatorGibbsOperator op = operator(indicator, prior, pOne);

        assertProposalMatchesExpected(indicator, prior, op, pOne);
    }

    @Test
    public void proposalFromOldOneReturnsExpectedGibbsHastingsTerm() {
        Randomizer.setSeed(202L);

        final IntegerParameter indicator = integerParameter("1");
        final RelaxedRatesPriorSVS prior = validPrior(indicator);
        final double pOne = 0.5;
        final IndicatorGibbsOperator op = operator(indicator, prior, pOne);

        assertProposalMatchesExpected(indicator, prior, op, pOne);
    }

    @Test
    public void asymmetricPOneIsIncludedInExpectedHastingsTerm() {
        Randomizer.setSeed(303L);

        final IntegerParameter indicator = integerParameter("0");
        final RelaxedRatesPriorSVS prior = validPrior(indicator);
        final double pOne = 0.73;
        final IndicatorGibbsOperator op = operator(indicator, prior, pOne);

        assertProposalMatchesExpected(indicator, prior, op, pOne);
    }

    @Test
    public void proposalMutatesOnlyIndicatorNotRatesOrHyperparameters() {
        Randomizer.setSeed(404L);

        final IntegerParameter indicator = integerParameter("0");
        final RealParameter rates = realParameter(POSITIVE_RATES);
        final RealParameter ucldStdev = realParameter("0.5");
        final RealParameter rootLogRate = realParameter("0.0");
        final RealParameter sigma2 = realParameter("0.2");
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), rates, indicator, ucldStdev, rootLogRate, sigma2);
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        final double[] ratesBefore = copyValues(rates);
        final double ucldBefore = ucldStdev.getValue(0);
        final double rootBefore = rootLogRate.getValue(0);
        final double sigmaBefore = sigma2.getValue(0);

        final double hr = op.proposal();

        assertTrue(Double.isFinite(hr));
        assertArrayEquals(ratesBefore, copyValues(rates), EPS);
        assertEquals(ucldBefore, ucldStdev.getValue(0), EPS);
        assertEquals(rootBefore, rootLogRate.getValue(0), EPS);
        assertEquals(sigmaBefore, sigma2.getValue(0), EPS);
        assertBinary(indicator.getValue(0));
    }

    @Test
    public void repeatedProposalsKeepIndicatorBinaryAndReturnFiniteHastingsTerms() {
        Randomizer.setSeed(505L);

        final IntegerParameter indicator = integerParameter("0");
        final RelaxedRatesPriorSVS prior = validPrior(indicator);
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        for (int i = 0; i < 100; i++) {
            final double hr = op.proposal();
            assertTrue(Double.isFinite(hr));
            assertBinary(indicator.getValue(0));
        }
    }

    private static void assertProposalMatchesExpected(final IntegerParameter indicator,
                                                      final RelaxedRatesPriorSVS prior,
                                                      final IndicatorGibbsOperator op,
                                                      final double pOne) {
        final int oldK = indicator.getValue(0);
        final Conditional conditional = conditional(prior, pOne);

        final double hr = op.proposal();
        final int newK = indicator.getValue(0);

        assertBinary(newK);
        assertEquals(conditional.expectedLogHastings(oldK, newK), hr, EPS);
    }

    private static Conditional conditional(final RelaxedRatesPriorSVS prior, final double pOne) {
        final double logW0 = Math.log(1.0 - pOne) + prior.logPriorUCOnly();
        final double logW1 = Math.log(pOne) + prior.logPriorACOnly();
        final double logDen = logSumExp(logW0, logW1);
        return new Conditional(logW0 - logDen, logW1 - logDen);
    }

    private static double logSumExp(final double a, final double b) {
        if (a > b) {
            return a + Math.log1p(Math.exp(b - a));
        }
        return b + Math.log1p(Math.exp(a - b));
    }

    private static final class Conditional {
        private final double logP0;
        private final double logP1;

        private Conditional(final double logP0, final double logP1) {
            this.logP0 = logP0;
            this.logP1 = logP1;
        }

        private double expectedLogHastings(final int oldK, final int newK) {
            final double logPold = oldK == 0 ? logP0 : logP1;
            final double logPnew = newK == 0 ? logP0 : logP1;
            return logPold - logPnew;
        }
    }

    private static IndicatorGibbsOperator operator(final IntegerParameter indicator,
                                                  final RelaxedRatesPriorSVS prior,
                                                  final double pOne) {
        final IndicatorGibbsOperator op = new IndicatorGibbsOperator();
        op.indicatorInput.setValue(indicator, op);
        op.priorInput.setValue(prior, op);
        op.pOneInput.setValue(pOne, op);
        op.initAndValidate();
        return op;
    }

    private static RelaxedRatesPriorSVS validPrior(final IntegerParameter indicator) {
        return prior(fixedTree(), realParameter(POSITIVE_RATES), indicator,
                realParameter("0.5"), realParameter("0.0"), realParameter("0.2"));
    }

    private static RelaxedRatesPriorSVS prior(final Tree tree,
                                              final RealParameter rates,
                                              final IntegerParameter indicator,
                                              final RealParameter ucldStdev,
                                              final RealParameter rootLogRate,
                                              final RealParameter sigma2) {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.initByName(
                "tree", tree,
                "rates", rates,
                "indicator", indicator,
                "ucldStdev", ucldStdev,
                "rootLogRate", rootLogRate,
                "sigma2", sigma2,
                "minBranchLength", 1.0e-12
        );
        return prior;
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealParameter realParameter(final String value) {
        return new RealParameter(value);
    }

    private static IntegerParameter integerParameter(final String value) {
        return new IntegerParameter(value);
    }

    private static double[] copyValues(final RealParameter parameter) {
        final double[] values = new double[parameter.getDimension()];
        for (int i = 0; i < values.length; i++) {
            values[i] = parameter.getValue(i);
        }
        return values;
    }

    private static void assertArrayEquals(final double[] expected, final double[] observed, final double eps) {
        assertEquals(expected.length, observed.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], observed[i], eps);
        }
    }

    private static void assertBinary(final int value) {
        assertTrue(value == 0 || value == 1);
    }
}
