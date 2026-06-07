package mixture.beast.evolution.operator;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Int;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import mixture.beast.evolution.mixture.RelaxedRatesPriorSVS;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class IndicatorGibbsOperatorTypedInputTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";
    private static final String POSITIVE_RATES = "0.8 1.1 1.4 0.9";

    @Test
    public void typedDefaultPOneIsOneHalfAndProposalUsesIt() {
        Randomizer.setSeed(12345L);

        final IntScalarParam<Int> indicator = indicator(0);
        final RelaxedRatesPriorSVS prior = typedPrior(indicator);
        final IndicatorGibbsOperator op = operatorWithDefaultPOne(indicator, prior);

        assertProposalMatchesExpected(indicator, prior, op, 0.5);
    }

    @Test
    public void typedProposalFromOldZeroReturnsExpectedGibbsHastingsTerm() {
        Randomizer.setSeed(101L);

        final IntScalarParam<Int> indicator = indicator(0);
        final RelaxedRatesPriorSVS prior = typedPrior(indicator);
        final double pOne = 0.5;
        final IndicatorGibbsOperator op = operator(indicator, prior, pOne);

        assertProposalMatchesExpected(indicator, prior, op, pOne);
    }

    @Test
    public void typedProposalFromOldOneReturnsExpectedGibbsHastingsTerm() {
        Randomizer.setSeed(202L);

        final IntScalarParam<Int> indicator = indicator(1);
        final RelaxedRatesPriorSVS prior = typedPrior(indicator);
        final double pOne = 0.5;
        final IndicatorGibbsOperator op = operator(indicator, prior, pOne);

        assertProposalMatchesExpected(indicator, prior, op, pOne);
    }

    @Test
    public void typedAsymmetricPOneIsIncludedInExpectedHastingsTerm() {
        Randomizer.setSeed(303L);

        final IntScalarParam<Int> indicator = indicator(0);
        final RelaxedRatesPriorSVS prior = typedPrior(indicator);
        final double pOne = 0.73;
        final IndicatorGibbsOperator op = operator(indicator, prior, pOne);

        assertProposalMatchesExpected(indicator, prior, op, pOne);
    }

    @Test
    public void typedInvalidOldIndicatorReturnsNegativeInfinityAndLeavesIndicatorUnchanged() {
        final IntScalarParam<Int> indicator = indicator(2);
        final RelaxedRatesPriorSVS prior = typedPrior(indicator);
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertEquals(2, indicator.get());
    }

    @Test
    public void typedBothPriorWeightsNegativeInfinityReturnsNegativeInfinityAndLeavesIndicatorUnchanged() {
        final IntScalarParam<Int> indicator = indicator(0);
        final RelaxedRatesPriorSVS prior = priorWithTypedIndicatorAndLegacyRates(indicator, "0.0 0.0 0.0 0.0");
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertEquals(0, indicator.get());
    }

    @Test
    public void typedProposalMutatesOnlyIndicatorNotRatesOrHyperparameters() {
        Randomizer.setSeed(404L);

        final IntScalarParam<Int> indicator = indicator(0);
        final RealVectorParam<PositiveReal> rates = rates(0.8, 1.1, 1.4, 0.9);
        final RealScalarParam<PositiveReal> ucldStdev = positiveScalar(0.5);
        final RealScalarParam<Real> rootLogRate = realScalar(0.0);
        final RealScalarParam<PositiveReal> sigma2 = positiveScalar(0.2);
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates, indicator, ucldStdev,
                rootLogRate, sigma2);
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        final double[] ratesBefore = copyValues(rates);
        final double ucldBefore = ucldStdev.get();
        final double rootBefore = rootLogRate.get();
        final double sigmaBefore = sigma2.get();

        final double hr = op.proposal();

        assertTrue(Double.isFinite(hr));
        assertArrayEquals(ratesBefore, copyValues(rates), EPS);
        assertEquals(ucldBefore, ucldStdev.get(), EPS);
        assertEquals(rootBefore, rootLogRate.get(), EPS);
        assertEquals(sigmaBefore, sigma2.get(), EPS);
        assertBinary(indicator.get());
    }

    @Test
    public void typedRepeatedProposalsKeepIndicatorBinaryAndReturnFiniteHastingsTerms() {
        Randomizer.setSeed(505L);

        final IntScalarParam<Int> indicator = indicator(0);
        final RelaxedRatesPriorSVS prior = typedPrior(indicator);
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        for (int i = 0; i < 100; i++) {
            final double hr = op.proposal();
            assertTrue(Double.isFinite(hr));
            assertBinary(indicator.get());
        }
    }

    @Test
    public void typedOperatorAndTypedPriorShareIndicatorState() {
        Randomizer.setSeed(606L);

        final IntScalarParam<Int> indicator = indicator(0);
        final RelaxedRatesPriorSVS prior = typedPrior(indicator);
        final IndicatorGibbsOperator op = operator(indicator, prior, 0.5);

        op.proposal();
        final int newK = indicator.get();

        if (newK == 0) {
            assertEquals(prior.logPriorUCOnly(), prior.calculateLogP(), EPS);
        } else {
            assertEquals(prior.logPriorACOnly(), prior.calculateLogP(), EPS);
        }
    }

    @Test
    public void typedAndLegacyOperatorsAreEquivalentForSameRandomDraw() {
        final long seed = 707L;
        final double pOne = 0.61;

        final IntegerParameter legacyIndicator = new IntegerParameter("0");
        final RelaxedRatesPriorSVS legacyPrior = legacyPrior(legacyIndicator);
        final IndicatorGibbsOperator legacyOperator = legacyOperator(legacyIndicator, legacyPrior, pOne);

        Randomizer.setSeed(seed);
        final double legacyHr = legacyOperator.proposal();
        final int legacyNewK = legacyIndicator.getValue(0);

        final IntScalarParam<Int> typedIndicator = indicator(0);
        final RelaxedRatesPriorSVS typedPrior = typedPrior(typedIndicator);
        final IndicatorGibbsOperator typedOperator = operator(typedIndicator, typedPrior, pOne);

        Randomizer.setSeed(seed);
        final double typedHr = typedOperator.proposal();
        final int typedNewK = typedIndicator.get();

        assertEquals(legacyNewK, typedNewK);
        assertEquals(legacyHr, typedHr, EPS);
    }

    @Test
    public void legacyAndTypedIndicatorCannotBothBeSpecified() {
        final IndicatorGibbsOperator op = new IndicatorGibbsOperator();
        op.indicatorInput.setValue(new IntegerParameter("0"), op);
        op.indicatorScalarInput.setValue(indicator(0), op);
        op.priorInput.setValue(typedPrior(indicator(0)), op);

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void missingIndicatorIsRejected() {
        final IndicatorGibbsOperator op = new IndicatorGibbsOperator();
        op.priorInput.setValue(typedPrior(indicator(0)), op);

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void typedInvalidPOneValidationIsPreserved() {
        assertThrows(IllegalArgumentException.class,
                () -> operator(indicator(0), typedPrior(indicator(0)), 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> operator(indicator(0), typedPrior(indicator(0)), 1.0));
        assertThrows(IllegalArgumentException.class,
                () -> operator(indicator(0), typedPrior(indicator(0)), -0.1));
        assertThrows(IllegalArgumentException.class,
                () -> operator(indicator(0), typedPrior(indicator(0)), 1.1));
    }

    private static void assertProposalMatchesExpected(final IntScalarParam<Int> indicator,
                                                      final RelaxedRatesPriorSVS prior,
                                                      final IndicatorGibbsOperator op,
                                                      final double pOne) {
        final int oldK = indicator.get();
        final Conditional conditional = conditional(prior, pOne);

        final double hr = op.proposal();
        final int newK = indicator.get();

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

    private static IndicatorGibbsOperator operatorWithDefaultPOne(final IntScalarParam<Int> indicator,
                                                                  final RelaxedRatesPriorSVS prior) {
        final IndicatorGibbsOperator op = new IndicatorGibbsOperator();
        op.indicatorScalarInput.setValue(indicator, op);
        op.priorInput.setValue(prior, op);
        op.initAndValidate();
        return op;
    }

    private static IndicatorGibbsOperator operator(final IntScalarParam<Int> indicator,
                                                  final RelaxedRatesPriorSVS prior,
                                                  final double pOne) {
        final IndicatorGibbsOperator op = new IndicatorGibbsOperator();
        op.indicatorScalarInput.setValue(indicator, op);
        op.priorInput.setValue(prior, op);
        op.pOneInput.setValue(pOne, op);
        op.initAndValidate();
        return op;
    }

    private static IndicatorGibbsOperator legacyOperator(final IntegerParameter indicator,
                                                        final RelaxedRatesPriorSVS prior,
                                                        final double pOne) {
        final IndicatorGibbsOperator op = new IndicatorGibbsOperator();
        op.indicatorInput.setValue(indicator, op);
        op.priorInput.setValue(prior, op);
        op.pOneInput.setValue(pOne, op);
        op.initAndValidate();
        return op;
    }

    private static RelaxedRatesPriorSVS typedPrior(final IntScalarParam<Int> indicator) {
        return typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9), indicator,
                positiveScalar(0.5), realScalar(0.0), positiveScalar(0.2));
    }

    private static RelaxedRatesPriorSVS typedPrior(final Tree tree,
                                                  final RealVectorParam<PositiveReal> rates,
                                                  final IntScalarParam<Int> indicator,
                                                  final RealScalarParam<PositiveReal> ucldStdev,
                                                  final RealScalarParam<Real> rootLogRate,
                                                  final RealScalarParam<PositiveReal> sigma2) {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.initByName(
                "tree", tree,
                "ratesVector", rates,
                "indicatorScalar", indicator,
                "ucldStdevScalar", ucldStdev,
                "rootLogRateScalar", rootLogRate,
                "sigma2Scalar", sigma2,
                "minBranchLength", 1.0e-12
        );
        return prior;
    }

    private static RelaxedRatesPriorSVS priorWithTypedIndicatorAndLegacyRates(final IntScalarParam<Int> indicator,
                                                                              final String rates) {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.initByName(
                "tree", fixedTree(),
                "rates", new RealParameter(rates),
                "indicatorScalar", indicator,
                "ucldStdevScalar", positiveScalar(0.5),
                "rootLogRateScalar", realScalar(0.0),
                "sigma2Scalar", positiveScalar(0.2),
                "minBranchLength", 1.0e-12
        );
        return prior;
    }

    private static RelaxedRatesPriorSVS legacyPrior(final IntegerParameter indicator) {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.initByName(
                "tree", fixedTree(),
                "rates", new RealParameter(POSITIVE_RATES),
                "indicator", indicator,
                "ucldStdev", new RealParameter("0.5"),
                "rootLogRate", new RealParameter("0.0"),
                "sigma2", new RealParameter("0.2"),
                "minBranchLength", 1.0e-12
        );
        return prior;
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealVectorParam<PositiveReal> rates(final double... values) {
        return new RealVectorParam<>(values, PositiveReal.INSTANCE);
    }

    private static IntScalarParam<Int> indicator(final int value) {
        return new IntScalarParam<>(value, Int.INSTANCE);
    }

    private static RealScalarParam<PositiveReal> positiveScalar(final double value) {
        return new RealScalarParam<>(value, PositiveReal.INSTANCE);
    }

    private static RealScalarParam<Real> realScalar(final double value) {
        return new RealScalarParam<>(value, Real.INSTANCE);
    }

    private static double[] copyValues(final RealVectorParam<PositiveReal> parameter) {
        final double[] values = new double[parameter.size()];
        for (int i = 0; i < values.length; i++) {
            values[i] = parameter.get(i);
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
