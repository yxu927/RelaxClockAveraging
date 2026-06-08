package mixture.lphybeast.tobeast.generators;

import beast.base.evolution.tree.TreeParser;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import mixture.beast.evolution.mixture.RelaxedRatesPriorSVS;
import mixture.beast.evolution.operator.ACSigma2NonCenteredOperator;
import mixture.beast.evolution.operator.ACSubtreeUIncrementOperator;
import mixture.beast.evolution.operator.IndicatorGibbsOperator;
import mixture.beast.evolution.operator.SingleRateScaleOperator;
import mixture.beast.evolution.operator.SubtreeRateScaleOperator;
import mixture.beast.evolution.operator.UCACSwitchBridgeOperator;
import mixture.beast.evolution.operator.UCLDStdevNonCenteredOperator;
import mixture.beast.evolution.util.BranchRateIndexHelper;
import mixture.lphybeast.tobeast.TypedParameterUtils;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class SVSRawBranchRatesToBEASTOperatorTest {

    @Test
    public void randomIndicatorAddsFullSVSOperatorSchedule() {
        final Fixture fixture = fixture(0);

        final List<Operator> operators = operatorSchedule(fixture, true);

        assertHasOperator(operators, SingleRateScaleOperator.class);
        assertHasOperator(operators, SubtreeRateScaleOperator.class);
        assertHasOperator(operators, IndicatorGibbsOperator.class);
        assertHasOperator(operators, ACSubtreeUIncrementOperator.class);
        assertHasOperator(operators, UCLDStdevNonCenteredOperator.class);
        assertHasOperator(operators, ACSigma2NonCenteredOperator.class);
        assertHasOperator(operators, UCACSwitchBridgeOperator.class);
    }

    @Test
    public void fixedUCIndicatorDoesNotAddIndicatorSwitchOperators() {
        final Fixture fixture = fixture(0);

        final List<Operator> operators = operatorSchedule(fixture, false);

        assertHasOperator(operators, SingleRateScaleOperator.class);
        assertHasOperator(operators, SubtreeRateScaleOperator.class);
        assertHasOperator(operators, UCLDStdevNonCenteredOperator.class);
        assertNoOperator(operators, IndicatorGibbsOperator.class);
        assertNoOperator(operators, ACSubtreeUIncrementOperator.class);
        assertNoOperator(operators, ACSigma2NonCenteredOperator.class);
        assertNoOperator(operators, UCACSwitchBridgeOperator.class);
    }

    @Test
    public void fixedACIndicatorDoesNotAddIndicatorSwitchOperators() {
        final Fixture fixture = fixture(1);

        final List<Operator> operators = operatorSchedule(fixture, false);

        assertHasOperator(operators, SingleRateScaleOperator.class);
        assertHasOperator(operators, SubtreeRateScaleOperator.class);
        assertHasOperator(operators, ACSubtreeUIncrementOperator.class);
        assertHasOperator(operators, ACSigma2NonCenteredOperator.class);
        assertNoOperator(operators, IndicatorGibbsOperator.class);
        assertNoOperator(operators, UCLDStdevNonCenteredOperator.class);
        assertNoOperator(operators, UCACSwitchBridgeOperator.class);
    }

    @Test
    public void nodeIndexedLPhyRawRatesDropRootSlotBeforeTypedRates() {
        final TreeParser tree = beastTree();
        final int nodeCount = tree.getNodeCount();
        final double[] nodeIndexed = new double[nodeCount];
        for (int i = 0; i < nodeIndexed.length; i++) {
            nodeIndexed[i] = 1.0 + 0.1 * i;
        }
        nodeIndexed[tree.getRoot().getNr()] = 0.0;

        final RealParameter lphyRawRates = new RealParameter(join(nodeIndexed));
        lphyRawRates.setID("rawRates.nodeIndexed");

        final RealVectorParam<?> rates = SVSRawBranchRatesToBEAST.positiveRatesFromLPhyValue(
                "fallbackRates", lphyRawRates, tree, nodeCount - 1, true);

        final BranchRateIndexHelper.Mapping mapping = BranchRateIndexHelper.buildDeterministic(tree);
        assertEquals(nodeCount - 1, rates.size());
        for (int nr = 0; nr < nodeCount; nr++) {
            final int rateIndex = mapping.idxForNodeNr(nr);
            if (rateIndex >= 0) {
                assertEquals(nodeIndexed[nr], rates.get(rateIndex), 1.0e-12);
            }
        }
    }

    private static List<Operator> operatorSchedule(final Fixture fixture, final boolean indicatorIsRandom) {
        return SVSRawBranchRatesToBEAST.createOperatorSchedule(
                "fixture",
                fixture.tree,
                fixture.rates,
                fixture.indicator,
                indicatorIsRandom,
                fixture.prior,
                fixture.ucldStdev,
                fixture.sigma2,
                fixture.rootLogRate
        );
    }

    private static Fixture fixture(final int indicatorValue) {
        final TreeParser tree = beastTree();
        final RealVectorParam<?> rates = TypedParameterUtils.positiveRealVector(
                "rawRates.fixture", new double[]{1.0, 1.0, 1.0, 1.0}, true);
        final IntScalarParam<?> indicator = TypedParameterUtils.intScalar("indicator.fixture", indicatorValue, true);
        final RealParameter ucldStdev = realParameter("ucldStdev", 0.4);
        final RealParameter sigma2 = realParameter("sigma2", 0.2);
        final RealParameter rootLogRate = realParameter("rootLogRate", 0.0);
        final RelaxedRatesPriorSVS prior = prior(tree, rates, indicator, ucldStdev, sigma2, rootLogRate);
        return new Fixture(tree, rates, indicator, ucldStdev, sigma2, rootLogRate, prior);
    }

    private static TreeParser beastTree() {
        return new TreeParser("((A:2.0,B:2.0):3.0,C:5.0);", false, true, true, 1);
    }

    private static RealParameter realParameter(final String id, final double value) {
        final RealParameter parameter = new RealParameter(Double.toString(value));
        parameter.setID(id);
        parameter.initAndValidate();
        return parameter;
    }

    private static RelaxedRatesPriorSVS prior(final TreeParser tree,
                                              final RealVectorParam<?> rates,
                                              final IntScalarParam<?> indicator,
                                              final RealParameter ucldStdev,
                                              final RealParameter sigma2,
                                              final RealParameter rootLogRate) {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.setID("SVSRelaxedClockPrior.fixture");
        prior.setInputValue("tree", tree);
        prior.setInputValue("ratesVector", rates);
        prior.setInputValue("indicatorScalar", indicator);
        prior.setInputValue("ucldStdev", ucldStdev);
        prior.setInputValue("sigma2", sigma2);
        prior.setInputValue("rootLogRate", rootLogRate);
        prior.initAndValidate();
        return prior;
    }

    private static void assertHasOperator(final List<Operator> operators,
                                          final Class<? extends Operator> operatorClass) {
        assertTrue(operatorClass.getSimpleName() + " should be generated",
                operators.stream().anyMatch(operatorClass::isInstance));
    }

    private static void assertNoOperator(final List<Operator> operators,
                                         final Class<? extends Operator> operatorClass) {
        assertFalse(operatorClass.getSimpleName() + " should not be generated",
                operators.stream().anyMatch(operatorClass::isInstance));
    }

    private static String join(final double[] values) {
        final StringBuilder builder = new StringBuilder();
        for (final double value : values) {
            if (!builder.isEmpty()) {
                builder.append(' ');
            }
            builder.append(value);
        }
        return builder.toString();
    }

    private record Fixture(TreeParser tree,
                           RealVectorParam<?> rates,
                           IntScalarParam<?> indicator,
                           RealParameter ucldStdev,
                           RealParameter sigma2,
                           RealParameter rootLogRate,
                           RelaxedRatesPriorSVS prior) {
    }
}
