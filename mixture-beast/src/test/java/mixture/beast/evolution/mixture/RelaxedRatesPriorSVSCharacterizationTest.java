package mixture.beast.evolution.mixture;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import mixture.beast.evolution.util.BranchRateIndexHelper;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

public class RelaxedRatesPriorSVSCharacterizationTest {

    private static final double EPS = 1.0e-10;
    private static final double LOG_2PI = Math.log(2.0 * Math.PI);
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";
    private static final String POSITIVE_RATES = "0.8 1.1 1.4 0.9";

    @Test
    public void indicatorZeroCalculateLogPMatchesUCOnly() {
        final RealParameter rates = realParameter(POSITIVE_RATES);
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), rates, integerParameter("0"),
                realParameter("0.5"), null, realParameter("0.2"), 1.0e-12);

        assertEquals(prior.logPriorUCOnly(), prior.calculateLogP(), EPS);
        assertEquals(expectedUCLogPrior(rates, 0.5), prior.calculateLogP(), EPS);
    }

    @Test
    public void indicatorOneCalculateLogPMatchesACOnly() {
        final Tree tree = fixedTree();
        final RealParameter rates = realParameter(POSITIVE_RATES);
        final RealParameter sigma2 = realParameter("0.2");
        final RealParameter rootLogRate = realParameter("0.0");
        final RelaxedRatesPriorSVS prior = prior(tree, rates, integerParameter("1"),
                realParameter("0.5"), rootLogRate, sigma2, 1.0e-12);

        assertEquals(prior.logPriorACOnly(), prior.calculateLogP(), EPS);
        assertEquals(expectedACLogPrior(tree, rates, sigma2.getValue(0), rootLogRate.getValue(0), 1.0e-12),
                prior.calculateLogP(), EPS);
    }

    @Test
    public void invalidIndicatorReturnsNegativeInfinity() {
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("2"), realParameter("0.5"), null, realParameter("0.2"), 1.0e-12);

        assertEquals(Double.NEGATIVE_INFINITY, prior.calculateLogP(), 0.0);
    }

    @Test
    public void ucPriorRejectsNonPositiveStdev() {
        final RelaxedRatesPriorSVS zero = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("0"), realParameter("0.0"), null, realParameter("0.2"), 1.0e-12);
        final RelaxedRatesPriorSVS negative = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("0"), realParameter("-0.1"), null, realParameter("0.2"), 1.0e-12);

        assertEquals(Double.NEGATIVE_INFINITY, zero.logPriorUCOnly(), 0.0);
        assertEquals(Double.NEGATIVE_INFINITY, negative.logPriorUCOnly(), 0.0);
    }

    @Test
    public void acPriorRejectsNonPositiveSigma2() {
        final RelaxedRatesPriorSVS zero = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), null, realParameter("0.0"), 1.0e-12);
        final RelaxedRatesPriorSVS negative = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), null, realParameter("-0.1"), 1.0e-12);

        assertEquals(Double.NEGATIVE_INFINITY, zero.logPriorACOnly(), 0.0);
        assertEquals(Double.NEGATIVE_INFINITY, negative.logPriorACOnly(), 0.0);
    }

    @Test
    public void acPriorRejectsNonPositiveMinBranchLength() {
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), null, realParameter("0.2"), 0.0);

        assertEquals(Double.NEGATIVE_INFINITY, prior.logPriorACOnly(), 0.0);
    }

    @Test
    public void acPriorRejectsBranchesAtOrBelowMinBranchLength() {
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), null, realParameter("0.2"), 5.0);

        assertEquals(Double.NEGATIVE_INFINITY, prior.logPriorACOnly(), 0.0);
    }

    @Test
    public void priorsRejectNonPositiveRates() {
        final RelaxedRatesPriorSVS uc = prior(fixedTree(), realParameter("0.8 0.0 1.4 0.9"),
                integerParameter("0"), realParameter("0.5"), null, realParameter("0.2"), 1.0e-12);
        final RelaxedRatesPriorSVS ac = prior(fixedTree(), realParameter("0.8 0.0 1.4 0.9"),
                integerParameter("1"), realParameter("0.5"), null, realParameter("0.2"), 1.0e-12);

        assertEquals(Double.NEGATIVE_INFINITY, uc.logPriorUCOnly(), 0.0);
        assertEquals(Double.NEGATIVE_INFINITY, ac.logPriorACOnly(), 0.0);
    }

    @Test
    public void missingRootLogRateEqualsExplicitZeroRootLogRate() {
        final Tree tree = fixedTree();
        final RelaxedRatesPriorSVS missing = prior(tree, realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), null, realParameter("0.2"), 1.0e-12);
        final RelaxedRatesPriorSVS explicitZero = prior(tree, realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), realParameter("0.0"),
                realParameter("0.2"), 1.0e-12);

        assertEquals(missing.logPriorACOnly(), explicitZero.logPriorACOnly(), EPS);
    }

    @Test
    public void rootLogRateChangesACPriorWhenRootChildrenExist() {
        final Tree tree = fixedTree();
        final RelaxedRatesPriorSVS zero = prior(tree, realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), realParameter("0.0"),
                realParameter("0.2"), 1.0e-12);
        final RelaxedRatesPriorSVS shifted = prior(tree, realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), realParameter("0.7"),
                realParameter("0.2"), 1.0e-12);

        assertNotEquals(zero.logPriorACOnly(), shifted.logPriorACOnly(), EPS);
    }

    @Test
    public void scalarRatesExpandToOneRatePerNonRootBranch() {
        final Tree tree = fixedTree();
        final RealParameter rates = realParameter("1.25");
        final RelaxedRatesPriorSVS prior = prior(tree, rates, integerParameter("0"),
                realParameter("0.5"), null, realParameter("0.2"), 1.0e-12);

        assertEquals(tree.getNodeCount() - 1, rates.getDimension());
        assertTrue(Double.isFinite(prior.calculateLogP()));
    }

    @Test
    public void requiresRecalculationRespondsToDirtyRates() {
        final RealParameter rates = realParameter(POSITIVE_RATES);
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), rates, integerParameter("0"),
                realParameter("0.5"), realParameter("0.0"), realParameter("0.2"), 1.0e-12);

        rates.setValue(0, 0.95);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void requiresRecalculationRespondsToDirtyIndicator() {
        final IntegerParameter indicator = integerParameter("0");
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter(POSITIVE_RATES), indicator,
                realParameter("0.5"), realParameter("0.0"), realParameter("0.2"), 1.0e-12);

        indicator.setValue(0, 1);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void requiresRecalculationRespondsToDirtyUcldStdev() {
        final RealParameter ucldStdev = realParameter("0.5");
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("0"), ucldStdev, realParameter("0.0"), realParameter("0.2"), 1.0e-12);

        ucldStdev.setValue(0, 0.6);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void requiresRecalculationRespondsToDirtySigma2() {
        final RealParameter sigma2 = realParameter("0.2");
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), realParameter("0.0"), sigma2, 1.0e-12);

        sigma2.setValue(0, 0.3);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void requiresRecalculationRespondsToDirtyRootLogRate() {
        final RealParameter rootLogRate = realParameter("0.0");
        final RelaxedRatesPriorSVS prior = prior(fixedTree(), realParameter(POSITIVE_RATES),
                integerParameter("1"), realParameter("0.5"), rootLogRate, realParameter("0.2"), 1.0e-12);

        rootLogRate.setValue(0, 0.1);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void getArgumentsAndConditionsRecordCurrentMetadataBehaviour() {
        final Tree tree = fixedTree();
        tree.setID("tree");
        final RealParameter rates = realParameter(POSITIVE_RATES);
        rates.setID("rates");
        final IntegerParameter indicator = integerParameter("0");
        indicator.setID("indicator");
        final RealParameter ucldStdev = realParameter("0.5");
        ucldStdev.setID("ucldStdev");
        final RealParameter rootLogRate = realParameter("0.0");
        rootLogRate.setID("rootLogRate");
        final RealParameter sigma2 = realParameter("0.2");
        sigma2.setID("sigma2");

        final RelaxedRatesPriorSVS prior = prior(tree, rates, indicator, ucldStdev,
                rootLogRate, sigma2, 1.0e-12);

        final List<String> arguments = prior.getArguments();
        final List<String> conditions = prior.getConditions();

        assertTrue(arguments.contains("rates"));
        assertTrue(conditions.contains("tree"));
        assertTrue(conditions.contains("indicator"));
        assertTrue(conditions.contains("ucldStdev"));
        assertTrue(conditions.contains("rootLogRate"));
        assertTrue(conditions.contains("sigma2"));
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

    private static RelaxedRatesPriorSVS prior(final Tree tree,
                                              final RealParameter rates,
                                              final IntegerParameter indicator,
                                              final RealParameter ucldStdev,
                                              final RealParameter rootLogRate,
                                              final RealParameter sigma2,
                                              final double minBranchLength) {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        if (rootLogRate == null) {
            prior.initByName(
                    "tree", tree,
                    "rates", rates,
                    "indicator", indicator,
                    "ucldStdev", ucldStdev,
                    "sigma2", sigma2,
                    "minBranchLength", minBranchLength
            );
        } else {
            prior.initByName(
                    "tree", tree,
                    "rates", rates,
                    "indicator", indicator,
                    "ucldStdev", ucldStdev,
                    "rootLogRate", rootLogRate,
                    "sigma2", sigma2,
                    "minBranchLength", minBranchLength
            );
        }
        return prior;
    }

    private static double expectedUCLogPrior(final RealParameter rates, final double s) {
        final double var = s * s;
        final double meanLog = -0.5 * var;
        double lp = 0.0;
        for (int i = 0; i < rates.getDimension(); i++) {
            lp += logLogNormalRate(rates.getValue(i), meanLog, var);
        }
        return lp;
    }

    private static double expectedACLogPrior(final Tree tree,
                                             final RealParameter rates,
                                             final double sigma2,
                                             final double rootLogRate,
                                             final double minBranchLength) {
        final BranchRateIndexHelper.Mapping mapping = BranchRateIndexHelper.buildDeterministic(tree);
        double lp = 0.0;
        for (int i = 0; i < mapping.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            if (node.isRoot()) {
                continue;
            }

            final double dt = node.getLength();
            if (!(dt > minBranchLength)) {
                return Double.NEGATIVE_INFINITY;
            }

            final double var = sigma2 * dt;
            final int idxChild = mapping.idxForNode(node);
            final double childRate = rates.getValue(idxChild);

            final Node parent = node.getParent();
            final double parentLogRate;
            if (parent.isRoot()) {
                parentLogRate = rootLogRate;
            } else {
                parentLogRate = Math.log(rates.getValue(mapping.idxForNode(parent)));
            }

            lp += logLogNormalRate(childRate, parentLogRate - 0.5 * var, var);
        }
        return lp;
    }

    private static double logLogNormalRate(final double r, final double meanLog, final double var) {
        if (!(r > 0.0) || !(var > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }
        final double x = Math.log(r);
        final double z = x - meanLog;
        return -x - 0.5 * (LOG_2PI + Math.log(var) + (z * z) / var);
    }
}
