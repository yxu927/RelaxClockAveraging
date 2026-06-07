package mixture.beast.evolution.mixture;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import mixture.beast.evolution.util.BranchRateIndexHelper;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

public class SharedRatesClockModelCharacterizationTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";
    private static final String FOUR_BRANCH_RATES = "0.5 1.0 2.0 4.0";

    @Test
    public void rootRateIsCurrentSentinelValue() {
        final Tree tree = fixedTree();
        final SharedRatesClockModel clock = clock(tree, realParameter(FOUR_BRANCH_RATES),
                false, realParameter("3.0"));

        assertEquals(1.0, clock.getRateForBranch(tree.getRoot()), EPS);
    }

    @Test
    public void normalizeFalseReturnsRawRatesTimesMeanRateForNonRootBranches() {
        final Tree tree = fixedTree();
        final RealParameter rates = realParameter(FOUR_BRANCH_RATES);
        final RealParameter meanRate = realParameter("3.0");
        final SharedRatesClockModel clock = clock(tree, rates, false, meanRate);
        final BranchRateIndexHelper.Mapping mapping = BranchRateIndexHelper.buildDeterministic(tree);

        for (int i = 0; i < tree.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            if (!node.isRoot()) {
                final int rateIndex = mapping.idxForNode(node);
                assertEquals(rates.getValue(rateIndex) * meanRate.getValue(),
                        clock.getRateForBranch(node), EPS);
            }
        }
    }

    @Test
    public void normalizeTrueMakesTimeWeightedMeanRateEqualMeanRate() {
        final Tree tree = fixedTree();
        final RealParameter meanRate = realParameter("3.0");
        final SharedRatesClockModel clock = clock(tree, realParameter(FOUR_BRANCH_RATES),
                true, meanRate);

        assertEquals(meanRate.getValue(), timeWeightedMeanRate(tree, clock), EPS);
    }

    @Test
    public void normalizeTrueRespondsToRawRateChangesAfterRecalculation() {
        final Tree tree = fixedTree();
        final RealParameter rates = realParameter(FOUR_BRANCH_RATES);
        final RealParameter meanRate = realParameter("3.0");
        final SharedRatesClockModel clock = clock(tree, rates, true, meanRate);
        final Node firstMappedBranch = nodeForRateIndex(tree, 0);
        final double branchRateBefore = clock.getRateForBranch(firstMappedBranch);

        rates.setValue(0, 10.0);
        rates.setValue(1, 0.25);

        assertTrue(clock.requiresRecalculation());
        assertNotEquals(branchRateBefore, clock.getRateForBranch(firstMappedBranch), EPS);
        assertEquals(meanRate.getValue(), timeWeightedMeanRate(tree, clock), EPS);
    }

    @Test
    public void normalizeTrueRespondsToMeanRateChangesAfterRecalculation() {
        final Tree tree = fixedTree();
        final RealParameter meanRate = realParameter("3.0");
        final SharedRatesClockModel clock = clock(tree, realParameter(FOUR_BRANCH_RATES),
                true, meanRate);

        meanRate.setValue(0, 5.0);

        assertTrue(clock.requiresRecalculation());
        assertEquals(5.0, timeWeightedMeanRate(tree, clock), EPS);
    }

    @Test
    public void scalarRatesExpandToOneRatePerNonRootBranch() {
        final Tree tree = fixedTree();
        final RealParameter rates = realParameter("2.5");
        final RealParameter meanRate = realParameter("2.0");
        final SharedRatesClockModel clock = clock(tree, rates, false, meanRate);

        assertEquals(tree.getNodeCount() - 1, rates.getDimension());
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            if (!node.isRoot()) {
                assertEquals(5.0, clock.getRateForBranch(node), EPS);
            }
        }
    }

    @Test
    public void storeRestoreWithRestoredInputsPreservesObservableRates() {
        final Tree tree = fixedTree();
        final RealParameter rates = realParameter(FOUR_BRANCH_RATES);
        final RealParameter meanRate = realParameter("3.0");
        final SharedRatesClockModel clock = clock(tree, rates, true, meanRate);
        final Node firstMappedBranch = nodeForRateIndex(tree, 0);
        final double branchRateBefore = clock.getRateForBranch(firstMappedBranch);
        final double weightedMeanBefore = timeWeightedMeanRate(tree, clock);

        clock.store();
        rates.setValue(0, 20.0);
        assertTrue(clock.requiresRecalculation());
        assertNotEquals(branchRateBefore, clock.getRateForBranch(firstMappedBranch), EPS);

        rates.restore();
        clock.restore();

        assertEquals(branchRateBefore, clock.getRateForBranch(firstMappedBranch), EPS);
        assertEquals(weightedMeanBefore, timeWeightedMeanRate(tree, clock), EPS);
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealParameter realParameter(final String value) {
        return new RealParameter(value);
    }

    private static SharedRatesClockModel clock(final Tree tree,
                                               final RealParameter rates,
                                               final boolean normalize,
                                               final RealParameter meanRate) {
        final SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.initByName("tree", tree, "rates", rates, "normalize", normalize, "meanRate", meanRate);
        return clock;
    }

    private static double timeWeightedMeanRate(final Tree tree, final SharedRatesClockModel clock) {
        double sumRateTime = 0.0;
        double sumTime = 0.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            if (node.isRoot()) {
                continue;
            }
            final double branchLength = node.getLength();
            if (branchLength > 0.0) {
                sumRateTime += clock.getRateForBranch(node) * branchLength;
                sumTime += branchLength;
            }
        }
        return sumRateTime / sumTime;
    }

    private static Node nodeForRateIndex(final Tree tree, final int rateIndex) {
        final BranchRateIndexHelper.Mapping mapping = BranchRateIndexHelper.buildDeterministic(tree);
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            if (!node.isRoot() && mapping.idxForNode(node) == rateIndex) {
                return node;
            }
        }
        throw new IllegalArgumentException("No node maps to rate index " + rateIndex);
    }
}
