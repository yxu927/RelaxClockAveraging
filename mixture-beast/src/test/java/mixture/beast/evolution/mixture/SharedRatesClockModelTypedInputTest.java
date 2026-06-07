package mixture.beast.evolution.mixture;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import mixture.beast.evolution.util.BranchRateIndexHelper;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class SharedRatesClockModelTypedInputTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";

    @Test
    public void typedNormalizeFalseReturnsRawRatesTimesMeanRateForNonRootBranches() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = realVectorParam(0.5, 1.0, 2.0, 4.0);
        final RealScalarParam<PositiveReal> meanRate = realScalarParam(3.0);
        final SharedRatesClockModel clock = typedClock(tree, rates, false, meanRate);
        final BranchRateIndexHelper.Mapping mapping = BranchRateIndexHelper.buildDeterministic(tree);

        for (int i = 0; i < tree.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            if (!node.isRoot()) {
                final int rateIndex = mapping.idxForNode(node);
                assertEquals(rates.get(rateIndex) * meanRate.get(),
                        clock.getRateForBranch(node), EPS);
            }
        }
    }

    @Test
    public void typedNormalizeTrueMakesTimeWeightedMeanRateEqualMeanRate() {
        final Tree tree = fixedTree();
        final RealScalarParam<PositiveReal> meanRate = realScalarParam(3.0);
        final SharedRatesClockModel clock = typedClock(tree,
                realVectorParam(0.5, 1.0, 2.0, 4.0), true, meanRate);

        assertEquals(meanRate.get(), timeWeightedMeanRate(tree, clock), EPS);
    }

    @Test
    public void typedRawRateMutationTriggersRecalculationAndKeepsWeightedMean() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = realVectorParam(0.5, 1.0, 2.0, 4.0);
        final RealScalarParam<PositiveReal> meanRate = realScalarParam(3.0);
        final SharedRatesClockModel clock = typedClock(tree, rates, true, meanRate);

        rates.set(0, 10.0);
        rates.set(1, 0.25);

        assertTrue(clock.requiresRecalculation());
        assertEquals(meanRate.get(), timeWeightedMeanRate(tree, clock), EPS);
    }

    @Test
    public void typedMeanRateMutationTriggersRecalculationAndUpdatesWeightedMean() {
        final Tree tree = fixedTree();
        final RealScalarParam<PositiveReal> meanRate = realScalarParam(3.0);
        final SharedRatesClockModel clock = typedClock(tree,
                realVectorParam(0.5, 1.0, 2.0, 4.0), true, meanRate);

        meanRate.set(5.0);

        assertTrue(clock.requiresRecalculation());
        assertEquals(5.0, timeWeightedMeanRate(tree, clock), EPS);
    }

    @Test
    public void typedScalarRatesExpandOrBroadcastToNonRootBranches() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = realVectorParam(2.5);
        final RealScalarParam<PositiveReal> meanRate = realScalarParam(2.0);
        final SharedRatesClockModel clock = typedClock(tree, rates, false, meanRate);

        assertTrue(rates.size() == tree.getNodeCount() - 1 || rates.size() == 1);
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            if (!node.isRoot()) {
                assertEquals(5.0, clock.getRateForBranch(node), EPS);
            }
        }
    }

    @Test
    public void typedAndLegacyInputsCannotBothBeSpecified() {
        final SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.treeInput.setValue(fixedTree(), clock);
        clock.ratesInput.setValue(new RealParameter("1.0 1.0 1.0 1.0"), clock);
        clock.ratesVectorInput.setValue(realVectorParam(1.0, 1.0, 1.0, 1.0), clock);

        assertThrows(IllegalArgumentException.class, clock::initAndValidate);
    }

    @Test
    public void noRatesInputIsRejected() {
        final SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.treeInput.setValue(fixedTree(), clock);

        assertThrows(IllegalArgumentException.class, clock::initAndValidate);
    }

    @Test
    public void typedMeanRateAndLegacyMeanRateCannotBothBeSpecified() {
        final SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.treeInput.setValue(fixedTree(), clock);
        clock.ratesVectorInput.setValue(realVectorParam(1.0, 1.0, 1.0, 1.0), clock);
        clock.meanRateInput.setValue(new RealParameter("1.0"), clock);
        clock.meanRateScalarInput.setValue(realScalarParam(1.0), clock);

        assertThrows(IllegalArgumentException.class, clock::initAndValidate);
    }

    @Test
    public void missingMeanRateDefaultsToOne() {
        final Tree constantTree = fixedTree();
        final SharedRatesClockModel constantClock = typedClockWithoutMeanRate(constantTree,
                realVectorParam(1.0, 1.0, 1.0, 1.0), false);

        for (int i = 0; i < constantTree.getNodeCount(); i++) {
            final Node node = constantTree.getNode(i);
            if (!node.isRoot()) {
                assertEquals(1.0, constantClock.getRateForBranch(node), EPS);
            }
        }

        final Tree normalizedTree = fixedTree();
        final SharedRatesClockModel normalizedClock = typedClockWithoutMeanRate(normalizedTree,
                realVectorParam(0.5, 1.0, 2.0, 4.0), true);

        assertEquals(1.0, timeWeightedMeanRate(normalizedTree, normalizedClock), EPS);
    }

    @Test
    public void existingLegacyPathStillWorks() {
        final Tree tree = fixedTree();
        final RealParameter rates = new RealParameter("0.5 1.0 2.0 4.0");
        final RealParameter meanRate = new RealParameter("3.0");
        final SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.initByName("tree", tree, "rates", rates, "normalize", true, "meanRate", meanRate);

        assertEquals(meanRate.getValue(), timeWeightedMeanRate(tree, clock), EPS);
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealVectorParam<PositiveReal> realVectorParam(final double... values) {
        return new RealVectorParam<>(values, PositiveReal.INSTANCE);
    }

    private static RealScalarParam<PositiveReal> realScalarParam(final double value) {
        return new RealScalarParam<>(value, PositiveReal.INSTANCE);
    }

    private static SharedRatesClockModel typedClock(final Tree tree,
                                                    final RealVectorParam<PositiveReal> rates,
                                                    final boolean normalize,
                                                    final RealScalarParam<PositiveReal> meanRate) {
        final SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.initByName("tree", tree, "ratesVector", rates, "normalize", normalize,
                "meanRateScalar", meanRate);
        return clock;
    }

    private static SharedRatesClockModel typedClockWithoutMeanRate(final Tree tree,
                                                                   final RealVectorParam<PositiveReal> rates,
                                                                   final boolean normalize) {
        final SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.initByName("tree", tree, "ratesVector", rates, "normalize", normalize);
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
}
