package mixture.beast.evolution.operator;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import mixture.beast.evolution.util.BranchRateIndexHelper;
import org.junit.Test;

import java.util.ArrayDeque;
import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class SubtreeRateScaleOperatorCharacterizationTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,(C:1.5,D:2.5):3.5);";

    @Test
    public void scalarRatesExpandToOneRatePerNonRootBranch() {
        final Tree tree = fixedTree();
        final RealParameter rates = realParameter("2.0");

        final SubtreeRateScaleOperator op = operator(tree, rates, 0.25);

        assertTrue(op != null);
        assertEquals(tree.getNodeCount() - 1, rates.getDimension());
    }

    @Test
    public void invalidRatesDimensionIsRejectedAtInit() {
        final Tree tree = fixedTree();
        final RealParameter rates = realParameter("1.0 2.0");

        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(tree, op);
        op.ratesInput.setValue(rates, op);

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void invalidWindowReturnsNegativeInfinityAndDoesNotMutateRates() {
        final Tree tree = fixedTree();
        final RealParameter rates = fullRates(tree, 1.0);
        final double[] before = copyValues(rates);
        final SubtreeRateScaleOperator op = operator(tree, rates, 0.0);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void negativeWindowReturnsNegativeInfinityAndDoesNotMutateRates() {
        final Tree tree = fixedTree();
        final RealParameter rates = fullRates(tree, 1.0);
        final double[] before = copyValues(rates);
        final SubtreeRateScaleOperator op = operator(tree, rates, -0.5);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void allNonPositiveRatesReturnNegativeInfinityAndDoNotMutateRates() {
        Randomizer.setSeed(11L);

        final Tree tree = fixedTree();
        final RealParameter rates = fullRates(tree, 0.0);
        final double[] before = copyValues(rates);
        final SubtreeRateScaleOperator op = operator(tree, rates, 0.5);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void validProposalScalesExactlyOneSubtreeAndReturnsCountTimesLogScaleFactor() {
        Randomizer.setSeed(101L);

        final Tree tree = fixedTree();
        final RealParameter rates = increasingRates(tree);
        final double[] before = copyValues(rates);
        final double window = 0.31;
        final SubtreeRateScaleOperator op = operator(tree, rates, window);

        final double hr = op.proposal();

        assertTrue(Double.isFinite(hr));

        final double[] after = copyValues(rates);
        final BranchRateIndexHelper.Mapping mapping = BranchRateIndexHelper.buildDeterministic(tree);
        final Set<Integer> changedIndices = changedIndices(before, after);
        assertTrue(changedIndices.size() >= 1);

        final double eps = commonLogRatio(before, after, changedIndices);
        assertTrue(Math.abs(eps) <= window + EPS);
        assertEquals(changedIndices.size() * eps, hr, EPS);

        assertChangedIndicesFormOneSubtree(tree, mapping, changedIndices);

        for (int i = 0; i < before.length; i++) {
            if (!changedIndices.contains(i)) {
                assertEquals(before[i], after[i], EPS);
            }
        }
    }

    @Test
    public void defaultWindowIsPointFiveByObservedSubtreeScaleBound() {
        Randomizer.setSeed(202L);

        final Tree tree = fixedTree();
        final RealParameter rates = increasingRates(tree);
        final double[] before = copyValues(rates);
        final SubtreeRateScaleOperator op = operatorWithDefaultWindow(tree, rates);

        final double hr = op.proposal();

        assertTrue(Double.isFinite(hr));
        final Set<Integer> changedIndices = changedIndices(before, copyValues(rates));
        assertTrue(changedIndices.size() >= 1);
        final double eps = hr / changedIndices.size();
        assertTrue(Math.abs(eps) <= 0.5 + EPS);
    }

    @Test
    public void repeatedValidProposalsKeepRatesPositiveAndReturnFiniteHastingsTerms() {
        Randomizer.setSeed(303L);

        final Tree tree = fixedTree();
        final RealParameter rates = increasingRates(tree);
        final SubtreeRateScaleOperator op = operator(tree, rates, 0.2);

        for (int i = 0; i < 100; i++) {
            final double hr = op.proposal();
            assertTrue(Double.isFinite(hr));
            for (int j = 0; j < rates.getDimension(); j++) {
                assertTrue(rates.getValue(j) > 0.0);
            }
        }
    }

    private static SubtreeRateScaleOperator operator(final Tree tree,
                                                    final RealParameter rates,
                                                    final double window) {
        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(tree, op);
        op.ratesInput.setValue(rates, op);
        op.windowInput.setValue(window, op);
        op.initAndValidate();
        return op;
    }

    private static SubtreeRateScaleOperator operatorWithDefaultWindow(final Tree tree,
                                                                      final RealParameter rates) {
        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(tree, op);
        op.ratesInput.setValue(rates, op);
        op.initAndValidate();
        return op;
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealParameter realParameter(final String value) {
        return new RealParameter(value);
    }

    private static RealParameter fullRates(final Tree tree, final double value) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tree.getNodeCount() - 1; i++) {
            if (i > 0) {
                sb.append(' ');
            }
            sb.append(value);
        }
        return new RealParameter(sb.toString());
    }

    private static RealParameter increasingRates(final Tree tree) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tree.getNodeCount() - 1; i++) {
            if (i > 0) {
                sb.append(' ');
            }
            sb.append(1.0 + i * 0.25);
        }
        return new RealParameter(sb.toString());
    }

    private static double[] copyValues(final RealParameter parameter) {
        final double[] values = new double[parameter.getDimension()];
        for (int i = 0; i < values.length; i++) {
            values[i] = parameter.getValue(i);
        }
        return values;
    }

    private static Set<Integer> changedIndices(final double[] before, final double[] after) {
        final Set<Integer> indices = new HashSet<>();
        for (int i = 0; i < before.length; i++) {
            if (Math.abs(before[i] - after[i]) > EPS) {
                indices.add(i);
            }
        }
        return indices;
    }

    private static double commonLogRatio(final double[] before,
                                         final double[] after,
                                         final Set<Integer> changedIndices) {
        Double eps = null;
        for (final int idx : changedIndices) {
            final double thisEps = Math.log(after[idx] / before[idx]);
            if (eps == null) {
                eps = thisEps;
            } else {
                assertEquals(eps, thisEps, EPS);
            }
        }
        if (eps == null) {
            throw new IllegalArgumentException("No changed indices");
        }
        return eps;
    }

    private static void assertChangedIndicesFormOneSubtree(final Tree tree,
                                                           final BranchRateIndexHelper.Mapping mapping,
                                                           final Set<Integer> changedIndices) {
        final Set<Node> changedNodes = new HashSet<>();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            if (!node.isRoot() && changedIndices.contains(mapping.idxForNode(node))) {
                changedNodes.add(node);
            }
        }

        Node top = null;
        for (final Node node : changedNodes) {
            final Node parent = node.getParent();
            if (parent == null || parent.isRoot() || !changedNodes.contains(parent)) {
                if (top != null) {
                    throw new AssertionError("Changed nodes have more than one top subtree node");
                }
                top = node;
            }
        }

        if (top == null) {
            throw new AssertionError("No top subtree node identified");
        }

        final Set<Node> expected = nonRootDescendantsIncludingSelf(top);
        assertEquals(expected, changedNodes);
    }

    private static Set<Node> nonRootDescendantsIncludingSelf(final Node top) {
        final Set<Node> nodes = new HashSet<>();
        final ArrayDeque<Node> stack = new ArrayDeque<>();
        stack.push(top);
        while (!stack.isEmpty()) {
            final Node node = stack.pop();
            if (!node.isRoot()) {
                nodes.add(node);
            }
            for (int i = 0; i < node.getChildCount(); i++) {
                stack.push(node.getChild(i));
            }
        }
        return nodes;
    }

    private static void assertArrayEquals(final double[] expected, final double[] observed, final double eps) {
        assertEquals(expected.length, observed.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], observed[i], eps);
        }
    }
}
