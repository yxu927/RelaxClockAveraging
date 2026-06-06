package mixture.beast.evolution.operator;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import mixture.beast.evolution.util.BranchRateIndexHelper;
import org.junit.Test;

import java.util.ArrayDeque;
import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class SubtreeRateScaleOperatorTypedInputTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,(C:1.5,D:2.5):3.5);";

    @Test
    public void typedScalarRatesExpandToOneRatePerNonRootBranch() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = positiveRates(2.0);

        final SubtreeRateScaleOperator op = typedOperator(tree, rates, 0.25);

        assertTrue(op != null);
        assertEquals(tree.getNodeCount() - 1, rates.size());
        for (int i = 0; i < rates.size(); i++) {
            assertEquals(2.0, rates.get(i), EPS);
        }
    }

    @Test
    public void typedInvalidRatesDimensionIsRejectedAtInit() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = positiveRates(1.0, 2.0);

        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(tree, op);
        op.ratesVectorInput.setValue(rates, op);

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void typedInvalidWindowReturnsNegativeInfinityAndDoesNotMutateRates() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = fullTypedRates(tree, 1.0);
        final double[] before = copyValues(rates);
        final SubtreeRateScaleOperator op = typedOperator(tree, rates, 0.0);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void typedNegativeWindowReturnsNegativeInfinityAndDoesNotMutateRates() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = fullTypedRates(tree, 1.0);
        final double[] before = copyValues(rates);
        final SubtreeRateScaleOperator op = typedOperator(tree, rates, -0.5);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void typedAllNonPositiveRatesReturnNegativeInfinityAndDoNotMutateRates() {
        Randomizer.setSeed(11L);

        final Tree tree = fixedTree();
        final RealVectorParam<Real> rates = fullRealTypedRates(tree, 0.0);
        final double[] before = copyValues(rates);
        final SubtreeRateScaleOperator op = typedOperator(tree, rates, 0.5);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void typedValidProposalScalesExactlyOneSubtreeAndReturnsCountTimesLogScaleFactor() {
        Randomizer.setSeed(101L);

        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = increasingTypedRates(tree);
        final double[] before = copyValues(rates);
        final double window = 0.31;
        final SubtreeRateScaleOperator op = typedOperator(tree, rates, window);

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
    public void typedDefaultWindowIsPointFiveByObservedSubtreeScaleBound() {
        Randomizer.setSeed(202L);

        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = increasingTypedRates(tree);
        final double[] before = copyValues(rates);
        final SubtreeRateScaleOperator op = typedOperatorWithDefaultWindow(tree, rates);

        final double hr = op.proposal();

        assertTrue(Double.isFinite(hr));
        final Set<Integer> changedIndices = changedIndices(before, copyValues(rates));
        assertTrue(changedIndices.size() >= 1);
        final double eps = hr / changedIndices.size();
        assertTrue(Math.abs(eps) <= 0.5 + EPS);
    }

    @Test
    public void typedRepeatedValidProposalsKeepRatesPositiveAndReturnFiniteHastingsTerms() {
        Randomizer.setSeed(303L);

        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = increasingTypedRates(tree);
        final SubtreeRateScaleOperator op = typedOperator(tree, rates, 0.2);

        for (int i = 0; i < 100; i++) {
            final double hr = op.proposal();
            assertTrue(Double.isFinite(hr));
            for (int j = 0; j < rates.size(); j++) {
                assertTrue(rates.get(j) > 0.0);
            }
        }
    }

    @Test
    public void legacyAndTypedRatesCannotBothBeSpecified() {
        final Tree tree = fixedTree();
        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(tree, op);
        op.ratesInput.setValue(fullLegacyRates(tree, 1.0), op);
        op.ratesVectorInput.setValue(fullTypedRates(tree, 1.0), op);

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void missingRatesIsRejected() {
        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(fixedTree(), op);

        assertThrows(IllegalArgumentException.class, op::initAndValidate);
    }

    @Test
    public void typedAndLegacyOperatorsEquivalentForSameRandomDraw() {
        final long seed = 808L;
        final double window = 0.31;

        final Tree legacyTree = fixedTree();
        final RealParameter legacyRates = increasingLegacyRates(legacyTree);
        final SubtreeRateScaleOperator legacyOperator = legacyOperator(legacyTree, legacyRates, window);
        final double[] legacyBefore = copyValues(legacyRates);

        Randomizer.setSeed(seed);
        final double legacyHr = legacyOperator.proposal();

        final Tree typedTree = fixedTree();
        final RealVectorParam<PositiveReal> typedRates = increasingTypedRates(typedTree);
        final SubtreeRateScaleOperator typedOperator = typedOperator(typedTree, typedRates, window);
        final double[] typedBefore = copyValues(typedRates);

        Randomizer.setSeed(seed);
        final double typedHr = typedOperator.proposal();

        assertArrayEquals(legacyBefore, typedBefore, EPS);
        assertEquals(legacyHr, typedHr, EPS);
        assertArrayEquals(copyValues(legacyRates), copyValues(typedRates), EPS);
        assertEquals(changedIndices(legacyBefore, copyValues(legacyRates)),
                changedIndices(typedBefore, copyValues(typedRates)));
    }

    private static SubtreeRateScaleOperator typedOperator(final Tree tree,
                                                         final RealVectorParam<?> rates,
                                                         final double window) {
        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(tree, op);
        op.ratesVectorInput.setValue(rates, op);
        op.windowInput.setValue(window, op);
        op.initAndValidate();
        return op;
    }

    private static SubtreeRateScaleOperator typedOperatorWithDefaultWindow(final Tree tree,
                                                                          final RealVectorParam<?> rates) {
        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(tree, op);
        op.ratesVectorInput.setValue(rates, op);
        op.initAndValidate();
        return op;
    }

    private static SubtreeRateScaleOperator legacyOperator(final Tree tree,
                                                          final RealParameter rates,
                                                          final double window) {
        final SubtreeRateScaleOperator op = new SubtreeRateScaleOperator();
        op.treeInput.setValue(tree, op);
        op.ratesInput.setValue(rates, op);
        op.windowInput.setValue(window, op);
        op.initAndValidate();
        return op;
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealVectorParam<PositiveReal> positiveRates(final double... values) {
        return new RealVectorParam<>(values, PositiveReal.INSTANCE);
    }

    private static RealVectorParam<Real> realRates(final double... values) {
        return new RealVectorParam<>(values, Real.INSTANCE);
    }

    private static RealVectorParam<PositiveReal> fullTypedRates(final Tree tree, final double value) {
        final double[] values = fullValues(tree, value);
        return positiveRates(values);
    }

    private static RealVectorParam<Real> fullRealTypedRates(final Tree tree, final double value) {
        final double[] values = fullValues(tree, value);
        return realRates(values);
    }

    private static RealVectorParam<PositiveReal> increasingTypedRates(final Tree tree) {
        return positiveRates(increasingValues(tree));
    }

    private static RealParameter fullLegacyRates(final Tree tree, final double value) {
        return new RealParameter(join(fullValues(tree, value)));
    }

    private static RealParameter increasingLegacyRates(final Tree tree) {
        return new RealParameter(join(increasingValues(tree)));
    }

    private static double[] fullValues(final Tree tree, final double value) {
        final double[] values = new double[tree.getNodeCount() - 1];
        for (int i = 0; i < values.length; i++) {
            values[i] = value;
        }
        return values;
    }

    private static double[] increasingValues(final Tree tree) {
        final double[] values = new double[tree.getNodeCount() - 1];
        for (int i = 0; i < values.length; i++) {
            values[i] = 1.0 + i * 0.25;
        }
        return values;
    }

    private static String join(final double[] values) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < values.length; i++) {
            if (i > 0) {
                sb.append(' ');
            }
            sb.append(values[i]);
        }
        return sb.toString();
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
