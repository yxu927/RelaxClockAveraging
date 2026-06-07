package mixture.beast.evolution.operator;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import mixture.beast.evolution.util.BranchRateIndexHelper;
import org.junit.Test;

import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class RemainingOperatorsCharacterizationTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,(C:1.5,D:2.5):3.5);";

    @Test
    public void ucldNonCenteredPreservesLatentZAndReturnsJacobianHastings() {
        Randomizer.setSeed(101L);

        final RealParameter rates = realParameter(0.8, 1.1, 1.4, 0.9);
        final IntegerParameter indicator = integerParameter(0);
        final RealParameter stdev = realParameter(0.5);
        final UCLDStdevNonCenteredOperator op = legacyUcld(rates, indicator, stdev, 0.2);

        final double oldS = stdev.getValue(0);
        final double[] before = copyValues(rates);
        final double[] zBefore = ucldZ(before, oldS);

        final double hr = op.proposal();

        final double newS = stdev.getValue(0);
        final double[] after = copyValues(rates);
        final double eps = Math.log(newS / oldS);
        final double expected = sumLogRatio(before, after) + (before.length + 1.0) * eps;

        assertTrue(Double.isFinite(hr));
        assertEquals(expected, hr, EPS);
        assertArrayEquals(zBefore, ucldZ(after, newS), EPS);
        assertAllPositive(after);
    }

    @Test
    public void ucldNonCenteredRejectsWhenIndicatorIsACAndDoesNotMutate() {
        final RealParameter rates = realParameter(0.8, 1.1, 1.4, 0.9);
        final IntegerParameter indicator = integerParameter(1);
        final RealParameter stdev = realParameter(0.5);
        final UCLDStdevNonCenteredOperator op = legacyUcld(rates, indicator, stdev, 0.2);

        final double[] before = copyValues(rates);
        final double oldS = stdev.getValue(0);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
        assertEquals(oldS, stdev.getValue(0), EPS);
    }

    @Test
    public void acSigma2NonCenteredPreservesACIncrementsAndReturnsJacobianHastings() {
        Randomizer.setSeed(202L);

        final Tree tree = fixedTree();
        final RealParameter rates = increasingLegacyRates(tree);
        final IntegerParameter indicator = integerParameter(1);
        final RealParameter sigma2 = realParameter(0.25);
        final RealParameter rootLogRate = realParameter(0.0);
        final ACSigma2NonCenteredOperator op = legacyAcSigma(tree, rates, indicator, sigma2, rootLogRate, 0.15);

        final double oldSigma2 = sigma2.getValue(0);
        final double[] before = copyValues(rates);
        final double[] uBefore = acU(tree, before, oldSigma2, rootLogRate.getValue(0), 1e-12);

        final double hr = op.proposal();

        final double newSigma2 = sigma2.getValue(0);
        final double[] after = copyValues(rates);
        final double eps = Math.log(newSigma2 / oldSigma2);
        final double expected = sumLogRatio(before, after) + (0.5 * before.length + 1.0) * eps;

        assertTrue(Double.isFinite(hr));
        assertEquals(expected, hr, EPS);
        assertArrayEquals(uBefore, acU(tree, after, newSigma2, rootLogRate.getValue(0), 1e-12), EPS);
        assertAllPositive(after);
    }

    @Test
    public void acSigma2NonCenteredRejectsWhenIndicatorIsUCAndDoesNotMutate() {
        final Tree tree = fixedTree();
        final RealParameter rates = increasingLegacyRates(tree);
        final IntegerParameter indicator = integerParameter(0);
        final RealParameter sigma2 = realParameter(0.25);
        final ACSigma2NonCenteredOperator op = legacyAcSigma(tree, rates, indicator, sigma2, null, 0.15);

        final double[] before = copyValues(rates);
        final double oldSigma2 = sigma2.getValue(0);

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
        assertEquals(oldSigma2, sigma2.getValue(0), EPS);
    }

    @Test
    public void acSubtreeUIncrementReturnsSubtreeLogJacobianAndKeepsRatesPositive() {
        Randomizer.setSeed(303L);

        final Tree tree = fixedTree();
        final RealParameter rates = increasingLegacyRates(tree);
        final IntegerParameter indicator = integerParameter(1);
        final RealParameter sigma2 = realParameter(0.25);
        final ACSubtreeUIncrementOperator op = legacyAcSubtree(tree, rates, indicator, sigma2, null, 0.25);

        final double[] before = copyValues(rates);

        final double hr = op.proposal();

        final double[] after = copyValues(rates);
        final Set<Integer> changed = changedIndices(before, after);

        assertTrue(Double.isFinite(hr));
        assertTrue(changed.size() >= 1);
        assertEquals(sumLogRatio(before, after, changed), hr, EPS);
        assertChangedIndicesFormOneSubtree(tree, BranchRateIndexHelper.buildDeterministic(tree), changed);
        assertAllPositive(after);
        assertEquals(0.25, sigma2.getValue(0), EPS);
    }

    @Test
    public void ucacBridgeRoundTripRestoresRatesIndicatorAndHastingsTermsCancel() {
        final Tree tree = fixedTree();
        final RealParameter rates = increasingLegacyRates(tree);
        final IntegerParameter indicator = integerParameter(0);
        final RealParameter stdev = realParameter(0.5);
        final RealParameter sigma2 = realParameter(0.25);
        final RealParameter rootLogRate = realParameter(0.0);
        final UCACSwitchBridgeOperator op = legacyBridge(tree, rates, indicator, stdev, sigma2, rootLogRate);

        final double[] before = copyValues(rates);

        final double ucToAc = op.proposal();
        assertEquals(1, indicator.getValue(0).intValue());

        final double acToUc = op.proposal();
        assertEquals(0, indicator.getValue(0).intValue());

        assertTrue(Double.isFinite(ucToAc));
        assertTrue(Double.isFinite(acToUc));
        assertEquals(0.0, ucToAc + acToUc, EPS);
        assertArrayEquals(before, copyValues(rates), EPS);
    }

    @Test
    public void alphaAnnealingMovesLinearlyToEndAndThenStops() {
        final RealParameter alpha = realParameter(1.0);
        final AlphaAnnealingOperator op = legacyAlpha(alpha, 1.0, 0.0, 4);

        assertEquals(0.0, op.proposal(), EPS);
        assertEquals(0.75, alpha.getValue(0), EPS);
        assertEquals(0.0, op.proposal(), EPS);
        assertEquals(0.5, alpha.getValue(0), EPS);
        assertEquals(0.0, op.proposal(), EPS);
        assertEquals(0.25, alpha.getValue(0), EPS);
        assertEquals(0.0, op.proposal(), EPS);
        assertEquals(0.0, alpha.getValue(0), EPS);
        assertEquals(0.0, op.proposal(), EPS);
        assertEquals(0.0, alpha.getValue(0), EPS);
    }

    private static UCLDStdevNonCenteredOperator legacyUcld(final RealParameter rates,
                                                           final IntegerParameter indicator,
                                                           final RealParameter stdev,
                                                           final double window) {
        final UCLDStdevNonCenteredOperator op = new UCLDStdevNonCenteredOperator();
        op.ratesInput.setValue(rates, op);
        op.indicatorInput.setValue(indicator, op);
        op.ucldStdevInput.setValue(stdev, op);
        op.windowInput.setValue(window, op);
        op.initAndValidate();
        return op;
    }

    private static ACSigma2NonCenteredOperator legacyAcSigma(final Tree tree,
                                                             final RealParameter rates,
                                                             final IntegerParameter indicator,
                                                             final RealParameter sigma2,
                                                             final RealParameter rootLogRate,
                                                             final double window) {
        final ACSigma2NonCenteredOperator op = new ACSigma2NonCenteredOperator();
        op.treeInput.setValue(tree, op);
        op.ratesInput.setValue(rates, op);
        op.indicatorInput.setValue(indicator, op);
        op.sigma2Input.setValue(sigma2, op);
        if (rootLogRate != null) {
            op.rootLogRateInput.setValue(rootLogRate, op);
        }
        op.windowInput.setValue(window, op);
        op.initAndValidate();
        return op;
    }

    private static ACSubtreeUIncrementOperator legacyAcSubtree(final Tree tree,
                                                               final RealParameter rates,
                                                               final IntegerParameter indicator,
                                                               final RealParameter sigma2,
                                                               final RealParameter rootLogRate,
                                                               final double delta) {
        final ACSubtreeUIncrementOperator op = new ACSubtreeUIncrementOperator();
        op.treeInput.setValue(tree, op);
        op.ratesInput.setValue(rates, op);
        op.indicatorInput.setValue(indicator, op);
        op.sigma2Input.setValue(sigma2, op);
        if (rootLogRate != null) {
            op.rootLogRateInput.setValue(rootLogRate, op);
        }
        op.deltaInput.setValue(delta, op);
        op.initAndValidate();
        return op;
    }

    private static UCACSwitchBridgeOperator legacyBridge(final Tree tree,
                                                         final RealParameter rates,
                                                         final IntegerParameter indicator,
                                                         final RealParameter stdev,
                                                         final RealParameter sigma2,
                                                         final RealParameter rootLogRate) {
        final UCACSwitchBridgeOperator op = new UCACSwitchBridgeOperator();
        op.treeInput.setValue(tree, op);
        op.ratesInput.setValue(rates, op);
        op.indicatorInput.setValue(indicator, op);
        op.ucldStdevInput.setValue(stdev, op);
        op.sigma2Input.setValue(sigma2, op);
        if (rootLogRate != null) {
            op.rootLogRateInput.setValue(rootLogRate, op);
        }
        op.initAndValidate();
        return op;
    }

    private static AlphaAnnealingOperator legacyAlpha(final RealParameter alpha,
                                                      final double start,
                                                      final double end,
                                                      final int steps) {
        final AlphaAnnealingOperator op = new AlphaAnnealingOperator();
        op.alphaInput.setValue(alpha, op);
        op.alphaStartInput.setValue(start, op);
        op.alphaEndInput.setValue(end, op);
        op.alphaStepsInput.setValue(steps, op);
        op.initAndValidate();
        return op;
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealParameter realParameter(final double... values) {
        return new RealParameter(join(values));
    }

    private static IntegerParameter integerParameter(final int value) {
        return new IntegerParameter(Integer.toString(value));
    }

    private static RealParameter increasingLegacyRates(final Tree tree) {
        final double[] values = new double[tree.getNodeCount() - 1];
        for (int i = 0; i < values.length; i++) {
            values[i] = 0.8 + i * 0.15;
        }
        return realParameter(values);
    }

    private static double[] ucldZ(final double[] rates, final double stdev) {
        final double[] z = new double[rates.length];
        for (int i = 0; i < rates.length; i++) {
            z[i] = (Math.log(rates[i]) + 0.5 * stdev * stdev) / stdev;
        }
        return z;
    }

    private static double[] acU(final Tree tree,
                                final double[] rates,
                                final double sigma2,
                                final double rootLogRate,
                                final double minDt) {
        final double[] u = new double[rates.length];
        fillAcU(tree.getRoot(), rootLogRate, rates, u, sigma2, minDt,
                BranchRateIndexHelper.buildDeterministic(tree));
        return u;
    }

    private static void fillAcU(final Node parent,
                                final double logParent,
                                final double[] rates,
                                final double[] u,
                                final double sigma2,
                                final double minDt,
                                final BranchRateIndexHelper.Mapping mapping) {
        for (int i = 0; i < parent.getChildCount(); i++) {
            final Node child = parent.getChild(i);
            final int idx = mapping.idxForNode(child);
            final double var = sigma2 * child.getLength();
            if (!(child.getLength() > minDt) || !(var > 0.0)) {
                throw new AssertionError("Invalid AC branch in test tree");
            }
            final double x = Math.log(rates[idx]);
            final double mean = logParent - 0.5 * var;
            u[idx] = (x - mean) / Math.sqrt(var);
            fillAcU(child, x, rates, u, sigma2, minDt, mapping);
        }
    }

    private static double[] copyValues(final RealParameter parameter) {
        final double[] values = new double[parameter.getDimension()];
        for (int i = 0; i < values.length; i++) {
            values[i] = parameter.getValue(i);
        }
        return values;
    }

    private static double sumLogRatio(final double[] before, final double[] after) {
        double sum = 0.0;
        for (int i = 0; i < before.length; i++) {
            sum += Math.log(after[i] / before[i]);
        }
        return sum;
    }

    private static double sumLogRatio(final double[] before,
                                      final double[] after,
                                      final Set<Integer> indices) {
        double sum = 0.0;
        for (final int i : indices) {
            sum += Math.log(after[i] / before[i]);
        }
        return sum;
    }

    private static Set<Integer> changedIndices(final double[] before, final double[] after) {
        final Set<Integer> changed = new HashSet<>();
        for (int i = 0; i < before.length; i++) {
            if (Math.abs(before[i] - after[i]) > EPS) {
                changed.add(i);
            }
        }
        return changed;
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

        final Set<Node> expected = new HashSet<>();
        collectNonRootDescendants(top, expected);
        assertEquals(expected, changedNodes);
    }

    private static void collectNonRootDescendants(final Node node, final Set<Node> out) {
        if (!node.isRoot()) {
            out.add(node);
        }
        for (int i = 0; i < node.getChildCount(); i++) {
            collectNonRootDescendants(node.getChild(i), out);
        }
    }

    private static void assertAllPositive(final double[] values) {
        for (final double value : values) {
            assertTrue(value > 0.0);
        }
    }

    private static void assertArrayEquals(final double[] expected, final double[] observed, final double eps) {
        assertEquals(expected.length, observed.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], observed[i], eps);
        }
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
}
