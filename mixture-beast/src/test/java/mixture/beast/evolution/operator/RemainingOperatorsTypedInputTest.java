package mixture.beast.evolution.operator;

import beast.base.evolution.tree.Node;
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
import mixture.beast.evolution.util.BranchRateIndexHelper;
import org.junit.Test;

import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class RemainingOperatorsTypedInputTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,(C:1.5,D:2.5):3.5);";

    @Test
    public void typedUcldNonCenteredPreservesLatentZAndReturnsJacobianHastings() {
        final RealVectorParam<PositiveReal> typedRates = positiveRates(0.8, 1.1, 1.4, 0.9);
        final IntScalarParam<Int> typedIndicator = intScalar(0);
        final RealScalarParam<PositiveReal> typedStdev = positiveScalar(0.5);
        final UCLDStdevNonCenteredOperator typed = typedUcld(typedRates, typedIndicator, typedStdev, 0.2);

        final double oldS = typedStdev.get();
        final double[] before = copyValues(typedRates);
        final double[] zBefore = ucldZ(before, oldS);

        Randomizer.setSeed(101L);
        final double typedHr = typed.proposal();

        final double newS = typedStdev.get();
        final double[] after = copyValues(typedRates);
        final double eps = Math.log(newS / oldS);
        final double expected = sumLogRatio(before, after) + (before.length + 1.0) * eps;

        assertTrue(Double.isFinite(typedHr));
        assertEquals(expected, typedHr, EPS);
        assertArrayEquals(zBefore, ucldZ(after, newS), EPS);
        assertAllPositive(after);
        assertEquals(0, typedIndicator.get());
    }

    @Test
    public void typedAcSigma2NonCenteredPreservesACIncrementsAndReturnsJacobianHastings() {
        final Tree typedTree = fixedTree();
        final RealVectorParam<PositiveReal> typedRates = increasingTypedRates(typedTree);
        final IntScalarParam<Int> typedIndicator = intScalar(1);
        final RealScalarParam<PositiveReal> typedSigma2 = positiveScalar(0.25);
        final RealScalarParam<Real> typedRootLogRate = realScalar(0.0);
        final ACSigma2NonCenteredOperator typed = typedAcSigma(typedTree, typedRates, typedIndicator,
                typedSigma2, typedRootLogRate, 0.15);

        final double oldSigma2 = typedSigma2.get();
        final double[] before = copyValues(typedRates);
        final double[] uBefore = acU(typedTree, before, oldSigma2, typedRootLogRate.get(), 1e-12);

        Randomizer.setSeed(202L);
        final double typedHr = typed.proposal();

        final double newSigma2 = typedSigma2.get();
        final double[] after = copyValues(typedRates);
        final double eps = Math.log(newSigma2 / oldSigma2);
        final double expected = sumLogRatio(before, after) + (0.5 * before.length + 1.0) * eps;

        assertTrue(Double.isFinite(typedHr));
        assertEquals(expected, typedHr, EPS);
        assertArrayEquals(uBefore, acU(typedTree, after, newSigma2, typedRootLogRate.get(), 1e-12), EPS);
        assertAllPositive(after);
        assertEquals(1, typedIndicator.get());
        assertEquals(0.0, typedRootLogRate.get(), EPS);
    }

    @Test
    public void typedAcSubtreeUIncrementReturnsSubtreeLogJacobianAndKeepsRatesPositive() {
        final Tree typedTree = fixedTree();
        final RealVectorParam<PositiveReal> typedRates = increasingTypedRates(typedTree);
        final IntScalarParam<Int> typedIndicator = intScalar(1);
        final RealScalarParam<PositiveReal> typedSigma2 = positiveScalar(0.25);
        final ACSubtreeUIncrementOperator typed = typedAcSubtree(typedTree, typedRates, typedIndicator,
                typedSigma2, null, 0.25);

        final double[] before = copyValues(typedRates);

        Randomizer.setSeed(303L);
        final double typedHr = typed.proposal();

        final double[] after = copyValues(typedRates);
        final Set<Integer> changed = changedIndices(before, after);

        assertTrue(Double.isFinite(typedHr));
        assertTrue(changed.size() >= 1);
        assertEquals(sumLogRatio(before, after, changed), typedHr, EPS);
        assertChangedIndicesFormOneSubtree(typedTree, BranchRateIndexHelper.buildDeterministic(typedTree), changed);
        assertAllPositive(after);
        assertEquals(1, typedIndicator.get());
        assertEquals(0.25, typedSigma2.get(), EPS);
    }

    @Test
    public void ucacBridgeLegacyAndTypedAreEquivalentForBothDirections() {
        final Tree legacyTree = fixedTree();
        final RealParameter legacyRates = increasingLegacyRates(legacyTree);
        final IntegerParameter legacyIndicator = integerParameter(0);
        final RealParameter legacyStdev = realParameter(0.5);
        final RealParameter legacySigma2 = realParameter(0.25);
        final RealParameter legacyRootLogRate = realParameter(0.0);
        final UCACSwitchBridgeOperator legacy = legacyBridge(legacyTree, legacyRates, legacyIndicator,
                legacyStdev, legacySigma2, legacyRootLogRate);
        final double[] initialRates = copyValues(legacyRates);

        final double legacyForwardHr = legacy.proposal();
        final double[] legacyForwardRates = copyValues(legacyRates);
        final int legacyForwardIndicator = legacyIndicator.getValue(0);
        final double legacyReverseHr = legacy.proposal();

        final Tree typedTree = legacyTree;
        final RealVectorParam<PositiveReal> typedRates = positiveRates(initialRates);
        final IntScalarParam<Int> typedIndicator = intScalar(0);
        final RealScalarParam<PositiveReal> typedStdev = positiveScalar(0.5);
        final RealScalarParam<PositiveReal> typedSigma2 = positiveScalar(0.25);
        final RealScalarParam<Real> typedRootLogRate = realScalar(0.0);
        final UCACSwitchBridgeOperator typed = typedBridge(typedTree, typedRates, typedIndicator,
                typedStdev, typedSigma2, typedRootLogRate);

        final double typedForwardHr = typed.proposal();
        final double[] typedForwardRates = copyValues(typedRates);
        final int typedForwardIndicator = typedIndicator.get();
        final double typedReverseHr = typed.proposal();

        assertEquals(legacyForwardHr, typedForwardHr, EPS);
        assertArrayEquals(legacyForwardRates, typedForwardRates, EPS);
        assertEquals(legacyForwardIndicator, typedForwardIndicator);
        assertEquals(legacyReverseHr, typedReverseHr, EPS);
        assertArrayEquals(copyValues(legacyRates), copyValues(typedRates), EPS);
        assertEquals(legacyIndicator.getValue(0).intValue(), typedIndicator.get());
    }

    @Test
    public void alphaAnnealingLegacyAndTypedAreEquivalent() {
        final RealParameter legacyAlpha = realParameter(1.0);
        final AlphaAnnealingOperator legacy = legacyAlpha(legacyAlpha, 1.0, 0.0, 4);

        final RealScalarParam<Real> typedAlpha = realScalar(1.0);
        final AlphaAnnealingOperator typed = typedAlpha(typedAlpha, 1.0, 0.0, 4);

        for (int i = 0; i < 5; i++) {
            assertEquals(legacy.proposal(), typed.proposal(), EPS);
            assertEquals(legacyAlpha.getValue(0), typedAlpha.get(), EPS);
        }
    }

    @Test
    public void typedScalarRatesExpandForTreeOperators() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = positiveRates(1.5);

        final ACSigma2NonCenteredOperator op = typedAcSigma(tree, rates, intScalar(1),
                positiveScalar(0.25), null, 0.15);

        assertEquals(tree.getNodeCount() - 1, rates.size());
        for (int i = 0; i < rates.size(); i++) {
            assertEquals(1.5, rates.get(i), EPS);
        }
        assertEquals(op, op);
    }

    @Test
    public void typedUcldRejectsWhenIndicatorIsACAndDoesNotMutate() {
        final RealVectorParam<PositiveReal> rates = positiveRates(0.8, 1.1, 1.4, 0.9);
        final IntScalarParam<Int> indicator = intScalar(1);
        final RealScalarParam<PositiveReal> stdev = positiveScalar(0.5);
        final UCLDStdevNonCenteredOperator op = typedUcld(rates, indicator, stdev, 0.2);

        final double[] before = copyValues(rates);
        final double stdevBefore = stdev.get();

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
        assertEquals(1, indicator.get());
        assertEquals(stdevBefore, stdev.get(), EPS);
    }

    @Test
    public void typedAcSigmaRejectsWhenIndicatorIsUCAndDoesNotMutate() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = increasingTypedRates(tree);
        final IntScalarParam<Int> indicator = intScalar(0);
        final RealScalarParam<PositiveReal> sigma2 = positiveScalar(0.25);
        final RealScalarParam<Real> rootLogRate = realScalar(0.0);
        final ACSigma2NonCenteredOperator op = typedAcSigma(tree, rates, indicator, sigma2, rootLogRate, 0.15);

        final double[] before = copyValues(rates);
        final double sigma2Before = sigma2.get();
        final double rootBefore = rootLogRate.get();

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
        assertEquals(0, indicator.get());
        assertEquals(sigma2Before, sigma2.get(), EPS);
        assertEquals(rootBefore, rootLogRate.get(), EPS);
    }

    @Test
    public void typedAcSubtreeRejectsWhenIndicatorIsUCAndDoesNotMutate() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> rates = increasingTypedRates(tree);
        final IntScalarParam<Int> indicator = intScalar(0);
        final RealScalarParam<PositiveReal> sigma2 = positiveScalar(0.25);
        final ACSubtreeUIncrementOperator op = typedAcSubtree(tree, rates, indicator, sigma2, null, 0.25);

        final double[] before = copyValues(rates);
        final double sigma2Before = sigma2.get();

        final double hr = op.proposal();

        assertEquals(Double.NEGATIVE_INFINITY, hr, 0.0);
        assertArrayEquals(before, copyValues(rates), EPS);
        assertEquals(0, indicator.get());
        assertEquals(sigma2Before, sigma2.get(), EPS);
    }

    @Test
    public void remainingOperatorsRejectBothLegacyAndTypedInputs() {
        assertThrows(IllegalArgumentException.class, () -> {
            final UCLDStdevNonCenteredOperator op = new UCLDStdevNonCenteredOperator();
            op.ratesInput.setValue(realParameter(1.0), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorScalarInput.setValue(intScalar(0), op);
            op.ucldStdevScalarInput.setValue(positiveScalar(0.5), op);
            op.initAndValidate();
        });

        assertThrows(IllegalArgumentException.class, () -> {
            final ACSigma2NonCenteredOperator op = new ACSigma2NonCenteredOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorInput.setValue(integerParameter(1), op);
            op.indicatorScalarInput.setValue(intScalar(1), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });

        assertThrows(IllegalArgumentException.class, () -> {
            final ACSubtreeUIncrementOperator op = new ACSubtreeUIncrementOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorScalarInput.setValue(intScalar(1), op);
            op.sigma2Input.setValue(realParameter(0.25), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });

        assertThrows(IllegalArgumentException.class, () -> {
            final UCACSwitchBridgeOperator op = new UCACSwitchBridgeOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorScalarInput.setValue(intScalar(0), op);
            op.ucldStdevInput.setValue(realParameter(0.5), op);
            op.ucldStdevScalarInput.setValue(positiveScalar(0.5), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });

        assertThrows(IllegalArgumentException.class, () -> {
            final AlphaAnnealingOperator op = new AlphaAnnealingOperator();
            op.alphaInput.setValue(realParameter(1.0), op);
            op.alphaScalarInput.setValue(realScalar(1.0), op);
            op.initAndValidate();
        });
    }

    @Test
    public void remainingOperatorsRejectMissingMandatoryInputs() {
        assertThrows(IllegalArgumentException.class, () -> {
            final UCLDStdevNonCenteredOperator op = new UCLDStdevNonCenteredOperator();
            op.indicatorScalarInput.setValue(intScalar(0), op);
            op.ucldStdevScalarInput.setValue(positiveScalar(0.5), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final UCLDStdevNonCenteredOperator op = new UCLDStdevNonCenteredOperator();
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.ucldStdevScalarInput.setValue(positiveScalar(0.5), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final UCLDStdevNonCenteredOperator op = new UCLDStdevNonCenteredOperator();
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorScalarInput.setValue(intScalar(0), op);
            op.initAndValidate();
        });

        assertThrows(IllegalArgumentException.class, () -> {
            final ACSigma2NonCenteredOperator op = new ACSigma2NonCenteredOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.indicatorScalarInput.setValue(intScalar(1), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final ACSigma2NonCenteredOperator op = new ACSigma2NonCenteredOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final ACSigma2NonCenteredOperator op = new ACSigma2NonCenteredOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorScalarInput.setValue(intScalar(1), op);
            op.initAndValidate();
        });

        assertThrows(IllegalArgumentException.class, () -> {
            final ACSubtreeUIncrementOperator op = new ACSubtreeUIncrementOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.indicatorScalarInput.setValue(intScalar(1), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final ACSubtreeUIncrementOperator op = new ACSubtreeUIncrementOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final ACSubtreeUIncrementOperator op = new ACSubtreeUIncrementOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorScalarInput.setValue(intScalar(1), op);
            op.initAndValidate();
        });

        assertThrows(IllegalArgumentException.class, () -> {
            final UCACSwitchBridgeOperator op = new UCACSwitchBridgeOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.indicatorScalarInput.setValue(intScalar(0), op);
            op.ucldStdevScalarInput.setValue(positiveScalar(0.5), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final UCACSwitchBridgeOperator op = new UCACSwitchBridgeOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.ucldStdevScalarInput.setValue(positiveScalar(0.5), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final UCACSwitchBridgeOperator op = new UCACSwitchBridgeOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorScalarInput.setValue(intScalar(0), op);
            op.sigma2ScalarInput.setValue(positiveScalar(0.25), op);
            op.initAndValidate();
        });
        assertThrows(IllegalArgumentException.class, () -> {
            final UCACSwitchBridgeOperator op = new UCACSwitchBridgeOperator();
            op.treeInput.setValue(fixedTree(), op);
            op.ratesVectorInput.setValue(positiveRates(1.0), op);
            op.indicatorScalarInput.setValue(intScalar(0), op);
            op.ucldStdevScalarInput.setValue(positiveScalar(0.5), op);
            op.initAndValidate();
        });

        assertThrows(IllegalArgumentException.class, () -> {
            final AlphaAnnealingOperator op = new AlphaAnnealingOperator();
            op.initAndValidate();
        });
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

    private static UCLDStdevNonCenteredOperator typedUcld(final RealVectorParam<?> rates,
                                                          final IntScalarParam<?> indicator,
                                                          final RealScalarParam<?> stdev,
                                                          final double window) {
        final UCLDStdevNonCenteredOperator op = new UCLDStdevNonCenteredOperator();
        op.ratesVectorInput.setValue(rates, op);
        op.indicatorScalarInput.setValue(indicator, op);
        op.ucldStdevScalarInput.setValue(stdev, op);
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

    private static ACSigma2NonCenteredOperator typedAcSigma(final Tree tree,
                                                            final RealVectorParam<?> rates,
                                                            final IntScalarParam<?> indicator,
                                                            final RealScalarParam<?> sigma2,
                                                            final RealScalarParam<?> rootLogRate,
                                                            final double window) {
        final ACSigma2NonCenteredOperator op = new ACSigma2NonCenteredOperator();
        op.treeInput.setValue(tree, op);
        op.ratesVectorInput.setValue(rates, op);
        op.indicatorScalarInput.setValue(indicator, op);
        op.sigma2ScalarInput.setValue(sigma2, op);
        if (rootLogRate != null) {
            op.rootLogRateScalarInput.setValue(rootLogRate, op);
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

    private static ACSubtreeUIncrementOperator typedAcSubtree(final Tree tree,
                                                              final RealVectorParam<?> rates,
                                                              final IntScalarParam<?> indicator,
                                                              final RealScalarParam<?> sigma2,
                                                              final RealScalarParam<?> rootLogRate,
                                                              final double delta) {
        final ACSubtreeUIncrementOperator op = new ACSubtreeUIncrementOperator();
        op.treeInput.setValue(tree, op);
        op.ratesVectorInput.setValue(rates, op);
        op.indicatorScalarInput.setValue(indicator, op);
        op.sigma2ScalarInput.setValue(sigma2, op);
        if (rootLogRate != null) {
            op.rootLogRateScalarInput.setValue(rootLogRate, op);
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

    private static UCACSwitchBridgeOperator typedBridge(final Tree tree,
                                                        final RealVectorParam<?> rates,
                                                        final IntScalarParam<?> indicator,
                                                        final RealScalarParam<?> stdev,
                                                        final RealScalarParam<?> sigma2,
                                                        final RealScalarParam<?> rootLogRate) {
        final UCACSwitchBridgeOperator op = new UCACSwitchBridgeOperator();
        op.treeInput.setValue(tree, op);
        op.ratesVectorInput.setValue(rates, op);
        op.indicatorScalarInput.setValue(indicator, op);
        op.ucldStdevScalarInput.setValue(stdev, op);
        op.sigma2ScalarInput.setValue(sigma2, op);
        if (rootLogRate != null) {
            op.rootLogRateScalarInput.setValue(rootLogRate, op);
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

    private static AlphaAnnealingOperator typedAlpha(final RealScalarParam<?> alpha,
                                                     final double start,
                                                     final double end,
                                                     final int steps) {
        final AlphaAnnealingOperator op = new AlphaAnnealingOperator();
        op.alphaScalarInput.setValue(alpha, op);
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

    private static RealVectorParam<PositiveReal> positiveRates(final double... values) {
        return new RealVectorParam<>(values, PositiveReal.INSTANCE);
    }

    private static IntScalarParam<Int> intScalar(final int value) {
        return new IntScalarParam<>(value, Int.INSTANCE);
    }

    private static RealScalarParam<PositiveReal> positiveScalar(final double value) {
        return new RealScalarParam<>(value, PositiveReal.INSTANCE);
    }

    private static RealScalarParam<Real> realScalar(final double value) {
        return new RealScalarParam<>(value, Real.INSTANCE);
    }

    private static RealParameter increasingLegacyRates(final Tree tree) {
        return realParameter(increasingValues(tree));
    }

    private static RealVectorParam<PositiveReal> increasingTypedRates(final Tree tree) {
        return positiveRates(increasingValues(tree));
    }

    private static double[] increasingValues(final Tree tree) {
        final double[] values = new double[tree.getNodeCount() - 1];
        for (int i = 0; i < values.length; i++) {
            values[i] = 0.8 + i * 0.15;
        }
        return values;
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
