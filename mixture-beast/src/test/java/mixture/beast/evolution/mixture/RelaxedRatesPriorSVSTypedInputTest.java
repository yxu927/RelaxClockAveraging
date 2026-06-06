package mixture.beast.evolution.mixture;

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
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class RelaxedRatesPriorSVSTypedInputTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";

    @Test
    public void typedIndicatorZeroCalculateLogPMatchesUCOnly() {
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                indicator(0), positiveScalar(0.5), null, positiveScalar(0.2), 1.0e-12);

        assertEquals(prior.logPriorUCOnly(), prior.calculateLogP(), EPS);
        assertTrue(Double.isFinite(prior.calculateLogP()));
    }

    @Test
    public void typedIndicatorOneCalculateLogPMatchesACOnly() {
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                indicator(1), positiveScalar(0.5), realScalar(0.0), positiveScalar(0.2), 1.0e-12);

        assertEquals(prior.logPriorACOnly(), prior.calculateLogP(), EPS);
        assertTrue(Double.isFinite(prior.calculateLogP()));
    }

    @Test
    public void typedInvalidIndicatorReturnsNegativeInfinity() {
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                indicator(2), positiveScalar(0.5), null, positiveScalar(0.2), 1.0e-12);

        assertEquals(Double.NEGATIVE_INFINITY, prior.calculateLogP(), 0.0);
    }

    @Test
    public void typedMissingRootLogRateEqualsExplicitZeroRootLogRate() {
        final Tree tree = fixedTree();
        final RelaxedRatesPriorSVS missing = typedPrior(tree, rates(0.8, 1.1, 1.4, 0.9),
                indicator(1), positiveScalar(0.5), null, positiveScalar(0.2), 1.0e-12);
        final RelaxedRatesPriorSVS explicitZero = typedPrior(tree, rates(0.8, 1.1, 1.4, 0.9),
                indicator(1), positiveScalar(0.5), realScalar(0.0), positiveScalar(0.2), 1.0e-12);

        assertEquals(missing.logPriorACOnly(), explicitZero.logPriorACOnly(), EPS);
    }

    @Test
    public void typedRootLogRateChangesACPrior() {
        final Tree tree = fixedTree();
        final RelaxedRatesPriorSVS zero = typedPrior(tree, rates(0.8, 1.1, 1.4, 0.9),
                indicator(1), positiveScalar(0.5), realScalar(0.0), positiveScalar(0.2), 1.0e-12);
        final RelaxedRatesPriorSVS shifted = typedPrior(tree, rates(0.8, 1.1, 1.4, 0.9),
                indicator(1), positiveScalar(0.5), realScalar(0.7), positiveScalar(0.2), 1.0e-12);

        assertNotEquals(zero.logPriorACOnly(), shifted.logPriorACOnly(), EPS);
    }

    @Test
    public void typedUCNonPositiveStdevRejectedByDomainOrPrior() {
        try {
            final RealScalarParam<PositiveReal> stdev = positiveScalar(0.0);
            final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                    indicator(0), stdev, null, positiveScalar(0.2), 1.0e-12);
            assertEquals(Double.NEGATIVE_INFINITY, prior.logPriorUCOnly(), 0.0);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("not valid"));
        }
    }

    @Test
    public void typedACNonPositiveSigma2RejectedByDomainOrPrior() {
        try {
            final RealScalarParam<PositiveReal> sigma2 = positiveScalar(0.0);
            final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                    indicator(1), positiveScalar(0.5), null, sigma2, 1.0e-12);
            assertEquals(Double.NEGATIVE_INFINITY, prior.logPriorACOnly(), 0.0);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("not valid"));
        }
    }

    @Test
    public void typedPriorsRejectNonPositiveRatesByDomainOrPrior() {
        try {
            final RealVectorParam<PositiveReal> invalidRates = rates(0.8, 0.0, 1.4, 0.9);
            final RelaxedRatesPriorSVS uc = typedPrior(fixedTree(), invalidRates, indicator(0),
                    positiveScalar(0.5), null, positiveScalar(0.2), 1.0e-12);
            final RelaxedRatesPriorSVS ac = typedPrior(fixedTree(), invalidRates, indicator(1),
                    positiveScalar(0.5), null, positiveScalar(0.2), 1.0e-12);
            assertEquals(Double.NEGATIVE_INFINITY, uc.logPriorUCOnly(), 0.0);
            assertEquals(Double.NEGATIVE_INFINITY, ac.logPriorACOnly(), 0.0);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("not valid"));
        }
    }

    @Test
    public void typedScalarRatesExpandOrBroadcast() {
        final Tree tree = fixedTree();
        final RealVectorParam<PositiveReal> scalarRates = rates(1.25);
        final RelaxedRatesPriorSVS prior = typedPrior(tree, scalarRates, indicator(0),
                positiveScalar(0.5), null, positiveScalar(0.2), 1.0e-12);

        assertTrue(scalarRates.size() == tree.getNodeCount() - 1 || scalarRates.size() == 1);
        assertTrue(Double.isFinite(prior.calculateLogP()));
        assertTrue(Double.isFinite(prior.logPriorUCOnly()));
    }

    @Test
    public void typedDirtyRatesTriggerRequiresRecalculation() {
        final RealVectorParam<PositiveReal> rates = rates(0.8, 1.1, 1.4, 0.9);
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates, indicator(0),
                positiveScalar(0.5), realScalar(0.0), positiveScalar(0.2), 1.0e-12);

        rates.set(0, 0.95);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void typedDirtyIndicatorTriggerRequiresRecalculation() {
        final IntScalarParam<Int> indicator = indicator(0);
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                indicator, positiveScalar(0.5), realScalar(0.0), positiveScalar(0.2), 1.0e-12);

        indicator.set(1);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void typedDirtyUcldStdevTriggerRequiresRecalculation() {
        final RealScalarParam<PositiveReal> ucldStdev = positiveScalar(0.5);
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                indicator(0), ucldStdev, realScalar(0.0), positiveScalar(0.2), 1.0e-12);

        ucldStdev.set(0.6);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void typedDirtySigma2TriggerRequiresRecalculation() {
        final RealScalarParam<PositiveReal> sigma2 = positiveScalar(0.2);
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                indicator(1), positiveScalar(0.5), realScalar(0.0), sigma2, 1.0e-12);

        sigma2.set(0.3);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void typedDirtyRootLogRateTriggerRequiresRecalculation() {
        final RealScalarParam<Real> rootLogRate = realScalar(0.0);
        final RelaxedRatesPriorSVS prior = typedPrior(fixedTree(), rates(0.8, 1.1, 1.4, 0.9),
                indicator(1), positiveScalar(0.5), rootLogRate, positiveScalar(0.2), 1.0e-12);

        rootLogRate.set(0.1);

        assertTrue(prior.requiresRecalculation());
    }

    @Test
    public void typedGetArgumentsAndConditionsUseTypedIDs() {
        final Tree tree = fixedTree();
        tree.setID("tree");
        final RealVectorParam<PositiveReal> rates = rates(0.8, 1.1, 1.4, 0.9);
        rates.setID("ratesVector");
        final IntScalarParam<Int> indicator = indicator(0);
        indicator.setID("indicatorScalar");
        final RealScalarParam<PositiveReal> ucldStdev = positiveScalar(0.5);
        ucldStdev.setID("ucldStdevScalar");
        final RealScalarParam<Real> rootLogRate = realScalar(0.0);
        rootLogRate.setID("rootLogRateScalar");
        final RealScalarParam<PositiveReal> sigma2 = positiveScalar(0.2);
        sigma2.setID("sigma2Scalar");

        final RelaxedRatesPriorSVS prior = typedPrior(tree, rates, indicator, ucldStdev,
                rootLogRate, sigma2, 1.0e-12);
        final List<String> arguments = prior.getArguments();
        final List<String> conditions = prior.getConditions();

        assertTrue(arguments.contains("ratesVector"));
        assertTrue(conditions.contains("tree"));
        assertTrue(conditions.contains("indicatorScalar"));
        assertTrue(conditions.contains("ucldStdevScalar"));
        assertTrue(conditions.contains("rootLogRateScalar"));
        assertTrue(conditions.contains("sigma2Scalar"));
    }

    @Test
    public void legacyAndTypedRatesCannotBothBeSpecified() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.ratesInput.setValue(new RealParameter("0.8 1.1 1.4 0.9"), prior);
        prior.ratesVectorInput.setValue(rates(0.8, 1.1, 1.4, 0.9), prior);

        assertThrows(IllegalArgumentException.class, prior::initAndValidate);
    }

    @Test
    public void legacyAndTypedIndicatorCannotBothBeSpecified() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.indicatorInput.setValue(new IntegerParameter("0"), prior);
        prior.indicatorScalarInput.setValue(indicator(0), prior);

        assertThrows(IllegalArgumentException.class, prior::initAndValidate);
    }

    @Test
    public void legacyAndTypedUcldStdevCannotBothBeSpecified() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.ucldStdevInput.setValue(new RealParameter("0.5"), prior);
        prior.ucldStdevScalarInput.setValue(positiveScalar(0.5), prior);

        assertThrows(IllegalArgumentException.class, prior::initAndValidate);
    }

    @Test
    public void legacyAndTypedRootLogRateCannotBothBeSpecified() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.rootLogRateInput.setValue(new RealParameter("0.0"), prior);
        prior.rootLogRateScalarInput.setValue(realScalar(0.0), prior);

        assertThrows(IllegalArgumentException.class, prior::initAndValidate);
    }

    @Test
    public void legacyAndTypedSigma2CannotBothBeSpecified() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.sigma2Input.setValue(new RealParameter("0.2"), prior);
        prior.sigma2ScalarInput.setValue(positiveScalar(0.2), prior);

        assertThrows(IllegalArgumentException.class, prior::initAndValidate);
    }

    @Test
    public void missingMandatoryInputsAreRejected() {
        assertThrows(IllegalArgumentException.class, missingRatesPrior()::initAndValidate);
        assertThrows(IllegalArgumentException.class, missingIndicatorPrior()::initAndValidate);
        assertThrows(IllegalArgumentException.class, missingUcldStdevPrior()::initAndValidate);
        assertThrows(IllegalArgumentException.class, missingSigma2Prior()::initAndValidate);
    }

    @Test
    public void existingLegacyPathStillWorks() {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.initByName(
                "tree", fixedTree(),
                "rates", new RealParameter("0.8 1.1 1.4 0.9"),
                "indicator", new IntegerParameter("1"),
                "ucldStdev", new RealParameter("0.5"),
                "rootLogRate", new RealParameter("0.0"),
                "sigma2", new RealParameter("0.2"),
                "minBranchLength", 1.0e-12
        );

        assertEquals(prior.logPriorACOnly(), prior.calculateLogP(), EPS);
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

    private static RelaxedRatesPriorSVS typedPrior(final Tree tree,
                                                   final RealVectorParam<PositiveReal> rates,
                                                   final IntScalarParam<Int> indicator,
                                                   final RealScalarParam<PositiveReal> ucldStdev,
                                                   final RealScalarParam<Real> rootLogRate,
                                                   final RealScalarParam<PositiveReal> sigma2,
                                                   final double minBranchLength) {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        if (rootLogRate == null) {
            prior.initByName(
                    "tree", tree,
                    "ratesVector", rates,
                    "indicatorScalar", indicator,
                    "ucldStdevScalar", ucldStdev,
                    "sigma2Scalar", sigma2,
                    "minBranchLength", minBranchLength
            );
        } else {
            prior.initByName(
                    "tree", tree,
                    "ratesVector", rates,
                    "indicatorScalar", indicator,
                    "ucldStdevScalar", ucldStdev,
                    "rootLogRateScalar", rootLogRate,
                    "sigma2Scalar", sigma2,
                    "minBranchLength", minBranchLength
            );
        }
        return prior;
    }

    private static RelaxedRatesPriorSVS requiredTypedPriorShell() {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.treeInput.setValue(fixedTree(), prior);
        prior.ratesVectorInput.setValue(rates(0.8, 1.1, 1.4, 0.9), prior);
        prior.indicatorScalarInput.setValue(indicator(0), prior);
        prior.ucldStdevScalarInput.setValue(positiveScalar(0.5), prior);
        prior.sigma2ScalarInput.setValue(positiveScalar(0.2), prior);
        return prior;
    }

    private static RelaxedRatesPriorSVS missingRatesPrior() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.ratesVectorInput.setValue(null, prior);
        return prior;
    }

    private static RelaxedRatesPriorSVS missingIndicatorPrior() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.indicatorScalarInput.setValue(null, prior);
        return prior;
    }

    private static RelaxedRatesPriorSVS missingUcldStdevPrior() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.ucldStdevScalarInput.setValue(null, prior);
        return prior;
    }

    private static RelaxedRatesPriorSVS missingSigma2Prior() {
        final RelaxedRatesPriorSVS prior = requiredTypedPriorShell();
        prior.sigma2ScalarInput.setValue(null, prior);
        return prior;
    }
}
