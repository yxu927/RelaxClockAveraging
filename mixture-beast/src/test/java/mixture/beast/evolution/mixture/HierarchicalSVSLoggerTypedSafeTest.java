package mixture.beast.evolution.mixture;

import beast.base.evolution.likelihood.GenericTreeLikelihood;
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

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class HierarchicalSVSLoggerTypedSafeTest {

    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";

    @Test
    public void typedTopMixtureDefaultThreeModelWeightsMatchLegacy() {
        final HierarchicalSVSLogger legacy = logger(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                "hier"
        );
        final HierarchicalSVSLogger typed = logger(
                typedTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                "hier"
        );

        assertLoggerOutputsEqual(legacy, typed);
    }

    @Test
    public void typedSvsPriorDefaultThreeModelWeightsMatchLegacy() {
        final HierarchicalSVSLogger legacy = logger(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                "hier"
        );
        final HierarchicalSVSLogger typed = logger(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                typedSvsPrior(),
                "hier"
        );

        assertLoggerOutputsEqual(legacy, typed);
    }

    @Test
    public void typedTopAndTypedSvsTogetherMatchLegacy() {
        final HierarchicalSVSLogger legacy = logger(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                "hier"
        );
        final HierarchicalSVSLogger typed = logger(
                typedTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                typedSvsPrior(),
                "hier"
        );

        assertLoggerOutputsEqual(legacy, typed);
    }

    @Test
    public void typedInnerWeightsVectorMatchesLegacyInnerWeights() {
        final HierarchicalSVSLogger legacy = loggerWithLegacyInnerWeights(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                new RealParameter("0.25 0.75"),
                "hier"
        );
        final HierarchicalSVSLogger typed = loggerWithTypedInnerWeights(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                new RealVectorParam<>(new double[]{0.25, 0.75}, Real.INSTANCE),
                "hier"
        );

        assertLoggerOutputsEqual(legacy, typed);
    }

    @Test
    public void typedPrintTopMatchesLegacy() {
        final HierarchicalSVSLogger legacy = logger(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                "hier"
        );
        final HierarchicalSVSLogger typed = logger(
                typedTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                "hier"
        );
        configureTopOnly(legacy);
        configureTopOnly(typed);

        assertLoggerOutputsEqual(legacy, typed);
    }

    @Test
    public void typedPrintInnerMatchesLegacy() {
        final HierarchicalSVSLogger legacy = logger(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                "hier"
        );
        final HierarchicalSVSLogger typed = logger(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                typedSvsPrior(),
                "hier"
        );
        configureInnerOnly(legacy);
        configureInnerOnly(typed);

        assertLoggerOutputsEqual(legacy, typed);
    }

    @Test
    public void typedRbStatsMatchLegacy() {
        final HierarchicalSVSLogger legacy = loggerWithLegacyInnerWeights(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                new RealParameter("0.25 0.75"),
                "hier"
        );
        final HierarchicalSVSLogger typed = loggerWithTypedInnerWeights(
                typedTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                typedSvsPrior(),
                new RealVectorParam<>(new double[]{0.25, 0.75}, Real.INSTANCE),
                "hier"
        );
        configureRbOnly(legacy);
        configureRbOnly(typed);

        assertLoggerOutputsEqual(legacy, typed);
    }

    @Test
    public void typedXmlLikeRbStatsDefaultAlsoPrintsThreeModelWeights() {
        final HierarchicalSVSLogger legacy = logger(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                "hier"
        );
        final HierarchicalSVSLogger typed = logger(
                typedTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                typedSvsPrior(),
                "hier"
        );
        legacy.printRBStatsInput.setValue(true, legacy);
        typed.printRBStatsInput.setValue(true, typed);
        legacy.initAndValidate();
        typed.initAndValidate();

        final String[] legacyHeader = tokens(captureInit(legacy));
        final String[] typedHeader = tokens(captureInit(typed));
        final String[] legacyValues = tokens(captureLog(legacy));
        final String[] typedValues = tokens(captureLog(typed));

        assertEquals(11, legacyHeader.length);
        assertArrayEquals(legacyHeader, typedHeader);
        assertArrayEquals(legacyValues, typedValues);
        for (final String token : typedValues) {
            assertTrue(Double.isFinite(parse(token)));
        }
    }

    @Test
    public void innerWeightsAndInnerWeightsVectorCannotBothBeSpecified() {
        final HierarchicalSVSLogger logger = loggerWithLegacyInnerWeights(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                new RealParameter("0.25 0.75"),
                "hier"
        );
        logger.innerWeightsVectorInput.setValue(new RealVectorParam<>(new double[]{0.25, 0.75}, Real.INSTANCE),
                logger);

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
    }

    @Test
    public void invalidTypedInnerWeightsDimensionRejected() {
        final HierarchicalSVSLogger logger = loggerWithTypedInnerWeights(
                legacyTopMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0),
                legacySvsPrior(),
                new RealVectorParam<>(new double[]{1.0}, Real.INSTANCE),
                "hier"
        );

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
    }

    @Test
    public void typedTopWithoutWeightsRejectedClearly() {
        final MixtureTreeLikelihood top = rawTopMixture(-10.0, -12.0);
        final HierarchicalSVSLogger logger = logger(top, legacySvsPrior(), "hier");

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
    }

    private static HierarchicalSVSLogger logger(final MixtureTreeLikelihood top,
                                                final RelaxedRatesPriorSVS svs,
                                                final String id) {
        final HierarchicalSVSLogger logger = new HierarchicalSVSLogger();
        logger.setID(id);
        logger.topMixtureInput.setValue(top, logger);
        logger.svsPriorInput.setValue(svs, logger);
        return logger;
    }

    private static HierarchicalSVSLogger loggerWithLegacyInnerWeights(final MixtureTreeLikelihood top,
                                                                      final RelaxedRatesPriorSVS svs,
                                                                      final RealParameter innerWeights,
                                                                      final String id) {
        final HierarchicalSVSLogger logger = logger(top, svs, id);
        logger.innerWeightsInput.setValue(innerWeights, logger);
        return logger;
    }

    private static HierarchicalSVSLogger loggerWithTypedInnerWeights(final MixtureTreeLikelihood top,
                                                                     final RelaxedRatesPriorSVS svs,
                                                                     final RealVectorParam<Real> innerWeights,
                                                                     final String id) {
        final HierarchicalSVSLogger logger = logger(top, svs, id);
        logger.innerWeightsVectorInput.setValue(innerWeights, logger);
        return logger;
    }

    private static void configureTopOnly(final HierarchicalSVSLogger logger) {
        logger.printTopInput.setValue(true, logger);
        logger.printThreeModelWeightsInput.setValue(false, logger);
    }

    private static void configureInnerOnly(final HierarchicalSVSLogger logger) {
        logger.printInnerInput.setValue(true, logger);
        logger.printThreeModelWeightsInput.setValue(false, logger);
    }

    private static void configureRbOnly(final HierarchicalSVSLogger logger) {
        logger.printThreeModelWeightsInput.setValue(false, logger);
        logger.printRBStatsInput.setValue(true, logger);
    }

    private static MixtureTreeLikelihood legacyTopMixture(final double[] weights, final double... logL) {
        final MixtureTreeLikelihood mix = rawTopMixture(logL);
        mix.weightsInput.setValue(realParameter(weights), mix);
        mix.initAndValidate();
        return mix;
    }

    private static MixtureTreeLikelihood typedTopMixture(final double[] weights, final double... logL) {
        final MixtureTreeLikelihood mix = rawTopMixture(logL);
        mix.weightsVectorInput.setValue(new RealVectorParam<>(weights, Real.INSTANCE), mix);
        mix.initAndValidate();
        return mix;
    }

    private static MixtureTreeLikelihood rawTopMixture(final double... logL) {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        for (final double value : logL) {
            mix.subLikelihoodsInput.get().add(new ConstantLikelihood(value));
        }
        return mix;
    }

    private static RelaxedRatesPriorSVS legacySvsPrior() {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.initByName(
                "tree", fixedTree(),
                "rates", new RealParameter("0.8 1.1 1.4 0.9"),
                "indicator", new IntegerParameter("0"),
                "ucldStdev", new RealParameter("0.5"),
                "rootLogRate", new RealParameter("0.0"),
                "sigma2", new RealParameter("0.2"),
                "minBranchLength", 1.0e-12
        );
        return prior;
    }

    private static RelaxedRatesPriorSVS typedSvsPrior() {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.initByName(
                "tree", fixedTree(),
                "ratesVector", new RealVectorParam<>(new double[]{0.8, 1.1, 1.4, 0.9}, PositiveReal.INSTANCE),
                "indicatorScalar", new IntScalarParam<>(0, Int.INSTANCE),
                "ucldStdevScalar", new RealScalarParam<>(0.5, PositiveReal.INSTANCE),
                "rootLogRateScalar", new RealScalarParam<>(0.0, Real.INSTANCE),
                "sigma2Scalar", new RealScalarParam<>(0.2, PositiveReal.INSTANCE),
                "minBranchLength", 1.0e-12
        );
        return prior;
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealParameter realParameter(final double[] values) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < values.length; i++) {
            if (i > 0) {
                sb.append(' ');
            }
            sb.append(values[i]);
        }
        return new RealParameter(sb.toString());
    }

    private static void assertLoggerOutputsEqual(final HierarchicalSVSLogger legacyLogger,
                                                 final HierarchicalSVSLogger typedLogger) {
        legacyLogger.initAndValidate();
        typedLogger.initAndValidate();
        assertArrayEquals(tokens(captureInit(legacyLogger)), tokens(captureInit(typedLogger)));
        assertArrayEquals(tokens(captureLog(legacyLogger)), tokens(captureLog(typedLogger)));
    }

    private static String captureInit(final HierarchicalSVSLogger logger) {
        return capture(logger::init);
    }

    private static String captureLog(final HierarchicalSVSLogger logger) {
        return capture(out -> logger.log(0L, out));
    }

    private static String capture(final StreamAction action) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (PrintStream out = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            action.write(out);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }

    private static String[] tokens(final String output) {
        String normalized = output;
        while (normalized.endsWith("\t")) {
            normalized = normalized.substring(0, normalized.length() - 1);
        }
        if (normalized.isEmpty()) {
            return new String[0];
        }
        return normalized.split("\t");
    }

    private static double parse(final String token) {
        if ("NaN".equals(token)) {
            return Double.NaN;
        }
        if ("Inf".equals(token)) {
            return Double.POSITIVE_INFINITY;
        }
        if ("-Inf".equals(token)) {
            return Double.NEGATIVE_INFINITY;
        }
        return Double.parseDouble(token);
    }

    private interface StreamAction {
        void write(PrintStream out);
    }

    public static final class ConstantLikelihood extends GenericTreeLikelihood {
        private final double value;

        private ConstantLikelihood(final double value) {
            this.value = value;
        }

        @Override
        public void initAndValidate() {
            // Do not require data/tree/siteModel for this constant test stub.
        }

        @Override
        public double calculateLogP() {
            logP = value;
            return logP;
        }
    }
}
