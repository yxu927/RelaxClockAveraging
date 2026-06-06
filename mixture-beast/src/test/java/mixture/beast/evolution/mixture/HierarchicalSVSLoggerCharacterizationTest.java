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

public class HierarchicalSVSLoggerCharacterizationTest {

    private static final double EPS = 1.0e-6;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";

    @Test
    public void defaultLoggerPrintsThreeModelWeightsOnly() {
        final MixtureTreeLikelihood top = legacyMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final RelaxedRatesPriorSVS svs = legacySvsPrior(new IntegerParameter("0"));
        final HierarchicalSVSLogger logger = logger(top, svs, null, "hier");
        logger.initAndValidate();

        assertArrayEquals(new String[]{"hier.pStrict", "hier.pRelaxUC", "hier.pRelaxAC"}, tokens(captureInit(logger)));

        final String[] values = tokens(captureLog(logger));
        final TopExpected topExpected = topExpected(new double[]{1.0 / 3.0, 2.0 / 3.0}, new double[]{-10.0, -12.0});
        final InnerExpected innerExpected = innerExpected(svs, new double[]{0.5, 0.5});

        assertEquals(topExpected.rStrict, parse(values[0]), EPS);
        assertEquals(topExpected.rRelax * innerExpected.rUC, parse(values[1]), EPS);
        assertEquals(topExpected.rRelax * innerExpected.rAC, parse(values[2]), EPS);
    }

    @Test
    public void printTopHeaderAndValuesMatchSource() {
        final MixtureTreeLikelihood top = legacyMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final RelaxedRatesPriorSVS svs = legacySvsPrior(new IntegerParameter("0"));
        final HierarchicalSVSLogger logger = logger(top, svs, null, "hier");
        logger.printTopInput.setValue(true, logger);
        logger.printThreeModelWeightsInput.setValue(false, logger);
        logger.initAndValidate();

        assertArrayEquals(new String[]{
                "hier.topWStrict",
                "hier.topWRelax",
                "hier.topLogLStrict",
                "hier.topLogLRelax",
                "hier.rStrict",
                "hier.rRelax"
        }, tokens(captureInit(logger)));

        final String[] values = tokens(captureLog(logger));
        final TopExpected expected = topExpected(new double[]{1.0 / 3.0, 2.0 / 3.0}, new double[]{-10.0, -12.0});

        assertEquals(1.0 / 3.0, parse(values[0]), EPS);
        assertEquals(2.0 / 3.0, parse(values[1]), EPS);
        assertEquals(-10.0, parse(values[2]), EPS);
        assertEquals(-12.0, parse(values[3]), EPS);
        assertEquals(expected.rStrict, parse(values[4]), EPS);
        assertEquals(expected.rRelax, parse(values[5]), EPS);
    }

    @Test
    public void printInnerWithoutInnerWeightsUsesDefaultHalfHalfWeights() {
        final MixtureTreeLikelihood top = legacyMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final RelaxedRatesPriorSVS svs = legacySvsPrior(new IntegerParameter("0"));
        final HierarchicalSVSLogger logger = logger(top, svs, null, "hier");
        logger.printInnerInput.setValue(true, logger);
        logger.printThreeModelWeightsInput.setValue(false, logger);
        logger.initAndValidate();

        assertArrayEquals(new String[]{"hier.logPUC", "hier.logPAC"}, tokens(captureInit(logger)));

        final String[] values = tokens(captureLog(logger));

        assertEquals(svs.logPriorUCOnly(), parse(values[0]), EPS);
        assertEquals(svs.logPriorACOnly(), parse(values[1]), EPS);
    }

    @Test
    public void explicitInnerWeightsAffectThreeModelWeights() {
        final MixtureTreeLikelihood top = legacyMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final RelaxedRatesPriorSVS svs = legacySvsPrior(new IntegerParameter("0"));
        final RealParameter innerWeights = new RealParameter("0.25 0.75");
        final HierarchicalSVSLogger logger = logger(top, svs, innerWeights, "hier");
        logger.initAndValidate();

        final String[] values = tokens(captureLog(logger));
        final TopExpected topExpected = topExpected(new double[]{1.0 / 3.0, 2.0 / 3.0}, new double[]{-10.0, -12.0});
        final InnerExpected innerExpected = innerExpected(svs, new double[]{0.25, 0.75});

        assertEquals(topExpected.rStrict, parse(values[0]), EPS);
        assertEquals(topExpected.rRelax * innerExpected.rUC, parse(values[1]), EPS);
        assertEquals(topExpected.rRelax * innerExpected.rAC, parse(values[2]), EPS);
    }

    @Test
    public void printRbStatsHeaderAndFiniteOutputMatchSource() {
        final MixtureTreeLikelihood top = legacyMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final RelaxedRatesPriorSVS svs = legacySvsPrior(new IntegerParameter("0"));
        final HierarchicalSVSLogger logger = logger(top, svs, new RealParameter("0.25 0.75"), "hier");
        logger.printThreeModelWeightsInput.setValue(false, logger);
        logger.printRBStatsInput.setValue(true, logger);
        logger.initAndValidate();

        assertArrayEquals(new String[]{
                "hier.rb_ucld_num",
                "hier.rb_ucld_den",
                "hier.rb_sigma2_num",
                "hier.rb_sigma2_den",
                "hier.rb_relax_ucld_num",
                "hier.rb_relax_ucld_den",
                "hier.rb_relax_sigma2_num",
                "hier.rb_relax_sigma2_den"
        }, tokens(captureInit(logger)));

        final String[] values = tokens(captureLog(logger));
        final TopExpected topExpected = topExpected(new double[]{1.0 / 3.0, 2.0 / 3.0}, new double[]{-10.0, -12.0});
        final InnerExpected innerExpected = innerExpected(svs, new double[]{0.25, 0.75});
        final double pRelaxUC = topExpected.rRelax * innerExpected.rUC;
        final double pRelaxAC = topExpected.rRelax * innerExpected.rAC;
        final double ucldStdev = 0.5;
        final double sigma2 = 0.2;

        assertEquals(innerExpected.rUC * ucldStdev, parse(values[0]), EPS);
        assertEquals(innerExpected.rUC, parse(values[1]), EPS);
        assertEquals(innerExpected.rAC * sigma2, parse(values[2]), EPS);
        assertEquals(innerExpected.rAC, parse(values[3]), EPS);
        assertEquals(pRelaxUC * ucldStdev, parse(values[4]), EPS);
        assertEquals(pRelaxUC, parse(values[5]), EPS);
        assertEquals(pRelaxAC * sigma2, parse(values[6]), EPS);
        assertEquals(pRelaxAC, parse(values[7]), EPS);
    }

    @Test
    public void xmlLikeRbStatsDefaultAlsoPrintsThreeModelWeights() {
        final MixtureTreeLikelihood top = legacyMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final RelaxedRatesPriorSVS svs = legacySvsPrior(new IntegerParameter("0"));
        final HierarchicalSVSLogger logger = logger(top, svs, null, "hier");
        logger.printRBStatsInput.setValue(true, logger);
        logger.initAndValidate();

        final String[] header = tokens(captureInit(logger));
        final String[] values = tokens(captureLog(logger));

        assertEquals(11, header.length);
        assertEquals(header.length, values.length);
        for (final String token : values) {
            assertTrue(Double.isFinite(parse(token)));
        }
    }

    @Test
    public void typedTopMixtureNowMatchesEquivalentLegacyTopMixture() {
        final MixtureTreeLikelihood legacyTop = legacyMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final MixtureTreeLikelihood typedTop = typedMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final HierarchicalSVSLogger legacyLogger = logger(legacyTop, legacySvsPrior(new IntegerParameter("0")), null,
                "hier");
        final HierarchicalSVSLogger typedLogger = logger(typedTop, legacySvsPrior(new IntegerParameter("0")), null,
                "hier");
        legacyLogger.initAndValidate();
        typedLogger.initAndValidate();

        assertArrayEquals(tokens(captureInit(legacyLogger)), tokens(captureInit(typedLogger)));
        assertArrayEquals(tokens(captureLog(legacyLogger)), tokens(captureLog(typedLogger)));
    }

    @Test
    public void typedSvsPriorNowMatchesEquivalentLegacySvsPrior() {
        final MixtureTreeLikelihood top = legacyMixture(new double[]{1.0 / 3.0, 2.0 / 3.0}, -10.0, -12.0);
        final HierarchicalSVSLogger legacyLogger = logger(top, legacySvsPrior(new IntegerParameter("0")), null,
                "hier");
        final HierarchicalSVSLogger typedLogger = logger(top, typedSvsPrior(), null, "hier");
        legacyLogger.initAndValidate();
        typedLogger.initAndValidate();

        assertArrayEquals(tokens(captureInit(legacyLogger)), tokens(captureInit(typedLogger)));
        assertArrayEquals(tokens(captureLog(legacyLogger)), tokens(captureLog(typedLogger)));
    }

    private static HierarchicalSVSLogger logger(final MixtureTreeLikelihood top,
                                                final RelaxedRatesPriorSVS svs,
                                                final RealParameter innerWeights,
                                                final String id) {
        final HierarchicalSVSLogger logger = new HierarchicalSVSLogger();
        logger.setID(id);
        logger.topMixtureInput.setValue(top, logger);
        logger.svsPriorInput.setValue(svs, logger);
        if (innerWeights != null) {
            logger.innerWeightsInput.setValue(innerWeights, logger);
        }
        return logger;
    }

    private static MixtureTreeLikelihood legacyMixture(final double[] weights, final double... logL) {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        for (final double value : logL) {
            mix.subLikelihoodsInput.get().add(new ConstantLikelihood(value));
        }
        mix.weightsInput.setValue(realParameter(weights), mix);
        mix.initAndValidate();
        return mix;
    }

    private static MixtureTreeLikelihood typedMixture(final double[] weights, final double... logL) {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        for (final double value : logL) {
            mix.subLikelihoodsInput.get().add(new ConstantLikelihood(value));
        }
        mix.weightsVectorInput.setValue(new RealVectorParam<>(weights, Real.INSTANCE), mix);
        mix.initAndValidate();
        return mix;
    }

    private static RelaxedRatesPriorSVS legacySvsPrior(final IntegerParameter indicator) {
        final RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.initByName(
                "tree", fixedTree(),
                "rates", new RealParameter("0.8 1.1 1.4 0.9"),
                "indicator", indicator,
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

    private static TopExpected topExpected(final double[] weights, final double[] logL) {
        final double sStrict = Math.log(weights[0]) + logL[0];
        final double sRelax = Math.log(weights[1]) + logL[1];
        final double m = Math.max(sStrict, sRelax);
        final double denom = Math.exp(sStrict - m) + Math.exp(sRelax - m);
        final double rStrict = Math.exp(sStrict - m) / denom;
        return new TopExpected(rStrict, 1.0 - rStrict);
    }

    private static InnerExpected innerExpected(final RelaxedRatesPriorSVS svs, final double[] weights) {
        final double logPUC = svs.logPriorUCOnly();
        final double logPAC = svs.logPriorACOnly();
        final double tUC = Math.log(weights[0]) + logPUC;
        final double tAC = Math.log(weights[1]) + logPAC;
        final double m = Math.max(tUC, tAC);
        final double denom = Math.exp(tUC - m) + Math.exp(tAC - m);
        final double rUC = Math.exp(tUC - m) / denom;
        return new InnerExpected(rUC, 1.0 - rUC);
    }

    private interface StreamAction {
        void write(PrintStream out);
    }

    private static final class TopExpected {
        private final double rStrict;
        private final double rRelax;

        private TopExpected(final double rStrict, final double rRelax) {
            this.rStrict = rStrict;
            this.rRelax = rRelax;
        }
    }

    private static final class InnerExpected {
        private final double rUC;
        private final double rAC;

        private InnerExpected(final double rUC, final double rAC) {
            this.rUC = rUC;
            this.rAC = rAC;
        }
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
