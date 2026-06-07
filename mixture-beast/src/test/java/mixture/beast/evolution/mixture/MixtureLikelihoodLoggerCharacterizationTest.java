package mixture.beast.evolution.mixture;

import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

public class MixtureLikelihoodLoggerCharacterizationTest {

    private static final double EPS = 1.0e-8;

    @Test
    public void initRejectsTopMixtureWithFewerThanTwoComponents() {
        final MixtureTreeLikelihood mix = rawMixture(new double[]{1.0}, null, -1.0);
        final MixtureLikelihoodLogger logger = logger(mix, "mixLog");

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
    }

    @Test
    public void initRejectsWeightsDimensionMismatch() {
        final MixtureTreeLikelihood mix = rawMixture(new double[]{1.0}, null, -1.0, -2.0);
        final MixtureLikelihoodLogger logger = logger(mix, "mixLog");

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
    }

    @Test
    public void initRejectsAlphaDimensionMismatch() {
        final MixtureTreeLikelihood mix = rawMixture(new double[]{0.5, 0.5}, "0.0 1.0", -1.0, -2.0);
        final MixtureLikelihoodLogger logger = logger(mix, "mixLog");

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
    }

    @Test
    public void defaultLoggerPrintsResponsibilitiesOnly() {
        final MixtureTreeLikelihood mix = legacyMixture(new double[]{0.25, 0.75}, null, -10.0, -12.0);
        final MixtureLikelihoodLogger logger = logger(mix, "mixLog");
        logger.initAndValidate();

        assertArrayEquals(new String[]{"mixLog.r[1]", "mixLog.r[2]"}, tokens(captureInit(logger)));

        final String[] values = tokens(captureLog(logger));
        final ExpectedMixture expected = expected(new double[]{0.25, 0.75}, new double[]{-10.0, -12.0}, 0.0);

        assertEquals(expected.resp[0], parse(values[0]), EPS);
        assertEquals(expected.resp[1], parse(values[1]), EPS);
    }

    @Test
    public void allPrintFlagsProduceExpectedHeaders() {
        final MixtureTreeLikelihood mix = legacyMixture(new double[]{0.25, 0.75}, null, -10.0, -12.0);
        final MixtureLikelihoodLogger logger = loggerWithAllFlags(mix, "mixLog");
        logger.initAndValidate();

        assertArrayEquals(new String[]{
                "mixLog.totalLogP",
                "mixLog.logMix",
                "mixLog.alpha",
                "mixLog.couplingTerm",
                "mixLog.sumW",
                "mixLog.w[1]",
                "mixLog.w[2]",
                "mixLog.logL[1]",
                "mixLog.logL[2]",
                "mixLog.logS[1]",
                "mixLog.logS[2]",
                "mixLog.r[1]",
                "mixLog.r[2]",
                "mixLog.maxLogL",
                "mixLog.maxLogS",
                "mixLog.logMixMinusMaxLogL",
                "mixLog.logMixMinusMaxLogS",
                "mixLog.violationMixGTMaxLogL"
        }, tokens(captureInit(logger)));
    }

    @Test
    public void allPrintFlagsProduceExpectedValuesWithoutAlpha() {
        final double[] weights = {0.25, 0.75};
        final double[] logL = {-10.0, -12.0};
        final MixtureTreeLikelihood mix = legacyMixture(weights, null, logL);
        final MixtureLikelihoodLogger logger = loggerWithAllFlags(mix, "mixLog");
        logger.initAndValidate();

        final String[] values = tokens(captureLog(logger));
        final ExpectedMixture expected = expected(weights, logL, 0.0);

        assertEquals(expected.totalLogP, parse(values[0]), EPS);
        assertEquals(expected.logMix, parse(values[1]), EPS);
        assertEquals(0.0, parse(values[2]), EPS);
        assertEquals(0.0, parse(values[3]), EPS);
        assertEquals(1.0, parse(values[4]), EPS);
        assertEquals(weights[0], parse(values[5]), EPS);
        assertEquals(weights[1], parse(values[6]), EPS);
        assertEquals(logL[0], parse(values[7]), EPS);
        assertEquals(logL[1], parse(values[8]), EPS);
        assertEquals(expected.logS[0], parse(values[9]), EPS);
        assertEquals(expected.logS[1], parse(values[10]), EPS);
        assertEquals(expected.resp[0], parse(values[11]), EPS);
        assertEquals(expected.resp[1], parse(values[12]), EPS);
        assertEquals(expected.maxLogL, parse(values[13]), EPS);
        assertEquals(expected.maxLogS, parse(values[14]), EPS);
        assertEquals(expected.logMix - expected.maxLogL, parse(values[15]), EPS);
        assertEquals(expected.logMix - expected.maxLogS, parse(values[16]), EPS);
        assertEquals(0.0, parse(values[17]), EPS);
    }

    @Test
    public void positiveAlphaProducesCouplingTermAndTotal() {
        final double[] weights = {0.25, 0.75};
        final double[] logL = {-10.0, -12.0};
        final double alpha = 0.2;
        final MixtureTreeLikelihood mix = legacyMixture(weights, alpha, logL);
        final MixtureLikelihoodLogger logger = loggerWithAllFlags(mix, "mixLog");
        logger.initAndValidate();

        final String[] values = tokens(captureLog(logger));
        final ExpectedMixture expected = expected(weights, logL, alpha);

        assertEquals(expected.totalLogP, parse(values[0]), EPS);
        assertEquals(expected.logMix, parse(values[1]), EPS);
        assertEquals(alpha, parse(values[2]), EPS);
        assertEquals(expected.couplingTerm, parse(values[3]), EPS);
    }

    @Test
    public void zeroWeightComponentHasZeroResponsibilityEvenWhenComponentLogLIsNan() {
        final MixtureTreeLikelihood mix = legacyMixture(new double[]{0.0, 1.0}, null, Double.NaN, -3.0);
        final MixtureLikelihoodLogger logger = logger(mix, "mixLog");
        logger.initAndValidate();

        final String[] values = tokens(captureLog(logger));

        assertEquals(0.0, parse(values[0]), EPS);
        assertEquals(1.0, parse(values[1]), EPS);
    }

    @Test
    public void nonFinitePositiveWeightComponentMatchesSourceResponsibilityHandling() {
        final MixtureTreeLikelihood mix = legacyMixture(
                new double[]{0.5, 0.5},
                null,
                Double.NEGATIVE_INFINITY,
                -2.0
        );
        final MixtureLikelihoodLogger logger = loggerWithAllFlags(mix, "mixLog");
        logger.initAndValidate();

        final String[] values = tokens(captureLog(logger));

        assertEquals(Double.NEGATIVE_INFINITY, parse(values[9]), 0.0);
        assertEquals(Math.log(0.5) - 2.0, parse(values[10]), EPS);
        assertEquals(0.0, parse(values[11]), EPS);
        assertEquals(1.0, parse(values[12]), EPS);
    }

    @Test
    public void closeIsNoOp() {
        final MixtureLikelihoodLogger logger = logger(legacyMixture(new double[]{0.5, 0.5}, null, -1.0, -2.0),
                "mixLog");
        logger.initAndValidate();

        assertEquals("", captureClose(logger));
    }

    @Test
    public void typedMixtureNowMatchesEquivalentLegacyMixture() {
        final MixtureTreeLikelihood legacy = legacyMixture(new double[]{0.25, 0.75}, 0.2, -10.0, -12.0);
        final MixtureTreeLikelihood typed = typedMixture(new double[]{0.25, 0.75}, 0.2, -10.0, -12.0);
        final MixtureLikelihoodLogger legacyLogger = loggerWithAllFlags(legacy, "mixLog");
        final MixtureLikelihoodLogger typedLogger = loggerWithAllFlags(typed, "mixLog");
        legacyLogger.initAndValidate();
        typedLogger.initAndValidate();

        assertArrayEquals(tokens(captureInit(legacyLogger)), tokens(captureInit(typedLogger)));
        assertArrayEquals(tokens(captureLog(legacyLogger)), tokens(captureLog(typedLogger)));
    }

    private static MixtureLikelihoodLogger logger(final MixtureTreeLikelihood mix, final String id) {
        final MixtureLikelihoodLogger logger = new MixtureLikelihoodLogger();
        logger.setID(id);
        logger.mixtureInput.setValue(mix, logger);
        return logger;
    }

    private static MixtureLikelihoodLogger loggerWithAllFlags(final MixtureTreeLikelihood mix, final String id) {
        final MixtureLikelihoodLogger logger = logger(mix, id);
        logger.printTotalLogPInput.setValue(true, logger);
        logger.printLogMixInput.setValue(true, logger);
        logger.printAlphaInput.setValue(true, logger);
        logger.printCouplingTermInput.setValue(true, logger);
        logger.printSumWeightsInput.setValue(true, logger);
        logger.printWeightsInput.setValue(true, logger);
        logger.printLogLInput.setValue(true, logger);
        logger.printLogSInput.setValue(true, logger);
        logger.printRespInput.setValue(true, logger);
        logger.printMaxLogLInput.setValue(true, logger);
        logger.printMaxLogSInput.setValue(true, logger);
        logger.printMixMinusMaxLogLInput.setValue(true, logger);
        logger.printMixMinusMaxLogSInput.setValue(true, logger);
        logger.printViolationFlagInput.setValue(true, logger);
        return logger;
    }

    private static MixtureTreeLikelihood legacyMixture(final double[] weights,
                                                       final Double alpha,
                                                       final double... logL) {
        final MixtureTreeLikelihood mix = rawMixture(weights, alpha == null ? null : Double.toString(alpha), logL);
        mix.initAndValidate();
        return mix;
    }

    private static MixtureTreeLikelihood rawMixture(final double[] weights,
                                                   final String alpha,
                                                   final double... logL) {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        for (final double value : logL) {
            mix.subLikelihoodsInput.get().add(new ConstantLikelihood(value));
        }
        mix.weightsInput.setValue(realParameter(weights), mix);
        if (alpha != null) {
            mix.alphaInput.setValue(new RealParameter(alpha), mix);
        }
        return mix;
    }

    private static MixtureTreeLikelihood typedMixture(final double[] weights,
                                                      final Double alpha,
                                                      final double... logL) {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        for (final double value : logL) {
            mix.subLikelihoodsInput.get().add(new ConstantLikelihood(value));
        }
        mix.weightsVectorInput.setValue(new RealVectorParam<>(weights, Real.INSTANCE), mix);
        if (alpha != null) {
            mix.alphaScalarInput.setValue(new RealScalarParam<>(alpha, Real.INSTANCE), mix);
        }
        mix.initAndValidate();
        return mix;
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

    private static String captureInit(final MixtureLikelihoodLogger logger) {
        return capture(logger::init);
    }

    private static String captureLog(final MixtureLikelihoodLogger logger) {
        return capture(out -> logger.log(0L, out));
    }

    private static String captureClose(final MixtureLikelihoodLogger logger) {
        return capture(logger::close);
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

    private static ExpectedMixture expected(final double[] weights, final double[] logL, final double alpha) {
        final double[] logS = new double[weights.length];
        double maxLogL = Double.NEGATIVE_INFINITY;
        double maxLogS = Double.NEGATIVE_INFINITY;
        double m = Double.NEGATIVE_INFINITY;
        double sumLogL = 0.0;
        for (int k = 0; k < weights.length; k++) {
            if (weights[k] > 0.0 && Double.isFinite(logL[k])) {
                maxLogL = Math.max(maxLogL, logL[k]);
                logS[k] = Math.log(weights[k]) + logL[k];
                maxLogS = Math.max(maxLogS, logS[k]);
                m = Math.max(m, logS[k]);
            } else {
                logS[k] = Double.NEGATIVE_INFINITY;
            }
            if (alpha != 0.0) {
                sumLogL = Double.isFinite(sumLogL) && Double.isFinite(logL[k])
                        ? sumLogL + logL[k]
                        : Double.NEGATIVE_INFINITY;
            }
        }

        double denomExp = 0.0;
        if (Double.isFinite(m)) {
            for (final double v : logS) {
                if (Double.isFinite(v)) {
                    denomExp += Math.exp(v - m);
                }
            }
        }
        final double logMix = denomExp > 0.0 ? m + Math.log(denomExp) : Double.NEGATIVE_INFINITY;
        final double[] resp = new double[weights.length];
        for (int k = 0; k < weights.length; k++) {
            resp[k] = denomExp > 0.0 && Double.isFinite(logS[k]) ? Math.exp(logS[k] - m) / denomExp : 0.0;
        }
        final double couplingTerm = alpha == 0.0 ? 0.0
                : (Double.isFinite(sumLogL) ? alpha * sumLogL : Double.NEGATIVE_INFINITY);
        final double totalLogP = Double.isFinite(logMix) && Double.isFinite(couplingTerm)
                ? logMix + couplingTerm
                : Double.NEGATIVE_INFINITY;
        return new ExpectedMixture(logMix, couplingTerm, totalLogP, logS, resp, maxLogL, maxLogS);
    }

    private interface StreamAction {
        void write(PrintStream out);
    }

    private static final class ExpectedMixture {
        private final double logMix;
        private final double couplingTerm;
        private final double totalLogP;
        private final double[] logS;
        private final double[] resp;
        private final double maxLogL;
        private final double maxLogS;

        private ExpectedMixture(final double logMix,
                                final double couplingTerm,
                                final double totalLogP,
                                final double[] logS,
                                final double[] resp,
                                final double maxLogL,
                                final double maxLogS) {
            this.logMix = logMix;
            this.couplingTerm = couplingTerm;
            this.totalLogP = totalLogP;
            this.logS = logS;
            this.resp = resp;
            this.maxLogL = maxLogL;
            this.maxLogS = maxLogS;
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
