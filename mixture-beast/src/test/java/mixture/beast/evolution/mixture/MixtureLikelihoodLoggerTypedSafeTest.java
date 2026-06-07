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
import static org.junit.Assert.assertThrows;

public class MixtureLikelihoodLoggerTypedSafeTest {

    @Test
    public void typedMixtureDefaultLoggerMatchesLegacyResponsibilities() {
        final MixtureTreeLikelihood legacy = legacyMixture(new double[]{0.25, 0.75}, null, -10.0, -12.0);
        final MixtureTreeLikelihood typed = typedMixture(new double[]{0.25, 0.75}, null, -10.0, -12.0);

        assertLoggerOutputsEqual(logger(legacy, "mixLog"), logger(typed, "mixLog"));
    }

    @Test
    public void typedMixtureAllPrintFlagsMatchLegacyWithoutAlpha() {
        final MixtureTreeLikelihood legacy = legacyMixture(new double[]{0.25, 0.75}, null, -10.0, -12.0);
        final MixtureTreeLikelihood typed = typedMixture(new double[]{0.25, 0.75}, null, -10.0, -12.0);

        assertLoggerOutputsEqual(loggerWithAllFlags(legacy, "mixLog"), loggerWithAllFlags(typed, "mixLog"));
    }

    @Test
    public void typedMixtureAllPrintFlagsMatchLegacyWithPositiveAlpha() {
        final MixtureTreeLikelihood legacy = legacyMixture(new double[]{0.25, 0.75}, 0.2, -10.0, -12.0);
        final MixtureTreeLikelihood typed = typedMixture(new double[]{0.25, 0.75}, 0.2, -10.0, -12.0);

        assertLoggerOutputsEqual(loggerWithAllFlags(legacy, "mixLog"), loggerWithAllFlags(typed, "mixLog"));
    }

    @Test
    public void typedZeroWeightAndNonFiniteHandlingMatchesLegacy() {
        final MixtureTreeLikelihood legacy = legacyMixture(new double[]{0.0, 1.0}, null, Double.NaN, -3.0);
        final MixtureTreeLikelihood typed = typedMixture(new double[]{0.0, 1.0}, null, Double.NaN, -3.0);

        assertLoggerOutputsEqual(logger(legacy, "mixLog"), logger(typed, "mixLog"));
    }

    @Test
    public void typedMixtureLoggerRejectsMissingWeights() {
        final MixtureTreeLikelihood mix = rawMixture(-1.0, -2.0);
        final MixtureLikelihoodLogger logger = logger(mix, "mixLog");

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
    }

    @Test
    public void typedMixtureLoggerRejectsBothLegacyAndTypedWeights() {
        final MixtureTreeLikelihood mix = rawMixture(-1.0, -2.0);
        mix.weightsInput.setValue(realParameter(new double[]{0.25, 0.75}), mix);
        mix.weightsVectorInput.setValue(new RealVectorParam<>(new double[]{0.25, 0.75}, Real.INSTANCE), mix);
        final MixtureLikelihoodLogger logger = logger(mix, "mixLog");

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
    }

    @Test
    public void typedMixtureLoggerRejectsBothLegacyAndTypedAlpha() {
        final MixtureTreeLikelihood mix = rawMixture(-1.0, -2.0);
        mix.weightsVectorInput.setValue(new RealVectorParam<>(new double[]{0.25, 0.75}, Real.INSTANCE), mix);
        mix.alphaInput.setValue(new RealParameter("0.2"), mix);
        mix.alphaScalarInput.setValue(new RealScalarParam<>(0.2, Real.INSTANCE), mix);
        final MixtureLikelihoodLogger logger = logger(mix, "mixLog");

        assertThrows(IllegalArgumentException.class, logger::initAndValidate);
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
        final MixtureTreeLikelihood mix = rawMixture(logL);
        mix.weightsInput.setValue(realParameter(weights), mix);
        if (alpha != null) {
            mix.alphaInput.setValue(new RealParameter(Double.toString(alpha)), mix);
        }
        mix.initAndValidate();
        return mix;
    }

    private static MixtureTreeLikelihood typedMixture(final double[] weights,
                                                      final Double alpha,
                                                      final double... logL) {
        final MixtureTreeLikelihood mix = rawMixture(logL);
        mix.weightsVectorInput.setValue(new RealVectorParam<>(weights, Real.INSTANCE), mix);
        if (alpha != null) {
            mix.alphaScalarInput.setValue(new RealScalarParam<>(alpha, Real.INSTANCE), mix);
        }
        mix.initAndValidate();
        return mix;
    }

    private static MixtureTreeLikelihood rawMixture(final double... logL) {
        final MixtureTreeLikelihood mix = new MixtureTreeLikelihood();
        for (final double value : logL) {
            mix.subLikelihoodsInput.get().add(new ConstantLikelihood(value));
        }
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

    private static void assertLoggerOutputsEqual(final MixtureLikelihoodLogger legacyLogger,
                                                 final MixtureLikelihoodLogger typedLogger) {
        legacyLogger.initAndValidate();
        typedLogger.initAndValidate();
        assertArrayEquals(tokens(captureInit(legacyLogger)), tokens(captureInit(typedLogger)));
        assertArrayEquals(tokens(captureLog(legacyLogger)), tokens(captureLog(typedLogger)));
    }

    private static String captureInit(final MixtureLikelihoodLogger logger) {
        return capture(logger::init);
    }

    private static String captureLog(final MixtureLikelihoodLogger logger) {
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
