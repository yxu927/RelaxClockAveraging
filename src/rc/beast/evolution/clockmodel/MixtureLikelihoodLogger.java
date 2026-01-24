package rc.beast.evolution.clockmodel;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.parameter.RealParameter;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class MixtureLikelihoodLogger extends BEASTObject implements Loggable {

    public final Input<MixtureTreeLikelihood> mixtureInput =
            new Input<>("mixture", "MixtureTreeLikelihood to log", Validate.REQUIRED);

    /** Whether to print the total mixture log-likelihood column. */
    public final Input<Boolean> printMixLogPInput =
            new Input<>("printMixLogP", "print mixture log-likelihood mixLogP", true);

    /** Whether to print mixture weights w_k. */
    public final Input<Boolean> printWeightsInput =
            new Input<>("printWeights", "print mixture weights w_k", true);

    /** Whether to print component log-likelihoods logL_k. */
    public final Input<Boolean> printLogLInput =
            new Input<>("printLogL", "print component log-likelihoods logL_k", true);

    /** Whether to print posterior responsibilities r_k. */
    public final Input<Boolean> printRespInput =
            new Input<>("printResp", "print posterior responsibilities r_k", true);

    private int K;
    private List<GenericTreeLikelihood> subLiks;
    private RealParameter weights;

    // Cached flags (avoid Input.get() calls in the hot path)
    private boolean printMixLogP;
    private boolean printWeights;
    private boolean printLogL;
    private boolean printResp;

    // Reused buffers (avoid allocating every log() call)
    private double[] logL;
    private double[] wk;
    private double[] logS;

    private static final ThreadLocal<DecimalFormat> DF = ThreadLocal.withInitial(() -> {
        DecimalFormatSymbols sym = DecimalFormatSymbols.getInstance(Locale.US);
        DecimalFormat df = new DecimalFormat("0.000000", sym);
        df.setGroupingUsed(false);
        return df;
    });

    @Override
    public void initAndValidate() {
        final MixtureTreeLikelihood mix = mixtureInput.get();

        this.subLiks = new ArrayList<>(mix.subLikelihoodsInput.get());
        this.weights = mix.weightsInput.get();
        this.K = subLiks.size();

        if (K < 2) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": need K >= 2");
        }
        if (weights.getDimension() != K) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": weights dimension != K");
        }

        this.printMixLogP = printMixLogPInput.get();
        this.printWeights = printWeightsInput.get();
        this.printLogL = printLogLInput.get();
        this.printResp = printRespInput.get();

        this.logL = new double[K];
        this.wk = new double[K];
        this.logS = new double[K];
    }

    @Override
    public void init(final PrintStream out) {
        final String prefix = (getID() == null) ? "" : (getID() + ".");

        if (printMixLogP) {
            out.print(prefix + "mixLogP\t");
        }
        if (printWeights) {
            for (int k = 0; k < K; k++) {
                out.print(prefix + "w[" + (k + 1) + "]\t");
            }
        }
        if (printLogL) {
            for (int k = 0; k < K; k++) {
                out.print(prefix + "logL[" + (k + 1) + "]\t");
            }
        }
        if (printResp) {
            for (int k = 0; k < K; k++) {
                out.print(prefix + "r[" + (k + 1) + "]\t");
            }
        }
    }

    @Override
    public void log(final long sample, final PrintStream out) {
        // Only compute what is needed (big speed-up if you disable logL/resp/mixLogP).
        final boolean needWeights = printWeights || printMixLogP || printResp;
        final boolean needLogL    = printLogL || printMixLogP || printResp;

        if (needWeights) {
            for (int k = 0; k < K; k++) {
                wk[k] = weights.getArrayValue(k);
            }
        }

        if (needLogL) {
            for (int k = 0; k < K; k++) {

                logL[k] = subLiks.get(k).calculateLogP();
            }
        }

        double m = Double.NEGATIVE_INFINITY;
        double denomExp = 0.0;
        double mixLogP = Double.NEGATIVE_INFINITY;

        if (printMixLogP || printResp) {
            // logS[k] = log(wk[k]) + logL[k] for valid components, else -Inf.
            for (int k = 0; k < K; k++) {
                final double w = wk[k];
                final double li = logL[k];

                if (w > 0.0 && Double.isFinite(li)) {
                    final double v = Math.log(w) + li;
                    logS[k] = v;
                    if (v > m) m = v;
                } else {
                    logS[k] = Double.NEGATIVE_INFINITY;
                }
            }

            if (Double.isFinite(m)) {
                for (int k = 0; k < K; k++) {
                    final double v = logS[k];
                    if (Double.isFinite(v)) {
                        denomExp += Math.exp(v - m);
                    }
                }
                if (denomExp > 0.0) {
                    mixLogP = m + Math.log(denomExp);
                }
            }
        }

        if (printMixLogP) {
            out.print(fmt(mixLogP));
            out.print('\t');
        }

        if (printWeights) {
            for (int k = 0; k < K; k++) {
                out.print(fmt(wk[k]));
                out.print('\t');
            }
        }

        if (printLogL) {
            for (int k = 0; k < K; k++) {
                out.print(fmt(logL[k]));
                out.print('\t');
            }
        }

        if (printResp) {
            if (denomExp > 0.0 && Double.isFinite(m)) {
                for (int k = 0; k < K; k++) {
                    final double v = logS[k];
                    final double rk = Double.isFinite(v) ? (Math.exp(v - m) / denomExp) : 0.0;
                    out.print(fmt(rk));
                    out.print('\t');
                }
            } else {
                for (int k = 0; k < K; k++) {
                    out.print("NaN\t");
                }
            }
        }
    }

    @Override
    public void close(final PrintStream out) {
    }

    private static String fmt(final double x) {
        if (Double.isNaN(x)) return "NaN";
        if (x == Double.POSITIVE_INFINITY) return "Inf";
        if (x == Double.NEGATIVE_INFINITY) return "-Inf";
        return DF.get().format(x);
    }
}
