package mixture.beast.evolution.mixture;

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

    public final Input<Boolean> printTotalLogPInput =
            new Input<>("printTotalLogP", "print total logP (logMix + couplingTerm)", false);

    public final Input<Boolean> printLogMixInput =
            new Input<>("printLogMix", "print logMix only (log-sum-exp of w_k * exp(logL_k))", false);

    public final Input<Boolean> printAlphaInput =
            new Input<>("printAlpha", "print alpha (coupling exponent)", false);

    public final Input<Boolean> printCouplingTermInput =
            new Input<>("printCouplingTerm", "print couplingTerm = alpha * sum_k logL_k", false);

    public final Input<Boolean> printSumWeightsInput =
            new Input<>("printSumWeights", "print sum_k w_k", false);

    public final Input<Boolean> printWeightsInput =
            new Input<>("printWeights", "print mixture weights w_k", false);

    public final Input<Boolean> printLogLInput =
            new Input<>("printLogL", "print component log-likelihoods logL_k", false);

    public final Input<Boolean> printLogSInput =
            new Input<>("printLogS", "print logS_k = log(w_k) + logL_k", false);

    public final Input<Boolean> printRespInput =
            new Input<>("printResp", "print responsibilities r_k (from logS)", true);

    public final Input<Boolean> printMaxLogLInput =
            new Input<>("printMaxLogL", "print max_k logL_k over positive-weight components", false);

    public final Input<Boolean> printMaxLogSInput =
            new Input<>("printMaxLogS", "print max_k logS_k over positive-weight components", false);

    public final Input<Boolean> printMixMinusMaxLogLInput =
            new Input<>("printMixMinusMaxLogL", "print (logMix - maxLogL); should be <= 0 for a proper mixture", false);

    public final Input<Boolean> printMixMinusMaxLogSInput =
            new Input<>("printMixMinusMaxLogS", "print (logMix - maxLogS); should be >= 0 for a proper mixture", false);

    public final Input<Boolean> printViolationFlagInput =
            new Input<>("printViolationFlag", "print 1 if logMix > maxLogL + 1e-10, else 0", false);

    private int K;
    private List<GenericTreeLikelihood> subLiks;
    private RealParameter weights;
    private RealParameter alpha;

    private boolean printTotalLogP;
    private boolean printLogMix;
    private boolean printAlpha;
    private boolean printCouplingTerm;
    private boolean printSumWeights;
    private boolean printWeights;
    private boolean printLogL;
    private boolean printLogS;
    private boolean printResp;
    private boolean printMaxLogL;
    private boolean printMaxLogS;
    private boolean printMixMinusMaxLogL;
    private boolean printMixMinusMaxLogS;
    private boolean printViolationFlag;

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
        this.alpha = mix.alphaInput.get(); // may be null
        this.K = subLiks.size();

        if (K < 2) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": need K >= 2");
        }
        if (weights.getDimension() != K) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": weights dimension != K");
        }
        if (alpha != null && alpha.getDimension() != 1) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": alpha must have dimension=1");
        }

        this.printTotalLogP = printTotalLogPInput.get();
        this.printLogMix = printLogMixInput.get();
        this.printAlpha = printAlphaInput.get();
        this.printCouplingTerm = printCouplingTermInput.get();
        this.printSumWeights = printSumWeightsInput.get();
        this.printWeights = printWeightsInput.get();
        this.printLogL = printLogLInput.get();
        this.printLogS = printLogSInput.get();
        this.printResp = printRespInput.get();
        this.printMaxLogL = printMaxLogLInput.get();
        this.printMaxLogS = printMaxLogSInput.get();
        this.printMixMinusMaxLogL = printMixMinusMaxLogLInput.get();
        this.printMixMinusMaxLogS = printMixMinusMaxLogSInput.get();
        this.printViolationFlag = printViolationFlagInput.get();

        this.logL = new double[K];
        this.wk = new double[K];
        this.logS = new double[K];
    }

    @Override
    public void init(final PrintStream out) {
        final String prefix = (getID() == null) ? "" : (getID() + ".");

        if (printTotalLogP) out.print(prefix + "totalLogP\t");
        if (printLogMix) out.print(prefix + "logMix\t");
        if (printAlpha) out.print(prefix + "alpha\t");
        if (printCouplingTerm) out.print(prefix + "couplingTerm\t");
        if (printSumWeights) out.print(prefix + "sumW\t");

        if (printWeights) {
            for (int k = 0; k < K; k++) out.print(prefix + "w[" + (k + 1) + "]\t");
        }
        if (printLogL) {
            for (int k = 0; k < K; k++) out.print(prefix + "logL[" + (k + 1) + "]\t");
        }
        if (printLogS) {
            for (int k = 0; k < K; k++) out.print(prefix + "logS[" + (k + 1) + "]\t");
        }
        if (printResp) {
            for (int k = 0; k < K; k++) out.print(prefix + "r[" + (k + 1) + "]\t");
        }

        if (printMaxLogL) out.print(prefix + "maxLogL\t");
        if (printMaxLogS) out.print(prefix + "maxLogS\t");
        if (printMixMinusMaxLogL) out.print(prefix + "logMixMinusMaxLogL\t");
        if (printMixMinusMaxLogS) out.print(prefix + "logMixMinusMaxLogS\t");
        if (printViolationFlag) out.print(prefix + "violationMixGTMaxLogL\t");
    }

    @Override
    public void log(final long sample, final PrintStream out) {
        final double a = (alpha == null ? 0.0 : alpha.getArrayValue(0));

        final boolean needWeights = printWeights || printTotalLogP || printLogMix || printResp || printLogS
                || printSumWeights || printMaxLogS || printMixMinusMaxLogS;
        final boolean needLogL = printLogL || printTotalLogP || printLogMix || printResp || printLogS
                || printMaxLogL || printMaxLogS || printMixMinusMaxLogL || printMixMinusMaxLogS || printViolationFlag
                || (a != 0.0) || printCouplingTerm;

        if (needWeights) {
            for (int k = 0; k < K; k++) wk[k] = weights.getArrayValue(k);
        }

        if (needLogL) {
            for (int k = 0; k < K; k++) logL[k] = subLiks.get(k).calculateLogP();
        }

        double sumW = 0.0;
        if (printSumWeights) {
            for (int k = 0; k < K; k++) sumW += wk[k];
        }

        double maxLogL = Double.NEGATIVE_INFINITY;
        double maxLogS = Double.NEGATIVE_INFINITY;

        double m = Double.NEGATIVE_INFINITY;
        double denomExp = 0.0;

        double sumLogL = 0.0;
        if (a != 0.0) {
            for (int k = 0; k < K; k++) {
                final double li = logL[k];
                if (!Double.isFinite(li)) {
                    sumLogL = Double.NEGATIVE_INFINITY;
                    break;
                }
                sumLogL += li;
            }
        }

        for (int k = 0; k < K; k++) {
            final double w = wk[k];
            final double li = logL[k];

            if (w > 0.0 && Double.isFinite(li)) {
                if (li > maxLogL) maxLogL = li;
                final double v = Math.log(w) + li;
                logS[k] = v;
                if (v > maxLogS) maxLogS = v;
                if (v > m) m = v;
            } else {
                logS[k] = Double.NEGATIVE_INFINITY;
            }
        }

        double logMix = Double.NEGATIVE_INFINITY;
        if (Double.isFinite(m)) {
            for (int k = 0; k < K; k++) {
                final double v = logS[k];
                if (Double.isFinite(v)) denomExp += Math.exp(v - m);
            }
            if (denomExp > 0.0) logMix = m + Math.log(denomExp);
        }

        final double couplingTerm = (a == 0.0 ? 0.0 : (Double.isFinite(sumLogL) ? a * sumLogL : Double.NEGATIVE_INFINITY));
        final double totalLogP = (Double.isFinite(logMix) && Double.isFinite(couplingTerm)) ? (logMix + couplingTerm) : Double.NEGATIVE_INFINITY;

        if (printTotalLogP) { out.print(fmt(totalLogP)); out.print('\t'); }
        if (printLogMix) { out.print(fmt(logMix)); out.print('\t'); }
        if (printAlpha) { out.print(fmt(a)); out.print('\t'); }
        if (printCouplingTerm) { out.print(fmt(couplingTerm)); out.print('\t'); }
        if (printSumWeights) { out.print(fmt(sumW)); out.print('\t'); }

        if (printWeights) {
            for (int k = 0; k < K; k++) { out.print(fmt(wk[k])); out.print('\t'); }
        }
        if (printLogL) {
            for (int k = 0; k < K; k++) { out.print(fmt(logL[k])); out.print('\t'); }
        }
        if (printLogS) {
            for (int k = 0; k < K; k++) { out.print(fmt(logS[k])); out.print('\t'); }
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
                for (int k = 0; k < K; k++) out.print("NaN\t");
            }
        }

        if (printMaxLogL) { out.print(fmt(maxLogL)); out.print('\t'); }
        if (printMaxLogS) { out.print(fmt(maxLogS)); out.print('\t'); }

        if (printMixMinusMaxLogL) {
            out.print(fmt(Double.isFinite(logMix) && Double.isFinite(maxLogL) ? (logMix - maxLogL) : Double.NaN));
            out.print('\t');
        }
        if (printMixMinusMaxLogS) {
            out.print(fmt(Double.isFinite(logMix) && Double.isFinite(maxLogS) ? (logMix - maxLogS) : Double.NaN));
            out.print('\t');
        }
        if (printViolationFlag) {
            final double eps = 1e-10;
            final double flag = (Double.isFinite(logMix) && Double.isFinite(maxLogL) && (logMix > maxLogL + eps)) ? 1.0 : 0.0;
            out.print(fmt(flag));
            out.print('\t');
        }
    }

    @Override
    public void close(final PrintStream out) { }

    private static String fmt(final double x) {
        if (Double.isNaN(x)) return "NaN";
        if (x == Double.POSITIVE_INFINITY) return "Inf";
        if (x == Double.NEGATIVE_INFINITY) return "-Inf";
        return DF.get().format(x);
    }
}
