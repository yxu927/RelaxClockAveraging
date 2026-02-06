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

public class HierarchicalSVSLogger extends BEASTObject implements Loggable {

    public final Input<MixtureTreeLikelihood> topMixtureInput =
            new Input<>("topMixture", "Top-level mixture (strict vs relaxed).", Validate.REQUIRED);

    public final Input<RelaxedRatesPriorSVS> svsPriorInput =
            new Input<>("svsPrior", "SVS prior inside relaxed clock (UC vs AC).", Validate.REQUIRED);

    public final Input<RealParameter> innerWeightsInput =
            new Input<>("innerWeights", "Inner weights for UC vs AC (dim=2). Default is 0.5 0.5.", (RealParameter) null);

    public final Input<Long> startSampleForSummaryInput =
            new Input<>("startSampleForSummary",
                    "Only accumulate conditional summaries when sample >= this value (use to exclude burn-in).",
                    0L);

    public final Input<Boolean> printThreeModelWeightsInput =
            new Input<>("printThreeModelWeights", "Print pStrict, pRelaxUC, pRelaxAC.", true);

    public final Input<Boolean> printConditionalMeansInput =
            new Input<>("printConditionalMeans",
                    "Print conditional means E[ucldStdev|relaxed&UC] and E[sigma2|relaxed&AC] estimated by ratio-of-sums.",
                    true);

    public final Input<Boolean> printInnerLogPInput =
            new Input<>("printInnerLogP", "Print logPUC and logPAC.", false);

    public final Input<Boolean> printInnerRespInput =
            new Input<>("printInnerResp", "Print rUC and rAC.", false);

    public final Input<Boolean> printFinalSummaryToStderrInput =
            new Input<>("printFinalSummaryToStderr", "Print final conditional summaries to stderr in close().", true);

    private List<GenericTreeLikelihood> topSubLiks;
    private RealParameter topWeights;
    private RealParameter innerWeights;
    private RelaxedRatesPriorSVS svs;

    private long startSampleForSummary;

    private boolean printThree;
    private boolean printCondMeans;
    private boolean printInnerLogP;
    private boolean printInnerResp;
    private boolean printFinalSummaryToStderr;

    // Accumulators for conditional means (after burn-in threshold).
    // E[ucldStdev | relaxed&UC] ≈ sum(pRelaxUC * ucldStdev) / sum(pRelaxUC)
    // E[sigma2    | relaxed&AC] ≈ sum(pRelaxAC * sigma2)    / sum(pRelaxAC)
    private double sumRelaxUC = 0.0;
    private double sumRelaxAC = 0.0;
    private double sumRelaxUC_ucld = 0.0;
    private double sumRelaxAC_sigma2 = 0.0;

    private static final ThreadLocal<DecimalFormat> DF = ThreadLocal.withInitial(() -> {
        DecimalFormatSymbols sym = DecimalFormatSymbols.getInstance(Locale.US);
        DecimalFormat df = new DecimalFormat("0.000000", sym);
        df.setGroupingUsed(false);
        return df;
    });

    @Override
    public void initAndValidate() {
        MixtureTreeLikelihood top = topMixtureInput.get();
        this.topSubLiks = new ArrayList<>(top.subLikelihoodsInput.get());
        this.topWeights = top.weightsInput.get();

        if (topSubLiks.size() != 2) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": expects K=2 at top level.");
        }
        if (topWeights.getDimension() != 2) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": top weights must have dimension 2.");
        }

        this.svs = svsPriorInput.get();

        this.innerWeights = innerWeightsInput.get();
        if (innerWeights == null) {
            innerWeights = new RealParameter("0.5 0.5");
        }
        if (innerWeights.getDimension() != 2) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": innerWeights must have dimension 2.");
        }

        this.startSampleForSummary = startSampleForSummaryInput.get();

        this.printThree = printThreeModelWeightsInput.get();
        this.printCondMeans = printConditionalMeansInput.get();
        this.printInnerLogP = printInnerLogPInput.get();
        this.printInnerResp = printInnerRespInput.get();
        this.printFinalSummaryToStderr = printFinalSummaryToStderrInput.get();
    }

    @Override
    public void init(final PrintStream out) {
        final String p = (getID() == null) ? "" : (getID() + ".");

        if (printInnerLogP) {
            out.print(p + "logPUC\t");
            out.print(p + "logPAC\t");
        }
        if (printInnerResp) {
            out.print(p + "rUC\t");
            out.print(p + "rAC\t");
        }
        if (printThree) {
            out.print(p + "pStrict\t");
            out.print(p + "pRelaxUC\t");
            out.print(p + "pRelaxAC\t");
        }
        if (printCondMeans) {
            out.print(p + "E_ucldStdev_relaxUC\t");
            out.print(p + "E_sigma2_relaxAC\t");
        }
    }

    @Override
    public void log(final long sample, final PrintStream out) {

        // ----- Top level responsibilities (strict vs relaxed) -----
        final double wStrict = topWeights.getArrayValue(0);
        final double wRelax  = topWeights.getArrayValue(1);

        final double logLStrict = topSubLiks.get(0).calculateLogP();
        final double logLRelax  = topSubLiks.get(1).calculateLogP();

        final double sStrict = (wStrict > 0.0 && Double.isFinite(logLStrict))
                ? (Math.log(wStrict) + logLStrict) : Double.NEGATIVE_INFINITY;

        final double sRelax = (wRelax > 0.0 && Double.isFinite(logLRelax))
                ? (Math.log(wRelax) + logLRelax) : Double.NEGATIVE_INFINITY;

        final double mTop = Math.max(sStrict, sRelax);
        final double denomTop = (Double.isFinite(mTop))
                ? (Math.exp(sStrict - mTop) + Math.exp(sRelax - mTop)) : 0.0;

        final double rStrict = (denomTop > 0.0 && Double.isFinite(sStrict))
                ? (Math.exp(sStrict - mTop) / denomTop) : Double.NaN;

        final double rRelaxed = Double.isFinite(rStrict) ? (1.0 - rStrict) : Double.NaN;

        // ----- Inner responsibilities (UC vs AC) based on SVS prior densities -----
        final double wu = innerWeights.getArrayValue(0);
        final double wa = innerWeights.getArrayValue(1);

        final double logPUC = svs.logPriorUCOnly();
        final double logPAC = svs.logPriorACOnly();

        final double tUC = (wu > 0.0 && Double.isFinite(logPUC))
                ? (Math.log(wu) + logPUC) : Double.NEGATIVE_INFINITY;

        final double tAC = (wa > 0.0 && Double.isFinite(logPAC))
                ? (Math.log(wa) + logPAC) : Double.NEGATIVE_INFINITY;

        final double mInner = Math.max(tUC, tAC);
        final double denomInner = (Double.isFinite(mInner))
                ? (Math.exp(tUC - mInner) + Math.exp(tAC - mInner)) : 0.0;

        final double rUC = (denomInner > 0.0 && Double.isFinite(tUC))
                ? (Math.exp(tUC - mInner) / denomInner) : Double.NaN;

        final double rAC = Double.isFinite(rUC) ? (1.0 - rUC) : Double.NaN;

        // ----- Three-model weights (Rao-Blackwell over both latent indicators) -----
        final double pStrict  = rStrict;
        final double pRelaxUC = (Double.isFinite(rRelaxed) && Double.isFinite(rUC)) ? (rRelaxed * rUC) : Double.NaN;
        final double pRelaxAC = (Double.isFinite(rRelaxed) && Double.isFinite(rAC)) ? (rRelaxed * rAC) : Double.NaN;

        // ----- Parameters to summarize -----
        final double ucldStdev = svs.ucldStdevInput.get().getValue();
        final double sigma2 = svs.sigma2Input.get().getValue();

        // ----- Accumulate conditional means after burn-in threshold -----
        if (sample >= startSampleForSummary) {
            if (Double.isFinite(pRelaxUC) && Double.isFinite(ucldStdev) && pRelaxUC > 0.0) {
                sumRelaxUC += pRelaxUC;
                sumRelaxUC_ucld += pRelaxUC * ucldStdev;
            }
            if (Double.isFinite(pRelaxAC) && Double.isFinite(sigma2) && pRelaxAC > 0.0) {
                sumRelaxAC += pRelaxAC;
                sumRelaxAC_sigma2 += pRelaxAC * sigma2;
            }
        }

        final double estUcld = (sumRelaxUC > 0.0) ? (sumRelaxUC_ucld / sumRelaxUC) : Double.NaN;
        final double estSigma2 = (sumRelaxAC > 0.0) ? (sumRelaxAC_sigma2 / sumRelaxAC) : Double.NaN;

        // ----- Output -----
        if (printInnerLogP) {
            out.print(fmt(logPUC)); out.print('\t');
            out.print(fmt(logPAC)); out.print('\t');
        }
        if (printInnerResp) {
            out.print(fmt(rUC)); out.print('\t');
            out.print(fmt(rAC)); out.print('\t');
        }
        if (printThree) {
            out.print(fmt(pStrict)); out.print('\t');
            out.print(fmt(pRelaxUC)); out.print('\t');
            out.print(fmt(pRelaxAC)); out.print('\t');
        }
        if (printCondMeans) {
            out.print(fmt(estUcld)); out.print('\t');
            out.print(fmt(estSigma2)); out.print('\t');
        }
    }

    @Override
    public void close(final PrintStream out) {
        if (printFinalSummaryToStderr) {
            final double estUcld = (sumRelaxUC > 0.0) ? (sumRelaxUC_ucld / sumRelaxUC) : Double.NaN;
            final double estSigma2 = (sumRelaxAC > 0.0) ? (sumRelaxAC_sigma2 / sumRelaxAC) : Double.NaN;

            System.err.println("[HierarchicalSVSLogger] startSampleForSummary=" + startSampleForSummary);
            System.err.println("[HierarchicalSVSLogger] E[ucldStdev | relaxed&UC] ≈ " + estUcld
                    + "  (den=" + sumRelaxUC + ")");
            System.err.println("[HierarchicalSVSLogger] E[sigma2 | relaxed&AC] ≈ " + estSigma2
                    + "  (den=" + sumRelaxAC + ")");
        }
    }

    private static String fmt(final double x) {
        if (Double.isNaN(x)) return "NaN";
        if (x == Double.POSITIVE_INFINITY) return "Inf";
        if (x == Double.NEGATIVE_INFINITY) return "-Inf";
        return DF.get().format(x);
    }
}
