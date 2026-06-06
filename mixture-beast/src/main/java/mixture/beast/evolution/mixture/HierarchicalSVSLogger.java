package mixture.beast.evolution.mixture;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class HierarchicalSVSLogger extends BEASTObject implements Loggable {

    public final Input<MixtureTreeLikelihood> topMixtureInput =
            new Input<>("topMixture", "Top-level mixture (strict vs relaxed)", Validate.REQUIRED);

    public final Input<RelaxedRatesPriorSVS> svsPriorInput =
            new Input<>("svsPrior", "SVS prior used inside relaxed clock (UC vs AC)", Validate.REQUIRED);

    public final Input<RealParameter> innerWeightsInput =
            new Input<>("innerWeights", "Inner weights for UC vs AC (dim=2). Default is 0.5 0.5.", (RealParameter) null);

    public final Input<RealVector> innerWeightsVectorInput =
            new Input<>("innerWeightsVector",
                    "BEAST3 typed inner weights for UC vs AC (dim=2). If absent, defaults to 0.5 0.5.",
                    Validate.OPTIONAL);

    public final Input<Boolean> printTopInput =
            new Input<>("printTop", "print top-level rStrict/rRelax and top weights/logL", false);

    public final Input<Boolean> printInnerInput =
            new Input<>("printInner", "print inner UC/AC logP and responsibilities", false);

    public final Input<Boolean> printThreeModelWeightsInput =
            new Input<>("printThreeModelWeights", "print pStrict, pRelaxUC, pRelaxAC", true);

    public final Input<Boolean> printRBStatsInput =
            new Input<>("printRBStats", "print Rao-Blackwell numerator/denominator columns for ucldStdev and sigma2", false);

    private List<GenericTreeLikelihood> topSubLiks;
    private RealParameter legacyTopWeights;
    private RealVector typedTopWeights;
    private RealParameter legacyInnerWeights;
    private RealVector typedInnerWeights;
    private RelaxedRatesPriorSVS svs;

    private boolean printTop;
    private boolean printInner;
    private boolean printThree;
    private boolean printRB;

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
        this.legacyTopWeights = top.weightsInput.get();
        this.typedTopWeights = top.weightsVectorInput.get();

        if (topSubLiks.size() != 2) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": currently expects K=2 at the top level.");
        }
        if (legacyTopWeights == null && typedTopWeights == null) {
            throw new IllegalArgumentException("HierarchicalSVSLogger: top mixture must provide weights or weightsVector.");
        }
        if (legacyTopWeights != null && typedTopWeights != null) {
            throw new IllegalArgumentException("HierarchicalSVSLogger: top mixture specifies both weights and weightsVector.");
        }
        if (topWeightsDimension() != 2) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": top weights must have dimension 2.");
        }

        this.svs = svsPriorInput.get();

        this.legacyInnerWeights = innerWeightsInput.get();
        this.typedInnerWeights = innerWeightsVectorInput.get();
        if (legacyInnerWeights != null && typedInnerWeights != null) {
            throw new IllegalArgumentException("HierarchicalSVSLogger: specify only one of innerWeights or innerWeightsVector.");
        }
        if (innerWeightsDimension() != 2) {
            throw new IllegalArgumentException(getClass().getSimpleName() + ": innerWeights must have dimension 2.");
        }

        this.printTop = printTopInput.get();
        this.printInner = printInnerInput.get();
        this.printThree = printThreeModelWeightsInput.get();
        this.printRB = printRBStatsInput.get();
    }

    @Override
    public void init(final PrintStream out) {
        final String p = (getID() == null) ? "" : (getID() + ".");

        if (printTop) {
            out.print(p + "topWStrict\t");
            out.print(p + "topWRelax\t");
            out.print(p + "topLogLStrict\t");
            out.print(p + "topLogLRelax\t");
            out.print(p + "rStrict\t");
            out.print(p + "rRelax\t");
        }

        if (printInner) {

            out.print(p + "logPUC\t");
            out.print(p + "logPAC\t");

        }

        if (printThree) {
            out.print(p + "pStrict\t");
            out.print(p + "pRelaxUC\t");
            out.print(p + "pRelaxAC\t");
        }

        if (printRB) {
            out.print(p + "rb_ucld_num\t");
            out.print(p + "rb_ucld_den\t");
            out.print(p + "rb_sigma2_num\t");
            out.print(p + "rb_sigma2_den\t");
            out.print(p + "rb_relax_ucld_num\t");
            out.print(p + "rb_relax_ucld_den\t");
            out.print(p + "rb_relax_sigma2_num\t");
            out.print(p + "rb_relax_sigma2_den\t");
        }
    }

    @Override
    public void log(final long sample, final PrintStream out) {
        final double wStrict = topWeightValue(0);
        final double wRelax = topWeightValue(1);

        final double logLStrict = topSubLiks.get(0).calculateLogP();
        final double logLRelax = topSubLiks.get(1).calculateLogP();

        final double sStrict = (wStrict > 0.0 && Double.isFinite(logLStrict)) ? (Math.log(wStrict) + logLStrict) : Double.NEGATIVE_INFINITY;
        final double sRelax  = (wRelax  > 0.0 && Double.isFinite(logLRelax))  ? (Math.log(wRelax)  + logLRelax)  : Double.NEGATIVE_INFINITY;

        final double m = Math.max(sStrict, sRelax);
        final double denom = (Double.isFinite(m)) ? (Math.exp(sStrict - m) + Math.exp(sRelax - m)) : 0.0;

        final double rStrict = (denom > 0.0 && Double.isFinite(sStrict)) ? (Math.exp(sStrict - m) / denom) : Double.NaN;
        final double rRelaxed = (Double.isFinite(rStrict)) ? (1.0 - rStrict) : Double.NaN;

        final double wu = innerWeightValue(0);
        final double wa = innerWeightValue(1);

        final double logPUC = svs.logPriorUCOnly();
        final double logPAC = svs.logPriorACOnly();

        final double tUC = (wu > 0.0 && Double.isFinite(logPUC)) ? (Math.log(wu) + logPUC) : Double.NEGATIVE_INFINITY;
        final double tAC = (wa > 0.0 && Double.isFinite(logPAC)) ? (Math.log(wa) + logPAC) : Double.NEGATIVE_INFINITY;

        final double mm = Math.max(tUC, tAC);
        final double denom2 = (Double.isFinite(mm)) ? (Math.exp(tUC - mm) + Math.exp(tAC - mm)) : 0.0;

        final double rUC = (denom2 > 0.0 && Double.isFinite(tUC)) ? (Math.exp(tUC - mm) / denom2) : Double.NaN;
        final double rAC = (Double.isFinite(rUC)) ? (1.0 - rUC) : Double.NaN;

        final double ucldStdev = svsUcldStdevValue();
        final double sigma2 = svsSigma2Value();

        final double pStrict = rStrict;
        final double pRelaxUC = (Double.isFinite(rRelaxed) && Double.isFinite(rUC)) ? (rRelaxed * rUC) : Double.NaN;
        final double pRelaxAC = (Double.isFinite(rRelaxed) && Double.isFinite(rAC)) ? (rRelaxed * rAC) : Double.NaN;

        if (printTop) {
            out.print(fmt(wStrict)); out.print('\t');
            out.print(fmt(wRelax)); out.print('\t');
            out.print(fmt(logLStrict)); out.print('\t');
            out.print(fmt(logLRelax)); out.print('\t');
            out.print(fmt(rStrict)); out.print('\t');
            out.print(fmt(rRelaxed)); out.print('\t');
        }

        if (printInner) {
            out.print(fmt(logPUC)); out.print('\t');
            out.print(fmt(logPAC)); out.print('\t');
        }

        if (printThree) {
            out.print(fmt(pStrict)); out.print('\t');
            out.print(fmt(pRelaxUC)); out.print('\t');
            out.print(fmt(pRelaxAC)); out.print('\t');
        }

        if (printRB) {
            final double rbU_num = (Double.isFinite(rUC)) ? (rUC * ucldStdev) : Double.NaN;
            final double rbU_den = rUC;
            final double rbA_num = (Double.isFinite(rAC)) ? (rAC * sigma2) : Double.NaN;
            final double rbA_den = rAC;

            final double rbRU_num = (Double.isFinite(pRelaxUC)) ? (pRelaxUC * ucldStdev) : Double.NaN;
            final double rbRU_den = pRelaxUC;
            final double rbRA_num = (Double.isFinite(pRelaxAC)) ? (pRelaxAC * sigma2) : Double.NaN;
            final double rbRA_den = pRelaxAC;

            out.print(fmt(rbU_num)); out.print('\t');
            out.print(fmt(rbU_den)); out.print('\t');
            out.print(fmt(rbA_num)); out.print('\t');
            out.print(fmt(rbA_den)); out.print('\t');

            out.print(fmt(rbRU_num)); out.print('\t');
            out.print(fmt(rbRU_den)); out.print('\t');
            out.print(fmt(rbRA_num)); out.print('\t');
            out.print(fmt(rbRA_den)); out.print('\t');
        }
    }

    @Override
    public void close(final PrintStream out) { }

    private int topWeightsDimension() {
        return legacyTopWeights != null ? legacyTopWeights.getDimension() : typedTopWeights.size();
    }

    private double topWeightValue(final int k) {
        return legacyTopWeights != null ? legacyTopWeights.getArrayValue(k) : typedTopWeights.get(k);
    }

    private boolean hasExplicitInnerWeights() {
        return legacyInnerWeights != null || typedInnerWeights != null;
    }

    private int innerWeightsDimension() {
        if (legacyInnerWeights != null) {
            return legacyInnerWeights.getDimension();
        }
        if (typedInnerWeights != null) {
            return typedInnerWeights.size();
        }
        return 2;
    }

    private double innerWeightValue(final int k) {
        if (legacyInnerWeights != null) {
            return legacyInnerWeights.getArrayValue(k);
        }
        if (typedInnerWeights != null) {
            return typedInnerWeights.get(k);
        }
        return 0.5;
    }

    private double svsUcldStdevValue() {
        final RealParameter legacy = svs.ucldStdevInput.get();
        if (legacy != null) {
            return legacy.getArrayValue(0);
        }
        final RealScalar typed = svs.ucldStdevScalarInput.get();
        if (typed != null) {
            return typed.get();
        }
        return Double.NaN;
    }

    private double svsSigma2Value() {
        final RealParameter legacy = svs.sigma2Input.get();
        if (legacy != null) {
            return legacy.getArrayValue(0);
        }
        final RealScalar typed = svs.sigma2ScalarInput.get();
        if (typed != null) {
            return typed.get();
        }
        return Double.NaN;
    }

    private static String fmt(final double x) {
        if (Double.isNaN(x)) return "NaN";
        if (x == Double.POSITIVE_INFINITY) return "Inf";
        if (x == Double.NEGATIVE_INFINITY) return "-Inf";
        return DF.get().format(x);
    }

    private static String fmt(final int x) {
        return Integer.toString(x);
    }
}
