package mixture.beast.evolution.mixture;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;

public class MixtureTreeLikelihood extends Distribution {

    public final Input<List<GenericTreeLikelihood>> subLikelihoodsInput = new Input<>(
            "subLikelihood",
            "Component likelihoods L_i = P(D|T,theta_i). Each can be TreeLikelihood or CompoundDistribution.",
            new ArrayList<>());

    public final Input<RealParameter> weightsInput = new Input<>(
            "weights",
            "Legacy mixture weights w (dimension K). Usually sum to 1.",
            Validate.OPTIONAL);

    public final Input<RealVector> weightsVectorInput = new Input<>(
            "weightsVector",
            "BEAST3 typed mixture weights w (dimension K). Usually sum to 1.",
            Validate.OPTIONAL);

    public final Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Legacy coupling exponent alpha (dimension 1). If provided, total logP = logMix + alpha * sum_i logL_i.",
            (RealParameter) null);

    public final Input<RealScalar> alphaScalarInput = new Input<>(
            "alphaScalar",
            "BEAST3 typed coupling exponent alpha. If provided, total logP = logMix + alpha * sum_i logL_i.",
            Validate.OPTIONAL);

    private List<GenericTreeLikelihood> subLikelihoods;
    private RealParameter legacyWeights;
    private RealVector typedWeights;
    private RealParameter legacyAlpha;
    private RealScalar typedAlpha;
    private int K;

    @Override
    public void initAndValidate() {
        subLikelihoods = subLikelihoodsInput.get();
        legacyWeights = weightsInput.get();
        typedWeights = weightsVectorInput.get();
        legacyAlpha = alphaInput.get();
        typedAlpha = alphaScalarInput.get();
        K = subLikelihoods.size();

        if (K < 2) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: need at least two subLikelihoods.");
        }
        if (legacyWeights == null && typedWeights == null) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: either weights or weightsVector must be specified.");
        }
        if (legacyWeights != null && typedWeights != null) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: specify only one of weights or weightsVector.");
        }
        if (legacyAlpha != null && typedAlpha != null) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: specify only one of alpha or alphaScalar.");
        }
        if (weightsDimension() != K) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: weights dimension ("
                    + weightsDimension() + ") != number of subLikelihoods (" + K + ").");
        }
        if (legacyAlpha != null && legacyAlpha.getDimension() != 1) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: alpha must have dimension 1.");
        }

        double wsum = 0.0;
        for (int i = 0; i < K; i++) {
            double wi = weightValue(i);
            if (wi < 0.0) {
                throw new IllegalArgumentException("MixtureTreeLikelihood: negative weight w[" + i + "]=" + wi);
            }
            wsum += wi;
        }
        if (Math.abs(wsum - 1.0) > 1e-6) {
            Log.warning("MixtureTreeLikelihood: weights do not sum to 1 (sum=" + wsum
                    + "). If you expect a true mixture, make sure sum(weights)=1.");
        }
    }

    @Override
    public double calculateLogP() {

        final double a = alphaValue();

        double maxTerm = Double.NEGATIVE_INFINITY;
        double sumLogL = 0.0;
        boolean needCoupling = (a > 0.0);

        final double[] term = new double[K];

        for (int i = 0; i < K; i++) {

            final double wi = weightValue(i);
            final boolean needThisLi = needCoupling || (wi > 0.0);

            if (!needThisLi) {
                term[i] = Double.NEGATIVE_INFINITY;
                continue;
            }

            final double li = subLikelihoods.get(i).calculateLogP();


            if (needCoupling) {
                if (!Double.isFinite(li)) {
                    logP = Double.NEGATIVE_INFINITY;
                    Log.warning("MixtureTreeLikelihood: component " + i + " logL is -Inf/NaN while alpha>0 -> total -Inf");
                    return logP;
                }
                sumLogL += li;
            }

            // mixture term uses only wi>0
            if (wi <= 0.0 || !Double.isFinite(li)) {
                term[i] = Double.NEGATIVE_INFINITY;
                continue;
            }

            final double ti = Math.log(wi) + li;
            term[i] = ti;
            if (ti > maxTerm) {
                maxTerm = ti;
            }
        }

        if (!Double.isFinite(maxTerm)) {
            logP = Double.NEGATIVE_INFINITY;
            Log.warning("MixtureTreeLikelihood: all positive-weight components are -Inf/NaN -> returning -Inf");
            return logP;
        }

        // log-sum-exp
        double sumExp = 0.0;
        for (int i = 0; i < K; i++) {
            if (Double.isFinite(term[i])) {
                sumExp += Math.exp(term[i] - maxTerm);
            }
        }
        if (!(sumExp > 0.0)) {
            logP = Double.NEGATIVE_INFINITY;
            Log.warning("MixtureTreeLikelihood: terms underflowed/NaN -> returning -Inf");
            return logP;
        }

        final double logMix = Math.log(sumExp) + maxTerm;
        final double logCouple = (a == 0.0 ? 0.0 : a * sumLogL);

        logP = logMix + logCouple;
        return logP;
    }

    @Override
    public List<String> getArguments() { return List.of(); }

    @Override
    public List<String> getConditions() { return List.of(); }

    @Override
    public void sample(State state, Random random) { }

    @Override
    protected boolean requiresRecalculation() { return true; }

    private int weightsDimension() {
        return legacyWeights != null ? legacyWeights.getDimension() : typedWeights.size();
    }

    private double weightValue(final int i) {
        return legacyWeights != null ? legacyWeights.getArrayValue(i) : typedWeights.get(i);
    }

    private double alphaValue() {
        if (legacyAlpha != null) {
            return legacyAlpha.getArrayValue(0);
        }
        if (typedAlpha != null) {
            return typedAlpha.get();
        }
        return 0.0;
    }
}
