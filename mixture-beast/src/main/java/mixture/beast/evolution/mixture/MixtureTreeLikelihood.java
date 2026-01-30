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

public class MixtureTreeLikelihood extends Distribution {

    public final Input<List<GenericTreeLikelihood>> subLikelihoodsInput = new Input<>(
            "subLikelihood",
            "Component likelihoods L_i = P(D|T,theta_i). Each can be TreeLikelihood or CompoundDistribution.",
            new ArrayList<>());

    public final Input<RealParameter> weightsInput = new Input<>(
            "weights",
            "Mixture weights w (dimension K). Usually sum to 1.",
            Validate.REQUIRED);

    // alpha is OPTIONAL: if not set -> behave like pure mixture (alpha=0)
    public final Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Coupling exponent alpha (dimension 1). If provided, total logP = logMix + alpha * sum_i logL_i.",
            (RealParameter) null);

    private List<GenericTreeLikelihood> subLikelihoods;
    private RealParameter weights;
    private RealParameter alpha;
    private int K;

    @Override
    public void initAndValidate() {
        subLikelihoods = subLikelihoodsInput.get();
        weights = weightsInput.get();
        alpha = alphaInput.get();
        K = subLikelihoods.size();

        if (K < 2) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: need at least two subLikelihoods.");
        }
        if (weights.getDimension() != K) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: weights dimension ("
                    + weights.getDimension() + ") != number of subLikelihoods (" + K + ").");
        }
        if (alpha != null && alpha.getDimension() != 1) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: alpha must have dimension 1.");
        }

        double wsum = 0.0;
        for (int i = 0; i < K; i++) {
            double wi = weights.getArrayValue(i);
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

        final double a = (alpha == null ? 0.0 : alpha.getArrayValue(0));

        double maxTerm = Double.NEGATIVE_INFINITY;
        double sumLogL = 0.0;
        boolean needCoupling = (a > 0.0);

        final double[] term = new double[K];

        for (int i = 0; i < K; i++) {

            final double wi = weights.getArrayValue(i);
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
}
