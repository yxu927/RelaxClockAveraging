
package rc.beast.evolution.clockmodel;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Input;
import beast.base.core.Input.Validate;

import beast.base.core.Log;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;


public class MixtureTreeLikelihood extends Distribution {



    public final Input<List<TreeLikelihood>> subLikelihoodsInput = new Input<>("subLikelihood",
                    "List of component tree likelihoods to be mixed", new ArrayList<>());


    public final Input<RealParameter> weightsInput = new Input<>("weights",
            "Mixture weights (Dirichlet / Delta‑exchange parameter). " + "Dimension must equal number of subLikelihoods.",
                    Validate.REQUIRED);



    private List<TreeLikelihood> subLikelihoods;
    private RealParameter        weights;
    private int                  K;
    private beast.base.evolution.alignment.Alignment alignment;



    @Override
    public void initAndValidate() {

        subLikelihoods = subLikelihoodsInput.get();
        weights        = weightsInput.get();
        K              = subLikelihoods.size();

        if (K < 2) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: need at least two subLikelihoods.");
        }
        if (weights.getDimension() != K) {
            throw new IllegalArgumentException("MixtureTreeLikelihood: weights dimension ("
                    + weights.getDimension() + ") ≠ number of subLikelihoods (" + K + ").");
        }

        alignment = subLikelihoods.get(0).dataInput.get();
        for (TreeLikelihood tl : subLikelihoods) {
            if (tl.dataInput.get() != alignment) {
                throw new IllegalArgumentException(
                        "All subLikelihoods must reference the *same* Alignment instance.");
            }
        }


    }



    @Override
    public double calculateLogP() {

        double maxLog = Double.NEGATIVE_INFINITY;
        final double[] logLi = new double[K];

        for (int i = 0; i < K; i++) {
            final double wi = weights.getArrayValue(i);

            if (wi == 0.0) {
                logLi[i] = Double.NEGATIVE_INFINITY;
                continue;
            }

            logLi[i] = subLikelihoods.get(i).calculateLogP();

            if (logLi[i] > maxLog) {
                maxLog = logLi[i];
            }
        }

        double sum = 0.0;
        for (int i = 0; i < K; i++) {
            final double wi = weights.getArrayValue(i);
            if (wi == 0.0) {
                continue;
            }
            sum += wi * Math.exp(logLi[i] - maxLog);
        }

        if (sum <= 0.0 || Double.isNaN(sum)) {
            logP = Double.NEGATIVE_INFINITY;
            Log.warning("MixtureTreeLikelihood: weighted likelihood sum underflowed to zero – returning -Inf");
            return logP;
        }

        logP = Math.log(sum) + maxLog;
        return logP;
    }

    @Override
    public List<String> getArguments() {
        return List.of();
    }

    @Override
    public List<String> getConditions() {
        return List.of();
    }

    @Override
    public void sample(State state, Random random) {

    }


    @Override
    protected boolean requiresRecalculation() {
        return true;
    }


    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }


}
