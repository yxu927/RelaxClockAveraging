package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.inference.operator.kernel.BactrianRandomWalkOperator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;

import lphybeast.GeneratorToBEAST;
import mixture.beast.evolution.mixture.RelaxedRatesPriorSVS;
import mixture.beast.evolution.operator.IndicatorFlipOperator;
import mixture.lphy.evolution.auto.SVSRawBranchRates;

import java.util.StringJoiner;

public class SVSRawBranchRatesToBEAST implements GeneratorToBEAST<SVSRawBranchRates, RelaxedRatesPriorSVS> {

    @Override
    public RelaxedRatesPriorSVS generatorToBEAST(SVSRawBranchRates dist, BEASTInterface beastValue, BEASTContext context) {

        // Inputs from LPhy
        @SuppressWarnings("unchecked")
        Value<TimeTree> treeVal = (Value<TimeTree>) dist.getParams().get(SVSRawBranchRates.TREE);

        @SuppressWarnings("unchecked")
        Value<Integer> indVal = (Value<Integer>) dist.getParams().get(SVSRawBranchRates.INDICATOR);

        @SuppressWarnings("unchecked")
        Value<Double> ucldVal = (Value<Double>) dist.getParams().get(SVSRawBranchRates.UCLD_STDEV);

        @SuppressWarnings("unchecked")
        Value<Double> sigma2Val = (Value<Double>) dist.getParams().get(SVSRawBranchRates.SIGMA2);

        @SuppressWarnings("unchecked")
        Value<Double> rootVal = (Value<Double>) dist.getParams().get(SVSRawBranchRates.ROOT_LOG_RATE);

        TreeInterface beastTree = (TreeInterface) context.getBEASTObject(treeVal);

        Value<?> outVal = context.getOutput(dist);

        RealParameter ratesParam;
        if (beastValue instanceof RealParameter rp) {
            ratesParam = rp;
        } else {
            // Create from scratch if needed
            ratesParam = new RealParameter();
            ratesParam.setID("rawRates." + dist.getUniqueId());
        }

        // Ensure dimension (nodeCount - 1)
        int nNodes = ((beast.base.evolution.tree.Tree) beastTree).getNodeCount();
        int targetDim = nNodes - 1;

        if (ratesParam.getDimension() != targetDim) {
            // Replace mapping with a fresh parameter of correct size
            RealParameter newRates = new RealParameter();
            newRates.setID("rawRates." + dist.getUniqueId());

            StringJoiner sj = new StringJoiner(" ");
            for (int i = 0; i < targetDim; i++) sj.add("1.0");

            newRates.setInputValue("value", sj.toString());
            newRates.setInputValue("lower", "0.0");
            newRates.setInputValue("estimate", true);
            newRates.initAndValidate();

            // remove old object if it was registered
            if (beastValue != null) {
                context.removeBEASTObject(beastValue);
            }
            context.putBEASTObject(outVal, newRates);
            ratesParam = newRates;
        } else {
            // Make sure it's positive-constrained
            ratesParam.setLower(0.0);
        }

        // ---- Indicator parameter ----
        // Prefer IntegerParameter so the flip operator is clean.
        IntegerParameter indParam;
        BEASTInterface indObj = context.getBEASTObject(indVal);
        if (indObj instanceof IntegerParameter ip) {
            indParam = ip;
        } else {
            // Fallback: construct it (e.g. if the framework doesn't provide IntegerParameter)
            indParam = new IntegerParameter();
            indParam.setID("indicator." + dist.getUniqueId());
            indParam.setInputValue("value", "0");
            indParam.setInputValue("lower", 0);
            indParam.setInputValue("upper", 1);
            indParam.setInputValue("estimate", true);
            indParam.initAndValidate();
            context.putBEASTObject(indVal, indParam);
        }

        // Hyperparameters
        RealParameter ucldStdev = context.getAsRealParameter(ucldVal);
        RealParameter sigma2 = context.getAsRealParameter(sigma2Val);
        RealParameter rootLogRate = context.getAsRealParameter(rootVal);

        // ---- SVS prior distribution ----
        RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.setID("SVSRelaxedClockPrior." + dist.getUniqueId());

        prior.setInputValue("tree", beastTree);
        prior.setInputValue("rates", ratesParam);
        prior.setInputValue("indicator", indParam);
        prior.setInputValue("ucldStdev", ucldStdev);
        prior.setInputValue("sigma2", sigma2);
        prior.setInputValue("rootLogRate", rootLogRate);

        prior.initAndValidate();
        context.addBEASTObject(prior, dist);

        // ---- Operators ----
        // 1) Rate vector move (placeholder: RW in rate-space; OK to start, but scale-type is usually better)
        BactrianRandomWalkOperator rw = new BactrianRandomWalkOperator();
        rw.setID("rawRates.rw." + dist.getUniqueId());
        rw.setInputValue("parameter", ratesParam);
        rw.setInputValue("windowSize", 0.10);
        rw.setInputValue("weight", 10.0);
        rw.initAndValidate();
        context.addExtraOperator(rw);

        // 2) Flip indicator 0<->1
        IndicatorFlipOperator flip = new IndicatorFlipOperator();
        flip.setID("indicator.flip." + dist.getUniqueId());
        flip.setInputValue("indicator", indParam);
        flip.setInputValue("weight", 2.0);
        flip.initAndValidate();
        context.addExtraOperator(flip);

        // Prevent default operators / loggers
        if (ratesParam instanceof StateNode) context.addSkipOperator((StateNode) ratesParam);
        if (indParam instanceof StateNode) context.addSkipOperator((StateNode) indParam);

        context.addSkipLoggable(ratesParam);

        return prior;
    }

    @Override public Class<SVSRawBranchRates> getGeneratorClass() { return SVSRawBranchRates.class; }
    @Override public Class<RelaxedRatesPriorSVS> getBEASTClass() { return RelaxedRatesPriorSVS.class; }
}
