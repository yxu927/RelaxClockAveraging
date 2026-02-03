package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import mixture.beast.evolution.mixture.RelaxedRatesPriorSVS;

import lphy.base.evolution.continuous.SVSRawBranchRates;

import mixture.beast.evolution.operator.IndicatorGibbsOperator;
import mixture.beast.evolution.operator.SingleRateScaleOperator;
import mixture.beast.evolution.operator.SubtreeRateScaleOperator;

import java.util.StringJoiner;

public class SVSRawBranchRatesToBEAST implements GeneratorToBEAST<SVSRawBranchRates, RelaxedRatesPriorSVS> {

    @Override
    public RelaxedRatesPriorSVS generatorToBEAST(SVSRawBranchRates dist, BEASTInterface beastValue, BEASTContext context) {

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

        // ---- rates parameter ----
        RealParameter ratesParam;
        if (beastValue instanceof RealParameter rp) {
            ratesParam = rp;
        } else {
            ratesParam = new RealParameter();
            ratesParam.setID("rawRates." + dist.getUniqueId());
        }

        int nNodes = ((beast.base.evolution.tree.Tree) beastTree).getNodeCount();
        int targetDim = nNodes - 1;

        if (ratesParam.getDimension() != targetDim) {
            RealParameter newRates = new RealParameter();
            newRates.setID("rawRates." + dist.getUniqueId());

            StringJoiner sj = new StringJoiner(" ");
            for (int i = 0; i < targetDim; i++) sj.add("1.0");

            newRates.setInputValue("value", sj.toString());
            newRates.setInputValue("lower", "0.0");
            newRates.setInputValue("estimate", true);
            newRates.initAndValidate();

            if (beastValue != null) {
                context.removeBEASTObject(beastValue);
            }
            context.putBEASTObject(outVal, newRates);
            ratesParam = newRates;
        } else {
            ratesParam.setLower(0.0);
        }

        // ---- indicator parameter ----
        IntegerParameter indParam;
        BEASTInterface indObj = context.getBEASTObject(indVal);
        if (indObj instanceof IntegerParameter ip) {
            indParam = ip;
        } else {
            indParam = new IntegerParameter();
            indParam.setID("indicator." + dist.getUniqueId());
            indParam.setInputValue("value", "0");
            indParam.setInputValue("lower", 0);
            indParam.setInputValue("upper", 1);
            indParam.setInputValue("estimate", true);
            indParam.initAndValidate();
            context.putBEASTObject(indVal, indParam);
        }

        // ---- hyperparameters ----
        RealParameter ucldStdev = context.getAsRealParameter(ucldVal);
        RealParameter sigma2 = context.getAsRealParameter(sigma2Val);
        RealParameter rootLogRate = context.getAsRealParameter(rootVal);

        // ---- SVS prior ----
        RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.setID("SVSRelaxedClockPrior." + dist.getUniqueId());

        prior.setInputValue("tree", beastTree);
        prior.setInputValue("rates", ratesParam);
        prior.setInputValue("indicator", indParam);
        prior.setInputValue("ucldStdev", ucldStdev);
        prior.setInputValue("sigma2", sigma2);
        prior.setInputValue("rootLogRate", rootLogRate);

        // Optional: if you want to be lenient with tiny branches
        // prior.setInputValue("minBranchLength", 1e-12);

        prior.initAndValidate();
        context.addBEASTObject(prior, dist);

        // ---- Operators ----

        // 1) Single-dimension multiplicative scale on positive rates
        SingleRateScaleOperator oneScale = new SingleRateScaleOperator();
        oneScale.setID("rawRates.oneScale." + dist.getUniqueId());
        oneScale.setInputValue("rates", ratesParam);
        oneScale.setInputValue("window", 0.3);
        oneScale.setInputValue("weight", 15.0);
        oneScale.initAndValidate();
        context.addExtraOperator(oneScale);

        // 2) Subtree correlated scale move (very useful for AC)
        SubtreeRateScaleOperator subScale = new SubtreeRateScaleOperator();
        subScale.setID("rawRates.subtreeScale." + dist.getUniqueId());
        subScale.setInputValue("tree", beastTree);
        subScale.setInputValue("rates", ratesParam);
        subScale.setInputValue("window", 0.25);
        subScale.setInputValue("weight", 10.0);
        subScale.initAndValidate();
        context.addExtraOperator(subScale);

        // 3) Gibbs-style indicator update (replaces flip)
        IndicatorGibbsOperator gibbs = new IndicatorGibbsOperator();
        gibbs.setID("indicator.gibbs." + dist.getUniqueId());
        gibbs.setInputValue("indicator", indParam);
        gibbs.setInputValue("prior", prior);
        gibbs.setInputValue("pOne", 0.5);   // should match your Categorical prior weights
        gibbs.setInputValue("weight", 2.0);
        gibbs.initAndValidate();
        context.addExtraOperator(gibbs);

        // Prevent default operators/loggers for these
        if (ratesParam instanceof StateNode) context.addSkipOperator((StateNode) ratesParam);
        if (indParam instanceof StateNode) context.addSkipOperator((StateNode) indParam);

        context.addSkipLoggable(ratesParam);

        return prior;
    }

    @Override public Class<SVSRawBranchRates> getGeneratorClass() { return SVSRawBranchRates.class; }
    @Override public Class<RelaxedRatesPriorSVS> getBEASTClass() { return RelaxedRatesPriorSVS.class; }
}
