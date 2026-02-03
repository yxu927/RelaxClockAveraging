package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import mixture.beast.evolution.mixture.SharedRatesClockModel;
import lphy.base.evolution.continuous.SharedRatesClock;

public class SharedRatesClockToBEAST implements GeneratorToBEAST<SharedRatesClock, SharedRatesClockModel> {

    @Override
    public SharedRatesClockModel generatorToBEAST(SharedRatesClock gen, BEASTInterface beastValue, BEASTContext context) {

        // Read LPhy inputs from params map (no need for getters)
        @SuppressWarnings("unchecked")
        Value<TimeTree> treeVal = (Value<TimeTree>) gen.getParams().get(SharedRatesClock.TREE);

        @SuppressWarnings("unchecked")
        Value<Double[]> ratesVal = (Value<Double[]>) gen.getParams().get(SharedRatesClock.RATES);

        @SuppressWarnings("unchecked")
        Value<Double> meanRateVal = (Value<Double>) gen.getParams().get(SharedRatesClock.MEAN_RATE);

        @SuppressWarnings("unchecked")
        Value<Boolean> normVal = (Value<Boolean>) gen.getParams().get(SharedRatesClock.NORMALIZE);

        TreeInterface beastTree = (TreeInterface) context.getBEASTObject(treeVal);

        // This should resolve to the RealParameter created for SVSRawBranchRates output
        RealParameter beastRates = context.getAsRealParameter(ratesVal);

        RealParameter beastMeanRate;
        if (meanRateVal != null) {
            beastMeanRate = context.getAsRealParameter(meanRateVal);
        } else {
            RealParameter fixed = new RealParameter("1.0");
            fixed.setID("meanRate.fixed." + gen.getUniqueId());
            fixed.initAndValidate();
            beastMeanRate = fixed;
        }

        boolean doNormalize = (normVal != null && Boolean.TRUE.equals(normVal.value()));

        SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.setID("SharedRatesClockModel." + gen.getUniqueId());

        clock.setInputValue("tree", beastTree);
        clock.setInputValue("rates", beastRates);
        clock.setInputValue("meanRate", beastMeanRate);
        clock.setInputValue("normalize", doNormalize);

        clock.initAndValidate();
        return clock;
    }

    @Override public Class<SharedRatesClock> getGeneratorClass() { return SharedRatesClock.class; }
    @Override public Class<SharedRatesClockModel> getBEASTClass() { return SharedRatesClockModel.class; }
}
