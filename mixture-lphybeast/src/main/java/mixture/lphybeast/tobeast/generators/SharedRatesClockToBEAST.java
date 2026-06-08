package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import mixture.beast.evolution.mixture.SharedRatesClockModel;
import mixture.lphy.evolution.auto.SharedRatesClock;
import mixture.lphybeast.tobeast.TypedParameterUtils;

public class SharedRatesClockToBEAST implements GeneratorToBEAST<SharedRatesClock, SharedRatesClockModel> {

    @Override
    public SharedRatesClockModel generatorToBEAST(SharedRatesClock gen, BEASTInterface beastValue, BEASTContext context) {

        @SuppressWarnings("unchecked")
        Value<TimeTree> treeVal = (Value<TimeTree>) gen.getParams().get(SharedRatesClock.TREE);

        @SuppressWarnings("unchecked")
        Value<Double[]> ratesVal = (Value<Double[]>) gen.getParams().get(SharedRatesClock.RATES);

        @SuppressWarnings("unchecked")
        Value<Double> meanRateVal = (Value<Double>) gen.getParams().get(SharedRatesClock.MEAN_RATE);

        @SuppressWarnings("unchecked")
        Value<Boolean> normVal = (Value<Boolean>) gen.getParams().get(SharedRatesClock.NORMALIZE);

        TreeInterface beastTree = (TreeInterface) context.getBEASTObject(treeVal);

        BEASTInterface beastRates = context.getBEASTObject(ratesVal);
        BEASTInterface beastMeanRate = meanRateVal == null ? null : context.getBEASTObject(meanRateVal);

        boolean doNormalize = (normVal != null && Boolean.TRUE.equals(normVal.value()));

        SharedRatesClockModel clock = new SharedRatesClockModel();
        clock.setID("SharedRatesClockModel." + gen.getUniqueId());

        clock.setInputValue("tree", beastTree);
        if (beastRates instanceof RealVector<?>) {
            clock.setInputValue("ratesVector", beastRates);
        } else {
            clock.setInputValue("rates", context.getAsRealParameter(ratesVal));
        }
        if (beastMeanRate instanceof RealScalar<?>) {
            clock.setInputValue("meanRateScalar", beastMeanRate);
        } else if (meanRateVal != null) {
            clock.setInputValue("meanRate", context.getAsRealParameter(meanRateVal));
        } else {
            clock.setInputValue("meanRateScalar",
                    TypedParameterUtils.positiveRealScalar("meanRate.fixed." + gen.getUniqueId(), 1.0, false));
        }
        clock.setInputValue("normalize", doNormalize);

        clock.initAndValidate();
        return clock;
    }

    @Override public Class<SharedRatesClock> getGeneratorClass() { return SharedRatesClock.class; }
    @Override public Class<SharedRatesClockModel> getBEASTClass() { return SharedRatesClockModel.class; }
}
