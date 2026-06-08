package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.type.RealVector;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import mixture.beast.evolution.mixture.SharedRatesClockModelSpec;
import mixture.lphy.evolution.auto.SharedRatesClock;
import mixture.lphybeast.tobeast.TypedParameterUtils;

public class SharedRatesClockToBEAST implements GeneratorToBEAST<SharedRatesClock, SharedRatesClockModelSpec> {

    @Override
    public SharedRatesClockModelSpec generatorToBEAST(SharedRatesClock gen, BEASTInterface beastValue, BEASTContext context) {

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

        boolean doNormalize = (normVal != null && Boolean.TRUE.equals(normVal.value()));

        SharedRatesClockModelSpec clock = new SharedRatesClockModelSpec();
        clock.setID("SharedRatesClockModelSpec." + gen.getUniqueId());

        clock.setInputValue("tree", beastTree);
        if (beastRates instanceof RealVector<?>) {
            clock.setInputValue("ratesVector", beastRates);
        } else {
            clock.setInputValue("ratesVector", context.getAsRealVector(ratesVal));
        }
        if (meanRateVal != null) {
            clock.setInputValue("meanRateScalar", context.getAsRealScalar(meanRateVal));
        } else {
            clock.setInputValue("meanRateScalar",
                    TypedParameterUtils.positiveRealScalar("meanRate.fixed." + gen.getUniqueId(), 1.0, false));
        }
        clock.setInputValue("normalize", doNormalize);

        clock.initAndValidate();
        return clock;
    }

    @Override public Class<SharedRatesClock> getGeneratorClass() { return SharedRatesClock.class; }
    @Override public Class<SharedRatesClockModelSpec> getBEASTClass() { return SharedRatesClockModelSpec.class; }
}
