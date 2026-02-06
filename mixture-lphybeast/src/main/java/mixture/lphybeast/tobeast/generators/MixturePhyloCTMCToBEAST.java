package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.ThreadedTreeLikelihood;
import beast.base.inference.Distribution;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.likelihood.AbstractPhyloCTMC;
import lphy.base.evolution.likelihood.MixturePhyloCTMC;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.generators.PhyloCTMCToBEAST;
import mixture.beast.evolution.mixture.MixtureTreeLikelihood;
import mutablealignment.MATreeLikelihood;
import mixture.beast.evolution.mixture.MixtureLikelihoodLogger;

import java.util.ArrayList;
import java.util.List;

public class MixturePhyloCTMCToBEAST implements
        GeneratorToBEAST<MixturePhyloCTMC, MixtureTreeLikelihood> {

    @Override
    public MixtureTreeLikelihood generatorToBEAST(
            MixturePhyloCTMC mix,
            BEASTInterface value,
            BEASTContext context) {

        if (!(value instanceof Alignment beastAlignment)) {
            throw new IllegalArgumentException("MixturePhyloCTMCToBEAST only supports BEAST Alignment data");
        }


        @SuppressWarnings("unchecked")
        Value<lphy.base.evolution.alignment.Alignment> comp1Val =
                (Value<lphy.base.evolution.alignment.Alignment>) mix.getParams().get(MixturePhyloCTMC.COMP1);
        @SuppressWarnings("unchecked")
        Value<lphy.base.evolution.alignment.Alignment> comp2Val =
                (Value<lphy.base.evolution.alignment.Alignment>) mix.getParams().get(MixturePhyloCTMC.COMP2);
        @SuppressWarnings("unchecked")
        Value<lphy.base.evolution.alignment.Alignment> comp3Val =
                (Value<lphy.base.evolution.alignment.Alignment>) mix.getParams().get(MixturePhyloCTMC.COMP3); // 可能为 null（2 分量时）

        List<Value<lphy.base.evolution.alignment.Alignment>> compVals = new ArrayList<>(3);
        compVals.add(comp1Val);
        compVals.add(comp2Val);
        if (comp3Val != null) compVals.add(comp3Val);


        for (Value<?> compVal : compVals) {
            BEASTInterface alg = context.getBEASTObject(compVal);
            if (alg instanceof Alignment && alg != beastAlignment) {

                context.removeBEASTObject(alg);
            }

            context.putBEASTObject(compVal, beastAlignment);
        }

        Value<?> mixOutput = context.getOutput(mix);
        boolean observed = mixOutput != null && context.isObserved(mixOutput);

        final int K = compVals.size();
        final List<GenericTreeLikelihood> subTL = new ArrayList<>(K);
        boolean addBranchOperatorsThisTime = true;

        for (int i = 0; i < K; i++) {
            AbstractPhyloCTMC comp = (AbstractPhyloCTMC) compVals.get(i).getGenerator();

            GenericTreeLikelihood tl;
            if (observed) {
                tl = new ThreadedTreeLikelihood();
                tl.setInputValue("useAmbiguities", true);
            } else {
                tl = new MATreeLikelihood();
                tl.setInputValue("useAmbiguities", false);
            }

            tl.setInputValue("data", beastAlignment);

            PhyloCTMCToBEAST.constructTreeAndBranchRate(
                    asPhyloCTMC(comp),
                    tl,
                    context,
                    /* skipBranchOperators = */ !addBranchOperatorsThisTime
            );
            addBranchOperatorsThisTime = false;

            // siteModel
            tl.setInputValue("siteModel", PhyloCTMCToBEAST.constructSiteModel(asPhyloCTMC(comp), context));

            tl.initAndValidate();
            tl.setID(beastAlignment.getID() + ".treeLikelihood.component" + i);

            context.addExtraLoggable(tl);
            subTL.add(tl);
        }

        MixtureTreeLikelihood beastMix = new MixtureTreeLikelihood();
        beastMix.setInputValue("subLikelihood", subTL);

        @SuppressWarnings("unchecked")
        Value<Double[]> wVal = (Value<Double[]>) mix.getParams().get(MixturePhyloCTMC.WEIGHTS);

        Double[] w = wVal.value();
        if (w == null || w.length != K) {
            throw new IllegalArgumentException("MixturePhyloCTMCToBEAST: weights length (" +
                    (w == null ? 0 : w.length) + ") must equal number of components (" + K + ").");
        }

        RealParameter wParam = context.getAsRealParameter(wVal);
        beastMix.setInputValue("weights", wParam);

        beastMix.initAndValidate();
        beastMix.setID(beastAlignment.getID() + ".mixtureTreeLikelihood");
        context.addExtraLoggable(beastMix);

        for (Value<?> compVal : compVals) {
            AbstractPhyloCTMC comp = (AbstractPhyloCTMC) compVal.getGenerator();
            BEASTInterface tlSingle = context.getBEASTObject(comp);
            if (tlSingle instanceof Distribution || tlSingle instanceof GenericTreeLikelihood) {
                context.removeBEASTObject(tlSingle);
            }
        }

        MixtureLikelihoodLogger allocLogger = new MixtureLikelihoodLogger();
        allocLogger.setInputValue("mixture", beastMix);
        allocLogger.setInputValue("printResp", true);

        allocLogger.initAndValidate();
        allocLogger.setID(beastAlignment.getID() + ".mixtureAlloc");
        context.addExtraLoggable(allocLogger);

        return beastMix;
    }

    private static PhyloCTMC asPhyloCTMC(AbstractPhyloCTMC ctmc) {
        if (ctmc instanceof PhyloCTMC) return (PhyloCTMC) ctmc;
        throw new IllegalArgumentException(
                "MixturePhyloCTMCToBEAST now only supports PhyloCTMC, got: " + ctmc.getClass().getName());
    }

    @Override
    public Class<MixturePhyloCTMC> getGeneratorClass() {
        return MixturePhyloCTMC.class;
    }

    @Override
    public Class<MixtureTreeLikelihood> getBEASTClass() {
        return MixtureTreeLikelihood.class;
    }
}
