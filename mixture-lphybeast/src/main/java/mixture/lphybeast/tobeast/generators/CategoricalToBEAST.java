package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;

import lphy.base.distribution.Categorical;
import lphy.base.distribution.DistributionConstants;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import mixture.beast.evolution.mixture.CategoricalDistribution;
import mixture.lphybeast.tobeast.TypedParameterUtils;

public class CategoricalToBEAST implements GeneratorToBEAST<Categorical, CategoricalDistribution> {

    @Override
    public CategoricalDistribution generatorToBEAST(Categorical generator, BEASTInterface value, BEASTContext context) {

        CategoricalDistribution cat = new CategoricalDistribution();


        Value<?> pVal = generator.getParams().get(DistributionConstants.pParamName);
        BEASTInterface pObject = context.getBEASTObject(pVal);
        if (!TypedParameterUtils.isRandom(pVal) && pVal.value() instanceof Double[] probabilities) {
            BEASTInterface typedP = TypedParameterUtils.realVectorFrom(
                    pObject,
                    pVal.getId(),
                    TypedParameterUtils.values(probabilities),
                    false
            );
            TypedParameterUtils.replaceInContext(context, pVal, pObject, typedP);
            cat.setInputValue("pVector", typedP);
        } else {
            cat.setInputValue("p", pObject);
        }

        Value<?> outVal = context.getOutput(generator);
        BEASTInterface parameter = TypedParameterUtils.intScalarFrom(
                value,
                outVal != null ? outVal.getId() : "categorical." + generator.getUniqueId(),
                outVal != null && outVal.value() instanceof Integer ? (Integer) outVal.value() : 0,
                outVal != null && outVal.isRandom()
        );
        TypedParameterUtils.replaceInContext(context, outVal, value, parameter);
        cat.setInputValue("parameterScalar", parameter);

        cat.initAndValidate();
        return cat;
    }

    @Override
    public Class<Categorical> getGeneratorClass() {
        return Categorical.class;
    }

    @Override
    public Class<CategoricalDistribution> getBEASTClass() {
        return CategoricalDistribution.class;
    }
}
