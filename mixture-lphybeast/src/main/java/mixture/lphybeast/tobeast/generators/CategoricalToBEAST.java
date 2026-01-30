package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;

import lphy.base.distribution.Categorical;
import lphy.base.distribution.DistributionConstants;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import mixture.beast.evolution.clockmodel.auto.CategoricalDistribution;

public class CategoricalToBEAST implements GeneratorToBEAST<Categorical, CategoricalDistribution> {

    @Override
    public CategoricalDistribution generatorToBEAST(Categorical generator, BEASTInterface value, BEASTContext context) {

        CategoricalDistribution cat = new CategoricalDistribution();


        Value<?> pVal = generator.getParams().get(DistributionConstants.pParamName);
        cat.setInputValue("p", context.getBEASTObject(pVal));


        cat.setInputValue("parameter", value);

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
