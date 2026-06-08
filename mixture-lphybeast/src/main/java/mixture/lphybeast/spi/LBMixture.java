package mixture.lphybeast.spi;

import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Generator;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTMapping;
import mixture.lphybeast.tobeast.generators.CategoricalToBEAST;
import mixture.lphybeast.tobeast.generators.MixturePhyloCTMCToBEAST;
import mixture.lphybeast.tobeast.generators.SVSRawBranchRatesToBEAST;
import mixture.lphybeast.tobeast.generators.SharedRatesClockToBEAST;

import java.util.List;
import java.util.Map;

/**
 * The "Container" provider class of SPI
 * which include a list of {@link ValueToBEAST},
 * {@link GeneratorToBEAST}, and {@link DataType}
 * to extend.
 *
 */
public class LBMixture implements LPhyBEASTMapping {

    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return List.of();
    }

    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return List.of(
                MixturePhyloCTMCToBEAST.class,
                CategoricalToBEAST.class,
                SharedRatesClockToBEAST.class,
                SVSRawBranchRatesToBEAST.class
        );
    }

    @Override
    public Map<SequenceType, beast.base.evolution.datatype.DataType> getDataTypeMap() {
        return Map.of();
    }

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return List.of();
    }

    @Override
    public List<Class> getExcludedValueType() {
        return List.of(TimeTreeNode.class);
    }

}
