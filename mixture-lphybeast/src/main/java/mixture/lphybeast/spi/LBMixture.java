package mixture.lphybeast.spi;



import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.datatype.DataType;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Generator;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;

import mixture.lphybeast.tobeast.generators.*;

import java.util.Arrays;
import java.util.List;
import java.util.Map;


/**
 * The "Container" provider class of SPI
 * which include a list of {@link ValueToBEAST},
 * {@link GeneratorToBEAST}, and {@link DataType}
 * to extend.
 *
 */
public class LBMixture implements LPhyBEASTExt {

    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return Arrays.asList(

        );
    }

    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return Arrays.asList(
                AutoCorrelatedClockToBEAST.class, AutoCorrelatedLogRatesToBEAST.class,
              MixturePhyloCTMCToBEAST.class,
                CategoricalToBEAST.class
        );
    }

    @Override
    public Map<SequenceType, beast.base.evolution.datatype.DataType> getDataTypeMap() {
        return Map.of();
    }

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return Arrays.asList(
        );
    }

    @Override
    public List<Class> getExcludedValueType() {
        return Arrays.asList(TimeTreeNode.class
        );
    }

}
