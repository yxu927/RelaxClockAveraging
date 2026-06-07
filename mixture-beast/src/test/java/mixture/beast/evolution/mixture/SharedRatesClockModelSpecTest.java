package mixture.beast.evolution.mixture;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class SharedRatesClockModelSpecTest {

    private static final double EPS = 1.0e-10;
    private static final String TREE_NEWICK = "((A:2.0,B:3.0):4.0,C:5.0);";

    @Test
    public void specWrapperIsAssignableToBeast3ClockModelBase() {
        assertTrue(new SharedRatesClockModelSpec() instanceof beast.base.spec.evolution.branchratemodel.Base);
    }

    @Test
    public void legacyInputsMatchExistingClockModelRates() {
        final Tree legacyTree = fixedTree();
        final Tree specTree = fixedTree();

        final SharedRatesClockModel legacy = new SharedRatesClockModel();
        legacy.initByName("tree", legacyTree,
                "rates", new RealParameter("0.5 1.0 2.0 4.0"),
                "normalize", true,
                "meanRate", new RealParameter("3.0"));

        final SharedRatesClockModelSpec spec = new SharedRatesClockModelSpec();
        spec.initByName("tree", specTree,
                "rates", new RealParameter("0.5 1.0 2.0 4.0"),
                "normalize", true,
                "meanRate", new RealParameter("3.0"));

        assertSameBranchRates(legacyTree, legacy, specTree, spec);
    }

    @Test
    public void typedInputsMatchExistingClockModelRates() {
        final Tree legacyTree = fixedTree();
        final Tree specTree = fixedTree();

        final SharedRatesClockModel legacy = new SharedRatesClockModel();
        legacy.initByName("tree", legacyTree,
                "ratesVector", realVectorParam(0.5, 1.0, 2.0, 4.0),
                "normalize", true,
                "meanRateScalar", realScalarParam(3.0));

        final SharedRatesClockModelSpec spec = new SharedRatesClockModelSpec();
        spec.initByName("tree", specTree,
                "ratesVector", realVectorParam(0.5, 1.0, 2.0, 4.0),
                "normalize", true,
                "meanRateScalar", realScalarParam(3.0));

        assertSameBranchRates(legacyTree, legacy, specTree, spec);
    }

    private static Tree fixedTree() {
        return new TreeParser(TREE_NEWICK, false, true, true, 1);
    }

    private static RealVectorParam<PositiveReal> realVectorParam(final double... values) {
        return new RealVectorParam<>(values, PositiveReal.INSTANCE);
    }

    private static RealScalarParam<PositiveReal> realScalarParam(final double value) {
        return new RealScalarParam<>(value, PositiveReal.INSTANCE);
    }

    private static void assertSameBranchRates(final Tree legacyTree,
                                              final SharedRatesClockModel legacy,
                                              final Tree specTree,
                                              final SharedRatesClockModelSpec spec) {
        for (int i = 0; i < legacyTree.getNodeCount(); i++) {
            final Node legacyNode = legacyTree.getNode(i);
            final Node specNode = specTree.getNode(i);
            assertEquals(legacy.getRateForBranch(legacyNode), spec.getRateForBranch(specNode), EPS);
        }
    }
}
