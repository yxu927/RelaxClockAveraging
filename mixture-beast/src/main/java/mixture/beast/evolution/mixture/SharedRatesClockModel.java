package mixture.beast.evolution.mixture;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import mixture.beast.evolution.util.BranchRateIndexHelper;

public class SharedRatesClockModel extends BranchRateModel.Base {

    public final Input<Tree> treeInput = new Input<>(
            "tree",
            "tree used by this clock",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "positive branch rates for NON-root nodes; dimension=(tree.nodeCount - 1)",
            Input.Validate.REQUIRED
    );

    public final Input<Boolean> normalizeInput = new Input<>(
            "normalize",
            "if true: rescale time-weighted mean rate across the tree to 1",
            false
    );

    public final Input<RealParameter> meanRateInput = new Input<>(
            "meanRate",
            "overall rate multiplier (dimension=1, default=1.0)."
    );

    private Tree tree;
    private RealParameter rates;
    private RealParameter meanRate;
    private boolean doNormalize;

    private BranchRateIndexHelper.Mapping mapping;

    private double scaleFactor = 1.0;
    private double storedScaleFactor = 1.0;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();
        doNormalize = normalizeInput.get();

        meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }

        if (meanRate.getDimension() != 1) {
            throw new IllegalArgumentException("meanRate must have dimension=1.");
        }

        BranchRateIndexHelper.validateRatesDimension(tree, rates, "SharedRatesClockModel");
        mapping = BranchRateIndexHelper.buildDeterministic(tree);

        if (doNormalize) {
            computeScaleFactor();
        } else {
            scaleFactor = 1.0;
        }

        Log.info.println("SharedRatesClockModel init: nNodes=" + mapping.getNodeCount()
                + ", normalize=" + doNormalize);
    }

    private void ensureMappingUpToDate() {
        mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, rates, "SharedRatesClockModel");
    }

    @Override
    public double getRateForBranch(final Node node) {
        if (node.isRoot()) {
            return 1.0;
        }

        ensureMappingUpToDate();

        final double mr = meanRate.getValue();
        if (!(mr > 0.0)) {
            return 0.0;
        }

        final int idx = mapping.idxForNode(node);
        if (idx < 0) {
            return 0.0;
        }

        final double r = rates.getValue(idx);
        if (!(r > 0.0)) {
            return 0.0;
        }

        return r * scaleFactor * mr;
    }

    private void computeScaleFactor() {
        ensureMappingUpToDate();

        double sumRateTime = 0.0;
        double sumTime = 0.0;

        final int nNodes = mapping.getNodeCount();
        for (int i = 0; i < nNodes; i++) {
            final Node node = tree.getNode(i);
            if (node.isRoot()) {
                continue;
            }

            final double dt = node.getLength();
            if (!(dt > 0.0)) {
                continue;
            }

            final int idx = mapping.idxForNode(node);
            final double r = rates.getValue(idx);
            if (!(r > 0.0)) {
                continue;
            }

            sumRateTime += r * dt;
            sumTime += dt;
        }

        if (!(sumRateTime > 0.0) || !(sumTime > 0.0)) {
            scaleFactor = 1.0;
        } else {
            scaleFactor = sumTime / sumRateTime;
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        boolean dirty = false;

        if (tree != null && ((StateNode) tree).somethingIsDirty()) {
            dirty = true;
            ensureMappingUpToDate();
        }
        if (rates != null && rates.somethingIsDirty()) {
            dirty = true;
        }
        if (meanRate != null && meanRate.somethingIsDirty()) {
            dirty = true;
        }

        if (dirty) {
            if (doNormalize) {
                computeScaleFactor();
            } else {
                scaleFactor = 1.0;
            }
        }

        return dirty;
    }

    @Override
    public void store() {
        storedScaleFactor = scaleFactor;
        super.store();
    }

    @Override
    public void restore() {
        scaleFactor = storedScaleFactor;
        super.restore();
    }
}
