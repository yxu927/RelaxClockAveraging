package mixture.beast.evolution.mixture;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;

@Description("Clock model that reads a shared vector of positive branch rates (non-root nodes). "
        + "Optional normalization rescales the time-weighted mean rate across the tree to 1. "
        + "A separate meanRate multiplier can set the absolute scale.")
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

    private int nNodes;
    private int rootNr;

    // Map nodeNr -> index in rates vector; root -> -1
    private int[] idxMap;

    // Derived normalization factor (only used if normalize=true)
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

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();

        if (rates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("rates must have dimension (nodeCount - 1). Found "
                    + rates.getDimension() + " vs " + (nNodes - 1));
        }

        buildIndexMap();

        // Compute initial scale factor
        if (doNormalize) {
            computeScaleFactor();
        } else {
            scaleFactor = 1.0;
        }

        Log.info.println("SharedRatesClockModel init: nNodes=" + nNodes
                + ", normalize=" + doNormalize);
    }

    private void buildIndexMap() {
        idxMap = new int[nNodes];
        Arrays.fill(idxMap, -2);

        int k = 0;
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            int nr = node.getNr();
            if (nr < 0 || nr >= nNodes) {
                throw new IllegalArgumentException("Node nr out of range: " + nr + " (nNodes=" + nNodes + ")");
            }
            if (nr == rootNr) {
                idxMap[nr] = -1;
            } else {
                idxMap[nr] = k++;
            }
        }

        if (k != nNodes - 1) {
            throw new IllegalStateException("Internal error: expected to map " + (nNodes - 1)
                    + " non-root nodes but mapped " + k);
        }
    }

    @Override
    public double getRateForBranch(Node node) {
        if (node.isRoot()) {
            return 1.0;
        }

        final double mr = meanRate.getValue();
        if (!(mr > 0.0)) {
            return 0.0;
        }

        final int nr = node.getNr();
        final int idx = idxMap[nr];
        if (idx < 0) {
            // should never happen for non-root nodes
            return 0.0;
        }

        final double r = rates.getValue(idx);
        if (!(r > 0.0)) {
            return 0.0;
        }

        return r * scaleFactor * mr;
    }

    /**
     * Computes a scale factor so that the time-weighted mean rate across the tree equals 1.
     * Only used when normalize=true.
     */
    private void computeScaleFactor() {
        double sumRateTime = 0.0;
        double sumTime = 0.0;

        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            if (node.isRoot()) continue;

            final double dt = node.getLength(); // branch length in time units
            if (!(dt > 0.0)) continue;

            final int idx = idxMap[node.getNr()];
            final double r = rates.getValue(idx);
            if (!(r > 0.0)) continue;

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
