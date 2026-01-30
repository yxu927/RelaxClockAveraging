package mixture.beast.evolution.clockmodel.auto;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;

@Description("Auto-correlated clock: branch rate = E(Z) from bridged log-rates; "
        + "optional normalization to time-weighted mean 1, plus optional global meanRate multiplier.")
public class AutoCorrelatedClockModel extends BranchRateModel.Base {

    public final Input<Tree> treeInput = new Input<>(
            "tree",
            "tree used by this clock",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> nodeLogRatesInput = new Input<>(
            "nodeRates",
            "log(rate) for NON-root nodes, dimension = (tree.nodeCount - 1)",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "root log(rate), dimension=1. If normalize=true, this must be fixed at 0.",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "Brownian variance parameter (sigma^2)",
            Input.Validate.REQUIRED
    );

    public final Input<Integer> taylorOrderInput = new Input<>(
            "taylorOrder",
            "Taylor expansion order for E(Z) approximation (default=10)",
            10
    );

    public final Input<Boolean> normalizeInput = new Input<>(
            "normalize",
            "if true: rescale time-weighted mean rate across tree to 1",
            false
    );

    public final Input<RealParameter> meanRateInput = new Input<>(
            "meanRate",
            "overall rate multiplier (default=1.0). In normalize=true this is the absolute scale."
    );

    protected Tree tree;
    protected RealParameter nodeLogRates;
    protected RealParameter rootLogRate;
    protected RealParameter sigma2;
    protected RealParameter meanRate;

    protected int taylorOrder;
    protected boolean doNormalize;
    protected double[] branchRates;
    protected double[] storedBranchRates;
    protected double scaleFactor = 1.0;
    protected double storedScaleFactor = 1.0;
    protected int nNodes;
    protected int rootNr;
    protected int[] nodeIndexMapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        nodeLogRates = nodeLogRatesInput.get();
        rootLogRate = rootLogRateInput.get();
        sigma2 = sigma2Input.get();

        taylorOrder = taylorOrderInput.get();
        doNormalize = normalizeInput.get();

        meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }

        if (rootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate must have dimension=1.");
        }
        if (meanRate.getDimension() != 1) {
            throw new IllegalArgumentException("meanRate must have dimension=1.");
        }

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();

        if (nodeLogRates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("nodeRates must have dimension (nodeCount - 1). Found "
                    + nodeLogRates.getDimension() + " vs " + (nNodes - 1));
        }

        // normalize=true: enforce rootLogRate fixed at 0 (log-rate=0 => rate=1)
        if (doNormalize) {
            double v = rootLogRate.getValue();
            if (Math.abs(v) > 1e-12) {
                throw new IllegalArgumentException("normalize=true requires rootLogRate fixed to 0.0. Found: " + v);
            }
        }

        branchRates = new double[nNodes];
        storedBranchRates = new double[nNodes];

        // Build mapping by node number (nr)
        nodeIndexMapping = new int[nNodes];
        Arrays.fill(nodeIndexMapping, -2);

        int idx = 0;
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            int nr = node.getNr();
            if (nr < 0 || nr >= nNodes) {
                throw new IllegalArgumentException("Node nr out of range: " + nr + " (nNodes=" + nNodes + ")");
            }
            if (nr == rootNr) {
                nodeIndexMapping[nr] = -1;
            } else {
                nodeIndexMapping[nr] = idx++;
            }
        }

        recalcAllBranchRates();
        if (doNormalize) {
            computeScaleFactor();
        } else {
            scaleFactor = 1.0;
        }

        Log.info.println("AutoCorrelatedClockModel init: nNodes=" + nNodes
                + ", normalize=" + doNormalize
                + ", taylorOrder=" + taylorOrder);
    }

    @Override
    public double getRateForBranch(Node node) {
        if (node.isRoot()) {
            return 1.0;
        }

        double mr = meanRate.getValue();
        if (mr <= 0.0) {
            return 0.0;
        }

        return branchRates[node.getNr()] * scaleFactor * mr;
    }

    private void recalcAllBranchRates() {
        Arrays.fill(branchRates, 1.0);

        final double phi = sigma2.getValue();
        if (phi <= 0.0) {

            throw new IllegalArgumentException("sigma2 must be > 0. Found: " + phi);
        }

        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            if (node.isRoot()) continue;

            final int nr = node.getNr();
            final Node parent = node.getParent();
            final int parentNr = parent.getNr();

            final double vPar = parent.isRoot()
                    ? rootLogRate.getValue()
                    : nodeLogRates.getValue(nodeIndexMapping[parentNr]);

            final double vChi = nodeLogRates.getValue(nodeIndexMapping[nr]);

            final double rp = Math.exp(vPar);
            final double rc = Math.exp(vChi);

            double dt = parent.getHeight() - node.getHeight();
            if (dt <= 0.0) {

                branchRates[nr] = rp;
                continue;
            }

            branchRates[nr] = MeanZCalculator.computeMeanZ(rp, rc, dt, phi, taylorOrder);
        }
    }

    private void computeScaleFactor() {

        double sumRateTime = 0.0;
        double sumTime = 0.0;

        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            if (node.isRoot()) continue;

            final int nr = node.getNr();
            double dt = node.getParent().getHeight() - node.getHeight();
            if (dt <= 0.0) continue;

            sumRateTime += branchRates[nr] * dt;
            sumTime += dt;
        }

        if (sumRateTime <= 0.0) {
            scaleFactor = 1.0;
        } else {
            scaleFactor = sumTime / sumRateTime;
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        boolean ratesDirty = false;
        boolean anyDirty = false;

        if (tree != null && ((StateNode) tree).somethingIsDirty()) {
            ratesDirty = true;
            anyDirty = true;
        }
        if (nodeLogRates != null && nodeLogRates.somethingIsDirty()) {
            ratesDirty = true;
            anyDirty = true;
        }
        if (rootLogRate != null && rootLogRate.somethingIsDirty()) {
            // normalize=true: rootLogRate must not be estimated/proposed
            if (doNormalize) {
                throw new IllegalStateException("normalize=true: rootLogRate must be fixed at 0 and must not be proposed.");
            }
            ratesDirty = true;
            anyDirty = true;
        }
        if (sigma2 != null && sigma2.somethingIsDirty()) {
            ratesDirty = true;
            anyDirty = true;
        }

        if (meanRate != null && meanRate.somethingIsDirty()) {
            anyDirty = true;
        }

        if (ratesDirty) {
            recalcAllBranchRates();
            if (doNormalize) {
                computeScaleFactor();
            } else {
                scaleFactor = 1.0;
            }
        }

        return anyDirty;
    }

    @Override
    public void store() {
        System.arraycopy(branchRates, 0, storedBranchRates, 0, nNodes);
        storedScaleFactor = scaleFactor;
        super.store();
    }

    @Override
    public void restore() {
        double[] tmp = branchRates;
        branchRates = storedBranchRates;
        storedBranchRates = tmp;

        scaleFactor = storedScaleFactor;
        super.restore();
    }


}
