package rc.beast.evolution.clockmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.UpDownOperator;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;


@Description("AutoCorrelatedClockModel using bridging from parent->child log-rates and optional normalization.")
public class AutoCorrelatedClockModel extends BranchRateModel.Base implements UpDownOp {

    // Inputs
    public final Input<Tree> treeInput = new Input<>(
            "tree",
            "the tree used by this auto-correlated clock",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> nodeLogRatesInput = new Input<>(
            "nodeRates",
            "log(rate) for NON-root nodes, dimension=(tree.nodeCount - 1)",
            Input.Validate.REQUIRED
    );



    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "the log(rate) for the root node (dimension=1)",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "the Brownian variance parameter (sigma^2)",
            Input.Validate.REQUIRED
    );

    public final Input<Integer> taylorOrderInput = new Input<>(
            "taylorOrder",
            "Taylor expansion order for bridging integral (default=10)",
            10
    );

    public final Input<Boolean> normalizeInput = new Input<>(
            "normalize",
            "whether to scale the average rate to 1.0 across the tree (default false)",
            false
    );

    // Optional global multiplier for all branch rates
    public final Input<RealParameter> meanRateInput = new Input<>(
            "myMeanRate",  // <-- renamed!
            "global rate multiplier (default=1.0)"
    );

    // Class fields
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

    protected boolean recompute = true;
    protected boolean renormalize = true;

    protected int nNodes;
    protected int rootNr;

    // nodeIndexMapping[i] = index in nodeLogRates for node i, or -1 if root
    protected int[] nodeIndexMapping;

    @Override
    public void initAndValidate() {
        // Retrieve inputs
        tree = treeInput.get();
        nodeLogRates = nodeLogRatesInput.get();
        rootLogRate = rootLogRateInput.get();
        sigma2 = sigma2Input.get();
        taylorOrder = taylorOrderInput.get();
        doNormalize = normalizeInput.get();

        // Handle optional meanRate
        meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }


        nNodes = tree.getNodeCount();
        Node root = tree.getRoot();
        rootNr = root.getNr();

        if (rootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate must have dimension=1");
        }
        if (nodeLogRates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException(
                    "nodeRates dimension must be (tree.nodeCount - 1). Found: " +
                            nodeLogRates.getDimension() + " vs " + (nNodes - 1)
            );
        }

        branchRates = new double[nNodes];
        storedBranchRates = new double[nNodes];

        nodeIndexMapping = new int[nNodes];
        int idx = 0;
        for (int i = 0; i < nNodes; i++) {
            if (i == rootNr) {
                nodeIndexMapping[i] = -1; // root
            } else {
                nodeIndexMapping[i] = idx;
                idx++;
            }
        }

        Log.info.println("AutoCorrelatedClockModel init: nNodes=" + nNodes +
                ", doNormalize=" + doNormalize + ", taylorOrder=" + taylorOrder);
    }

    @Override
    public double getRateForBranch(Node node) {
        if (node.isRoot()) {
            // root has no parent => no branch
            return 1.0;
        }
        synchronized (this) {
            if (recompute) {
                recalcAllBranchRates();
                recompute = false;
                renormalize = true;
            }
            if (renormalize && doNormalize) {
                computeScaleFactor();
                renormalize = false;
            }
        }

        return branchRates[node.getNr()] * scaleFactor;
    }

    private void recalcAllBranchRates() {
        Arrays.fill(branchRates, 1.0);

        double phi = sigma2.getValue();
        double globalRateVal = meanRate.getValue();

        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            if (!node.isRoot()) {
                Node parent = node.getParent();
                int parentNr = parent.getNr();

                double vPar = (parent.isRoot())
                        ? rootLogRate.getValue()
                        : nodeLogRates.getValue(nodeIndexMapping[parentNr]);

                double vChi = nodeLogRates.getValue(nodeIndexMapping[i]);

                double rp = Math.exp(vPar);
                double rc = Math.exp(vChi);

                double dt = parent.getHeight() - node.getHeight();
                if (dt < 0.0) {
                    dt = 0.0;
                }
                // eZ = mean rate along branch from bridging
                double eZ = MeanZCalculator.computeMeanZ(rp, rc, dt, phi, taylorOrder);

                branchRates[i] = eZ * globalRateVal;
            }
        }
    }

    /**
     * If doNormalize = true, we scale so that the average rate across the tree = 1.0
     */
    private void computeScaleFactor() {
        double sumRate = 0.0;
        double sumTime = 0.0;
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            if (!node.isRoot()) {
                double dt = node.getParent().getHeight() - node.getHeight();
                if (dt < 0) dt = 0;
                sumRate += branchRates[i] * dt;
                sumTime += dt;
            }
        }
        if (sumRate <= 0.0) {
            scaleFactor = 1.0;
        } else {
            scaleFactor = sumTime / sumRate;
        }
    }

    @Override
    protected boolean requiresRecalculation() {

        recompute = true;
        return true;
    }

    @Override
    public void store() {

        System.arraycopy(branchRates, 0, storedBranchRates, 0, nNodes);
        storedScaleFactor = scaleFactor;
        super.store();
    }

    @Override
    public void restore() {
        double[] temp = branchRates;
        branchRates = storedBranchRates;
        storedBranchRates = temp;
        scaleFactor = storedScaleFactor;

        recompute = false;
        renormalize = false;
        super.restore();
    }

    @Override
    public UpDownOperator getUpDownOperator1(Tree tree) {

        UpDownOperator upDownOperator = new UpDownOperator();
        String idStr = getID() + "Up" + tree.getID() + "DownOperator1";
        upDownOperator.setID(idStr);

        upDownOperator.setInputValue("scaleFactor", 0.75);
        upDownOperator.setInputValue("weight", 3.0);

        upDownOperator.setInputValue("down", sigma2Input.get());
        upDownOperator.setInputValue("up", tree);

        upDownOperator.initAndValidate();
        return upDownOperator;
    }

    @Override
    public UpDownOperator getUpDownOperator2(Tree tree) {
        UpDownOperator upDownOperator = new UpDownOperator();
        String idStr = getID() + "Up" + tree.getID() + "DownOperator2";
        upDownOperator.setID(idStr);

        upDownOperator.setInputValue("scaleFactor", 0.75);
        upDownOperator.setInputValue("weight", 3.0);

        upDownOperator.setInputValue("up", rootLogRateInput.get());
        upDownOperator.setInputValue("down", tree);

        upDownOperator.initAndValidate();
        return upDownOperator;
    }


}
