package rc.beast.evolution.clockmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.parameter.RealParameter;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * AutoCorrelatedPrior:
 *  - Imposes (vChild - vParent) ~ Normal(0, sigma^2 * dt).
 *  - Also imposes root log-rate prior: root ~ Normal(rootMean, rootStdev^2).
 */
@Description("Auto-correlated prior for node log-rates in a Brownian sense.")
public class AutoCorrelatedPrior extends Distribution {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "the phylogenetic tree",
            Input.Validate.REQUIRED
    );

    public Input<RealParameter> nodeLogRatesInput = new Input<>(
            "nodeRates",
            "log(rate) at each node",
            Input.Validate.REQUIRED
    );

    public Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "the Brownian variance param",
            Input.Validate.REQUIRED
    );

    final public Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "the log-rate at the root (dimension=1)",
            Input.Validate.REQUIRED
    );

    protected Tree tree;
    protected RealParameter nodeLogRates;
    protected RealParameter sigma2;
    protected RealParameter rootLogRate;
    protected int nNodes;                 // total number of nodes in the tree
    protected int rootNr;                 // index of root node
    protected int[] nodeIndexMapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        nodeLogRates = nodeLogRatesInput.get();
        rootLogRate = rootLogRateInput.get();
        sigma2 = sigma2Input.get();

        nNodes = tree.getNodeCount();
        Node rootNode = tree.getRoot();
        rootNr = rootNode.getNr();

        if (rootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate must have dimension=1");
        }

        if (nodeLogRates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException(
                    "nodeRates dimension must be (tree.nodeCount - 1). Found: "
                            + nodeLogRates.getDimension() + " vs " + (nNodes - 1)
            );
        }

        // Build an index mapping for nodeLogRates
        nodeIndexMapping = new int[nNodes];
        int idx = 0;
        for (int i = 0; i < nNodes; i++) {
            if (i == rootNr) {
                // Root node does not appear in nodeLogRates
                nodeIndexMapping[i] = -1;
            } else {
                nodeIndexMapping[i] = idx;
                idx++;
            }
        }
    }

    @Override
    public double calculateLogP() {
        // We reset logP
        logP = 0.0;

        // sigma^2 parameter
        double s2 = sigma2.getValue();

        // Loop over all nodes in the tree
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);

            // If this node is NOT the root, then we compute the increment from parent
            if (!node.isRoot()) {
                Node parent = node.getParent();
                int parentNr = parent.getNr();

                // vPar = parent's log-rate
                double vPar;
                if (parent.isRoot()) {
                    // If the parent is the actual root node
                    // then we use rootLogRate
                    vPar = rootLogRate.getValue();
                } else {
                    // Otherwise, the parent's log-rate is in nodeLogRates
                    int idxPar = nodeIndexMapping[parentNr];
                    vPar = nodeLogRates.getValue(idxPar);
                }

                // vChi = child's log-rate
                int idxChild = nodeIndexMapping[i];
                double vChi = nodeLogRates.getValue(idxChild);

                // dt = parent.height - child.height
                double dt = parent.getHeight() - node.getHeight();
                if (dt < 0.0) {
                    // If dt is negative, we consider that invalid
                    logP = Double.NEGATIVE_INFINITY;
                    return logP;
                }

                // var = sigma^2 * dt
                double var = s2 * dt;
                if (var <= 0.0) {
                    logP = Double.NEGATIVE_INFINITY;
                    return logP;
                }

                double diff = vChi - vPar;

                double thisLP = -0.5 * (
                        Math.log(2.0 * Math.PI * var) +
                                (diff * diff) / var
                );

                // Accumulate into logP
                logP += thisLP;
            }
        }

        return logP;
    }

    @Override
    public List<String> getArguments() {
        // these are the parameters that get log-density from this distribution
        ArrayList<String> args = new ArrayList<>();
        args.add(nodeLogRates.getID());
        args.add(rootLogRate.getID());
        args.add(sigma2.getID());
        return args;
    }

    @Override
    public List<String> getConditions() {
        // the tree is a condition
        ArrayList<String> conds = new ArrayList<>();
        conds.add(tree.getID());
        return conds;
    }

    @Override
    public void sample(beast.base.inference.State state, Random random) {
        // Typically we let MCMC operator do proposals, so no direct sampling needed.
    }
}
