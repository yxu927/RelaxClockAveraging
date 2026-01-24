package rc.beast.evolution.clockmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

@Description("Brownian increment prior on log-rates: (vChild - vParent) ~ N(0, sigma2*dt). "
        + "This distribution does NOT include a root prior. Add a separate root prior externally if needed.")
public class AutoCorrelatedPrior extends Distribution {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public final Input<RealParameter> nodeLogRatesInput = new Input<>(
            "nodeRates",
            "log(rate) for NON-root nodes, dimension = (tree.nodeCount - 1)",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "root log-rate, dimension=1 (use fixed 0.0 in normalized mode)",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "Brownian variance parameter (sigma^2)",
            Input.Validate.REQUIRED
    );

    protected Tree tree;
    protected RealParameter nodeLogRates;
    protected RealParameter rootLogRate;
    protected RealParameter sigma2;

    protected int nNodes;
    protected int rootNr;
    protected int[] nodeIndexMapping;

    @Override
    public void initAndValidate() {

        Log.warning.println("AutoCorrelatedPrior loaded from: "
                + getClass().getProtectionDomain().getCodeSource().getLocation());

        tree = treeInput.get();
        nodeLogRates = nodeLogRatesInput.get();
        rootLogRate = rootLogRateInput.get();
        sigma2 = sigma2Input.get();

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();

        if (rootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate must have dimension=1.");
        }
        if (nodeLogRates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("nodeRates must have dimension (nodeCount - 1). Found "
                    + nodeLogRates.getDimension() + " vs " + (nNodes - 1));
        }

        nodeIndexMapping = new int[nNodes];
        int idx = 0;
        for (int i = 0; i < nNodes; i++) {
            int nr = tree.getNode(i).getNr();
            if (nr == rootNr) {
                nodeIndexMapping[nr] = -1;
            } else {
                nodeIndexMapping[nr] = idx++;
            }
        }
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        final double s2 = sigma2.getValue();
        if (s2 <= 0.0) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        final double tiny = 1e-12;

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

            double dt = parent.getHeight() - node.getHeight();
            if (dt < 0.0) dt = 0.0;

            final double var = s2 * dt;

            if (var <= tiny) {
                // deterministic limit: vChi must equal vPar
                if (Math.abs(vChi - vPar) > 1e-9) {
                    logP = Double.NEGATIVE_INFINITY;
                    return logP;
                }
                continue;
            }

            final double diff = vChi - vPar;
            logP += -0.5 * (Math.log(2.0 * Math.PI * var) + (diff * diff) / var);
        }

        return logP;
    }

    @Override
    public List<String> getArguments() {
        ArrayList<String> args = new ArrayList<>();
        args.add(nodeLogRates.getID());
        return args;
    }

    @Override
    public List<String> getConditions() {
        ArrayList<String> conds = new ArrayList<>();
        conds.add(tree.getID());
        conds.add(rootLogRate.getID());
        conds.add(sigma2.getID());
        return conds;
    }

    @Override
    public void sample(beast.base.inference.State state, Random random) {
    }
}
