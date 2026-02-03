package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.util.ArrayDeque;
import java.util.Arrays;

@Description("Scale (multiply) all rates in a randomly chosen subtree by a common factor exp(eps). "
        + "Useful for mixing when rates are shared and AC prior induces correlations along the tree.")
public class SubtreeRateScaleOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "positive branch rates for NON-root nodes; dimension=(tree.nodeCount - 1)",
            Input.Validate.REQUIRED
    );

    public final Input<Double> windowInput = new Input<>(
            "window",
            "proposal window on log scale: eps ~ Uniform(-window, +window)",
            0.5
    );

    private Tree tree;
    private RealParameter rates;

    private int nNodes;
    private int rootNr;
    private int[] idxMap;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();

        if (rates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("SubtreeRateScaleOperator: rates must have dimension nodeCount-1.");
        }

        buildIndexMap();
    }

    private void buildIndexMap() {
        idxMap = new int[nNodes];
        Arrays.fill(idxMap, -2);

        int k = 0;
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            int nr = node.getNr();
            if (nr == rootNr) {
                idxMap[nr] = -1;
            } else {
                idxMap[nr] = k++;
            }
        }
    }

    @Override
    public double proposal() {
        final double window = windowInput.get();
        if (!(window > 0.0)) return Double.NEGATIVE_INFINITY;

        // pick a random non-root node as subtree root
        Node subRoot;
        do {
            subRoot = tree.getNode(Randomizer.nextInt(nNodes));
        } while (subRoot.isRoot());

        // log multiplier
        final double eps = (Randomizer.nextDouble() * 2.0 - 1.0) * window;
        final double m = Math.exp(eps);

        // traverse subtree, scale all non-root node rates
        int count = 0;
        ArrayDeque<Node> stack = new ArrayDeque<>();
        stack.push(subRoot);

        while (!stack.isEmpty()) {
            Node node = stack.pop();

            if (!node.isRoot()) {
                final int idx = idxMap[node.getNr()];
                if (idx < 0) return Double.NEGATIVE_INFINITY;

                final double r = rates.getValue(idx);
                final double rNew = r * m;
                if (!(rNew > 0.0) || Double.isInfinite(rNew) || Double.isNaN(rNew)) {
                    return Double.NEGATIVE_INFINITY;
                }
                rates.setValue(idx, rNew);
                count++;
            }

            for (int i = 0; i < node.getChildCount(); i++) {
                stack.push(node.getChild(i));
            }
        }

        // Hastings ratio for scaling count dimensions by multiplier m:
        // logHR = count * log(m) = count * eps
        return count * eps;
    }
}
