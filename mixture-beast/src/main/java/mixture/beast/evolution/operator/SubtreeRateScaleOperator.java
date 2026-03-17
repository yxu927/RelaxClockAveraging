package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import mixture.beast.evolution.util.BranchRateIndexHelper;

import java.util.ArrayDeque;

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
    private BranchRateIndexHelper.Mapping mapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();

        BranchRateIndexHelper.validateRatesDimension(tree, rates, "SubtreeRateScaleOperator");
        mapping = BranchRateIndexHelper.buildDeterministic(tree);
    }

    private void ensureMappingUpToDate() {
        mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, rates, "SubtreeRateScaleOperator");
    }

    @Override
    public double proposal() {
        ensureMappingUpToDate();

        final double window = windowInput.get();
        if (!(window > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final int nNodes = mapping.getNodeCount();

        Node subRoot;
        do {
            subRoot = tree.getNode(Randomizer.nextInt(nNodes));
        } while (subRoot.isRoot());

        final double eps = (Randomizer.nextDouble() * 2.0 - 1.0) * window;
        final double m = Math.exp(eps);

        int count = 0;
        final ArrayDeque<Node> stack = new ArrayDeque<>();
        stack.push(subRoot);

        while (!stack.isEmpty()) {
            final Node node = stack.pop();

            if (!node.isRoot()) {
                final int idx = mapping.idxForNode(node);
                if (idx < 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                final double r = rates.getValue(idx);
                final double rNew = r * m;
                if (!(rNew > 0.0) || Double.isInfinite(rNew) || Double.isNaN(rNew)) {
                    return Double.NEGATIVE_INFINITY;
                }
                rates.setValue(idx, rNew);
                count++;
            }

            final int childCount = node.getChildCount();
            for (int i = 0; i < childCount; i++) {
                stack.push(node.getChild(i));
            }
        }

        return count * eps;
    }
}
