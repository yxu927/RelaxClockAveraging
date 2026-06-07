package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import mixture.beast.evolution.util.BranchRateIndexHelper;

import java.util.ArrayDeque;

@Description("Scale (multiply) all rates in a randomly chosen subtree by a common factor exp(eps). "
        + "Useful for mixing when rates are shared and AC prior induces correlations along the tree.")
public class SubtreeRateScaleOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "Legacy positive branch rates for NON-root nodes; dimension=(tree.nodeCount - 1).",
            Input.Validate.OPTIONAL
    );

    public final Input<RealVectorParam<?>> ratesVectorInput = new Input<>(
            "ratesVector",
            "BEAST3 typed mutable positive branch rates for NON-root nodes; dimension=(tree.nodeCount - 1).",
            Input.Validate.OPTIONAL
    );

    public final Input<Double> windowInput = new Input<>(
            "window",
            "proposal window on log scale: eps ~ Uniform(-window, +window)",
            0.5
    );

    private Tree tree;
    private RealParameter legacyRates;
    private RealVectorParam<?> typedRates;
    private BranchRateIndexHelper.Mapping mapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        legacyRates = ratesInput.get();
        typedRates = ratesVectorInput.get();

        if (legacyRates == null && typedRates == null) {
            throw new IllegalArgumentException("SubtreeRateScaleOperator: either rates or ratesVector must be specified.");
        }
        if (legacyRates != null && typedRates != null) {
            throw new IllegalArgumentException("SubtreeRateScaleOperator: specify only one of rates or ratesVector.");
        }
        validateOrExpandRatesDimension();
        mapping = BranchRateIndexHelper.buildDeterministic(tree);
    }

    private void ensureMappingUpToDate() {
        validateOrExpandRatesDimension();
        if (mapping == null || !mapping.matches(tree)) {
            mapping = BranchRateIndexHelper.buildDeterministic(tree);
        }
    }

    private int expectedRateDimension() {
        return tree.getNodeCount() - 1;
    }

    private int rateDimension() {
        return legacyRates != null ? legacyRates.getDimension() : typedRates.size();
    }

    private double rateValue(final int i) {
        return legacyRates != null ? legacyRates.getValue(i) : typedRates.get(i);
    }

    private void setRateValue(final int i, final double value) {
        if (legacyRates != null) {
            legacyRates.setValue(i, value);
        } else {
            typedRates.set(i, value);
        }
    }

    private void validateOrExpandRatesDimension() {
        if (legacyRates != null) {
            BranchRateIndexHelper.validateRatesDimension(tree, legacyRates, "SubtreeRateScaleOperator");
            return;
        }

        final int expected = expectedRateDimension();
        if (expected < 1) {
            throw new IllegalArgumentException("SubtreeRateScaleOperator: tree must contain at least one non-root branch. "
                    + "Found nodeCount=" + tree.getNodeCount());
        }

        final int observed = rateDimension();
        if (observed == expected) {
            return;
        }
        if (observed == 1) {
            final double initialValue = typedRates.get(0);
            typedRates.setDimension(expected);
            for (int i = 0; i < expected; i++) {
                typedRates.set(i, initialValue);
            }
            return;
        }

        throw new IllegalArgumentException("SubtreeRateScaleOperator: ratesVector must have dimension (nodeCount - 1). Found "
                + observed + " vs " + expected
                + ". Initialise typed shared rates as a scalar so they can be expanded automatically.");
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

                final double r = rateValue(idx);
                final double rNew = r * m;
                if (!(rNew > 0.0) || Double.isInfinite(rNew) || Double.isNaN(rNew)) {
                    return Double.NEGATIVE_INFINITY;
                }
                setRateValue(idx, rNew);
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
