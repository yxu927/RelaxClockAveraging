package rc.beast.evolution.clockmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Coupled move: scale tree heights by s, shift all log-rates by -log(s). "
        + "Optionally scale sigma2 and lambda by 1/s to match time units.")
public class TreeLogRateScaleOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "the tree to scale", Input.Validate.REQUIRED);
    public final Input<RealParameter> rootLogRateInput = new Input<>("rootLogRate", "root log-rate (dim=1)", Input.Validate.REQUIRED);
    public final Input<RealParameter> nodeLogRatesInput = new Input<>("nodeRates", "log-rates for non-root nodes", Input.Validate.REQUIRED);

    public final Input<RealParameter> sigma2Input = new Input<>("sigma2", "optional; if scaled then sigma2 <- sigma2 / s", Input.Validate.OPTIONAL);
    public final Input<RealParameter> lambdaInput = new Input<>("lambda", "optional; if scaled then lambda <- lambda / s", Input.Validate.OPTIONAL);

    public final Input<Boolean> scaleSigma2Input = new Input<>("scaleSigma2", "if true, sigma2 <- sigma2 / s", true);
    public final Input<Boolean> scaleLambdaInput = new Input<>("scaleLambda", "if true, lambda <- lambda / s", true);

    public final Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor",
            "scaleFactor in (0,1). Proposal uses log-uniform scaling: s = exp(u), u ~ Uniform(-log(1/sf), +log(1/sf)).",
            0.75
    );

    private Tree tree;
    private RealParameter rootLogRate;
    private RealParameter nodeLogRates;
    private RealParameter sigma2;
    private RealParameter lambda;
    private boolean doScaleSigma2;
    private boolean doScaleLambda;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rootLogRate = rootLogRateInput.get();
        nodeLogRates = nodeLogRatesInput.get();
        sigma2 = sigma2Input.get();
        lambda = lambdaInput.get();

        doScaleSigma2 = scaleSigma2Input.get();
        doScaleLambda = scaleLambdaInput.get();

        if (rootLogRate.getDimension() != 1) throw new IllegalArgumentException("rootLogRate must have dimension=1.");
        if (sigma2 != null && sigma2.getDimension() != 1) throw new IllegalArgumentException("sigma2 must have dimension=1.");
        if (lambda != null && lambda.getDimension() != 1) throw new IllegalArgumentException("lambda must have dimension=1.");
    }

    @Override
    public double proposal() {
        final double sf = scaleFactorInput.get();
        if (sf <= 0.0 || sf >= 1.0) {
            throw new IllegalArgumentException("scaleFactor must be in (0,1). Found: " + sf);
        }

        final double logRange = Math.log(1.0 / sf);
        final double logS = (Randomizer.nextDouble() - 0.5) * 2.0 * logRange;
        final double s = Math.exp(logS);

        tree.startEditing(this);
        rootLogRate.startEditing(this);
        nodeLogRates.startEditing(this);
        if (sigma2 != null && doScaleSigma2) sigma2.startEditing(this);
        if (lambda != null && doScaleLambda) lambda.startEditing(this);


        int k = 0;
        final int n = tree.getNodeCount();
        for (int i = 0; i < n; i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) continue;
            node.setHeight(node.getHeight() * s);
            k++;
        }

        // Shift all log-rates by -log(s) so rates scale by 1/s
        rootLogRate.setValue(0, rootLogRate.getValue() - logS);
        final int dim = nodeLogRates.getDimension();
        for (int i = 0; i < dim; i++) {
            nodeLogRates.setValue(i, nodeLogRates.getValue(i) - logS);
        }

        double logHastings = k * logS;

        if (sigma2 != null && doScaleSigma2) {
            sigma2.setValue(0, sigma2.getValue() / s);
            logHastings -= logS;
        }
        if (lambda != null && doScaleLambda) {
            lambda.setValue(0, lambda.getValue() / s);
            logHastings -= logS;
        }

        return logHastings;
    }
}
