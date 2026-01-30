//package rc.beast.evolution.clockmodel;
//
//import beast.base.core.Input;
//import beast.base.evolution.tree.Node;
//import beast.base.evolution.tree.Tree;
//import beast.base.inference.Operator;
//import beast.base.inference.parameter.RealParameter;
//import beast.base.util.Randomizer;
//import java.util.ArrayDeque;
//
//public class TreeLogRateScaleOperator extends Operator {
//
//    public final Input<Tree> treeInput = new Input<>("tree", "the tree to scale", Input.Validate.REQUIRED);
//    public final Input<RealParameter> rootLogRateInput = new Input<>("rootLogRate", "root log-rate (dim=1)", Input.Validate.REQUIRED);
//    public final Input<RealParameter> nodeLogRatesInput = new Input<>("nodeRates", "log-rates for non-root nodes", Input.Validate.REQUIRED);
//
//    public final Input<RealParameter> meanRateInput = new Input<>(
//            "meanRate",
//            "optional global mean rate (dim=1). Used when scaleMeanRate=true.",
//            Input.Validate.OPTIONAL
//    );
//
//    public final Input<RealParameter> sigma2Input = new Input<>(
//            "sigma2",
//            "optional; per-time variance parameter (dim=1)",
//            Input.Validate.OPTIONAL
//    );
//
//    public final Input<RealParameter> lambdaInput = new Input<>(
//            "lambda",
//            "optional; per-time rate parameter (dim=1)",
//            Input.Validate.OPTIONAL
//    );
//
//    public final Input<Boolean> shiftLogRatesInput = new Input<>(
//            "shiftLogRates",
//            "if true, shift all log-rates by log(g) so rates multiply by g",
//            true
//    );
//
//    public final Input<Boolean> scaleMeanRateInput = new Input<>(
//            "scaleMeanRate",
//            "if true, meanRate <- meanRate * g (recommended when normalize=true)",
//            false
//    );
//
//    public final Input<Boolean> useTreeLengthInput = new Input<>(
//            "useTreeLength",
//            "if true, g = oldTreeLength/newTreeLength; else g = 1/s",
//            true
//    );
//
//    public final Input<Boolean> scaleSigma2Input = new Input<>("scaleSigma2", "if true, scale sigma2", true);
//    public final Input<Boolean> scaleLambdaInput = new Input<>("scaleLambda", "if true, scale lambda", true);
//
//    public final Input<Boolean> hyperUseTreeLengthInput = new Input<>(
//            "hyperUseTreeLength",
//            "if true, sigma2/lambda are scaled by oldTreeLength/newTreeLength; else by 1/s",
//            true
//    );
//
//    public final Input<Double> scaleFactorInput = new Input<>(
//            "scaleFactor",
//            "scaleFactor in (0,1). Proposal uses log-uniform scaling: s = exp(u), u ~ Uniform(-log(1/sf), +log(1/sf)).",
//            0.75
//    );
//
//    private Tree tree;
//    private RealParameter rootLogRate;
//    private RealParameter nodeLogRates;
//    private RealParameter meanRate;
//    private RealParameter sigma2;
//    private RealParameter lambda;
//
//    private boolean doShiftLogRates;
//    private boolean doScaleMeanRate;
//    private boolean useTreeLength;
//    private boolean doScaleSigma2;
//    private boolean doScaleLambda;
//    private boolean hyperUseTreeLength;
//
//    @Override
//    public void initAndValidate() {
//        tree = treeInput.get();
//        rootLogRate = rootLogRateInput.get();
//        nodeLogRates = nodeLogRatesInput.get();
//        meanRate = meanRateInput.get();
//        sigma2 = sigma2Input.get();
//        lambda = lambdaInput.get();
//
//        doShiftLogRates = shiftLogRatesInput.get();
//        doScaleMeanRate = scaleMeanRateInput.get();
//        useTreeLength = useTreeLengthInput.get();
//        doScaleSigma2 = scaleSigma2Input.get();
//        doScaleLambda = scaleLambdaInput.get();
//        hyperUseTreeLength = hyperUseTreeLengthInput.get();
//
//        final double sf = scaleFactorInput.get();
//        if (!(sf > 0.0) || !(sf < 1.0) || Double.isNaN(sf) || Double.isInfinite(sf)) {
//            throw new IllegalArgumentException("scaleFactor must be in (0,1) and finite. Found: " + sf);
//        }
//
//        if (rootLogRate.getDimension() != 1) {
//            throw new IllegalArgumentException("rootLogRate must have dimension=1.");
//        }
//        if (nodeLogRates.getDimension() < 1) {
//            throw new IllegalArgumentException("nodeRates must have dimension >= 1.");
//        }
//        if (meanRate != null && meanRate.getDimension() != 1) {
//            throw new IllegalArgumentException("meanRate must have dimension=1 if provided.");
//        }
//        if (sigma2 != null && sigma2.getDimension() != 1) {
//            throw new IllegalArgumentException("sigma2 must have dimension=1 if provided.");
//        }
//        if (lambda != null && lambda.getDimension() != 1) {
//            throw new IllegalArgumentException("lambda must have dimension=1 if provided.");
//        }
//        if (doScaleMeanRate && meanRate == null) {
//            throw new IllegalArgumentException("scaleMeanRate=true but meanRate was not provided.");
//        }
//    }
//
//    protected double drawLogS(final double logRange) {
//        return (Randomizer.nextDouble() - 0.5) * 2.0 * logRange;
//    }
//
//    @Override
//    public double proposal() {
//        final double sf = scaleFactorInput.get();
//        final double logRange = Math.log(1.0 / sf);
//        final double logS = drawLogS(logRange);
//
//        final double s = Math.exp(logS);
//
//        final int n = tree.getNodeCount();
//        final double[] oldIntervals = new double[n];
//
//        int k = 0;
//        for (int i = 0; i < n; i++) {
//            final Node node = tree.getNode(i);
//            if (node.isLeaf()) continue;
//
//            double maxChild = Double.NEGATIVE_INFINITY;
//            final int cc = node.getChildCount();
//            for (int c = 0; c < cc; c++) {
//                final double h = node.getChild(c).getHeight();
//                if (h > maxChild) maxChild = h;
//            }
//
//            double interval = node.getHeight() - maxChild;
//            if (interval < 0.0) {
//                if (interval < -1e-10) {
//                    return Double.NEGATIVE_INFINITY;
//                }
//                interval = 0.0;
//            }
//
//            oldIntervals[node.getNr()] = interval;
//            k++;
//        }
//
//        final double oldLen = treeLengthStrict(tree);
//        if (!(oldLen > 0.0)) return Double.NEGATIVE_INFINITY;
//
//        tree.startEditing(this);
//
//        if (doShiftLogRates) {
//            rootLogRate.startEditing(this);
//            nodeLogRates.startEditing(this);
//        }
//        if (doScaleMeanRate && meanRate != null) {
//            meanRate.startEditing(this);
//        }
//        if (sigma2 != null && doScaleSigma2) {
//            sigma2.startEditing(this);
//        }
//        if (lambda != null && doScaleLambda) {
//            lambda.startEditing(this);
//        }
//
//        final Node[] post = buildPostOrder(tree.getRoot(), n);
//        for (Node node : post) {
//            if (node.isLeaf()) continue;
//
//            double maxChild = Double.NEGATIVE_INFINITY;
//            final int cc = node.getChildCount();
//            for (int c = 0; c < cc; c++) {
//                final double h = node.getChild(c).getHeight();
//                if (h > maxChild) maxChild = h;
//            }
//
//            final double newH = maxChild + s * oldIntervals[node.getNr()];
//            node.setHeight(newH);
//        }
//
//        final double newLen = treeLengthStrict(tree);
//        if (!(newLen > 0.0)) return Double.NEGATIVE_INFINITY;
//
//        final double g = useTreeLength ? (oldLen / newLen) : (1.0 / s);
//        if (!(g > 0.0) || Double.isNaN(g) || Double.isInfinite(g)) return Double.NEGATIVE_INFINITY;
//        final double logG = Math.log(g);
//
//        final double hyperScale = hyperUseTreeLength ? (oldLen / newLen) : (1.0 / s);
//        if (!(hyperScale > 0.0) || Double.isNaN(hyperScale) || Double.isInfinite(hyperScale)) {
//            return Double.NEGATIVE_INFINITY;
//        }
//        final double logHyper = Math.log(hyperScale);
//
//        double logHastings = k * logS;
//
//        if (doShiftLogRates) {
//            rootLogRate.setValue(0, rootLogRate.getValue(0) + logG);
//
//            final int dim = nodeLogRates.getDimension();
//            for (int i = 0; i < dim; i++) {
//                nodeLogRates.setValue(i, nodeLogRates.getValue(i) + logG);
//            }
//        }
//
//        if (doScaleMeanRate && meanRate != null) {
//            meanRate.setValue(0, meanRate.getValue(0) * g);
//            logHastings += logG;
//        }
//
//        if (sigma2 != null && doScaleSigma2) {
//            sigma2.setValue(0, sigma2.getValue(0) * hyperScale);
//            logHastings += logHyper;
//        }
//        if (lambda != null && doScaleLambda) {
//            lambda.setValue(0, lambda.getValue(0) * hyperScale);
//            logHastings += logHyper;
//        }
//
//        return logHastings;
//    }
//
//    private static double treeLengthStrict(final Tree tree) {
//        double sum = 0.0;
//        final int n = tree.getNodeCount();
//        for (int i = 0; i < n; i++) {
//            final Node node = tree.getNode(i);
//            if (node.isRoot()) continue;
//            final Node parent = node.getParent();
//            final double bl = parent.getHeight() - node.getHeight();
//            if (bl < -1e-10) return Double.NaN;
//            if (bl > 0.0) sum += bl;
//        }
//        return sum;
//    }
//
//    private static Node[] buildPostOrder(final Node root, final int nodeCount) {
//        final ArrayDeque<Node> s1 = new ArrayDeque<>();
//        final ArrayDeque<Node> s2 = new ArrayDeque<>();
//        s1.push(root);
//
//        while (!s1.isEmpty()) {
//            final Node node = s1.pop();
//            s2.push(node);
//            final int cc = node.getChildCount();
//            for (int c = 0; c < cc; c++) {
//                s1.push(node.getChild(c));
//            }
//        }
//
//        final Node[] out = new Node[nodeCount];
//        int idx = 0;
//        while (!s2.isEmpty() && idx < nodeCount) {
//            out[idx++] = s2.pop();
//        }
//        return out;
//    }
//}
