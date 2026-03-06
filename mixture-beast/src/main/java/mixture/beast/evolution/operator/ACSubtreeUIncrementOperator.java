package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Description("AC-mode non-centered subtree block update: random-walk on u-increments within a subtree, "
        + "then reconstruct log-rates and rates. Helps mixing on large trees under AC.")
public class ACSubtreeUIncrementOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);
    public final Input<RealParameter> ratesInput = new Input<>("rates", "shared positive branch rates (non-root nodes)", Input.Validate.REQUIRED);
    public final Input<IntegerParameter> indicatorInput = new Input<>("indicator", "0=UC, 1=AC", Input.Validate.REQUIRED);

    public final Input<RealParameter> sigma2Input = new Input<>("sigma2", "AC Brownian variance per unit time", Input.Validate.REQUIRED);
    public final Input<RealParameter> rootLogRateInput = new Input<>("rootLogRate", "optional root log-rate anchor (default 0)", Input.Validate.OPTIONAL);

    public final Input<Double> minBranchLengthInput = new Input<>("minBranchLength", "min dt allowed", 1e-12);

    public final Input<Double> deltaInput = new Input<>("delta", "stdev of Normal perturbation on u increments", 0.25);

    public final Input<Boolean> internalOnlyInput = new Input<>("internalOnly", "choose only internal nodes as subtree roots", true);
    public final Input<Boolean> allowRootInput = new Input<>("allowRoot", "allow choosing the tree root as subtree root (global block move)", false);

    public final Input<Integer> maxSubtreeEdgesInput = new Input<>("maxSubtreeEdges", "max number of edges (nodes) updated in a move; <=0 means no limit", 80);
    public final Input<Integer> maxTriesInput = new Input<>("maxTries", "tries to find a subtree satisfying constraints", 40);

    /** If true (default), operator returns -inf when indicator!=1. For mixture runs, set autoOptimize=false in XML. */
    public final Input<Boolean> rejectIfNotACInput = new Input<>("rejectIfNotAC", "reject move when indicator!=1", true);

    private Tree tree;
    private RealParameter rates;
    private IntegerParameter indicator;
    private RealParameter sigma2;
    private RealParameter rootLogRate;

    private int nNodes;
    private int rootNr;
    private int[] idxMap;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();
        indicator = indicatorInput.get();
        sigma2 = sigma2Input.get();
        rootLogRate = rootLogRateInput.get();

        if (indicator.getDimension() != 1) throw new IllegalArgumentException("indicator dimension must be 1");
        if (sigma2.getDimension() != 1) throw new IllegalArgumentException("sigma2 dimension must be 1");
        if (rootLogRate != null && rootLogRate.getDimension() != 1) throw new IllegalArgumentException("rootLogRate dimension must be 1");

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();

        if (rates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("rates must have dimension nodeCount-1");
        }

        rebuildIdxMap_likeYourPrior();
    }

    private void rebuildIdxMap_likeYourPrior() {
        idxMap = new int[nNodes];
        Arrays.fill(idxMap, -2);

        int k = 0;
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            int nr = node.getNr();
            if (nr < 0 || nr >= nNodes) throw new IllegalArgumentException("Node nr out of range: " + nr);
            if (nr == rootNr) idxMap[nr] = -1;
            else idxMap[nr] = k++;
        }

        if (k != nNodes - 1) {
            throw new IllegalStateException("Expected " + (nNodes - 1) + " non-root nodes, got " + k);
        }
    }

    private void ensureMapsUpToDate() {
        int newN = tree.getNodeCount();
        int newRoot = tree.getRoot().getNr();
        if (newN != nNodes || newRoot != rootNr) {
            nNodes = newN;
            rootNr = newRoot;
            if (rates.getDimension() != nNodes - 1) throw new IllegalStateException("rates dimension inconsistent with tree");
            rebuildIdxMap_likeYourPrior();
        }
    }

    private double rootLog() {
        return (rootLogRate == null ? 0.0 : rootLogRate.getValue(0));
    }

    /** Collect all non-root nodes (edges) in the subtree under node 'n'. */
    private void collectEdgesInSubtree(Node n, List<Node> out) {
        if (!n.isRoot()) out.add(n);
        final int cc = n.getChildCount();
        for (int i = 0; i < cc; i++) {
            collectEdgesInSubtree(n.getChild(i), out);
        }
    }

    private Node pickSubtreeRoot() {
        final boolean internalOnly = internalOnlyInput.get();
        final boolean allowRoot = allowRootInput.get();
        final int maxEdges = maxSubtreeEdgesInput.get();
        final int tries = maxTriesInput.get();

        for (int t = 0; t < tries; t++) {
            Node cand = tree.getNode(Randomizer.nextInt(nNodes));
            if (cand.isRoot() && !allowRoot) continue;
            if (!allowRoot && cand.isRoot()) continue;

            if (internalOnly && cand.getChildCount() == 0) continue;

            // Check subtree size constraint
            if (maxEdges > 0) {
                ArrayList<Node> tmp = new ArrayList<>();
                collectEdgesInSubtree(cand, tmp);
                if (tmp.size() > maxEdges) continue;
            }
            return cand;
        }

        // fallback: any non-root node
        Node cand;
        do {
            cand = tree.getNode(Randomizer.nextInt(nNodes));
        } while (cand.isRoot() && !allowRoot);
        return cand;
    }

    /** Compute old u for a subtree rooted at 'node' given parent log-rate 'logPar'. */
    private void computeUOld(Node node,
                             double logPar,
                             double sig2,
                             double minDt,
                             double[] xOldByNr,
                             double[] uOldByNr) {
        final int nr = node.getNr();
        final int idx = idxMap[nr];
        if (idx < 0) throw new ArithmeticException("Bad idxMap");

        final double dt = node.getLength();
        if (!(dt > minDt)) throw new ArithmeticException("dt too small");
        final double var = sig2 * dt;
        if (!(var > 0.0)) throw new ArithmeticException("var<=0");

        final double r = rates.getValue(idx);
        if (!(r > 0.0)) throw new ArithmeticException("rate<=0");
        final double x = Math.log(r);

        final double mean = logPar - 0.5 * var;
        final double u = (x - mean) / Math.sqrt(var);

        xOldByNr[nr] = x;
        uOldByNr[nr] = u;

        final int cc = node.getChildCount();
        for (int i = 0; i < cc; i++) {
            computeUOld(node.getChild(i), x, sig2, minDt, xOldByNr, uOldByNr);
        }
    }

    /** Build new x for subtree rooted at 'node' given parent new log-rate 'logParNew'. */
    private void buildXNew(Node node,
                           double logParNew,
                           double sig2,
                           double minDt,
                           double[] uNewByNr,
                           double[] xNewByNr) {
        final int nr = node.getNr();
        final int idx = idxMap[nr];
        if (idx < 0) throw new ArithmeticException("Bad idxMap");

        final double dt = node.getLength();
        if (!(dt > minDt)) throw new ArithmeticException("dt too small");
        final double var = sig2 * dt;
        if (!(var > 0.0)) throw new ArithmeticException("var<=0");

        final double u = uNewByNr[nr];
        final double x = logParNew - 0.5 * var + Math.sqrt(var) * u;

        xNewByNr[nr] = x;

        final int cc = node.getChildCount();
        for (int i = 0; i < cc; i++) {
            buildXNew(node.getChild(i), x, sig2, minDt, uNewByNr, xNewByNr);
        }
    }

    @Override
    public double proposal() {
        ensureMapsUpToDate();

        final int k = indicator.getValue(0);
        if (k != 1) {
            return rejectIfNotACInput.get() ? Double.NEGATIVE_INFINITY : 0.0;
        }

        final double sig2 = sigma2.getValue(0);
        final double minDt = minBranchLengthInput.get();
        final double delta = deltaInput.get();
        if (!(sig2 > 0.0) || !(minDt > 0.0) || !(delta > 0.0)) return Double.NEGATIVE_INFINITY;

        final Node subRoot = pickSubtreeRoot();

        // edges (non-root nodes) to update
        final ArrayList<Node> edges = new ArrayList<>();
        collectEdgesInSubtree(subRoot, edges);
        if (edges.isEmpty()) return Double.NEGATIVE_INFINITY;

        final double[] xOldByNr = new double[nNodes];
        final double[] uOldByNr = new double[nNodes];
        final double[] uNewByNr = new double[nNodes];
        final double[] xNewByNr = new double[nNodes];
        Arrays.fill(xOldByNr, Double.NaN);
        Arrays.fill(uOldByNr, Double.NaN);
        Arrays.fill(uNewByNr, Double.NaN);
        Arrays.fill(xNewByNr, Double.NaN);

        // boundary parent log-rate for the subtree root edge
        double boundaryLogPar;
        if (subRoot.isRoot()) {
            boundaryLogPar = rootLog();
            // root itself has no edge; compute for children
            final int cc = subRoot.getChildCount();
            try {
                for (int i = 0; i < cc; i++) {
                    computeUOld(subRoot.getChild(i), boundaryLogPar, sig2, minDt, xOldByNr, uOldByNr);
                }
            } catch (ArithmeticException bad) {
                return Double.NEGATIVE_INFINITY;
            }
        } else {
            final Node parent = subRoot.getParent();
            if (parent.isRoot()) {
                boundaryLogPar = rootLog();
            } else {
                final int idxPar = idxMap[parent.getNr()];
                if (idxPar < 0) return Double.NEGATIVE_INFINITY;
                final double rPar = rates.getValue(idxPar);
                if (!(rPar > 0.0)) return Double.NEGATIVE_INFINITY;
                boundaryLogPar = Math.log(rPar);
            }
            try {
                computeUOld(subRoot, boundaryLogPar, sig2, minDt, xOldByNr, uOldByNr);
            } catch (ArithmeticException bad) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Propose new u (RW) for all edges in the chosen subtree
        for (Node e : edges) {
            final int nr = e.getNr();
            final double uOld = uOldByNr[nr];
            if (!Double.isFinite(uOld)) return Double.NEGATIVE_INFINITY;
            uNewByNr[nr] = uOld + delta * Randomizer.nextGaussian();
        }

        // Rebuild new x from boundary using proposed u
        try {
            if (subRoot.isRoot()) {
                final int cc = subRoot.getChildCount();
                final double rootLog = rootLog();
                for (int i = 0; i < cc; i++) {
                    buildXNew(subRoot.getChild(i), rootLog, sig2, minDt, uNewByNr, xNewByNr);
                }
            } else {
                buildXNew(subRoot, boundaryLogPar, sig2, minDt, uNewByNr, xNewByNr);
            }
        } catch (ArithmeticException bad) {
            return Double.NEGATIVE_INFINITY;
        }

        // Hastings ratio in rate-space for symmetric RW in u:
        // logH = sum_{updated edges} (xNew - xOld)
        double sumDelta = 0.0;
        for (Node e : edges) {
            final int nr = e.getNr();
            final double xo = xOldByNr[nr];
            final double xn = xNewByNr[nr];
            if (!Double.isFinite(xo) || !Double.isFinite(xn)) return Double.NEGATIVE_INFINITY;
            sumDelta += (xn - xo);
        }

        // Commit new rates
        rates.startEditing(this);
        for (Node e : edges) {
            final int idx = idxMap[e.getNr()];
            rates.setValue(idx, Math.exp(xNewByNr[e.getNr()]));
        }

        return sumDelta;
    }
}
