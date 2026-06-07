package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import mixture.beast.evolution.util.BranchRateIndexHelper;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Description("AC-mode non-centered subtree block update: random-walk on u-increments within a subtree, "
        + "then reconstruct log-rates and rates. Helps mixing on large trees under AC.")
public class ACSubtreeUIncrementOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);
    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "Legacy shared positive branch rates (non-root nodes).",
            Input.Validate.OPTIONAL);
    public final Input<RealVectorParam<?>> ratesVectorInput = new Input<>(
            "ratesVector",
            "BEAST3 typed mutable shared positive branch rates (non-root nodes).",
            Input.Validate.OPTIONAL);
    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator",
            "Legacy scalar indicator; 0=UC, 1=AC.",
            Input.Validate.OPTIONAL);
    public final Input<IntScalarParam<?>> indicatorScalarInput = new Input<>(
            "indicatorScalar",
            "BEAST3 typed mutable scalar indicator; 0=UC, 1=AC.",
            Input.Validate.OPTIONAL);

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "Legacy AC Brownian variance per unit time.",
            Input.Validate.OPTIONAL);
    public final Input<RealScalarParam<?>> sigma2ScalarInput = new Input<>(
            "sigma2Scalar",
            "BEAST3 typed mutable AC Brownian variance per unit time.",
            Input.Validate.OPTIONAL);
    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "Legacy optional root log-rate anchor (default 0).",
            Input.Validate.OPTIONAL);
    public final Input<RealScalarParam<?>> rootLogRateScalarInput = new Input<>(
            "rootLogRateScalar",
            "BEAST3 typed optional root log-rate anchor (default 0).",
            Input.Validate.OPTIONAL);

    public final Input<Double> minBranchLengthInput = new Input<>("minBranchLength", "min dt allowed", 1e-12);

    public final Input<Double> deltaInput = new Input<>("delta", "stdev of Normal perturbation on u increments", 0.25);

    public final Input<Boolean> internalOnlyInput = new Input<>("internalOnly", "choose only internal nodes as subtree roots", true);
    public final Input<Boolean> allowRootInput = new Input<>("allowRoot", "allow choosing the tree root as subtree root (global block move)", false);

    public final Input<Integer> maxSubtreeEdgesInput = new Input<>("maxSubtreeEdges", "max number of edges (nodes) updated in a move; <=0 means no limit", 80);
    public final Input<Integer> maxTriesInput = new Input<>("maxTries", "tries to find a subtree satisfying constraints", 40);

    /** If true (default), operator returns -inf when indicator!=1. For mixture runs, set autoOptimize=false in XML. */
    public final Input<Boolean> rejectIfNotACInput = new Input<>("rejectIfNotAC", "reject move when indicator!=1", true);

    private Tree tree;
    private RealParameter legacyRates;
    private RealVectorParam<?> typedRates;
    private IntegerParameter legacyIndicator;
    private IntScalarParam<?> typedIndicator;
    private RealParameter legacySigma2;
    private RealScalarParam<?> typedSigma2;
    private RealParameter legacyRootLogRate;
    private RealScalarParam<?> typedRootLogRate;

    private BranchRateIndexHelper.Mapping mapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        legacyRates = ratesInput.get();
        typedRates = ratesVectorInput.get();
        legacyIndicator = indicatorInput.get();
        typedIndicator = indicatorScalarInput.get();
        legacySigma2 = sigma2Input.get();
        typedSigma2 = sigma2ScalarInput.get();
        legacyRootLogRate = rootLogRateInput.get();
        typedRootLogRate = rootLogRateScalarInput.get();

        requireExactlyOne(legacyRates, typedRates, "rates", "ratesVector");
        requireExactlyOne(legacyIndicator, typedIndicator, "indicator", "indicatorScalar");
        requireExactlyOne(legacySigma2, typedSigma2, "sigma2", "sigma2Scalar");
        requireAtMostOne(legacyRootLogRate, typedRootLogRate, "rootLogRate", "rootLogRateScalar");

        if (legacyIndicator != null && legacyIndicator.getDimension() != 1) {
            throw new IllegalArgumentException("ACSubtreeUIncrementOperator: indicator dimension must be 1");
        }
        if (legacySigma2 != null && legacySigma2.getDimension() != 1) {
            throw new IllegalArgumentException("ACSubtreeUIncrementOperator: sigma2 dimension must be 1");
        }
        if (legacyRootLogRate != null && legacyRootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("ACSubtreeUIncrementOperator: rootLogRate dimension must be 1");
        }

        validateOrExpandRatesDimension();
        mapping = BranchRateIndexHelper.buildDeterministic(tree);
    }

    private static void requireExactlyOne(final Object legacy,
                                          final Object typed,
                                          final String legacyName,
                                          final String typedName) {
        if (legacy == null && typed == null) {
            throw new IllegalArgumentException("ACSubtreeUIncrementOperator: either "
                    + legacyName + " or " + typedName + " must be specified.");
        }
        requireAtMostOne(legacy, typed, legacyName, typedName);
    }

    private static void requireAtMostOne(final Object legacy,
                                         final Object typed,
                                         final String legacyName,
                                         final String typedName) {
        if (legacy != null && typed != null) {
            throw new IllegalArgumentException("ACSubtreeUIncrementOperator: specify only one of "
                    + legacyName + " or " + typedName + ".");
        }
    }

    private void ensureMappingUpToDate() {
        if (legacyRates != null) {
            mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, legacyRates, "ACSubtreeUIncrementOperator");
        } else {
            mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, typedRates, "ACSubtreeUIncrementOperator");
        }
    }

    private void validateOrExpandRatesDimension() {
        if (legacyRates != null) {
            BranchRateIndexHelper.validateRatesDimension(tree, legacyRates, "ACSubtreeUIncrementOperator");
        } else {
            BranchRateIndexHelper.validateRatesDimension(tree, typedRates, "ACSubtreeUIncrementOperator");
        }
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

    private int indicatorValue() {
        return legacyIndicator != null ? legacyIndicator.getValue(0) : typedIndicator.get();
    }

    private double sigma2Value() {
        return legacySigma2 != null ? legacySigma2.getValue(0) : typedSigma2.get();
    }

    private double rootLog() {
        if (legacyRootLogRate != null) {
            return legacyRootLogRate.getValue(0);
        }
        if (typedRootLogRate != null) {
            return typedRootLogRate.get();
        }
        return 0.0;
    }

    /** Collect all non-root nodes (edges) in the subtree under node 'n'. */
    private void collectEdgesInSubtree(final Node n, final List<Node> out) {
        if (!n.isRoot()) {
            out.add(n);
        }
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
        final int nNodes = mapping.getNodeCount();

        for (int t = 0; t < tries; t++) {
            final Node cand = tree.getNode(Randomizer.nextInt(nNodes));
            if (cand.isRoot() && !allowRoot) {
                continue;
            }
            if (internalOnly && cand.getChildCount() == 0) {
                continue;
            }

            if (maxEdges > 0) {
                final ArrayList<Node> tmp = new ArrayList<>();
                collectEdgesInSubtree(cand, tmp);
                if (tmp.size() > maxEdges) {
                    continue;
                }
            }
            return cand;
        }

        Node cand;
        do {
            cand = tree.getNode(Randomizer.nextInt(nNodes));
        } while (cand.isRoot() && !allowRoot);
        return cand;
    }

    /** Compute old u for a subtree rooted at 'node' given parent log-rate 'logPar'. */
    private void computeUOld(final Node node,
                             final double logPar,
                             final double sig2Value,
                             final double minDt,
                             final double[] xOldByNr,
                             final double[] uOldByNr) {
        final int nr = node.getNr();
        final int idx = mapping.idxForNode(node);
        if (idx < 0) {
            throw new ArithmeticException("Bad mapping");
        }

        final double dt = node.getLength();
        if (!(dt > minDt)) {
            throw new ArithmeticException("dt too small");
        }
        final double var = sig2Value * dt;
        if (!(var > 0.0)) {
            throw new ArithmeticException("var<=0");
        }

        final double r = rateValue(idx);
        if (!(r > 0.0)) {
            throw new ArithmeticException("rate<=0");
        }
        final double x = Math.log(r);

        final double mean = logPar - 0.5 * var;
        final double u = (x - mean) / Math.sqrt(var);

        xOldByNr[nr] = x;
        uOldByNr[nr] = u;

        final int cc = node.getChildCount();
        for (int i = 0; i < cc; i++) {
            computeUOld(node.getChild(i), x, sig2Value, minDt, xOldByNr, uOldByNr);
        }
    }

    /** Build new x for subtree rooted at 'node' given parent new log-rate 'logParNew'. */
    private void buildXNew(final Node node,
                           final double logParNew,
                           final double sig2Value,
                           final double minDt,
                           final double[] uNewByNr,
                           final double[] xNewByNr) {
        final int nr = node.getNr();
        final int idx = mapping.idxForNode(node);
        if (idx < 0) {
            throw new ArithmeticException("Bad mapping");
        }

        final double dt = node.getLength();
        if (!(dt > minDt)) {
            throw new ArithmeticException("dt too small");
        }
        final double var = sig2Value * dt;
        if (!(var > 0.0)) {
            throw new ArithmeticException("var<=0");
        }

        final double u = uNewByNr[nr];
        final double x = logParNew - 0.5 * var + Math.sqrt(var) * u;

        xNewByNr[nr] = x;

        final int cc = node.getChildCount();
        for (int i = 0; i < cc; i++) {
            buildXNew(node.getChild(i), x, sig2Value, minDt, uNewByNr, xNewByNr);
        }
    }

    @Override
    public double proposal() {
        ensureMappingUpToDate();

        final int k = indicatorValue();
        if (k != 1) {
            return rejectIfNotACInput.get() ? Double.NEGATIVE_INFINITY : 0.0;
        }

        final double sig2Value = sigma2Value();
        final double minDt = minBranchLengthInput.get();
        final double delta = deltaInput.get();
        if (!(sig2Value > 0.0) || !(minDt > 0.0) || !(delta > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final Node subRoot = pickSubtreeRoot();

        final ArrayList<Node> edges = new ArrayList<>();
        collectEdgesInSubtree(subRoot, edges);
        if (edges.isEmpty()) {
            return Double.NEGATIVE_INFINITY;
        }

        final int nNodes = mapping.getNodeCount();
        final double[] xOldByNr = new double[nNodes];
        final double[] uOldByNr = new double[nNodes];
        final double[] uNewByNr = new double[nNodes];
        final double[] xNewByNr = new double[nNodes];
        Arrays.fill(xOldByNr, Double.NaN);
        Arrays.fill(uOldByNr, Double.NaN);
        Arrays.fill(uNewByNr, Double.NaN);
        Arrays.fill(xNewByNr, Double.NaN);

        final double boundaryLogPar;
        if (subRoot.isRoot()) {
            boundaryLogPar = rootLog();
            final int cc = subRoot.getChildCount();
            try {
                for (int i = 0; i < cc; i++) {
                    computeUOld(subRoot.getChild(i), boundaryLogPar, sig2Value, minDt, xOldByNr, uOldByNr);
                }
            } catch (ArithmeticException bad) {
                return Double.NEGATIVE_INFINITY;
            }
        } else {
            final Node parent = subRoot.getParent();
            if (parent.isRoot()) {
                boundaryLogPar = rootLog();
            } else {
                final int idxPar = mapping.idxForNode(parent);
                if (idxPar < 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                final double rPar = rateValue(idxPar);
                if (!(rPar > 0.0)) {
                    return Double.NEGATIVE_INFINITY;
                }
                boundaryLogPar = Math.log(rPar);
            }

            try {
                computeUOld(subRoot, boundaryLogPar, sig2Value, minDt, xOldByNr, uOldByNr);
            } catch (ArithmeticException bad) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        for (final Node e : edges) {
            final int nr = e.getNr();
            final double uOld = uOldByNr[nr];
            if (!Double.isFinite(uOld)) {
                return Double.NEGATIVE_INFINITY;
            }
            uNewByNr[nr] = uOld + delta * Randomizer.nextGaussian();
        }

        try {
            if (subRoot.isRoot()) {
                final int cc = subRoot.getChildCount();
                final double rootLog = rootLog();
                for (int i = 0; i < cc; i++) {
                    buildXNew(subRoot.getChild(i), rootLog, sig2Value, minDt, uNewByNr, xNewByNr);
                }
            } else {
                buildXNew(subRoot, boundaryLogPar, sig2Value, minDt, uNewByNr, xNewByNr);
            }
        } catch (ArithmeticException bad) {
            return Double.NEGATIVE_INFINITY;
        }

        double sumDelta = 0.0;
        for (final Node e : edges) {
            final int nr = e.getNr();
            final double xo = xOldByNr[nr];
            final double xn = xNewByNr[nr];
            if (!Double.isFinite(xo) || !Double.isFinite(xn)) {
                return Double.NEGATIVE_INFINITY;
            }
            sumDelta += (xn - xo);
        }

        if (legacyRates != null) {
            legacyRates.startEditing(this);
        }
        for (final Node e : edges) {
            final int idx = mapping.idxForNode(e);
            setRateValue(idx, Math.exp(xNewByNr[e.getNr()]));
        }

        return sumDelta;
    }
}
