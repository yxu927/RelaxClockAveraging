package mixture.lphy.evolution.auto;

import lphy.base.distribution.ParametricDistribution;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.*;


public class SVSRawBranchRates extends ParametricDistribution<Double[]> {

    public static final String TREE = "tree";
    public static final String INDICATOR = "indicator";
    public static final String UCLD_STDEV = "ucldStdev";
    public static final String SIGMA2 = "sigma2";
    public static final String ROOT_LOG_RATE = "rootLogRate";

    private Value<TimeTree> tree;
    private Value<Integer> indicator;
    private Value<Double> ucldStdev;
    private Value<Double> sigma2;
    private Value<Double> rootLogRate;

    private int rootIndex;
    private List<TimeTreeNode> allNodes;
    private int maxIndex;

    public SVSRawBranchRates(
            @ParameterInfo(name = TREE, description = "TimeTree") Value<TimeTree> tree,
            @ParameterInfo(name = INDICATOR, description = "0=UC, 1=AC") Value<Integer> indicator,
            @ParameterInfo(name = UCLD_STDEV, description = "UC lognormal stdev (log scale)") Value<Double> ucldStdev,
            @ParameterInfo(name = SIGMA2, description = "AC Brownian variance sigma^2") Value<Double> sigma2,
            @ParameterInfo(name = ROOT_LOG_RATE, description = "AC root log-rate") Value<Double> rootLogRate
    ) {
        this.tree = tree;
        this.indicator = indicator;
        this.ucldStdev = ucldStdev;
        this.sigma2 = sigma2;
        this.rootLogRate = rootLogRate;
        constructDistribution(random);
    }

    @Override
    protected void constructDistribution(RandomGenerator random) {
        TimeTree t = tree.value();
        TimeTreeNode root = t.getRoot();
        rootIndex = root.getIndex();

        allNodes = new ArrayList<>();
        collectAllNodes(root, allNodes);

        maxIndex = -1;
        for (TimeTreeNode nd : allNodes) {
            maxIndex = Math.max(maxIndex, nd.getIndex());
        }
        if (maxIndex < 0) {
            throw new IllegalStateException("Tree has no nodes?");
        }
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = new TreeMap<>();
        map.put(TREE, tree);
        map.put(INDICATOR, indicator);
        map.put(UCLD_STDEV, ucldStdev);
        map.put(SIGMA2, sigma2);
        map.put(ROOT_LOG_RATE, rootLogRate);
        return map;
    }

    @SuppressWarnings("unchecked")
    @Override
    public void setParam(String paramName, Value<?> value) {
        switch (paramName) {
            case TREE -> tree = (Value<TimeTree>) value;
            case INDICATOR -> indicator = (Value<Integer>) value;
            case UCLD_STDEV -> ucldStdev = (Value<Double>) value;
            case SIGMA2 -> sigma2 = (Value<Double>) value;
            case ROOT_LOG_RATE -> rootLogRate = (Value<Double>) value;
            default -> throw new RuntimeException("Unrecognized param: " + paramName);
        }
        super.setParam(paramName, value);
    }

    @GeneratorInfo(
            name = "SVSRawBranchRates",
            verbClause = "are generated under stochastic variable selection for relaxed clocks",
            narrativeName = "SVS branch rates",
            category = GeneratorCategory.PHYLO_LIKELIHOOD,
            description = """
                    Generates raw branch rates under either an uncorrelated (UC) or autocorrelated (AC) model,
                    selected by an indicator. This is intended for the shared-rates SVS design where UC vs AC
                    differs only in the prior on the same rate vector.
                    """
    )
    @Override
    public RandomVariable<Double[]> sample() {
        int k = indicator.value();
        if (k != 0 && k != 1) {
            throw new IllegalArgumentException("indicator must be 0 or 1, got " + k);
        }

        Double[] r = new Double[maxIndex + 1];
        Arrays.fill(r, 0.0);
        r[rootIndex] = 0.0; // unused

        if (k == 0) {
            sampleUC(r);
        } else {
            sampleAC(r);
        }

        return new RandomVariable<>(null, r, this);
    }

    private void sampleUC(Double[] r) {
        double s = ucldStdev.value();
        if (!(s > 0.0)) {
            throw new IllegalArgumentException("ucldStdev must be > 0");
        }
        double mu = -0.5 * s * s;

        NormalDistribution nd = new NormalDistribution(
                random, mu, s,
                NormalDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY
        );

        for (TimeTreeNode node : allNodes) {
            if (node.isRoot()) continue;
            double logRate = nd.sample();
            double rate = Math.exp(logRate);
            r[node.getIndex()] = rate;
        }
    }

    private void sampleAC(Double[] r) {
        double s2 = sigma2.value();
        if (!(s2 > 0.0)) {
            throw new IllegalArgumentException("sigma2 must be > 0");
        }

        TimeTreeNode root = tree.value().getRoot();
        double vRoot = rootLogRate.value();

        // recurse from root, sampling log-rates and storing exp(log-rate) as branch rate on the child edge
        sampleACRec(root, vRoot, r);
    }

    private void sampleACRec(TimeTreeNode parent, double vPar, Double[] r) {
        double s2 = sigma2.value();
        for (TimeTreeNode child : parent.getChildren()) {
            double dt = parent.getAge() - child.getAge();
            if (dt < 0.0) dt = 0.0;

            double var = s2 * dt;
            double vChild;
            if (var <= 0.0) {
                vChild = vPar;
            } else {
                NormalDistribution nd = new NormalDistribution(
                        random, vPar, Math.sqrt(var),
                        NormalDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY
                );
                vChild = nd.sample();
            }

            r[child.getIndex()] = Math.exp(vChild);
            sampleACRec(child, vChild, r);
        }
    }

    @Override
    public double logDensity(Double[] x) {
        int k = indicator.value();
        if (k == 0) {
            return logDensityUC(x);
        } else if (k == 1) {
            return logDensityAC(x);
        } else {
            return Double.NEGATIVE_INFINITY;
        }
    }

    private double logDensityUC(Double[] x) {
        double s = ucldStdev.value();
        if (!(s > 0.0)) return Double.NEGATIVE_INFINITY;

        double s2 = s * s;
        double mu = -0.5 * s2;

        double logS = Math.log(s);
        double logSqrt2Pi = 0.5 * Math.log(2.0 * Math.PI);

        double lp = 0.0;

        for (TimeTreeNode node : allNodes) {
            if (node.isRoot()) continue;

            double r = x[node.getIndex()];
            if (!(r > 0.0)) return Double.NEGATIVE_INFINITY;

            double lr = Math.log(r);
            double z = lr - mu;

            // LogNormal on rate r
            lp += -lr - logS - logSqrt2Pi - 0.5 * (z * z) / s2;
        }

        return lp;
    }

    private double logDensityAC(Double[] x) {
        double s2 = sigma2.value();
        if (!(s2 > 0.0)) return Double.NEGATIVE_INFINITY;

        TimeTreeNode root = tree.value().getRoot();
        double vRoot = rootLogRate.value();

        return logDensityACRec(root, vRoot, x);
    }

    private double logDensityACRec(TimeTreeNode parent, double vPar, Double[] x) {
        double s2 = sigma2.value();
        double lp = 0.0;

        for (TimeTreeNode child : parent.getChildren()) {
            double rChild = x[child.getIndex()];
            if (!(rChild > 0.0)) return Double.NEGATIVE_INFINITY;

            double vChild = Math.log(rChild);

            double dt = parent.getAge() - child.getAge();
            if (dt < 0.0) dt = 0.0;

            double var = s2 * dt;
            if (!(var > 0.0)) {
                // deterministic limit
                if (Math.abs(vChild - vPar) > 1e-9) return Double.NEGATIVE_INFINITY;
            } else {
                double diff = vChild - vPar;
                lp += -0.5 * (Math.log(2.0 * Math.PI * var) + (diff * diff) / var);
            }

            // Jacobian because RV is in rate-space but density is on log-rate increments
            lp += -vChild;

            lp += logDensityACRec(child, vChild, x);
        }

        return lp;
    }

    private void collectAllNodes(TimeTreeNode node, List<TimeTreeNode> out) {
        out.add(node);
        for (TimeTreeNode c : node.getChildren()) {
            collectAllNodes(c, out);
        }
    }
}
