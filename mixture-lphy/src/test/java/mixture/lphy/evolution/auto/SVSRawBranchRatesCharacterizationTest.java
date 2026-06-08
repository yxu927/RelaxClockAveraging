package mixture.lphy.evolution.auto;

import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

public class SVSRawBranchRatesCharacterizationTest {

    private static final double EPS = 1.0e-12;
    private static final double LOG_2PI = Math.log(2.0 * Math.PI);

    @Test
    public void acLogDensityUsesMartingaleMeanCorrection() {
        final TreeFixture fixture = fixedTree();
        final double sigma2 = 0.2;
        final double rootLogRate = 0.3;
        final Double[] rates = rates(fixture, 0.10, 0.25, -0.20, -0.15);

        final SVSRawBranchRates distribution = new SVSRawBranchRates(
                new Value<>("tree", fixture.tree),
                new Value<>("indicator", 1),
                new Value<>("ucldStdev", 0.5),
                new Value<>("sigma2", sigma2),
                new Value<>("rootLogRate", rootLogRate)
        );

        final double corrected = correctedACLogDensity(fixture, sigma2, rootLogRate, rates);
        final double uncorrected = uncorrectedACLogDensity(fixture, sigma2, rootLogRate, rates);
        final double actual = distribution.logDensity(rates);

        assertEquals(corrected, actual, EPS);
        assertNotEquals(uncorrected, actual, 1.0e-8);
    }

    @Test
    public void ucLogDensityIsMeanOneLogNormalInRateSpace() {
        final TreeFixture fixture = fixedTree();
        final double stdev = 0.5;
        final Double[] rates = rates(fixture, 0.10, 0.25, -0.20, -0.15);

        final SVSRawBranchRates distribution = new SVSRawBranchRates(
                new Value<>("tree", fixture.tree),
                new Value<>("indicator", 0),
                new Value<>("ucldStdev", stdev),
                new Value<>("sigma2", 0.2),
                new Value<>("rootLogRate", 0.0)
        );

        final double var = stdev * stdev;
        final double meanLog = -0.5 * var;
        double expected = 0.0;
        for (final TimeTreeNode node : fixture.nonRootNodes()) {
            expected += logNormalRateDensity(logRate(rates, node), meanLog, var);
        }

        assertEquals(expected, distribution.logDensity(rates), EPS);
    }

    private static TreeFixture fixedTree() {
        final TimeTree tree = new TimeTree();

        final TimeTreeNode root = node(tree, "root", 5.0, 0);
        final TimeTreeNode internal = node(tree, "internal", 2.0, 1);
        final TimeTreeNode a = node(tree, "A", 0.0, 2);
        final TimeTreeNode b = node(tree, "B", 0.0, 3);
        final TimeTreeNode c = node(tree, "C", 0.0, 4);

        root.addChild(internal);
        root.addChild(c);
        internal.addChild(a);
        internal.addChild(b);
        tree.setRoot(root);

        return new TreeFixture(tree, root, internal, a, b, c);
    }

    private static TimeTreeNode node(final TimeTree tree, final String id, final double age, final int index) {
        final TimeTreeNode node = new TimeTreeNode(id, tree);
        node.setAge(age);
        node.setIndex(index);
        return node;
    }

    private static Double[] rates(final TreeFixture fixture,
                                  final double internalLogRate,
                                  final double aLogRate,
                                  final double bLogRate,
                                  final double cLogRate) {
        final Double[] rates = new Double[maxIndex(fixture) + 1];
        for (int i = 0; i < rates.length; i++) {
            rates[i] = 0.0;
        }
        rates[fixture.internal.getIndex()] = Math.exp(internalLogRate);
        rates[fixture.a.getIndex()] = Math.exp(aLogRate);
        rates[fixture.b.getIndex()] = Math.exp(bLogRate);
        rates[fixture.c.getIndex()] = Math.exp(cLogRate);
        return rates;
    }

    private static int maxIndex(final TreeFixture fixture) {
        int max = fixture.root.getIndex();
        for (final TimeTreeNode node : fixture.nonRootNodes()) {
            max = Math.max(max, node.getIndex());
        }
        return max;
    }

    private static double correctedACLogDensity(final TreeFixture fixture,
                                                final double sigma2,
                                                final double rootLogRate,
                                                final Double[] rates) {
        double lp = 0.0;
        lp += acEdgeLogDensity(logRate(rates, fixture.internal), rootLogRate,
                sigma2 * branchLength(fixture.root, fixture.internal), true);
        lp += acEdgeLogDensity(logRate(rates, fixture.c), rootLogRate,
                sigma2 * branchLength(fixture.root, fixture.c), true);
        lp += acEdgeLogDensity(logRate(rates, fixture.a), logRate(rates, fixture.internal),
                sigma2 * branchLength(fixture.internal, fixture.a), true);
        lp += acEdgeLogDensity(logRate(rates, fixture.b), logRate(rates, fixture.internal),
                sigma2 * branchLength(fixture.internal, fixture.b), true);
        return lp;
    }

    private static double uncorrectedACLogDensity(final TreeFixture fixture,
                                                  final double sigma2,
                                                  final double rootLogRate,
                                                  final Double[] rates) {
        double lp = 0.0;
        lp += acEdgeLogDensity(logRate(rates, fixture.internal), rootLogRate,
                sigma2 * branchLength(fixture.root, fixture.internal), false);
        lp += acEdgeLogDensity(logRate(rates, fixture.c), rootLogRate,
                sigma2 * branchLength(fixture.root, fixture.c), false);
        lp += acEdgeLogDensity(logRate(rates, fixture.a), logRate(rates, fixture.internal),
                sigma2 * branchLength(fixture.internal, fixture.a), false);
        lp += acEdgeLogDensity(logRate(rates, fixture.b), logRate(rates, fixture.internal),
                sigma2 * branchLength(fixture.internal, fixture.b), false);
        return lp;
    }

    private static double logRate(final Double[] rates, final TimeTreeNode node) {
        return Math.log(rates[node.getIndex()]);
    }

    private static double branchLength(final TimeTreeNode parent, final TimeTreeNode child) {
        return parent.getAge() - child.getAge();
    }

    private static double acEdgeLogDensity(final double childLogRate,
                                           final double parentLogRate,
                                           final double var,
                                           final boolean martingaleCorrected) {
        final double meanLog = martingaleCorrected ? parentLogRate - 0.5 * var : parentLogRate;
        return logNormalRateDensity(childLogRate, meanLog, var);
    }

    private static double logNormalRateDensity(final double logRate, final double meanLog, final double var) {
        final double z = logRate - meanLog;
        return -logRate - 0.5 * (LOG_2PI + Math.log(var) + (z * z) / var);
    }

    private record TreeFixture(TimeTree tree,
                               TimeTreeNode root,
                               TimeTreeNode internal,
                               TimeTreeNode a,
                               TimeTreeNode b,
                               TimeTreeNode c) {

        private TimeTreeNode[] nonRootNodes() {
            return new TimeTreeNode[]{internal, a, b, c};
        }
    }
}
