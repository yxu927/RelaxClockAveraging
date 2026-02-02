package mixture.lphy.evolution.auto;

import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.*;

public class SharedRatesClock extends DeterministicFunction<Double[]> {

    public static final String TREE = "tree";
    public static final String RATES = "rates";
    public static final String MEAN_RATE = "meanRate";
    public static final String NORMALIZE = "normalize";

    private Value<TimeTree> tree;
    private Value<Double[]> rates;
    private Value<Double> meanRate;
    private Value<Boolean> normalize;

    public SharedRatesClock(
            @ParameterInfo(name = TREE, description = "TimeTree") Value<TimeTree> tree,
            @ParameterInfo(name = RATES, description = "raw branch rates by node index") Value<Double[]> rates,
            @ParameterInfo(name = MEAN_RATE, description = "global mean rate multiplier", optional = true) Value<Double> meanRate,
            @ParameterInfo(name = NORMALIZE, description = "normalize time-weighted mean to 1", optional = true) Value<Boolean> normalize
    ) {
        this.tree = tree;
        this.rates = rates;
        this.meanRate = meanRate;
        this.normalize = normalize;
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = new TreeMap<>();
        map.put(TREE, tree);
        map.put(RATES, rates);
        if (meanRate != null) map.put(MEAN_RATE, meanRate);
        if (normalize != null) map.put(NORMALIZE, normalize);
        return map;
    }

    @GeneratorInfo(
            name = "SharedRatesClock",
            verbClause = "are mapped to CTMC branch rates using shared-rate normalization",
            narrativeName = "shared rates clock",
            category = GeneratorCategory.PHYLO_LIKELIHOOD,
            description = """
                    Deterministically maps a vector of raw branch rates to the branch rates used in the CTMC:
                    optional time-weighted mean normalization (to 1), then multiplication by meanRate.
                    """
    )
    @Override
    public Value<Double[]> apply() {
        TimeTree t = tree.value();
        Double[] raw = rates.value();

        double mr = (meanRate == null ? 1.0 : meanRate.value());
        boolean doNorm = (normalize == null ? false : normalize.value());

        List<TimeTreeNode> nodes = new ArrayList<>();
        collectAllNodes(t.getRoot(), nodes);

        double scale = 1.0;
        if (doNorm) {
            double sumRateTime = 0.0;
            double sumTime = 0.0;

            for (TimeTreeNode node : nodes) {
                if (node.isRoot()) continue;

                double dt = node.getParent().getAge() - node.getAge();
                if (dt < 0.0) dt = 0.0;
                if (dt == 0.0) continue;

                double r = raw[node.getIndex()];
                if (!(r > 0.0)) continue;

                sumRateTime += r * dt;
                sumTime += dt;
            }

            if (sumRateTime > 0.0 && sumTime > 0.0) {
                scale = sumTime / sumRateTime;
            }
        }

        Double[] out = new Double[raw.length];
        Arrays.fill(out, 0.0);

        for (TimeTreeNode node : nodes) {
            if (node.isRoot()) {
                out[node.getIndex()] = 0.0;
                continue;
            }
            double r = raw[node.getIndex()];
            out[node.getIndex()] = (r > 0.0 ? r * scale * mr : 0.0);
        }

        return new Value<>(out, this);
    }

    private void collectAllNodes(TimeTreeNode node, List<TimeTreeNode> out) {
        out.add(node);
        for (TimeTreeNode c : node.getChildren()) {
            collectAllNodes(c, out);
        }
    }
}
