package mixture.beast.evolution.mixture;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import mixture.beast.evolution.util.BranchRateIndexHelper;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

@Description("SVS-style relaxed-clock rate prior on a shared vector of positive branch rates. "
        + "indicator=0: UC i.i.d. LogNormal on rates (E[r]=1). "
        + "indicator=1: AC lognormal increments along the tree with mean-correction so E[r_child|r_parent]=r_parent.")
public class RelaxedRatesPriorSVS extends Distribution {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates", "positive branch rates for NON-root nodes; dimension=(tree.nodeCount - 1).",
            Input.Validate.REQUIRED);

    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator", "dimension=1; 0=uncorrelated, 1=autocorrelated",
            Input.Validate.REQUIRED);

    public final Input<RealParameter> ucldStdevInput = new Input<>(
            "ucldStdev", "lognormal stdev (sigma on log scale) for UC prior, dimension=1",
            Input.Validate.REQUIRED);

    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate", "OPTIONAL root log-rate for AC (if absent, root log-rate = 0 => root rate = 1).",
            Input.Validate.OPTIONAL);

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2", "Brownian variance parameter (sigma^2) per unit time for AC prior, dimension=1",
            Input.Validate.REQUIRED);

    public final Input<Double> minBranchLengthInput = new Input<>(
            "minBranchLength", "minimum branch length (time) allowed in AC; if shorter, logP=-inf.", 1e-12);

    private static final double LOG_2PI = Math.log(2.0 * Math.PI);

    private Tree tree;
    private RealParameter rates;
    private IntegerParameter indicator;
    private RealParameter ucldStdev;
    private RealParameter rootLogRate;
    private RealParameter sigma2;

    private BranchRateIndexHelper.Mapping mapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();
        indicator = indicatorInput.get();
        ucldStdev = ucldStdevInput.get();
        rootLogRate = rootLogRateInput.get();
        sigma2 = sigma2Input.get();

        if (indicator.getDimension() != 1) {
            throw new IllegalArgumentException("indicator must have dimension=1.");
        }
        if (ucldStdev.getDimension() != 1) {
            throw new IllegalArgumentException("ucldStdev must have dimension=1.");
        }
        if (sigma2.getDimension() != 1) {
            throw new IllegalArgumentException("sigma2 must have dimension=1.");
        }
        if (rootLogRate != null && rootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate must have dimension=1.");
        }

        BranchRateIndexHelper.validateRatesDimension(tree, rates, "RelaxedRatesPriorSVS");
        mapping = BranchRateIndexHelper.buildDeterministic(tree);
    }

    private void ensureMappingUpToDate() {
        mapping = BranchRateIndexHelper.ensureUpToDate(tree, mapping, rates, "RelaxedRatesPriorSVS");
    }

    /** log pdf of LogNormal in RATE space: r>0, log(r) ~ Normal(meanLog, var). */
    private static double logLogNormalRate(final double r, final double meanLog, final double var) {
        if (!(r > 0.0) || !(var > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double x = Math.log(r);
        final double z = x - meanLog;
        return -x - 0.5 * (LOG_2PI + Math.log(var) + (z * z) / var);
    }

    @Override
    public double calculateLogP() {
        ensureMappingUpToDate();

        final int k = indicator.getValue(0);
        if (k == 0) {
            logP = logPriorUCOnly();
        } else if (k == 1) {
            logP = logPriorACOnly();
        } else {
            logP = Double.NEGATIVE_INFINITY;
        }
        return logP;
    }

    /** UC: r_i iid LogNormal with E[r]=1. */
    public double logPriorUCOnly() {
        final double s = ucldStdev.getValue(0);
        if (!(s > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double var = s * s;
        final double meanLog = -0.5 * var;

        double lp = 0.0;
        for (int i = 0; i < rates.getDimension(); i++) {
            lp += logLogNormalRate(rates.getValue(i), meanLog, var);
            if (Double.isInfinite(lp)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return lp;
    }

    /** AC: per-branch lognormal increments with mean-correction (martingale): E[r_child|r_parent]=r_parent. */
    public double logPriorACOnly() {
        ensureMappingUpToDate();

        final double s2 = sigma2.getValue(0);
        if (!(s2 > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double minDt = minBranchLengthInput.get();
        if (!(minDt > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double rootLog = (rootLogRate == null ? 0.0 : rootLogRate.getValue(0));

        double lp = 0.0;
        final int nNodes = mapping.getNodeCount();

        for (int i = 0; i < nNodes; i++) {
            final Node node = tree.getNode(i);
            if (node.isRoot()) {
                continue;
            }

            final double dt = node.getLength();
            if (!(dt > minDt)) {
                return Double.NEGATIVE_INFINITY;
            }

            final double var = s2 * dt;

            final int idxChi = mapping.idxForNode(node);
            if (idxChi < 0) {
                return Double.NEGATIVE_INFINITY;
            }

            final double rChi = rates.getValue(idxChi);
            if (!(rChi > 0.0)) {
                return Double.NEGATIVE_INFINITY;
            }

            final Node parent = node.getParent();
            final double logPar;
            if (parent.isRoot()) {
                logPar = rootLog;
            } else {
                final int idxPar = mapping.idxForNode(parent);
                if (idxPar < 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                final double rPar = rates.getValue(idxPar);
                if (!(rPar > 0.0)) {
                    return Double.NEGATIVE_INFINITY;
                }

                logPar = Math.log(rPar);
            }

            final double meanLog = logPar - 0.5 * var;
            lp += logLogNormalRate(rChi, meanLog, var);
            if (Double.isInfinite(lp)) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        return lp;
    }

    @Override
    public List<String> getArguments() {
        return Collections.singletonList(rates.getID());
    }

    @Override
    public List<String> getConditions() {
        if (rootLogRate == null) {
            return Arrays.asList(tree.getID(), indicator.getID(), ucldStdev.getID(), sigma2.getID());
        }
        return Arrays.asList(tree.getID(), indicator.getID(), ucldStdev.getID(), rootLogRate.getID(), sigma2.getID());
    }

    @Override
    public void sample(final State state, final Random random) {
        // optional
    }

    @Override
    protected boolean requiresRecalculation() {
        boolean dirty = false;

        if (tree != null && ((StateNode) tree).somethingIsDirty()) {
            dirty = true;
            ensureMappingUpToDate();
        }
        if (rates != null && rates.somethingIsDirty()) {
            dirty = true;
        }
        if (indicator != null && indicator.somethingIsDirty()) {
            dirty = true;
        }
        if (ucldStdev != null && ucldStdev.somethingIsDirty()) {
            dirty = true;
        }
        if (sigma2 != null && sigma2.somethingIsDirty()) {
            dirty = true;
        }
        if (rootLogRate != null && rootLogRate.somethingIsDirty()) {
            dirty = true;
        }

        return dirty;
    }
}
