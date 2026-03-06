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

    private int nNodes;
    private int rootNr;
    private int[] idxMap;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();
        indicator = indicatorInput.get();
        ucldStdev = ucldStdevInput.get();
        rootLogRate = rootLogRateInput.get(); // may be null
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

        nNodes = tree.getNodeCount();
        rootNr = tree.getRoot().getNr();

        if (rates.getDimension() != nNodes - 1) {
            throw new IllegalArgumentException("rates must have dimension (nodeCount - 1). Found "
                    + rates.getDimension() + " vs " + (nNodes - 1));
        }

        buildIndexMapDeterministic();
    }


    private void buildIndexMapDeterministic() {
        idxMap = new int[nNodes];
        Arrays.fill(idxMap, -2);

        // sanity: collect nodes by nr to ensure 0..nNodes-1 all present exactly once
        Node[] byNr = new Node[nNodes];
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            int nr = node.getNr();
            if (nr < 0 || nr >= nNodes) {
                throw new IllegalArgumentException("Node nr out of range: " + nr + " (nNodes=" + nNodes + ")");
            }
            if (byNr[nr] != null) {
                throw new IllegalStateException("Duplicate node nr encountered: " + nr);
            }
            byNr[nr] = node;
        }
        for (int nr = 0; nr < nNodes; nr++) {
            if (byNr[nr] == null) {
                throw new IllegalStateException("Missing node for nr=" + nr);
            }
        }

        int k = 0;
        for (int nr = 0; nr < nNodes; nr++) {
            if (nr == rootNr) {
                idxMap[nr] = -1;
            } else {
                idxMap[nr] = k++;
            }
        }
        if (k != nNodes - 1) {
            throw new IllegalStateException("Internal error: expected to map " + (nNodes - 1)
                    + " non-root nodes but mapped " + k);
        }
    }

    /** log pdf of LogNormal in RATE space: r>0, log(r) ~ Normal(meanLog, var). */
    private static double logLogNormalRate(final double r, final double meanLog, final double var) {
        if (!(r > 0.0) || !(var > 0.0)) return Double.NEGATIVE_INFINITY;

        final double x = Math.log(r);     // log(r)
        final double z = x - meanLog;
        // log p(r) = log p(x) + log|dx/dr| = log Normal(x; meanLog,var) - log(r)
        return -x - 0.5 * (LOG_2PI + Math.log(var) + (z * z) / var);
    }

    @Override
    public double calculateLogP() {
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
        if (!(s > 0.0)) return Double.NEGATIVE_INFINITY;

        final double var = s * s;
        final double meanLog = -0.5 * var; // ensures E[r]=1

        double lp = 0.0;
        for (int i = 0; i < rates.getDimension(); i++) {
            lp += logLogNormalRate(rates.getValue(i), meanLog, var);
            if (Double.isInfinite(lp)) return Double.NEGATIVE_INFINITY;
        }
        return lp;
    }

    /** AC: per-branch lognormal increments with mean-correction (martingale): E[r_child|r_parent]=r_parent. */
    public double logPriorACOnly() {
        final double s2 = sigma2.getValue(0);
        if (!(s2 > 0.0)) return Double.NEGATIVE_INFINITY;

        final double minDt = minBranchLengthInput.get();
        if (!(minDt > 0.0)) return Double.NEGATIVE_INFINITY;

        final double rootLog = (rootLogRate == null ? 0.0 : rootLogRate.getValue(0));

        double lp = 0.0;

        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            if (node.isRoot()) continue;

            final double dt = node.getLength();
            if (!(dt > minDt)) return Double.NEGATIVE_INFINITY;

            final double var = s2 * dt;

            final int idxChi = idxMap[node.getNr()];
            if (idxChi < 0) return Double.NEGATIVE_INFINITY;

            final double rChi = rates.getValue(idxChi);
            if (!(rChi > 0.0)) return Double.NEGATIVE_INFINITY;

            final Node parent = node.getParent();
            final double logPar;
            if (parent.isRoot()) {
                logPar = rootLog;
            } else {
                final int idxPar = idxMap[parent.getNr()];
                if (idxPar < 0) return Double.NEGATIVE_INFINITY;

                final double rPar = rates.getValue(idxPar);
                if (!(rPar > 0.0)) return Double.NEGATIVE_INFINITY;

                logPar = Math.log(rPar);
            }
            // mean-corrected: log(r_child) ~ N(logPar - 0.5*var, var)
            final double meanLog = logPar - 0.5 * var;
            lp += logLogNormalRate(rChi, meanLog, var);
            if (Double.isInfinite(lp)) return Double.NEGATIVE_INFINITY;
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
    public void sample(State state, Random random) {
        // optional
    }

    @Override
    protected boolean requiresRecalculation() {
        boolean dirty = false;

        // If tree changed in a way that alters nodeCount/rootNr, rebuild map.
        if (tree != null && ((StateNode) tree).somethingIsDirty()) {
            dirty = true;
            final int newN = tree.getNodeCount();
            final int newRoot = tree.getRoot().getNr();
            if (newN != nNodes || newRoot != rootNr) {
                nNodes = newN;
                rootNr = newRoot;
                buildIndexMapDeterministic();
            }
        }

        if (rates != null && rates.somethingIsDirty()) dirty = true;
        if (indicator != null && indicator.somethingIsDirty()) dirty = true;
        if (ucldStdev != null && ucldStdev.somethingIsDirty()) dirty = true;
        if (sigma2 != null && sigma2.somethingIsDirty()) dirty = true;
        if (rootLogRate != null && rootLogRate.somethingIsDirty()) dirty = true;

        return dirty;
    }
}
