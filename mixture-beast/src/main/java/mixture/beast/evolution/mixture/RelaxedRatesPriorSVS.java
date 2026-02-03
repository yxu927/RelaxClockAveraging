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
import java.util.List;
import java.util.Random;

@Description("SVS-style relaxed-clock rate prior on a shared vector of positive branch rates. "
        + "indicator=0: UC i.i.d. LogNormal on rates (E[r]=1). "
        + "indicator=1: AC lognormal increments along the tree with mean-correction so E[r_child|r_parent]=r_parent. "
        + "AC is implemented as Normal on log(rates) + Jacobian term for rate-space parameterization.")
public class RelaxedRatesPriorSVS extends Distribution {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "positive branch rates for NON-root nodes; dimension=(tree.nodeCount - 1). "
                    + "Interpretation: each non-root node holds the rate for the branch leading to that node.",
            Input.Validate.REQUIRED
    );

    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator",
            "dimension=1; 0=uncorrelated, 1=autocorrelated",
            Input.Validate.REQUIRED
    );

    // UC parameter (lognormal sigma on log scale)
    public final Input<RealParameter> ucldStdevInput = new Input<>(
            "ucldStdev",
            "lognormal stdev (sigma on log scale) for UC prior, dimension=1",
            Input.Validate.REQUIRED
    );

    // AC parameters
    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "OPTIONAL root log-rate for AC.(relative root rate=1).",
            Input.Validate.OPTIONAL
    );

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "Brownian variance parameter (sigma^2) per unit time for AC prior, dimension=1",
            Input.Validate.REQUIRED
    );

    // Numerical settings
    public final Input<Double> minBranchLengthInput = new Input<>(
            "minBranchLength",
            "minimum branch length (time) allowed in AC; if a branch is shorter, logP=-inf. "
                    + "Set small (e.g. 1e-12) if you want to allow extremely short branches.",
            1e-12
    );

    private Tree tree;
    private RealParameter rates;
    private IntegerParameter indicator;
    private RealParameter ucldStdev;
    private RealParameter rootLogRate; // optional
    private RealParameter sigma2;

    private int nNodes;
    private int rootNr;

    // nodeNr -> index in rates vector; root -> -1
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

        buildIndexMap();
    }

    private void buildIndexMap() {
        idxMap = new int[nNodes];
        Arrays.fill(idxMap, -2);

        int k = 0;
        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            int nr = node.getNr();
            if (nr < 0 || nr >= nNodes) {
                throw new IllegalArgumentException("Node nr out of range: " + nr + " (nNodes=" + nNodes + ")");
            }
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

    /**
     * Exposed UC log-density for use by a Gibbs-style indicator operator.
     * UC: r_i ~ LogNormal(mu, s), i.i.d. in RATE space, with E[r]=1.
     */
    public double logPriorUCOnly() {
        final double s = ucldStdev.getValue();
        if (!(s > 0.0)) return Double.NEGATIVE_INFINITY;

        final double s2 = s * s;
        final double mu = -0.5 * s2; // ensures E[r]=1

        final double logSqrt2Pi = 0.5 * Math.log(2.0 * Math.PI);

        double lp = 0.0;

        // LogNormal pdf on rate r:
        // log f(r) = -log(r) -log(s) -0.5*log(2*pi) - (log(r)-mu)^2/(2*s^2)
        for (int i = 0; i < rates.getDimension(); i++) {
            final double r = rates.getValue(i);
            if (!(r > 0.0)) return Double.NEGATIVE_INFINITY;

            final double x = Math.log(r);
            final double z = (x - mu);
            lp += -x - Math.log(s) - logSqrt2Pi - 0.5 * (z * z) / s2;
        }

        return lp;
    }

    /**
     * Exposed AC log-density for use by a Gibbs-style indicator operator.
     *
     * AC increments (mean-corrected):
     * log(r_child) | log(r_parent) ~ Normal( log(r_parent) - 0.5*var, var ), var = sigma2 * dt
     * => r_child | r_parent ~ LogNormal( log(r_parent) - 0.5*var, sqrt(var) ) with mean r_parent
     *
     * Since the parameter is r (rate space), we include Jacobian term -log(r_child) for each child rate.
     * rootLogRate (if provided) is used as log(r_parent) for edges descending from the root; otherwise 0.
     */
    public double logPriorACOnly() {
        final double s2 = sigma2.getValue();
        if (!(s2 > 0.0)) return Double.NEGATIVE_INFINITY;

        final double minDt = minBranchLengthInput.get();
        if (!(minDt > 0.0)) return Double.NEGATIVE_INFINITY;

        final double log2Pi = Math.log(2.0 * Math.PI);
        final double rootLog = (rootLogRate == null ? 0.0 : rootLogRate.getValue());

        double lp = 0.0;

        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            if (node.isRoot()) continue;

            final double dt = node.getLength();
            if (!(dt > minDt)) {
                // If you want to allow zero-length branches, set minBranchLength very small (e.g. 1e-12),
                // or change this to: dt = minDt (but that changes the model).
                return Double.NEGATIVE_INFINITY;
            }

            final double var = s2 * dt;

            final int idxChi = idxMap[node.getNr()];
            if (idxChi < 0) return Double.NEGATIVE_INFINITY;

            final double rChi = rates.getValue(idxChi);
            if (!(rChi > 0.0)) return Double.NEGATIVE_INFINITY;

            final double logChi = Math.log(rChi);

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

            // Mean correction so that E[r_child | r_parent] = r_parent
            final double mean = logPar - 0.5 * var;

            final double z = logChi - mean;
            lp += -0.5 * (log2Pi + Math.log(var) + (z * z) / var);

            // Jacobian for x=log(r): log|dx/dr| = -log(r)
            lp += -logChi;
        }

        return lp;
    }

    @Override
    public List<String> getArguments() {
        return List.of(rates.getID());
    }

    @Override
    public List<String> getConditions() {
        // rootLogRate may be null
        if (rootLogRate == null) {
            return List.of(tree.getID(), indicator.getID(), ucldStdev.getID(), sigma2.getID());
        }
        return List.of(tree.getID(), indicator.getID(), ucldStdev.getID(), rootLogRate.getID(), sigma2.getID());
    }

    @Override
    public void sample(State state, Random random) {
        // optional: implement if you need direct sampling; otherwise leave empty
    }

    @Override
    protected boolean requiresRecalculation() {
        boolean dirty = false;

        if (tree != null && ((StateNode) tree).somethingIsDirty()) dirty = true;
        if (rates != null && rates.somethingIsDirty()) dirty = true;
        if (indicator != null && indicator.somethingIsDirty()) dirty = true;
        if (ucldStdev != null && ucldStdev.somethingIsDirty()) dirty = true;
        if (sigma2 != null && sigma2.somethingIsDirty()) dirty = true;
        if (rootLogRate != null && rootLogRate.somethingIsDirty()) dirty = true;

        return dirty;
    }
}
