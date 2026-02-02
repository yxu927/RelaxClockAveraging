package mixture.beast.evolution.mixture;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

@Description("SVS-style relaxed-clock rate prior on a shared vector of positive branch rates. "
        + "indicator=0: UC i.i.d. LogNormal on rates (mean=1). "
        + "indicator=1: AC Brownian increments on log(rates) along the tree with Jacobian for rates-space.")
public class RelaxedRatesPriorSVS extends Distribution {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "positive branch rates for NON-root nodes; dimension=(tree.nodeCount - 1)",
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
            "root log-rate (dimension=1). If you want normalized root fixed, set this to 0 and fix it in XML.",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "Brownian variance parameter (sigma^2) for AC prior, dimension=1",
            Input.Validate.REQUIRED
    );

    private Tree tree;
    private RealParameter rates;
    private IntegerParameter indicator;
    private RealParameter ucldStdev;
    private RealParameter rootLogRate;
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
        rootLogRate = rootLogRateInput.get();
        sigma2 = sigma2Input.get();

        if (indicator.getDimension() != 1) {
            throw new IllegalArgumentException("indicator must have dimension=1.");
        }
        if (ucldStdev.getDimension() != 1) {
            throw new IllegalArgumentException("ucldStdev must have dimension=1.");
        }
        if (rootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate must have dimension=1.");
        }
        if (sigma2.getDimension() != 1) {
            throw new IllegalArgumentException("sigma2 must have dimension=1.");
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
            logP = logPriorUC();
        } else if (k == 1) {
            logP = logPriorAC();
        } else {
            logP = Double.NEGATIVE_INFINITY;
        }
        return logP;
    }

    /**
     * Uncorrelated prior: r_i ~ LogNormal(mu, s) independently (in RATE space).
     * Choose mu=-0.5*s^2 so E[r]=1.
     */
    private double logPriorUC() {
        final double s = ucldStdev.getValue();
        if (!(s > 0.0)) return Double.NEGATIVE_INFINITY;

        final double s2 = s * s;
        final double mu = -0.5 * s2;

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
     * Autocorrelated prior on log(rates) along the tree:
     * log(r_child) - log(r_parent) ~ N(0, sigma2*dt).
     * Because the parameter is r (rate space), include Jacobian term -log(r_child) for each child rate.
     * rootLogRate is used as the parent's log-rate for edges descending from the root.
     */
    private double logPriorAC() {
        final double s2 = sigma2.getValue();
        if (!(s2 > 0.0)) return Double.NEGATIVE_INFINITY;

        final double tinyVar = 1e-12;
        final double log2Pi = Math.log(2.0 * Math.PI);

        double lp = 0.0;

        for (int i = 0; i < nNodes; i++) {
            Node node = tree.getNode(i);
            if (node.isRoot()) continue;

            final int nr = node.getNr();
            final int idxChi = idxMap[nr];
            if (idxChi < 0) return Double.NEGATIVE_INFINITY;

            final double rChi = rates.getValue(idxChi);
            if (!(rChi > 0.0)) return Double.NEGATIVE_INFINITY;

            final double logChi = Math.log(rChi);

            final Node parent = node.getParent();
            final double logPar;
            if (parent.isRoot()) {
                logPar = rootLogRate.getValue();
            } else {
                final int idxPar = idxMap[parent.getNr()];
                if (idxPar < 0) return Double.NEGATIVE_INFINITY;
                final double rPar = rates.getValue(idxPar);
                if (!(rPar > 0.0)) return Double.NEGATIVE_INFINITY;
                logPar = Math.log(rPar);
            }

            final double diff = logChi - logPar;

            final double dt = node.getLength();
            final double var = s2 * (dt > 0.0 ? dt : 0.0);

            if (var <= tinyVar) {
                // Deterministic limit: diff must be ~ 0
                if (Math.abs(diff) > 1e-9) return Double.NEGATIVE_INFINITY;
            } else {
                lp += -0.5 * (log2Pi + Math.log(var) + (diff * diff) / var);
            }

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
        return List.of(tree.getID(), indicator.getID(), ucldStdev.getID(), rootLogRate.getID(), sigma2.getID());
    }

    @Override
    public void sample(State state, Random random) {

    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
}
