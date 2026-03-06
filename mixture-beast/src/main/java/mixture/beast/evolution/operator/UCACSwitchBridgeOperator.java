package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.util.Arrays;

@Description("UC<->AC switch operator that deterministically maps the shared rate vector through "
        + "a latent standard-normal vector u (non-centered bridge). "
        + "Greatly improves mixing of the indicator on large trees.")
public class UCACSwitchBridgeOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);
    public final Input<RealParameter> ratesInput = new Input<>("rates", "shared positive branch rates (non-root nodes)", Input.Validate.REQUIRED);
    public final Input<IntegerParameter> indicatorInput = new Input<>("indicator", "0=UC, 1=AC", Input.Validate.REQUIRED);

    public final Input<RealParameter> ucldStdevInput = new Input<>("ucldStdev", "UC lognormal stdev (sigma on log scale)", Input.Validate.REQUIRED);
    public final Input<RealParameter> sigma2Input = new Input<>("sigma2", "AC Brownian variance per unit time", Input.Validate.REQUIRED);
    public final Input<RealParameter> rootLogRateInput = new Input<>("rootLogRate", "optional root log-rate anchor (default 0)", Input.Validate.OPTIONAL);

    public final Input<Double> minBranchLengthInput = new Input<>("minBranchLength", "min dt allowed for AC mapping", 1e-12);

    private Tree tree;
    private RealParameter rates;
    private IntegerParameter indicator;
    private RealParameter ucldStdev;
    private RealParameter sigma2;
    private RealParameter rootLogRate;

    private int nNodes;
    private int rootNr;
    private int[] idxMap; // nodeNr -> ratesIndex, root -> -1

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rates = ratesInput.get();
        indicator = indicatorInput.get();
        ucldStdev = ucldStdevInput.get();
        sigma2 = sigma2Input.get();
        rootLogRate = rootLogRateInput.get();

        if (indicator.getDimension() != 1) throw new IllegalArgumentException("indicator dimension must be 1");
        if (ucldStdev.getDimension() != 1) throw new IllegalArgumentException("ucldStdev dimension must be 1");
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
            if (nr < 0 || nr >= nNodes) {
                throw new IllegalArgumentException("Node nr out of range: " + nr);
            }
            if (nr == rootNr) {
                idxMap[nr] = -1;
            } else {
                idxMap[nr] = k++;
            }
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
            if (rates.getDimension() != nNodes - 1) {
                throw new IllegalStateException("rates dimension inconsistent with tree");
            }
            rebuildIdxMap_likeYourPrior();
        }
    }

    private double rootLog() {
        return (rootLogRate == null ? 0.0 : rootLogRate.getValue(0));
    }

    private static final class Sum {
        double sumLogVar = 0.0;
    }

    @Override
    public double proposal() {
        ensureMapsUpToDate();

        final int k = indicator.getValue(0);
        if (k == 0) {
            return proposeUCtoAC();
        } else if (k == 1) {
            return proposeACtoUC();
        }
        return Double.NEGATIVE_INFINITY;
    }

    private double proposeUCtoAC() {
        final double s = ucldStdev.getValue(0);
        final double sig2 = sigma2.getValue(0);
        final double minDt = minBranchLengthInput.get();
        if (!(s > 0.0) || !(sig2 > 0.0) || !(minDt > 0.0)) return Double.NEGATIVE_INFINITY;

        final int nEdges = rates.getDimension();
        final double muUC = -0.5 * s * s;

        final double[] xOld = new double[nEdges];
        final double[] u = new double[nEdges];

        for (int i = 0; i < nEdges; i++) {
            final double r = rates.getValue(i);
            if (!(r > 0.0)) return Double.NEGATIVE_INFINITY;
            final double x = Math.log(r);
            xOld[i] = x;
            u[i] = (x - muUC) / s;
        }

        final double[] xNew = new double[nEdges];
        final Sum sum = new Sum();

        try {
            fillLogRatesFromU(tree.getRoot(), rootLog(), u, xNew, sig2, minDt, sum);
        } catch (ArithmeticException bad) {
            return Double.NEGATIVE_INFINITY;
        }


        double sumDelta = 0.0;
        for (int i = 0; i < nEdges; i++) sumDelta += (xNew[i] - xOld[i]);
        final double logH = sumDelta + 0.5 * sum.sumLogVar - nEdges * Math.log(s);

        rates.startEditing(this);
        for (int i = 0; i < nEdges; i++) {
            rates.setValue(i, Math.exp(xNew[i]));
        }

        indicator.startEditing(this);
        indicator.setValue(0, 1);

        return logH;
    }


    private double proposeACtoUC() {
        final double s = ucldStdev.getValue(0);
        final double sig2 = sigma2.getValue(0);
        final double minDt = minBranchLengthInput.get();
        if (!(s > 0.0) || !(sig2 > 0.0) || !(minDt > 0.0)) return Double.NEGATIVE_INFINITY;

        final int nEdges = rates.getDimension();
        final double muUC = -0.5 * s * s;

        final double[] xOld = new double[nEdges];
        for (int i = 0; i < nEdges; i++) {
            final double r = rates.getValue(i);
            if (!(r > 0.0)) return Double.NEGATIVE_INFINITY;
            xOld[i] = Math.log(r);
        }

        final double[] u = new double[nEdges];
        final Sum sum = new Sum();

        try {
            fillUFromLogRatesAC(tree.getRoot(), rootLog(), xOld, u, sig2, minDt, sum);
        } catch (ArithmeticException bad) {
            return Double.NEGATIVE_INFINITY;
        }

        final double[] xNew = new double[nEdges];
        for (int i = 0; i < nEdges; i++) {
            xNew[i] = muUC + s * u[i];
        }


        double sumDelta = 0.0;
        for (int i = 0; i < nEdges; i++) sumDelta += (xNew[i] - xOld[i]);
        final double logH = sumDelta + nEdges * Math.log(s) - 0.5 * sum.sumLogVar;

        // Commit
        rates.startEditing(this);
        for (int i = 0; i < nEdges; i++) {
            rates.setValue(i, Math.exp(xNew[i]));
        }

        indicator.startEditing(this);
        indicator.setValue(0, 0);

        return logH;
    }

    private void fillLogRatesFromU(final Node parent,
                                   final double logPar,
                                   final double[] u,
                                   final double[] xOut,
                                   final double sigma2,
                                   final double minDt,
                                   final Sum sum) {
        final int cc = parent.getChildCount();
        for (int c = 0; c < cc; c++) {
            final Node child = parent.getChild(c);
            final int idx = idxMap[child.getNr()];
            if (idx < 0) throw new ArithmeticException("Bad idxMap for child");

            final double dt = child.getLength();
            if (!(dt > minDt)) throw new ArithmeticException("dt too small");
            final double var = sigma2 * dt;
            if (!(var > 0.0)) throw new ArithmeticException("var<=0");

            sum.sumLogVar += Math.log(var);

            final double logChild = logPar - 0.5 * var + Math.sqrt(var) * u[idx];
            xOut[idx] = logChild;

            fillLogRatesFromU(child, logChild, u, xOut, sigma2, minDt, sum);
        }
    }

    private void fillUFromLogRatesAC(final Node parent,
                                     final double logPar,
                                     final double[] x,
                                     final double[] uOut,
                                     final double sigma2,
                                     final double minDt,
                                     final Sum sum) {
        final int cc = parent.getChildCount();
        for (int c = 0; c < cc; c++) {
            final Node child = parent.getChild(c);
            final int idx = idxMap[child.getNr()];
            if (idx < 0) throw new ArithmeticException("Bad idxMap for child");

            final double dt = child.getLength();
            if (!(dt > minDt)) throw new ArithmeticException("dt too small");
            final double var = sigma2 * dt;
            if (!(var > 0.0)) throw new ArithmeticException("var<=0");

            sum.sumLogVar += Math.log(var);

            final double mean = logPar - 0.5 * var;
            uOut[idx] = (x[idx] - mean) / Math.sqrt(var);

            fillUFromLogRatesAC(child, x[idx], x, uOut, sigma2, minDt, sum);
        }
    }
}
