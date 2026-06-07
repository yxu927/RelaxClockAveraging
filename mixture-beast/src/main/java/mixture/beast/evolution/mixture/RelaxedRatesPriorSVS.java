package mixture.beast.evolution.mixture;

import beast.base.core.Description;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.IntScalar;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
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
            "rates",
            "Legacy BEAST parameter containing positive branch rates for NON-root nodes; "
                    + "dimension=(tree.nodeCount - 1).",
            Input.Validate.OPTIONAL);

    public final Input<RealVector> ratesVectorInput = new Input<>(
            "ratesVector",
            "BEAST3 typed real vector containing positive branch rates for NON-root nodes.",
            Input.Validate.OPTIONAL);

    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator",
            "Legacy integer parameter, dimension=1; 0=uncorrelated, 1=autocorrelated.",
            Input.Validate.OPTIONAL);

    public final Input<IntScalar> indicatorScalarInput = new Input<>(
            "indicatorScalar",
            "BEAST3 typed integer scalar; 0=uncorrelated, 1=autocorrelated.",
            Input.Validate.OPTIONAL);

    public final Input<RealParameter> ucldStdevInput = new Input<>(
            "ucldStdev",
            "Legacy lognormal stdev parameter for UC prior, dimension=1.",
            Input.Validate.OPTIONAL);

    public final Input<RealScalar> ucldStdevScalarInput = new Input<>(
            "ucldStdevScalar",
            "BEAST3 typed real scalar for lognormal stdev in UC prior.",
            Input.Validate.OPTIONAL);

    public final Input<RealParameter> rootLogRateInput = new Input<>(
            "rootLogRate",
            "Legacy optional root log-rate for AC; if absent, root log-rate = 0.",
            Input.Validate.OPTIONAL);

    public final Input<RealScalar> rootLogRateScalarInput = new Input<>(
            "rootLogRateScalar",
            "BEAST3 typed optional root log-rate for AC; if absent, root log-rate = 0.",
            Input.Validate.OPTIONAL);

    public final Input<RealParameter> sigma2Input = new Input<>(
            "sigma2",
            "Legacy Brownian variance parameter, dimension=1.",
            Input.Validate.OPTIONAL);

    public final Input<RealScalar> sigma2ScalarInput = new Input<>(
            "sigma2Scalar",
            "BEAST3 typed real scalar for Brownian variance parameter.",
            Input.Validate.OPTIONAL);

    public final Input<Double> minBranchLengthInput = new Input<>(
            "minBranchLength", "minimum branch length (time) allowed in AC; if shorter, logP=-inf.", 1e-12);

    private static final double LOG_2PI = Math.log(2.0 * Math.PI);

    private Tree tree;
    private RealParameter legacyRates;
    private RealVector typedRates;
    private IntegerParameter legacyIndicator;
    private IntScalar typedIndicator;
    private RealParameter legacyUcldStdev;
    private RealScalar typedUcldStdev;
    private RealParameter legacyRootLogRate;
    private RealScalar typedRootLogRate;
    private RealParameter legacySigma2;
    private RealScalar typedSigma2;

    private BranchRateIndexHelper.Mapping mapping;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        legacyRates = ratesInput.get();
        typedRates = ratesVectorInput.get();
        legacyIndicator = indicatorInput.get();
        typedIndicator = indicatorScalarInput.get();
        legacyUcldStdev = ucldStdevInput.get();
        typedUcldStdev = ucldStdevScalarInput.get();
        legacyRootLogRate = rootLogRateInput.get();
        typedRootLogRate = rootLogRateScalarInput.get();
        legacySigma2 = sigma2Input.get();
        typedSigma2 = sigma2ScalarInput.get();

        requireExactlyOne(legacyRates, typedRates, "rates", "ratesVector");
        requireExactlyOne(legacyIndicator, typedIndicator, "indicator", "indicatorScalar");
        requireExactlyOne(legacyUcldStdev, typedUcldStdev, "ucldStdev", "ucldStdevScalar");
        requireAtMostOne(legacyRootLogRate, typedRootLogRate, "rootLogRate", "rootLogRateScalar");
        requireExactlyOne(legacySigma2, typedSigma2, "sigma2", "sigma2Scalar");

        if (legacyIndicator != null && legacyIndicator.getDimension() != 1) {
            throw new IllegalArgumentException("indicator must have dimension=1.");
        }
        if (legacyUcldStdev != null && legacyUcldStdev.getDimension() != 1) {
            throw new IllegalArgumentException("ucldStdev must have dimension=1.");
        }
        if (legacySigma2 != null && legacySigma2.getDimension() != 1) {
            throw new IllegalArgumentException("sigma2 must have dimension=1.");
        }
        if (legacyRootLogRate != null && legacyRootLogRate.getDimension() != 1) {
            throw new IllegalArgumentException("rootLogRate must have dimension=1.");
        }

        validateOrExpandRatesDimension();
        mapping = BranchRateIndexHelper.buildDeterministic(tree);
    }

    private static void requireExactlyOne(final Object legacy,
                                          final Object typed,
                                          final String legacyName,
                                          final String typedName) {
        if (legacy == null && typed == null) {
            throw new IllegalArgumentException("Either " + legacyName + " or " + typedName + " must be specified.");
        }
        requireAtMostOne(legacy, typed, legacyName, typedName);
    }

    private static void requireAtMostOne(final Object legacy,
                                         final Object typed,
                                         final String legacyName,
                                         final String typedName) {
        if (legacy != null && typed != null) {
            throw new IllegalArgumentException("Specify only one of " + legacyName + " or " + typedName + ".");
        }
    }

    private void ensureMappingUpToDate() {
        validateOrExpandRatesDimension();
        if (mapping == null || !mapping.matches(tree)) {
            mapping = BranchRateIndexHelper.buildDeterministic(tree);
        }
    }

    private int expectedRateDimension() {
        return tree.getNodeCount() - 1;
    }

    private int rateDimension() {
        if (legacyRates != null) {
            return legacyRates.getDimension();
        }
        if (typedRates.size() == 1 && expectedRateDimension() > 1) {
            return expectedRateDimension();
        }
        return typedRates.size();
    }

    private double rateValue(final int i) {
        if (legacyRates != null) {
            return legacyRates.getValue(i);
        }
        if (typedRates.size() == 1) {
            return typedRates.get(0);
        }
        return typedRates.get(i);
    }

    private int indicatorValue() {
        if (legacyIndicator != null) {
            return legacyIndicator.getValue(0);
        }
        return typedIndicator.get();
    }

    private double ucldStdevValue() {
        if (legacyUcldStdev != null) {
            return legacyUcldStdev.getValue(0);
        }
        return typedUcldStdev.get();
    }

    private boolean hasRootLogRate() {
        return legacyRootLogRate != null || typedRootLogRate != null;
    }

    private double rootLogRateValue() {
        if (legacyRootLogRate != null) {
            return legacyRootLogRate.getValue(0);
        }
        if (typedRootLogRate != null) {
            return typedRootLogRate.get();
        }
        return 0.0;
    }

    private double sigma2Value() {
        if (legacySigma2 != null) {
            return legacySigma2.getValue(0);
        }
        return typedSigma2.get();
    }

    private boolean ratesDirty() {
        return legacyRates != null ? legacyRates.somethingIsDirty() : calculationNodeDirty(typedRates);
    }

    private boolean indicatorDirty() {
        return legacyIndicator != null ? legacyIndicator.somethingIsDirty() : calculationNodeDirty(typedIndicator);
    }

    private boolean ucldStdevDirty() {
        return legacyUcldStdev != null ? legacyUcldStdev.somethingIsDirty()
                : calculationNodeDirty(typedUcldStdev);
    }

    private boolean rootLogRateDirty() {
        if (legacyRootLogRate != null) {
            return legacyRootLogRate.somethingIsDirty();
        }
        return calculationNodeDirty(typedRootLogRate);
    }

    private boolean sigma2Dirty() {
        return legacySigma2 != null ? legacySigma2.somethingIsDirty() : calculationNodeDirty(typedSigma2);
    }

    private static boolean calculationNodeDirty(final Object object) {
        return object instanceof CalculationNode && ((CalculationNode) object).somethingIsDirty();
    }

    private static String idOf(final Object object) {
        if (object instanceof BEASTInterface) {
            return ((BEASTInterface) object).getID();
        }
        return null;
    }

    private Object ratesObject() {
        return legacyRates != null ? legacyRates : typedRates;
    }

    private Object indicatorObject() {
        return legacyIndicator != null ? legacyIndicator : typedIndicator;
    }

    private Object ucldStdevObject() {
        return legacyUcldStdev != null ? legacyUcldStdev : typedUcldStdev;
    }

    private Object rootLogRateObject() {
        return legacyRootLogRate != null ? legacyRootLogRate : typedRootLogRate;
    }

    private Object sigma2Object() {
        return legacySigma2 != null ? legacySigma2 : typedSigma2;
    }

    private void validateOrExpandRatesDimension() {
        if (legacyRates != null) {
            BranchRateIndexHelper.validateRatesDimension(tree, legacyRates, "RelaxedRatesPriorSVS");
            return;
        }

        final int expected = expectedRateDimension();
        if (expected < 1) {
            throw new IllegalArgumentException("RelaxedRatesPriorSVS: tree must contain at least one non-root branch. "
                    + "Found nodeCount=" + tree.getNodeCount());
        }

        final int observed = typedRates.size();
        if (observed == expected) {
            return;
        }

        if (observed == 1) {
            if (typedRates instanceof RealVectorParam<?>) {
                ((RealVectorParam<?>) typedRates).setDimension(expected);
            }
            return;
        }

        throw new IllegalArgumentException("RelaxedRatesPriorSVS: ratesVector must have dimension "
                + "(tree.nodeCount - 1). Found " + observed + " vs " + expected
                + ". Scalar typed rates can be expanded or broadcast across all non-root branches.");
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

        final int k = indicatorValue();
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
        final double s = ucldStdevValue();
        if (!(s > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double var = s * s;
        final double meanLog = -0.5 * var;

        double lp = 0.0;
        for (int i = 0; i < rateDimension(); i++) {
            lp += logLogNormalRate(rateValue(i), meanLog, var);
            if (Double.isInfinite(lp)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return lp;
    }

    /** AC: per-branch lognormal increments with mean-correction (martingale): E[r_child|r_parent]=r_parent. */
    public double logPriorACOnly() {
        ensureMappingUpToDate();

        final double s2 = sigma2Value();
        if (!(s2 > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double minDt = minBranchLengthInput.get();
        if (!(minDt > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double rootLog = rootLogRateValue();

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

            final double rChi = rateValue(idxChi);
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

                final double rPar = rateValue(idxPar);
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
        return Collections.singletonList(idOf(ratesObject()));
    }

    @Override
    public List<String> getConditions() {
        if (!hasRootLogRate()) {
            return Arrays.asList(tree.getID(), idOf(indicatorObject()), idOf(ucldStdevObject()), idOf(sigma2Object()));
        }
        return Arrays.asList(tree.getID(), idOf(indicatorObject()), idOf(ucldStdevObject()),
                idOf(rootLogRateObject()), idOf(sigma2Object()));
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
        if (ratesDirty()) {
            dirty = true;
        }
        if (indicatorDirty()) {
            dirty = true;
        }
        if (ucldStdevDirty()) {
            dirty = true;
        }
        if (sigma2Dirty()) {
            dirty = true;
        }
        if (rootLogRateDirty()) {
            dirty = true;
        }

        return dirty;
    }
}
