package mixture.beast.evolution.mixture;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import mixture.beast.evolution.util.BranchRateIndexHelper;

public class SharedRatesClockModelSpec extends Base {

    public final Input<Tree> treeInput = new Input<>(
            "tree",
            "tree used by this clock",
            Input.Validate.REQUIRED
    );

    public final Input<RealParameter> ratesInput = new Input<>(
            "rates",
            "Legacy BEAST parameter containing one rate multiplier per non-root branch.",
            Input.Validate.OPTIONAL
    );

    public final Input<RealVector> ratesVectorInput = new Input<>(
            "ratesVector",
            "BEAST3 typed real vector containing one rate multiplier per non-root branch.",
            Input.Validate.OPTIONAL
    );

    public final Input<Boolean> normalizeInput = new Input<>(
            "normalize",
            "if true: rescale time-weighted mean rate across the tree to 1",
            false
    );

    public final Input<RealParameter> legacyMeanRateInput = new Input<>(
            "meanRate",
            "Legacy BEAST parameter for the overall mean substitution rate.",
            Input.Validate.OPTIONAL
    );

    public final Input<RealScalar> meanRateScalarInput = new Input<>(
            "meanRateScalar",
            "BEAST3 typed real scalar for the overall mean substitution rate.",
            Input.Validate.OPTIONAL
    );

    private Tree tree;
    private RealParameter legacyRates;
    private RealVector typedRates;
    private RealParameter legacyMeanRate;
    private RealScalar typedMeanRate;
    private boolean doNormalize;

    private BranchRateIndexHelper.Mapping mapping;

    private double scaleFactor = 1.0;
    private double storedScaleFactor = 1.0;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        doNormalize = normalizeInput.get();

        legacyRates = ratesInput.get();
        typedRates = ratesVectorInput.get();
        if (legacyRates == null && typedRates == null) {
            throw new IllegalArgumentException("Either rates or ratesVector must be specified.");
        }
        if (legacyRates != null && typedRates != null) {
            throw new IllegalArgumentException("Specify only one of rates or ratesVector.");
        }

        legacyMeanRate = legacyMeanRateInput.get();
        typedMeanRate = meanRateScalarInput.get();
        if (legacyMeanRate != null && typedMeanRate != null) {
            throw new IllegalArgumentException("Specify only one of meanRate or meanRateScalar.");
        }
        if ((legacyMeanRate != null || typedMeanRate != null) && meanRateInput.get() != null) {
            throw new IllegalArgumentException("Specify only one of clock.rate, meanRate or meanRateScalar.");
        }
        if (legacyMeanRate != null && legacyMeanRate.getDimension() != 1) {
            throw new IllegalArgumentException("meanRate must have dimension=1.");
        }

        validateOrExpandRatesDimension();
        mapping = BranchRateIndexHelper.buildDeterministic(tree);

        if (doNormalize) {
            computeScaleFactor();
        } else {
            scaleFactor = 1.0;
        }

        Log.info.println("SharedRatesClockModelSpec init: nNodes=" + mapping.getNodeCount()
                + ", normalize=" + doNormalize);
    }

    private void ensureMappingUpToDate() {
        validateOrExpandRatesDimension();
        if (mapping == null || !mapping.matches(tree)) {
            mapping = BranchRateIndexHelper.buildDeterministic(tree);
        }
    }

    private int rateDimension() {
        if (legacyRates != null) {
            return legacyRates.getDimension();
        }
        return typedRates.size();
    }

    private double rawRateValue(final int i) {
        if (legacyRates != null) {
            return legacyRates.getValue(i);
        }
        if (typedRates.size() == 1) {
            return typedRates.get(0);
        }
        return typedRates.get(i);
    }

    private double meanRateValue() {
        if (legacyMeanRate != null) {
            return legacyMeanRate.getValue();
        }
        if (typedMeanRate != null) {
            return typedMeanRate.get();
        }
        final RealScalar<PositiveReal> clockRate = meanRateInput.get();
        if (clockRate != null) {
            return clockRate.get();
        }
        return 1.0;
    }

    private boolean ratesDirty() {
        if (legacyRates != null) {
            return legacyRates.somethingIsDirty();
        }
        return typedRates instanceof CalculationNode && ((CalculationNode) typedRates).somethingIsDirty();
    }

    private boolean meanRateDirty() {
        if (legacyMeanRate != null) {
            return legacyMeanRate.somethingIsDirty();
        }
        if (typedMeanRate instanceof CalculationNode) {
            return ((CalculationNode) typedMeanRate).somethingIsDirty();
        }
        final RealScalar<PositiveReal> clockRate = meanRateInput.get();
        return clockRate instanceof CalculationNode && ((CalculationNode) clockRate).somethingIsDirty();
    }

    private void validateOrExpandRatesDimension() {
        if (legacyRates != null) {
            BranchRateIndexHelper.validateRatesDimension(tree, legacyRates, "SharedRatesClockModelSpec");
            return;
        }

        final int expected = tree.getNodeCount() - 1;
        if (expected < 1) {
            throw new IllegalArgumentException("SharedRatesClockModelSpec: tree must contain at least one non-root branch. "
                    + "Found nodeCount=" + tree.getNodeCount());
        }

        final int observed = rateDimension();
        if (observed == expected) {
            return;
        }

        if (observed == 1) {
            if (typedRates instanceof RealVectorParam<?>) {
                ((RealVectorParam<?>) typedRates).setDimension(expected);
            }
            return;
        }

        throw new IllegalArgumentException("SharedRatesClockModelSpec: ratesVector must have dimension "
                + "(tree.nodeCount - 1). Found " + observed + " vs " + expected
                + ". Scalar typed rates can be expanded or broadcast across all non-root branches.");
    }

    @Override
    public double getRateForBranch(final Node node) {
        if (node.isRoot()) {
            return 1.0;
        }

        ensureMappingUpToDate();

        final double mr = meanRateValue();
        if (!(mr > 0.0)) {
            return 0.0;
        }

        final int idx = mapping.idxForNode(node);
        if (idx < 0) {
            return 0.0;
        }

        final double r = rawRateValue(idx);
        if (!(r > 0.0)) {
            return 0.0;
        }

        return r * scaleFactor * mr;
    }

    private void computeScaleFactor() {
        ensureMappingUpToDate();

        double sumRateTime = 0.0;
        double sumTime = 0.0;

        final int nNodes = mapping.getNodeCount();
        for (int i = 0; i < nNodes; i++) {
            final Node node = tree.getNode(i);
            if (node.isRoot()) {
                continue;
            }

            final double dt = node.getLength();
            if (!(dt > 0.0)) {
                continue;
            }

            final int idx = mapping.idxForNode(node);
            final double r = rawRateValue(idx);
            if (!(r > 0.0)) {
                continue;
            }

            sumRateTime += r * dt;
            sumTime += dt;
        }

        if (!(sumRateTime > 0.0) || !(sumTime > 0.0)) {
            scaleFactor = 1.0;
        } else {
            scaleFactor = sumTime / sumRateTime;
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        boolean dirty = false;

        if (tree != null && tree.somethingIsDirty()) {
            dirty = true;
            ensureMappingUpToDate();
        }

        if (ratesDirty() || meanRateDirty()) {
            dirty = true;
        }

        if (dirty && doNormalize) {
            computeScaleFactor();
        }

        return dirty;
    }

    @Override
    protected void store() {
        storedScaleFactor = scaleFactor;
        super.store();
    }

    @Override
    protected void restore() {
        scaleFactor = storedScaleFactor;
        super.restore();
    }
}
