package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import mixture.beast.evolution.mixture.RelaxedRatesPriorSVS;
import mixture.beast.evolution.util.BranchRateIndexHelper;

import mixture.lphy.evolution.auto.SVSRawBranchRates;

import mixture.beast.evolution.operator.ACSubtreeUIncrementOperator;
import mixture.beast.evolution.operator.ACSigma2NonCenteredOperator;
import mixture.beast.evolution.operator.IndicatorGibbsOperator;
import mixture.beast.evolution.operator.SingleRateScaleOperator;
import mixture.beast.evolution.operator.SubtreeRateScaleOperator;
import mixture.beast.evolution.operator.UCACSwitchBridgeOperator;
import mixture.beast.evolution.operator.UCLDStdevNonCenteredOperator;
import mixture.lphybeast.tobeast.TypedParameterUtils;

import java.util.ArrayList;
import java.util.List;

public class SVSRawBranchRatesToBEAST implements GeneratorToBEAST<SVSRawBranchRates, RelaxedRatesPriorSVS> {

    @Override
    public RelaxedRatesPriorSVS generatorToBEAST(SVSRawBranchRates dist, BEASTInterface beastValue, BEASTContext context) {

        @SuppressWarnings("unchecked")
        Value<TimeTree> treeVal = (Value<TimeTree>) dist.getParams().get(SVSRawBranchRates.TREE);

        @SuppressWarnings("unchecked")
        Value<Integer> indVal = (Value<Integer>) dist.getParams().get(SVSRawBranchRates.INDICATOR);

        @SuppressWarnings("unchecked")
        Value<Double> ucldVal = (Value<Double>) dist.getParams().get(SVSRawBranchRates.UCLD_STDEV);

        @SuppressWarnings("unchecked")
        Value<Double> sigma2Val = (Value<Double>) dist.getParams().get(SVSRawBranchRates.SIGMA2);

        @SuppressWarnings("unchecked")
        Value<Double> rootVal = (Value<Double>) dist.getParams().get(SVSRawBranchRates.ROOT_LOG_RATE);

        boolean indicatorIsRandom = indVal instanceof RandomVariable;
        boolean ucldStdevIsRandom = TypedParameterUtils.isRandom(ucldVal);
        boolean sigma2IsRandom = TypedParameterUtils.isRandom(sigma2Val);
        boolean rootLogRateIsRandom = TypedParameterUtils.isRandom(rootVal);

        TreeInterface beastTree = (TreeInterface) context.getBEASTObject(treeVal);
        Tree beastTreeImpl = (Tree) beastTree;
        Value<?> outVal = context.getOutput(dist);

        int nNodes = beastTreeImpl.getNodeCount();
        int targetDim = nNodes - 1;

        // ---- rates parameter ----
        RealVectorParam<?> ratesParam = positiveRatesFromLPhyValue(
                "rawRates." + dist.getUniqueId(),
                beastValue,
                beastTreeImpl,
                targetDim,
                true
        );
        TypedParameterUtils.replaceInContext(context, outVal, beastValue, ratesParam);

        // ---- indicator parameter ----
        BEASTInterface indObj = context.getBEASTObject(indVal);
        IntScalarParam<?> indParam = TypedParameterUtils.intScalarFrom(
                indObj,
                "indicator." + dist.getUniqueId(),
                indVal.value() != null ? indVal.value() : 0,
                indicatorIsRandom
        );
        TypedParameterUtils.replaceInContext(context, indVal, indObj, indParam);
        if (!indicatorIsRandom) {
            prepareFixedOperatorIndicator(context, indParam);
        }

        // ---- hyperparameters ----
        RealScalarParam<?> ucldStdev = realScalarParam(context, ucldVal, "ucldStdev");
        RealScalarParam<?> sigma2 = realScalarParam(context, sigma2Val, "sigma2");
        RealScalarParam<?> rootLogRate = realScalarParam(context, rootVal, "rootLogRate");
        prepareOperatorScalar(context, ucldStdev, ucldStdevIsRandom);
        prepareOperatorScalar(context, sigma2, sigma2IsRandom);
        prepareOperatorScalar(context, rootLogRate, rootLogRateIsRandom);

        // ---- SVS prior ----
        RelaxedRatesPriorSVS prior = new RelaxedRatesPriorSVS();
        prior.setID("SVSRelaxedClockPrior." + dist.getUniqueId());

        prior.setInputValue("tree", beastTree);
        prior.setInputValue("ratesVector", ratesParam);
        prior.setInputValue("indicatorScalar", indParam);
        prior.setInputValue("ucldStdevScalar", ucldStdev);
        prior.setInputValue("sigma2Scalar", sigma2);
        prior.setInputValue("rootLogRateScalar", rootLogRate);

        prior.initAndValidate();
        context.addBEASTObject(prior, dist);

        for (final Operator operator : createOperatorSchedule(dist.getUniqueId(), beastTree, ratesParam, indParam,
                indicatorIsRandom, prior, ucldStdev, ucldStdevIsRandom, sigma2, sigma2IsRandom, rootLogRate)) {
            context.addExtraOperator(operator);
        }

        context.addSkipOperator(ratesParam);
        if (indicatorIsRandom) {
            context.addSkipOperator(indParam);
        }
        context.addSkipLoggable(ratesParam);

        return prior;
    }

    private static double[] ones(final int dimension) {
        final double[] values = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            values[i] = 1.0;
        }
        return values;
    }

    static RealVectorParam<?> positiveRatesFromLPhyValue(final String fallbackId,
                                                         final BEASTInterface beastValue,
                                                         final Tree beastTree,
                                                         final int targetDimension,
                                                         final boolean estimate) {
        final double[] values = nonRootRateValues(beastValue, beastTree, targetDimension);
        return TypedParameterUtils.positiveRealVector(idOrFallback(beastValue, fallbackId), values, estimate);
    }

    private static double[] nonRootRateValues(final BEASTInterface beastValue,
                                              final Tree beastTree,
                                              final int targetDimension) {
        if (beastValue instanceof RealVectorParam<?> typed) {
            return nonRootRateValues(typed, beastTree, targetDimension);
        }
        if (beastValue instanceof RealParameter legacy) {
            return nonRootRateValues(legacy, beastTree, targetDimension);
        }
        return ones(targetDimension);
    }

    private static double[] nonRootRateValues(final RealVectorParam<?> typed,
                                              final Tree beastTree,
                                              final int targetDimension) {
        final int observed = typed.size();
        if (observed == 1) {
            return fill(targetDimension, typed.get(0));
        }
        if (observed == targetDimension) {
            final double[] values = new double[targetDimension];
            for (int i = 0; i < targetDimension; i++) {
                values[i] = typed.get(i);
            }
            return values;
        }
        if (observed == beastTree.getNodeCount()) {
            return stripRootSlot(beastTree, observed, typed::get);
        }
        return ones(targetDimension);
    }

    private static double[] nonRootRateValues(final RealParameter legacy,
                                              final Tree beastTree,
                                              final int targetDimension) {
        final int observed = legacy.getDimension();
        if (observed == 1) {
            return fill(targetDimension, legacy.getValue(0));
        }
        if (observed == targetDimension) {
            final double[] values = new double[targetDimension];
            for (int i = 0; i < targetDimension; i++) {
                values[i] = legacy.getValue(i);
            }
            return values;
        }
        if (observed == beastTree.getNodeCount()) {
            return stripRootSlot(beastTree, observed, legacy::getValue);
        }
        return ones(targetDimension);
    }

    private static double[] stripRootSlot(final Tree beastTree,
                                          final int observed,
                                          final IndexedDouble values) {
        final BranchRateIndexHelper.Mapping mapping = BranchRateIndexHelper.buildDeterministic(beastTree);
        final double[] nonRoot = new double[beastTree.getNodeCount() - 1];
        for (int i = 0; i < beastTree.getNodeCount(); i++) {
            final Node node = beastTree.getNode(i);
            final int nr = node.getNr();
            if (nr >= observed) {
                return ones(nonRoot.length);
            }
            final int rateIndex = mapping.idxForNodeNr(nr);
            if (rateIndex >= 0) {
                nonRoot[rateIndex] = values.get(nr);
            }
        }
        return nonRoot;
    }

    private static double[] fill(final int dimension, final double value) {
        final double[] values = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            values[i] = value;
        }
        return values;
    }

    private static String idOrFallback(final BEASTInterface beastValue, final String fallbackId) {
        return beastValue != null && beastValue.getID() != null ? beastValue.getID() : fallbackId;
    }

    private static RealScalarParam<?> realScalarParam(final BEASTContext context,
                                                      final Value<Double> value,
                                                      final String name) {
        final Object scalar = context.getAsRealScalar(value);
        if (scalar instanceof RealScalarParam<?> parameter) {
            return parameter;
        }
        throw new IllegalArgumentException("SVSRawBranchRatesToBEAST: " + name
                + " must map to RealScalarParam for operator scheduling, got "
                + scalar.getClass().getName());
    }

    private static void prepareOperatorScalar(final BEASTContext context,
                                              final RealScalarParam<?> scalar,
                                              final boolean estimate) {
        scalar.isEstimatedInput.setValue(estimate, scalar);
        if (estimate && !context.getState().contains(scalar)) {
            context.getState().add(scalar);
        }
        if (!estimate) {
            context.getState().remove(scalar);
            context.addSkipOperator(scalar);
            context.addSkipLoggable(scalar);
        }
    }

    private static void prepareFixedOperatorIndicator(final BEASTContext context,
                                                      final IntScalarParam<?> indicator) {
        indicator.isEstimatedInput.setValue(false, indicator);
        if (!context.getState().contains(indicator)) {
            context.getState().add(indicator);
        }
        context.addSkipOperator(indicator);
        context.addSkipLoggable(indicator);
    }

    @FunctionalInterface
    private interface IndexedDouble {
        double get(int index);
    }

    static List<Operator> createOperatorSchedule(final String uid,
                                                 final TreeInterface beastTree,
                                                 final RealVectorParam<?> ratesParam,
                                                 final IntScalarParam<?> indParam,
                                                 final boolean indicatorIsRandom,
                                                 final RelaxedRatesPriorSVS prior,
                                                 final RealScalarParam<?> ucldStdev,
                                                 final boolean ucldStdevIsRandom,
                                                 final RealScalarParam<?> sigma2,
                                                 final boolean sigma2IsRandom,
                                                 final RealScalarParam<?> rootLogRate) {
        final List<Operator> operators = new ArrayList<>();
        final int indicator = indParam.get();

        SingleRateScaleOperator oneScale = new SingleRateScaleOperator();
        oneScale.setID("rawRates.oneScale." + uid);
        oneScale.setInputValue("ratesVector", ratesParam);
        oneScale.setInputValue("window", 0.3);
        oneScale.setInputValue("weight", 15.0);
        oneScale.initAndValidate();
        operators.add(oneScale);

        SubtreeRateScaleOperator subScale = new SubtreeRateScaleOperator();
        subScale.setID("rawRates.subtreeScale." + uid);
        subScale.setInputValue("tree", beastTree);
        subScale.setInputValue("ratesVector", ratesParam);
        subScale.setInputValue("window", 0.25);
        subScale.setInputValue("weight", 10.0);
        subScale.initAndValidate();
        operators.add(subScale);

        if (indicatorIsRandom) {
            IndicatorGibbsOperator gibbs = new IndicatorGibbsOperator();
            gibbs.setID("indicator.gibbs." + uid);
            gibbs.setInputValue("indicatorScalar", indParam);
            gibbs.setInputValue("prior", prior);
            gibbs.setInputValue("pOne", 0.5);
            gibbs.setInputValue("weight", 2.0);
            gibbs.initAndValidate();
            operators.add(gibbs);
        }

        if (indicatorIsRandom || indicator == 1) {
            ACSubtreeUIncrementOperator acSubtree = new ACSubtreeUIncrementOperator();
            acSubtree.setID("acSubtreeU." + uid);
            acSubtree.setInputValue("tree", beastTree);
            acSubtree.setInputValue("ratesVector", ratesParam);
            acSubtree.setInputValue("indicatorScalar", indParam);
            acSubtree.setInputValue("sigma2Scalar", sigma2);
            acSubtree.setInputValue("rootLogRateScalar", rootLogRate);
            acSubtree.setInputValue("minBranchLength", 1e-12);
            acSubtree.setInputValue("delta", 0.25);
            acSubtree.setInputValue("internalOnly", true);
            acSubtree.setInputValue("allowRoot", false);
            acSubtree.setInputValue("maxSubtreeEdges", 80);
            acSubtree.setInputValue("rejectIfNotAC", true);
            acSubtree.setInputValue("weight", 15.0);
            acSubtree.initAndValidate();
            operators.add(acSubtree);

            if (sigma2IsRandom) {
                ACSigma2NonCenteredOperator sigma2NC = new ACSigma2NonCenteredOperator();
                sigma2NC.setID("acSigma2NC." + uid);
                sigma2NC.setInputValue("tree", beastTree);
                sigma2NC.setInputValue("ratesVector", ratesParam);
                sigma2NC.setInputValue("indicatorScalar", indParam);
                sigma2NC.setInputValue("sigma2Scalar", sigma2);
                sigma2NC.setInputValue("rootLogRateScalar", rootLogRate);
                sigma2NC.setInputValue("minBranchLength", 1e-12);
                sigma2NC.setInputValue("window", 0.15);
                sigma2NC.setInputValue("rejectIfNotAC", true);
                sigma2NC.setInputValue("weight", 15.0);
                sigma2NC.initAndValidate();
                operators.add(sigma2NC);
            }
        }

        if (ucldStdevIsRandom && (indicatorIsRandom || indicator == 0)) {
            UCLDStdevNonCenteredOperator ucldNC = new UCLDStdevNonCenteredOperator();
            ucldNC.setID("ucldSigmaNC." + uid);
            ucldNC.setInputValue("ratesVector", ratesParam);
            ucldNC.setInputValue("indicatorScalar", indParam);
            ucldNC.setInputValue("ucldStdevScalar", ucldStdev);
            ucldNC.setInputValue("window", 0.20);
            ucldNC.setInputValue("rejectIfNotUC", true);
            ucldNC.setInputValue("weight", 15.0);
            ucldNC.initAndValidate();
            operators.add(ucldNC);
        }

        if (indicatorIsRandom) {
            UCACSwitchBridgeOperator bridge = new UCACSwitchBridgeOperator();
            bridge.setID("ucacBridge." + uid);
            bridge.setInputValue("tree", beastTree);
            bridge.setInputValue("ratesVector", ratesParam);
            bridge.setInputValue("indicatorScalar", indParam);
            bridge.setInputValue("ucldStdevScalar", ucldStdev);
            bridge.setInputValue("sigma2Scalar", sigma2);
            bridge.setInputValue("rootLogRateScalar", rootLogRate);
            bridge.setInputValue("minBranchLength", 1e-12);
            bridge.setInputValue("weight", 10.0);
            bridge.initAndValidate();
            operators.add(bridge);
        }

        return operators;
    }

    @Override public Class<SVSRawBranchRates> getGeneratorClass() { return SVSRawBranchRates.class; }
    @Override public Class<RelaxedRatesPriorSVS> getBEASTClass() { return RelaxedRatesPriorSVS.class; }
}
