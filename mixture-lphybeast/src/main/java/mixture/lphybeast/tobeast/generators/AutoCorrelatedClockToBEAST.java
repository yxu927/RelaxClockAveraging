package mixture.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.continuous.AutoCorrelatedClock;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.loggers.MetaDataTreeLogger;
import mixture.beast.evolution.clockmodel.auto.AutoCorrelatedClockModel;
import mixture.beast.evolution.operator.GlobalLogRateShiftOperator;

import java.lang.reflect.Method;

public class AutoCorrelatedClockToBEAST implements GeneratorToBEAST<AutoCorrelatedClock, AutoCorrelatedClockModel> {

    @Override
    public AutoCorrelatedClockModel generatorToBEAST(AutoCorrelatedClock accLphy,
                                                     BEASTInterface beastValue,
                                                     BEASTContext context) {

        AutoCorrelatedClockModel accModel = new AutoCorrelatedClockModel();
        accModel.setID("AutoCorrClockModel." + accLphy.getUniqueId());

        Value<TimeTree> treeValue = accLphy.getTree();
        Value<Double[]> nodeLogRatesValue = accLphy.getNodeLogRates();
        Value<Double> rootLogRateValue = accLphy.getRootLogRate();
        Value<Double> sigma2Value = accLphy.getSigma2();
        Value<Double> meanRateValue = accLphy.getMeanRate();
        Value<Boolean> normalizeValue = accLphy.getNormalize();
        Value<Integer> taylorOrderValue = accLphy.getTaylorOrder();

        TreeInterface beastTree = (TreeInterface) context.getBEASTObject(treeValue);

        RealParameter beastNodeLogRates = context.getAsRealParameter(nodeLogRatesValue);
        RealParameter beastRootLogRate  = context.getAsRealParameter(rootLogRateValue);
        RealParameter beastSigma2       = context.getAsRealParameter(sigma2Value);

        final RealParameter beastMeanRate;
        final boolean meanRateProvided;
        if (meanRateValue != null) {
            beastMeanRate = context.getAsRealParameter(meanRateValue);
            meanRateProvided = true;
        } else {
            RealParameter fixed = new RealParameter();
            fixed.setID("meanRate.fixed." + accLphy.getUniqueId());
            fixed.setInputValue("value", "1.0");
            fixed.setInputValue("estimate", false);
            fixed.initAndValidate();
            beastMeanRate = fixed;
            meanRateProvided = false;
        }

        boolean doNormalize = false;
        if (normalizeValue != null) {
            doNormalize = normalizeValue.value();
            accModel.setInputValue("normalize", doNormalize);
        }

        // *** 关键：normalize=true 时，防止 rootLogRate 被任何默认 operator 动到 ***
        if (doNormalize && beastRootLogRate instanceof StateNode) {
            context.addSkipOperator((StateNode) beastRootLogRate);
        }

        accModel.setInputValue("tree", beastTree);
        accModel.setInputValue("nodeRates", beastNodeLogRates);
        accModel.setInputValue("rootLogRate", beastRootLogRate);
        accModel.setInputValue("sigma2", beastSigma2);
        accModel.setInputValue("meanRate", beastMeanRate);

        if (taylorOrderValue != null) {
            accModel.setInputValue("taylorOrder", taylorOrderValue.value());
        }

        accModel.initAndValidate();

        MetaDataTreeLogger treeLogger = new MetaDataTreeLogger(accModel, beastTree, context);
        context.addExtraLogger(treeLogger);

        addAutoCorrMixingOperators(accLphy, context, beastTree,
                beastNodeLogRates, beastRootLogRate, beastSigma2,
                beastMeanRate, doNormalize, meanRateProvided, treeValue);

        return accModel;
    }

    private void addAutoCorrMixingOperators(AutoCorrelatedClock accLphy,
                                            BEASTContext context,
                                            TreeInterface beastTree,
                                            RealParameter nodeLogRates,
                                            RealParameter rootLogRate,
                                            RealParameter sigma2,
                                            RealParameter meanRate,
                                            boolean doNormalize,
                                            boolean meanRateProvided,
                                            Value<TimeTree> treeValue) {

        // normalize=false 才加全局 log-rate shift
        if (!doNormalize) {
            GlobalLogRateShiftOperator shiftOp = new GlobalLogRateShiftOperator();
            shiftOp.setID("globalLogRateShift." + accLphy.getUniqueId());
            shiftOp.setInputValue("rootLogRate", rootLogRate);
            shiftOp.setInputValue("nodeRates", nodeLogRates);
            shiftOp.setInputValue("windowSize", 0.10);
            shiftOp.setInputValue("weight", 1.0);
            shiftOp.initAndValidate();
            context.addExtraOperator(shiftOp);
        }

        if (!(beastTree instanceof Tree)) {
            throw new IllegalArgumentException("Expected beastTree to be beast.base.evolution.tree.Tree, got: "
                    + beastTree.getClass().getName());
        }
        Tree tree = (Tree) beastTree;

//        TreeLogRateScaleOperator scaleOp = new TreeLogRateScaleOperator();
//        scaleOp.setID("treeLogRateScale." + accLphy.getUniqueId());
//
//        scaleOp.setInputValue("tree", tree);
//        scaleOp.setInputValue("rootLogRate", rootLogRate);
//        scaleOp.setInputValue("nodeRates", nodeLogRates);
//
//        // 通用：用 treeLength ratio（异时间采样下更稳）
//        scaleOp.setInputValue("useTreeLength", true);
//        scaleOp.setInputValue("hyperUseTreeLength", true);
//
//        // normalize=false：shift logRates
//        // normalize=true ：不 shift（不能碰 rootLogRate），改为 scale meanRate（绝对尺度）
//        scaleOp.setInputValue("shiftLogRates", !doNormalize);
//
//        if (doNormalize && meanRateProvided) {
//            scaleOp.setInputValue("meanRate", meanRate);
//            scaleOp.setInputValue("scaleMeanRate", true);
//        } else {
//            scaleOp.setInputValue("scaleMeanRate", false);
//            if (doNormalize && !meanRateProvided) {
//                Log.warning.println("normalize=true but meanRate not provided (fixed=1.0); TreeLogRateScaleOperator will not scale meanRate.");
//            }
//        }
//
//        scaleOp.setInputValue("sigma2", sigma2);
//        scaleOp.setInputValue("scaleSigma2", true);
//
//        RealParameter lambda = tryGetLambdaFromTreeGenerator(treeValue, context);
//        if (lambda != null) {
//            scaleOp.setInputValue("lambda", lambda);
//            scaleOp.setInputValue("scaleLambda", true);
//        } else {
//            scaleOp.setInputValue("scaleLambda", false);
//        }
//
//        scaleOp.setInputValue("scaleFactor", 0.75);
//        scaleOp.setInputValue("weight", 3.0);
//        scaleOp.initAndValidate();
//        context.addExtraOperator(scaleOp);
    }

    @SuppressWarnings("unchecked")
    private RealParameter tryGetLambdaFromTreeGenerator(Value<TimeTree> treeValue, BEASTContext context) {
        if (treeValue == null) return null;

        Object gen;
        try {
            Method getGen = treeValue.getClass().getMethod("getGenerator");
            gen = getGen.invoke(treeValue);
        } catch (Exception e) {
            return null;
        }
        if (gen == null) return null;

        String[] methodNames = new String[] { "getLambda", "getBirthRate", "getRate", "getSpeciationRate" };

        for (String mName : methodNames) {
            try {
                Method m = gen.getClass().getMethod(mName);
                Object out = m.invoke(gen);
                if (out instanceof Value) {
                    return context.getAsRealParameter((Value<Double>) out);
                }
            } catch (Exception ignored) {}
        }

        return null;
    }

    @Override
    public Class<AutoCorrelatedClock> getGeneratorClass() {
        return AutoCorrelatedClock.class;
    }

    @Override
    public Class<AutoCorrelatedClockModel> getBEASTClass() {
        return AutoCorrelatedClockModel.class;
    }
}
