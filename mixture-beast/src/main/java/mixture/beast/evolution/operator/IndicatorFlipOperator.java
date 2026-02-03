//package mixture.beast.evolution.operator;
//
//import beast.base.core.Description;
//import beast.base.core.Input;
//import beast.base.inference.Operator;
//import beast.base.inference.parameter.IntegerParameter;
//
//@Description("Flip a binary indicator parameter (0 <-> 1).")
//public class IndicatorFlipOperator extends Operator {
//
//    public final Input<IntegerParameter> indicatorInput = new Input<>(
//            "indicator",
//            "binary indicator parameter (dimension=1) to flip",
//            Input.Validate.REQUIRED
//    );
//
//    private IntegerParameter indicator;
//
//    @Override
//    public void initAndValidate() {
//        indicator = indicatorInput.get();
//        if (indicator.getDimension() != 1) {
//            throw new IllegalArgumentException("IndicatorFlipOperator: indicator must have dimension=1.");
//        }
//    }
//
//    @Override
//    public double proposal() {
//        final int v = indicator.getValue(0);
//        if (v == 0) {
//            indicator.setValue(0, 1);
//        } else if (v == 1) {
//            indicator.setValue(0, 0);
//        } else {
//            // If indicator is not binary, force rejection
//            return Double.NEGATIVE_INFINITY;
//        }
//        return 0.0; // symmetric proposal
//    }
//}
