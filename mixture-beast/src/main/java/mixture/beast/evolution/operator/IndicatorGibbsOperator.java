package mixture.beast.evolution.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.util.Randomizer;
import mixture.beast.evolution.mixture.RelaxedRatesPriorSVS;

@Description("Gibbs-style update for a binary indicator (0/1) used inside RelaxedRatesPriorSVS. "
        + "Proposes indicator from its conditional distribution given the shared rates vector "
        + "and returns a Hastings term that makes acceptance probability 1.")
public class IndicatorGibbsOperator extends Operator {

    public final Input<IntegerParameter> indicatorInput = new Input<>(
            "indicator",
            "binary indicator parameter (dimension=1)",
            Input.Validate.OPTIONAL
    );

    public final Input<IntScalarParam<?>> indicatorScalarInput = new Input<>(
            "indicatorScalar",
            "BEAST3 typed mutable integer scalar; 0=uncorrelated, 1=autocorrelated.",
            Input.Validate.OPTIONAL
    );

    public final Input<RelaxedRatesPriorSVS> priorInput = new Input<>(
            "prior",
            "the RelaxedRatesPriorSVS object providing UC/AC normalized log densities",
            Input.Validate.REQUIRED
    );

    public final Input<Double> pOneInput = new Input<>(
            "pOne",
            "prior probability P(indicator=1). Use 0.5 if you want equal prior model weights.",
            0.5
    );

    private IntegerParameter legacyIndicator;
    private IntScalarParam<?> typedIndicator;
    private RelaxedRatesPriorSVS prior;
    private double pOne;

    @Override
    public void initAndValidate() {
        legacyIndicator = indicatorInput.get();
        typedIndicator = indicatorScalarInput.get();
        prior = priorInput.get();
        pOne = pOneInput.get();

        if (legacyIndicator == null && typedIndicator == null) {
            throw new IllegalArgumentException("IndicatorGibbsOperator: either indicator or indicatorScalar must be specified.");
        }
        if (legacyIndicator != null && typedIndicator != null) {
            throw new IllegalArgumentException("IndicatorGibbsOperator: specify only one of indicator or indicatorScalar.");
        }
        if (legacyIndicator != null && legacyIndicator.getDimension() != 1) {
            throw new IllegalArgumentException("IndicatorGibbsOperator: indicator must have dimension=1.");
        }
        if (!(pOne > 0.0 && pOne < 1.0)) {
            throw new IllegalArgumentException("IndicatorGibbsOperator: pOne must be in (0,1).");
        }
    }

    private static double logSumExp(double a, double b) {
        if (a > b) {
            return a + Math.log1p(Math.exp(b - a));
        }
        return b + Math.log1p(Math.exp(a - b));
    }

    private int indicatorValue() {
        if (legacyIndicator != null) {
            return legacyIndicator.getValue(0);
        }
        return typedIndicator.get();
    }

    private void setIndicatorValue(final int value) {
        if (legacyIndicator != null) {
            legacyIndicator.setValue(0, value);
        } else {
            typedIndicator.set(value);
        }
    }

    @Override
    public double proposal() {
        final int oldK = indicatorValue();
        if (!(oldK == 0 || oldK == 1)) {
            return Double.NEGATIVE_INFINITY;
        }

        // Compute (unnormalized) log weights for k=0 and k=1
        final double logW0 = Math.log(1.0 - pOne) + prior.logPriorUCOnly();
        final double logW1 = Math.log(pOne) + prior.logPriorACOnly();

        if (Double.isInfinite(logW0) && Double.isInfinite(logW1)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double logDen = logSumExp(logW0, logW1);
        final double logP0 = logW0 - logDen;
        final double logP1 = logW1 - logDen;

        // Sample new indicator
        final double u = Randomizer.nextDouble();
        final int newK = (u < Math.exp(logP1)) ? 1 : 0;
        setIndicatorValue(newK);

        // Hastings ratio for Gibbs:
        // logHR = log q(old|new) - log q(new|old) = logP(old) - logP(new)
        final double logPold = (oldK == 0) ? logP0 : logP1;
        final double logPnew = (newK == 0) ? logP0 : logP1;

        return logPold - logPnew;
    }
}
