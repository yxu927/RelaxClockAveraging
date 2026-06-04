open module mixture.beast {
    requires beast.pkgmgmt;
    requires beast.base;

    exports mixture.beast.evolution.mixture;
    exports mixture.beast.evolution.operator;
    exports mixture.beast.evolution.util;

    provides beast.base.core.BEASTInterface with
        mixture.beast.evolution.mixture.CategoricalDistribution,
        mixture.beast.evolution.mixture.HierarchicalSVSLogger,
        mixture.beast.evolution.mixture.MixtureLikelihoodLogger,
        mixture.beast.evolution.mixture.MixtureTreeLikelihood,
        mixture.beast.evolution.mixture.RelaxedRatesPriorSVS,
        mixture.beast.evolution.mixture.SharedRatesClockModel,
        mixture.beast.evolution.operator.ACSigma2NonCenteredOperator,
        mixture.beast.evolution.operator.ACSubtreeUIncrementOperator,
        mixture.beast.evolution.operator.AlphaAnnealingOperator,
        mixture.beast.evolution.operator.IndicatorGibbsOperator,
        mixture.beast.evolution.operator.SingleRateScaleOperator,
        mixture.beast.evolution.operator.SubtreeRateScaleOperator,
        mixture.beast.evolution.operator.UCACSwitchBridgeOperator,
        mixture.beast.evolution.operator.UCLDStdevNonCenteredOperator;
}
