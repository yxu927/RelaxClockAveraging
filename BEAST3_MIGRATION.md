# BEAST3 migration status

## Summary

This branch now ships two BEAST XML examples:

- `examples/mixture.xml`
- `examples/mixture-typed.xml`

`examples/mixture.xml` is the legacy-compatible baseline. It keeps the original XML style while running under the BEAST3 migration branch.

`examples/mixture-typed.xml` is the BEAST3 typed/spec version of the same example.

Both examples validate and smoke-run under the local BEAST3 runner.

## Which example should I use?

| File | Purpose | Status |
| --- | --- | --- |
| `examples/mixture.xml` | Legacy-compatible baseline | Validate + smoke tested |
| `examples/mixture-typed.xml` | BEAST3 typed/spec example | Validate + smoke tested |

Use `examples/mixture-typed.xml` for BEAST3 typed/spec testing. Keep `examples/mixture.xml` for compatibility and comparison.

## Validation commands

```bash
mvn -q -pl mixture-beast test
mvn -q test
mvn -q -DskipTests package

scripts/beast3_validate_xml.sh examples/mixture.xml
scripts/beast3_validate_xml.sh examples/mixture-typed.xml
SMOKE_CHAIN_LENGTH=2000 scripts/validate_beast3_examples.sh
```

`scripts/validate_beast3_examples.sh` checks the final clock-mixing operator schedule, validates both XML files, and smoke-runs short-chain legacy and typed examples. The committed XML files keep their full `chainLength`.

## Typed XML changes

The typed example removes legacy `parameter.RealParameter` and `parameter.IntegerParameter` declarations from the XML. It also removes legacy `distribution.Prior` wrappers and legacy `distribution.LogNormalDistributionModel` declarations.

Typed parameters use BEAST3 `RealScalarParam`, `RealVectorParam`, and `IntScalarParam` specs. Typed prior wrappers use `beast.base.spec.inference.distribution.LogNormal`.

`FunctionOfTensor` wrappers are used where legacy BEAST components still require `Function` inputs.

| XML | legacy parameter specs | typed parameter specs | legacy Prior specs | typed LogNormal specs | custom mixture specs |
| --- | ---: | ---: | ---: | ---: | ---: |
| `examples/mixture.xml` | 17 | 3 | 4 | 1 | 13 |
| `examples/mixture-typed.xml` | 0 | 20 | 0 | 5 | 13 |

## Final Clock-Mixing Operators

Both examples use the same final relaxed-clock operator schedule. The legacy example uses legacy input names; the typed example uses mutable typed input names.

- `ACSubtreeUIncrementOperator`
- `UCLDStdevNonCenteredOperator`
- `ACSigma2NonCenteredOperator`
- `UCACSwitchBridgeOperator`

`AlphaAnnealingOperator` is typed-capable in Java but is intentionally not included in these examples. It requires an alpha-coupled `MixtureTreeLikelihood` setup, and adding it here would change the model meaning.

## Custom classes with typed support

These custom classes were migrated additively:

- `SharedRatesClockModel`
- `RelaxedRatesPriorSVS`
- `IndicatorGibbsOperator`
- `SingleRateScaleOperator`
- `SubtreeRateScaleOperator`
- `ACSubtreeUIncrementOperator`
- `UCLDStdevNonCenteredOperator`
- `ACSigma2NonCenteredOperator`
- `UCACSwitchBridgeOperator`
- `AlphaAnnealingOperator`
- `CategoricalDistribution`
- `MixtureTreeLikelihood`
- `MixtureLikelihoodLogger`
- `HierarchicalSVSLogger`

Legacy input names remain supported. Typed input names are used in `examples/mixture-typed.xml`.

## LPhy / LPhyBEAST status

The LPhy-side modules test and package successfully:

- `mixture-lphy`
- `mixture-lphybeast`
- `mixture-lphy-studio`

No safe, documented `.lphy -> XML` generation command was found in this repository. LPhyBEAST typed XML generation is therefore deferred until a stable command/test exists.

Likely future mapper files:

- `SVSRawBranchRatesToBEAST.java`
- `SharedRatesClockToBEAST.java`
- `MixturePhyloCTMCToBEAST.java`
- `CategoricalToBEAST.java`

## Remaining work

- Optionally add a stable `.lphy -> XML` LPhyBEAST command/test.
- Optionally update LPhyBEAST mappers to generate typed XML.
- Optional later work: add an alpha-coupled example if `AlphaAnnealingOperator` is needed for a separate annealed-likelihood workflow.

## Safety notes

`examples/mixture.xml` remains the legacy-compatible example. `examples/mixture-typed.xml` remains the separately validated BEAST3 typed/spec example.

Do not replace the legacy example until users and the LPhyBEAST generation path are ready.

This review branch intentionally excludes migration logs, intermediate XML candidates, and baseline snapshots.
