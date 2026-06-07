# RelaxClockAveraging

RelaxClockAveraging is a BEAST3 package for Bayesian model averaging over relaxed molecular clock models.

The package implements a mixture-clock analysis in which branch rates are shared across clock components and an indicator variable selects between relaxed-clock regimes during MCMC. The BEAST module contains the clock model, mixture prior/likelihood components, clock-switching operators, and loggers used by the example analyses.

## Project structure

RelaxClockAveraging contains four subprojects:

### mixture-beast

The BEAST3 package implementation.

- `SharedRatesClockModel` maps a shared branch-rate vector onto the tree.
- `RelaxedRatesPriorSVS` evaluates the strict / uncorrelated / autocorrelated relaxed-clock prior mixture.
- `CategoricalDistribution` provides the indicator prior.
- `MixtureTreeLikelihood` combines tree likelihood components by mixture weights.
- `IndicatorGibbsOperator` updates the clock indicator.
- The remaining clock operators update branch rates, hyperparameters, and UC/AC state switches.
- `MixtureLikelihoodLogger` and `HierarchicalSVSLogger` expose mixture diagnostics in the log file.

### mixture-lphy

LPhy model definitions and examples for the mixture-clock workflow.

### mixture-lphybeast

LPhyBEAST mapping code for converting LPhy objects into BEAST objects.

### mixture-lphy-studio

Studio-facing support for working with the LPhy examples.

## Examples

Two BEAST XML examples are provided:

- `examples/mixture.xml`
- `examples/mixture-typed.xml`

`examples/mixture.xml` is the legacy-style XML example. It is kept for comparison and for users who want an XML that still resembles the original BEAST2-era input style.

`examples/mixture-typed.xml` is the BEAST3 typed/spec version of the same analysis. It replaces legacy parameter and prior wrappers with BEAST3 typed parameters such as `RealScalarParam`, `RealVectorParam`, and `IntScalarParam`.

Both examples use the same final relaxed-clock operator schedule:

- `ACSubtreeUIncrementOperator`
- `UCLDStdevNonCenteredOperator`
- `ACSigma2NonCenteredOperator`
- `UCACSwitchBridgeOperator`

`AlphaAnnealingOperator` is implemented but is not used by these examples because it requires an alpha-coupled mixture likelihood setup.

## Build and test

Run from the repository root:

```bash
mvn -q -pl mixture-beast test
mvn -q test
mvn -q -DskipTests package
```

Validate the two XML examples:

```bash
scripts/beast3_validate_xml.sh examples/mixture.xml
scripts/beast3_validate_xml.sh examples/mixture-typed.xml
```

Run the full short-chain validation:

```bash
SMOKE_CHAIN_LENGTH=2000 scripts/validate_beast3_examples.sh
```

The validation script checks the final operator schedule, validates both XML files, and smoke-runs short-chain legacy-style and typed BEAST3 examples.

## Running BEAST3

The helper runner builds `mixture-beast`, constructs a path-safe module path, and launches the BEAST3 minimal runner with the package on `BEAST_PACKAGE_PATH`.

```bash
scripts/beast3_run.sh -overwrite examples/mixture-typed.xml
```

To run the legacy-style example under the same BEAST3 runtime:

```bash
scripts/beast3_run.sh -overwrite examples/mixture.xml
```

## BEAST2 and BEAST3 status

This repository now targets BEAST3. The legacy-style XML is preserved as a comparison example, but the package build and examples are intended to run with the BEAST3 Maven/module runtime.

The LPhy and LPhyBEAST modules build, but the repository does not yet provide a stable documented `.lphy -> XML` generation command for the typed BEAST3 example. The hand-maintained BEAST XML examples are therefore the current reference analyses.

See `BEAST3_MIGRATION.md` for the migration notes and typed XML details.
