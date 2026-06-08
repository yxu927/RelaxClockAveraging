# RelaxClockAveraging

RelaxClockAveraging is a BEAST package for Bayesian averaging over strict,
uncorrelated relaxed, and autocorrelated relaxed molecular-clock models.

The package uses a shared branch-rate vector and mixture likelihoods so that
clock-model uncertainty is sampled jointly with the tree and substitution
parameters. The current codebase supports direct BEAST XML use, BEAST3
typed/spec XML, a BEAUti clock-model template, and LPhy/LPhyBEAST XML
generation for the mixture relaxed-clock path.

## Repository Layout

- `mixture-beast`: BEAST package classes, operators, loggers, examples, and
  BEAUti template packaging.
- `mixture-lphy`: LPhy generators/functions for the mixture clock model.
- `mixture-lphybeast`: LPhyBEAST mappings from LPhy objects to BEAST XML.
- `mixture-lphybeast-launcher`: Maven launcher for real `.lphy -> XML`
  conversion.
- `mixture-lphy-studio`: LPhy Studio module.
- `examples`: hand-maintained BEAST XML examples.
- `scripts`: local BEAST3 validation and smoke-run helpers.

## Model Components

The main BEAST classes are:

- `MixtureTreeLikelihood`: log-sum-exp mixture over component tree likelihoods.
- `RelaxedRatesPriorSVS`: SVS prior over the shared branch-rate vector.
- `SharedRatesClockModel`: legacy-compatible shared branch-rate clock model.
- `SharedRatesClockModelSpec`: BEAST3 spec branch-rate model used by
  generated typed XML.
- `CategoricalDistribution`: categorical prior for clock-model indicators.
- `IndicatorGibbsOperator`: Gibbs update for the UC/AC indicator.
- `SingleRateScaleOperator` and `SubtreeRateScaleOperator`: shared-rate moves.
- `ACSubtreeUIncrementOperator`, `UCLDStdevNonCenteredOperator`,
  `ACSigma2NonCenteredOperator`, and `UCACSwitchBridgeOperator`: final
  relaxed-clock mixing operators.

`AlphaAnnealingOperator` is implemented, but it is not part of the shipped
examples. It needs an alpha-coupled mixture likelihood setup and should not be
added to the standard examples by default.

## Requirements

- JDK 25 for the BEAST3 typed modules.
- Maven.
- BEAST.base 2.8.0 beta5 or newer compatible BEAST3 artifacts.
- LPhy 1.8.0-beta1 and LPhyBEAST 2.0.0-SNAPSHOT for the LPhyBEAST conversion
  module.
- Optional: BEAGLE. If BEAGLE is not installed, BEAST falls back to the Java
  likelihood core.

If Maven cannot resolve `lphybeast:2.0.0-SNAPSHOT`, install the matching
LPhyBEAST 2.0 snapshot into your local Maven repository first.

## Build

From the repository root:

```bash
mvn -q test
mvn -q -DskipTests package
```

The BEAST package zip is produced at:

```text
mixture-beast/target/mixture-beast-0.1.0.zip
```

## Install as a BEAST Package

Build the package first:

```bash
mvn -q -DskipTests package
```

On macOS, install into the BEAST 2.8 package directory:

```bash
PKG_DIR="$HOME/Library/Application Support/BEAST/2.8/mixture"
rm -rf "$PKG_DIR"
mkdir -p "$PKG_DIR"
ditto -x -k mixture-beast/target/mixture-beast-0.1.0.zip "$PKG_DIR"
```

The installed package contains:

- `version.xml`
- `lib/mixture-beast-0.1.0.jar`
- `examples/`
- `fxtemplates/SVSRelaxedClockTemplate.xml`

## Run BEAST XML from the Checkout

The repository includes a path-safe BEAST3 runner that builds the local
`mixture-beast` module and runs BEAST with the local package classes.

Validate an XML without running MCMC:

```bash
scripts/beast3_validate_xml.sh examples/mixture.xml
scripts/beast3_validate_xml.sh examples/mixture-typed.xml
```

Run an XML:

```bash
scripts/beast3_run.sh -overwrite examples/mixture.xml
scripts/beast3_run.sh -overwrite examples/mixture-typed.xml
```

Run any XML by path:

```bash
scripts/beast3_run.sh -overwrite /path/to/analysis.xml
```

If you already installed the package into BEAST and have a BEAST executable on
your `PATH`, you can also run:

```bash
beast examples/mixture.xml
```

The local scripts are recommended for development because they always use the
classes from the current checkout.

## Shipped BEAST XML Examples

Two complete XML examples are maintained:

- `examples/mixture.xml`: legacy-compatible baseline XML.
- `examples/mixture-typed.xml`: BEAST3 typed/spec XML for the same analysis.

Both examples include the same final relaxed-clock operator schedule:

- `ACSubtreeUIncrementOperator`
- `UCLDStdevNonCenteredOperator`
- `ACSigma2NonCenteredOperator`
- `UCACSwitchBridgeOperator`

Both examples intentionally exclude `AlphaAnnealingOperator`.

Run the full local example check:

```bash
SMOKE_CHAIN_LENGTH=2000 scripts/validate_beast3_examples.sh
```

This script checks the operator schedule, runs Maven tests, validates both XML
files, creates short-chain copies under `target/beast3-smoke/`, and smoke-runs
both examples.

## BEAUti

The package ships a BEAUti clock-model template:

```text
fxtemplates/SVSRelaxedClockTemplate.xml
```

After installing the package zip into the BEAST 2.8 package directory, open
BEAUti 3 / BEAST.base 2.8+ and select:

```text
Clock Model -> SVSRelaxed Clock Mixture
```

The BEAUti template is legacy-compatible. It generates XML using ordinary BEAST
parameter objects, while the underlying Java classes also support BEAST3 typed
inputs. The template includes the final UC/AC clock-mixing operators and hides
structural inputs that users should not edit in BEAUti.

Recommended BEAUti workflow:

1. Open BEAUti 3.
2. Import an alignment.
3. Choose the site model and tree prior.
4. In the Clock Model tab, select `SVSRelaxed Clock Mixture`.
5. Set the visible clock parameters and priors.
6. Save the XML.
7. Validate the saved XML:

```bash
scripts/beast3_validate_xml.sh /path/to/saved.xml
```

## LPhy to BEAST XML

The current LPhyBEAST path can generate a BEAST3-valid XML from:

```text
mixture-lphy/examples/new.lphy
```

Generate XML:

```bash
mvn -q -pl mixture-lphybeast-launcher exec:exec \
  -Dlphybeast.args="convert --packagedir ../target/lphybeast-packages -o ../target/lphybeast-new.xml -l 2000 -le 100 ../mixture-lphy/examples/new.lphy"
```

The generated XML is written to:

```text
mixture-lphy/target/lphybeast-new.xml
```

The output path is relative to the `.lphy` file directory because LPhyBEAST
sets its working directory to the script location during conversion.

Validate and run the generated XML:

```bash
xmllint --noout mixture-lphy/target/lphybeast-new.xml
scripts/beast3_validate_xml.sh mixture-lphy/target/lphybeast-new.xml
scripts/beast3_run.sh -overwrite mixture-lphy/target/lphybeast-new.xml
```

The generated XML uses the BEAST3 typed/spec path for the custom mixture
components. It includes the final SVS clock operators and keeps simulation-only
top-level allocation variables out of the MCMC state.

## Regression Commands

Use these before opening or updating a PR:

```bash
mvn -q test
mvn -q -DskipTests package
SMOKE_CHAIN_LENGTH=2000 scripts/validate_beast3_examples.sh

mvn -q -pl mixture-lphybeast-launcher exec:exec \
  -Dlphybeast.args="convert --packagedir ../target/lphybeast-packages -o ../target/lphybeast-new.xml -l 2000 -le 100 ../mixture-lphy/examples/new.lphy"

xmllint --noout mixture-lphy/target/lphybeast-new.xml
scripts/beast3_validate_xml.sh mixture-lphy/target/lphybeast-new.xml
scripts/beast3_run.sh -overwrite mixture-lphy/target/lphybeast-new.xml
```

Warnings about missing BEAGLE JNI or BEAST operator tuning suggestions are not
test failures. BEAST will run with the Java likelihood core when BEAGLE is not
available.

## Citation

Citation information will be added when available.

## License

MIT - see `LICENSE.txt`.
