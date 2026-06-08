# RelaxClockAveraging

RelaxClockAveraging is a BEAST package for Bayesian averaging over strict, uncorrelated relaxed, and autocorrelated relaxed molecular-clock models. The model uses a shared branch-rate vector and mixture likelihoods so that clock-model uncertainty is sampled jointly with tree and substitution parameters.

The package supports direct BEAST XML use, BEAST3 typed/spec XML, a BEAUti clock-model template, and LPhy/LPhyBEAST XML generation.

## Repository Layout

- `mixture-beast`: BEAST classes, operators, loggers, examples, and BEAUti template packaging.
- `mixture-lphy`: LPhy generators and functions for the mixture clock model.
- `mixture-lphybeast`: LPhyBEAST mappings from LPhy objects to BEAST XML.
- `mixture-lphybeast-launcher`: Maven launcher for `.lphy -> XML` conversion.
- `mixture-lphy-studio`: optional LPhy Studio integration.
- `examples`: maintained BEAST XML examples.
- `scripts`: local BEAST3 validation and smoke-run helpers.

## Requirements

- JDK 25.
- Maven.
- BEAST.base 2.8.0 beta5 or a compatible BEAST3 build.
- LPhy 1.8.0-beta1 and LPhyBEAST 2.0.0-SNAPSHOT for LPhyBEAST conversion.
- Optional: BEAGLE. If BEAGLE is unavailable, BEAST runs with the Java likelihood core.

If Maven cannot resolve `lphybeast:2.0.0-SNAPSHOT`, install the matching LPhyBEAST 2.0 snapshot into the local Maven repository first.

## Build

Build and test from the repository root:

```bash
mvn -q test
mvn -q -DskipTests package
```

The BEAST package zip is written to:

```text
mixture-beast/target/mixture-beast-0.1.0.zip
```

## Install as a BEAST Package

On macOS, install the built package into the BEAST 2.8 package directory:

```bash
PKG_DIR="$HOME/Library/Application Support/BEAST/2.8/mixture"
rm -rf "$PKG_DIR"
mkdir -p "$PKG_DIR"
ditto -x -k mixture-beast/target/mixture-beast-0.1.0.zip "$PKG_DIR"
```

The installed package contains `version.xml`, `lib/mixture-beast-0.1.0.jar`, `examples/`, and `fxtemplates/SVSRelaxedClockTemplate.xml`.

## BEAST XML Examples and Validation

Two complete XML examples are maintained:

- `examples/mixture.xml`: legacy-compatible baseline example.
- `examples/mixture-typed.xml`: BEAST3 typed/spec example.

Both examples use the same relaxed-clock operator schedule: `ACSubtreeUIncrementOperator`, `UCLDStdevNonCenteredOperator`, `ACSigma2NonCenteredOperator`, and `UCACSwitchBridgeOperator`. `AlphaAnnealingOperator` is implemented but intentionally excluded from the shipped examples because it requires an alpha-coupled mixture-likelihood setup.

Validate the examples:

```bash
scripts/beast3_validate_xml.sh examples/mixture.xml
scripts/beast3_validate_xml.sh examples/mixture-typed.xml
```

Run the examples:

```bash
scripts/beast3_run.sh -overwrite examples/mixture.xml
scripts/beast3_run.sh -overwrite examples/mixture-typed.xml
```

Run any XML by path:

```bash
scripts/beast3_run.sh -overwrite /path/to/analysis.xml
```

Run the maintained short-chain example check:

```bash
SMOKE_CHAIN_LENGTH=2000 scripts/validate_beast3_examples.sh
```

The local scripts build the package from the checkout before invoking BEAST.

## BEAUti Template

The BEAUti clock-model template is:

```text
fxtemplates/SVSRelaxedClockTemplate.xml
```

After installing the package zip, open BEAUti 3 / BEAST.base 2.8+ and select:

```text
Clock Model -> SVSRelaxed Clock Mixture
```

The template remains legacy-compatible: it writes standard BEAST parameter objects, while the underlying Java classes also support BEAST3 typed inputs. It uses the same UC/AC clock-mixing operators as the maintained examples.

Validate XML saved from BEAUti:

```bash
scripts/beast3_validate_xml.sh /path/to/saved.xml
```

## LPhy to BEAST XML

The LPhyBEAST path generates BEAST3-valid XML from:

```text
mixture-lphy/examples/new.lphy
```

Generate XML from LPhy:

```bash
mvn -q -pl mixture-lphybeast-launcher exec:exec \
  -Dlphybeast.args="convert --packagedir ../target/lphybeast-packages -o ../target/lphybeast-new.xml -l 2000 -le 100 ../mixture-lphy/examples/new.lphy"
```

The generated XML is written to:

```text
mixture-lphy/target/lphybeast-new.xml
```

Validate and run the generated XML:

```bash
xmllint --noout mixture-lphy/target/lphybeast-new.xml
scripts/beast3_validate_xml.sh mixture-lphy/target/lphybeast-new.xml
scripts/beast3_run.sh -overwrite mixture-lphy/target/lphybeast-new.xml
```

The generated XML uses the BEAST3 typed/spec path for the custom mixture components, includes the same SVS relaxed-clock operators, and keeps simulation-only allocation variables out of the MCMC state.

## Regression Checks

Before updating the package or opening a pull request, run the build commands, the short-chain example check, and the LPhyBEAST conversion/validation commands above. Warnings about missing BEAGLE JNI and BEAST operator tuning suggestions are not test failures.

## Citation

Citation information will be added when available.

## License

MIT - see `LICENSE.txt`.
