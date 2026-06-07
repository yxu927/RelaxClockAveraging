# RelaxClockAveraging

A BEAST2/BEAST3 package for Bayesian averaging over strict and relaxed
molecular clock models using a mixture likelihood.

## Overview

RelaxClockAveraging implements a hierarchical mixture of three clock models:
a strict clock, an uncorrelated log-normal relaxed clock (UCLN), and an
autocorrelated relaxed clock (AC). All three models operate on a single
shared branch-rate vector. A binary SVS indicator, sampled at each MCMC
step, selects between the UCLN and AC parameterisations of the shared
rates. A top-level mixture combines the strict and relaxed likelihood
components, and posterior model support can be monitored through the
mixture loggers. Posterior clock-model probabilities are estimated jointly
with the phylogeny, propagating model uncertainty without trans-dimensional
jumps.

## Installation

This package is not currently listed in the BEAST2 package manager.
To install from a local build:

1. Build the package zip (see [Building from source](#building-from-source)).
2. Unpack `mixture-beast/target/mixture-beast-0.1.0.zip` into your BEAST2
   package directory. The zip contains `version.xml`, `lib/`, `examples/`,
   and `fxtemplates/`.

Requires BEAST.base >= 2.8.0.

## Usage

There is currently no BEAUti template. Analyses are configured directly
through XML. The key `spec` classes are:

- `mixture.beast.evolution.mixture.MixtureTreeLikelihood` - log-sum-exp
  mixture over component tree likelihoods
- `mixture.beast.evolution.mixture.RelaxedRatesPriorSVS` - SVS prior on
  the shared rate vector; indicator selects UC or AC
- `mixture.beast.evolution.mixture.SharedRatesClockModel` - maps the
  shared rate vector to per-branch CTMC rates
- `mixture.beast.evolution.mixture.CategoricalDistribution` - categorical
  prior on the UC/AC indicator
- `mixture.beast.evolution.operator.IndicatorGibbsOperator` - Gibbs
  update of the binary UC/AC indicator
- `mixture.beast.evolution.operator.UCACSwitchBridgeOperator`,
  `ACSigma2NonCenteredOperator`, `ACSubtreeUIncrementOperator`,
  `UCLDStdevNonCenteredOperator` - relaxed-clock mixing operators

See `examples/mixture.xml` for a complete annotated XML analysis.

## Example

```bash
beast examples/mixture.xml
```

Two example analyses are provided:

- `examples/mixture.xml` - legacy-compatible BEAST2 XML (31 taxa, 600
  sites, JC+Gamma site model, coalescent tree prior)
- `examples/mixture-typed.xml` - same analysis using BEAST3 typed/spec
  parameter inputs

## Citation

Citation information will be added when available.

## Building from source

Requires JDK 25 and Maven.

```bash
mvn -q test                   # compile and run tests
mvn -q -DskipTests package    # build JARs and package zip
```

## License

MIT - see [`LICENSE.txt`](LICENSE.txt).
