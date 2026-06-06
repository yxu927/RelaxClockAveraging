# BEAST3 typed example

This branch provides two BEAST XML examples:

| File | Purpose |
| --- | --- |
| `examples/mixture.xml` | Legacy-compatible baseline XML |
| `examples/mixture-typed.xml` | BEAST3 typed/spec XML example |

Use `examples/mixture-typed.xml` when testing the BEAST3 typed/spec migration. Keep `examples/mixture.xml` as a compatibility baseline.

## Validation

Run:

```bash
mvn -q -pl mixture-beast test
mvn -q test
mvn -q -DskipTests package

scripts/beast3_validate_xml.sh examples/mixture.xml
scripts/beast3_validate_xml.sh examples/mixture-typed.xml

SMOKE_CHAIN_LENGTH=2000 scripts/validate_beast3_examples.sh
```

The typed example keeps the full chain length from the original example. The validation script creates a short-chain temporary copy under `target/beast3-smoke/` for smoke testing.

## Typed XML status

`examples/mixture-typed.xml` removes legacy XML parameter declarations and legacy prior wrappers from the example path:

| XML | legacy parameter specs | typed parameter specs | legacy `distribution.Prior` specs | typed LogNormal specs |
| --- | ---: | ---: | ---: | ---: |
| `examples/mixture.xml` | 17 | 3 | 4 | 1 |
| `examples/mixture-typed.xml` | 0 | 20 | 0 | 5 |

The typed example uses BEAST3 typed/spec parameters such as:

- `RealScalarParam`
- `RealVectorParam`
- `IntScalarParam`

The custom model, prior, likelihood, operator, and logger classes used by the example support both the legacy-compatible XML path and the typed/spec XML path.

## LPhy / LPhyBEAST status

The LPhy-side modules build and package successfully. A stable, documented `.lphy -> BEAST XML` generation command was not found in this repository, so automatic LPhyBEAST generation of `examples/mixture-typed.xml` is deferred.

Future LPhyBEAST work should add a documented generation command/test before changing mapper output.
