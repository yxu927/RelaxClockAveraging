# Examples

- `mixture.xml`: legacy-compatible baseline.
- `mixture-typed.xml`: BEAST3 typed/spec version of the same example.

Both examples include the same final relaxed-clock mixing operators:

- `ACSubtreeUIncrementOperator`
- `UCLDStdevNonCenteredOperator`
- `ACSigma2NonCenteredOperator`
- `UCACSwitchBridgeOperator`

`AlphaAnnealingOperator` is not part of these examples because it requires an alpha-coupled mixture likelihood setup.

Validate:

```bash
scripts/beast3_validate_xml.sh mixture-beast/examples/mixture.xml
scripts/beast3_validate_xml.sh mixture-beast/examples/mixture-typed.xml
SMOKE_CHAIN_LENGTH=2000 scripts/validate_beast3_examples.sh
```
