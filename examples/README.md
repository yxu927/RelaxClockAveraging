# Examples

- `mixture.xml`: legacy-compatible baseline.
- `mixture-typed.xml`: BEAST3 typed/spec version of the same example.

Validate both examples:

```bash
scripts/beast3_validate_xml.sh examples/mixture.xml
scripts/beast3_validate_xml.sh examples/mixture-typed.xml
SMOKE_CHAIN_LENGTH=2000 scripts/validate_beast3_examples.sh
```
