#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

SMOKE_CHAIN_LENGTH="${SMOKE_CHAIN_LENGTH:-2000}"
SMOKE_DIR="$ROOT/target/beast3-smoke"
LEGACY_SMOKE_XML="$SMOKE_DIR/mixture-smoke.xml"
TYPED_SMOKE_XML="$SMOKE_DIR/mixture-typed-smoke.xml"
LEGACY_SMOKE_RUN_DIR="$SMOKE_DIR/legacy-run"
TYPED_SMOKE_RUN_DIR="$SMOKE_DIR/typed-run"

mkdir -p "$SMOKE_DIR" "$LEGACY_SMOKE_RUN_DIR" "$TYPED_SMOKE_RUN_DIR"

echo "== Check final operator schedule =="
python3 - <<'PY'
from pathlib import Path
import xml.etree.ElementTree as ET

targets = {
    "legacy": Path("examples/mixture.xml"),
    "typed": Path("examples/mixture-typed.xml"),
}
required = {
    "ACSubtreeUIncrementOperator",
    "UCLDStdevNonCenteredOperator",
    "ACSigma2NonCenteredOperator",
    "UCACSwitchBridgeOperator",
}

for label, path in targets.items():
    root = ET.parse(path).getroot()
    specs = [element.get("spec", "") for element in root.iter("operator")]
    print(f"{label}: {path}")
    for name in sorted(required):
        count = sum(1 for spec in specs if spec.endswith(name))
        print(f"  {name}: {count}")
        if count != 1:
            raise SystemExit(f"{path}: expected exactly one {name}, found {count}")
    alpha_count = sum(1 for spec in specs if spec.endswith("AlphaAnnealingOperator"))
    print(f"  AlphaAnnealingOperator: {alpha_count}")
    if alpha_count != 0:
        raise SystemExit(f"{path}: AlphaAnnealingOperator must not be in this example")
PY

echo "== Maven tests =="
mvn -q -pl mixture-beast -Dtest=RemainingOperatorsCharacterizationTest,RemainingOperatorsTypedInputTest test
mvn -q -pl mixture-beast test
mvn -q test

echo "== Package =="
mvn -q -DskipTests package

echo "== Validate legacy-compatible XML =="
scripts/beast3_validate_xml.sh examples/mixture.xml

echo "== Validate BEAST3 typed XML =="
scripts/beast3_validate_xml.sh examples/mixture-typed.xml

make_smoke_xml() {
  local src="$1"
  local dst="$2"

  python3 - "$src" "$dst" "$SMOKE_CHAIN_LENGTH" <<'PY'
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
chain = sys.argv[3]

parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
tree = ET.parse(src, parser=parser)
root = tree.getroot()

runs = [
    e for e in root.iter()
    if isinstance(e.tag, str) and e.tag == "run" and e.get("id") == "MCMC"
]
if len(runs) != 1:
    raise SystemExit(f"Expected exactly one MCMC run, found {len(runs)}")

runs[0].set("chainLength", chain)
dst.parent.mkdir(parents=True, exist_ok=True)
ET.indent(tree, space="    ")
tree.write(dst, encoding="utf-8", xml_declaration=True, short_empty_elements=True)
print(dst)
PY
}

echo "== Create short-chain smoke XML copies =="
make_smoke_xml examples/mixture.xml "$LEGACY_SMOKE_XML"
make_smoke_xml examples/mixture-typed.xml "$TYPED_SMOKE_XML"

echo "== Smoke run legacy-compatible XML =="
rm -rf "$LEGACY_SMOKE_RUN_DIR"
mkdir -p "$LEGACY_SMOKE_RUN_DIR"
(
  cd "$LEGACY_SMOKE_RUN_DIR"
  "$ROOT/scripts/beast3_run.sh" -overwrite -seed 1 "$LEGACY_SMOKE_XML"
)

echo "== Smoke run BEAST3 typed XML =="
rm -rf "$TYPED_SMOKE_RUN_DIR"
mkdir -p "$TYPED_SMOKE_RUN_DIR"
(
  cd "$TYPED_SMOKE_RUN_DIR"
  "$ROOT/scripts/beast3_run.sh" -overwrite -seed 1 "$TYPED_SMOKE_XML"
)

echo "BEAST3 examples validation passed."
