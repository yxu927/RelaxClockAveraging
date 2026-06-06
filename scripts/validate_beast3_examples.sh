#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

SMOKE_CHAIN_LENGTH="${SMOKE_CHAIN_LENGTH:-2000}"
SMOKE_DIR="$ROOT/target/beast3-smoke"
SMOKE_XML="$SMOKE_DIR/mixture-typed-smoke.xml"
SMOKE_RUN_DIR="$SMOKE_DIR/run"

mkdir -p "$SMOKE_DIR" "$SMOKE_RUN_DIR"

echo "== Maven tests =="
mvn -q -pl mixture-beast test
mvn -q test

echo "== Package =="
mvn -q -DskipTests package

echo "== Validate legacy-compatible XML =="
scripts/beast3_validate_xml.sh examples/mixture.xml

echo "== Validate BEAST3 typed XML =="
scripts/beast3_validate_xml.sh examples/mixture-typed.xml

echo "== Create short-chain typed XML smoke copy =="
python3 - examples/mixture-typed.xml "$SMOKE_XML" "$SMOKE_CHAIN_LENGTH" <<'PY'
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

echo "== Smoke run typed XML =="
rm -rf "$SMOKE_RUN_DIR"
mkdir -p "$SMOKE_RUN_DIR"
(
  cd "$SMOKE_RUN_DIR"
  "$ROOT/scripts/beast3_run.sh" -overwrite -seed 1 "$SMOKE_XML"
)

echo "BEAST3 examples validation passed."
