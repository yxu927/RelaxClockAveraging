#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 path/to/file.xml" >&2
  exit 2
fi

XML="$1"

if [ ! -f "$XML" ]; then
  echo "XML file not found: $XML" >&2
  exit 2
fi

echo "Validating XML with path-safe BEAST3 runner: $XML"
scripts/beast3_run.sh -validate "$XML"
