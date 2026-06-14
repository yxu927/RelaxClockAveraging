#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

OUT_DIR="${1:-/Users/xuyuan/workspace/b3simulation}"
COUNT="${2:-100}"
SEED_BASE="${3:-2026060900}"
LPHY_IN="../mixture-lphy/examples/simulation.lphy"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

mvn -q -DskipTests install

for ((i = 0; i < COUNT; i++)); do
  out="$OUT_DIR/mixture-${i}.xml"
  seed=$((SEED_BASE + i))
  log="$LOG_DIR/mixture-${i}.convert.log"
  echo "[$((i + 1))/$COUNT] $out seed=$seed"
  if ! mvn -q -pl mixture-lphybeast-launcher exec:exec \
    -Dlphybeast.args="convert --packagedir ../target/lphybeast-packages -o ${out} -l 2000 -le 100 -seed=${seed} ${LPHY_IN}" \
    >"$log" 2>&1; then
    echo "Failed to generate $out. Last log lines:"
    tail -80 "$log"
    exit 1
  fi
done

echo "Generated $COUNT XML files in $OUT_DIR"
