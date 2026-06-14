#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

OUT_DIR="${1:-/Users/xuyuan/workspace/b3simulation}"
COUNT="${2:-100}"
SEED_BASE="${3:-2026060900}"
CHAIN_LENGTH="${4:-2000}"

mkdir -p "$OUT_DIR"

mvn -q -DskipTests install

generate_group() {
  local label="$1"
  local lphy="$2"
  local seed_offset="$3"
  local group_dir="$OUT_DIR/$label"
  local log_dir="$group_dir/logs"

  mkdir -p "$group_dir" "$log_dir"

  for ((i = 0; i < COUNT; i++)); do
    local out="$group_dir/mixture-${i}.xml"
    local seed=$((SEED_BASE + seed_offset + i))
    local log="$log_dir/mixture-${i}.convert.log"
    echo "[$label $((i + 1))/$COUNT] $out seed=$seed"
    if ! mvn -q -pl mixture-lphybeast-launcher exec:exec \
      -Dlphybeast.args="convert --packagedir ../target/lphybeast-packages -o ${out} -l ${CHAIN_LENGTH} -le 100 -seed=${seed} ${lphy}" \
      >"$log" 2>&1; then
      echo "Failed to generate $out. Last log lines:"
      tail -80 "$log"
      exit 1
    fi
  done
}

generate_group strict "../mixture-lphy/examples/simulation-strict.lphy" 0
generate_group ucln "../mixture-lphy/examples/simulation-ucln.lphy" 100000
generate_group auto "../mixture-lphy/examples/simulation-auto.lphy" 200000

echo "Generated grouped simulation XML files in $OUT_DIR"
