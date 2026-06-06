#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

RUNTIME_DIR="$ROOT/target/beast3-runtime"
mkdir -p "$RUNTIME_DIR"

# Build mixture-beast first so target/classes and version.xml exist.
mvn -q -pl mixture-beast -DskipTests package

CP_FILE="$RUNTIME_DIR/mixture-beast-runtime.classpath"

mvn -q -pl mixture-beast dependency:build-classpath \
  -Dmdep.outputFile="$CP_FILE"

DEPS_CP="$(cat "$CP_FILE")"

# target/classes is an exploded JPMS module because it contains module-info.class.
# The whole module path is passed as one quoted argument so paths with spaces are safe.
MODULE_PATH="$ROOT/mixture-beast/target/classes:$DEPS_CP"
BEAST_PACKAGE_PATH_VALUE="$ROOT/mixture-beast/target/classes"

exec java \
  --module-path "$MODULE_PATH" \
  -DBEAST_PACKAGE_PATH="$BEAST_PACKAGE_PATH_VALUE" \
  -m beast.base/beast.base.minimal.BeastMain \
  "$@"
