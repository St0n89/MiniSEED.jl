#!/usr/bin/env bash
# Fetch the official libmseed reference test files used by test/validate.jl.
# Source: https://github.com/EarthScope/libmseed (Apache-2.0), test/data/.
set -euo pipefail

BASE="https://raw.githubusercontent.com/EarthScope/libmseed/main/test/data"
DIR="$(cd "$(dirname "$0")" && pwd)"

files=(
  reference-testdata-steim1.mseed2
  reference-testdata-steim1-LE.mseed2
  reference-testdata-steim2.mseed2
  reference-testdata-steim2-LE.mseed2
  reference-testdata-int16.mseed2
  reference-testdata-int32.mseed2
  reference-testdata-float32.mseed2
  reference-testdata-float64.mseed2
  reference-testdata-text.mseed2
  testdata-no-blockette1000-steim1.mseed2
)

for f in "${files[@]}"; do
  echo "fetching $f"
  curl -fsSL "$BASE/$f" -o "$DIR/$f"
done

echo "done: ${#files[@]} files in $DIR"
