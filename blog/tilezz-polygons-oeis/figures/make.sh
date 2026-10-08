#!/usr/bin/env bash
# Regenerate every figure in the post. Run from this directory; SVGs land in the
# parent (next to index.md). Each figure's parameters live either in the python
# script that draws it or in the invocation below -- nowhere else.
set -euo pipefail
cd "$(dirname "$0")"
OUT="${OUT:-..}"
REPO="${REPO:-..}"          # tilezz checkout, for the two Rust-backed families

usage() { echo "usage: $0 [places|automaton|periodic|coronas|all]"; exit 1; }
WHAT="${1:-all}"

# ---------------------------------------------------------------- pure python
# No dependencies beyond the standard library.
if [[ $WHAT == places || $WHAT == all ]]; then
  echo "== places (12-gon, its shadow, the zoom, the hexagon, the seesaw pair)"
  ( cd "$OUT" && python3 "$OLDPWD/fig_places.py" )
fi
if [[ $WHAT == automaton || $WHAT == all ]]; then
  echo "== automaton (walk over the cell grid, single cell)"
  ( cd "$OUT" && python3 "$OLDPWD/fig_automaton.py" )
fi

# ------------------------------------------------------------------ rust side
# fig_tilings.rs runs under rust-script, or as a cargo example. Cargo is used
# here because it reuses the checkout's target/ directory.
build_tilings() {
  # symlink, not a copy: figures/fig_tilings.rs stays the only source
  mkdir -p "$REPO/examples"
  ln -sfn ../figures/fig_tilings.rs "$REPO/examples/fig_tilings.rs"
  ( cd "$REPO" && cargo build --release --example fig_tilings --features cli )
}

if [[ $WHAT == periodic || $WHAT == all ]]; then
  echo "== periodic patches (one per detector)"
  build_tilings
  R="$REPO/target/release/examples/fig_tilings"
  # out                              word                                detector radius cap  colour
  "$R" "$OUT/fig-per-translation.svg" "-5,1,3,1,2,2,-2,-1,5,2,-2,-1,5,2"  bn      5      60   2AA8C0
  "$R" "$OUT/fig-per-conway.svg"      "-5,-3,-3,2,4,1,0,5,-4,3,1,3,3,5"   conway  6      40
  "$R" "$OUT/fig-per-isohedral.svg"   "-3,-2,4,2,3,-1,2,3,4"              iso     5      60
  "$R" "$OUT/fig-per-anisohedral.svg" "-1,4,1,3,1,4"                      aniso   6      40
  "$R" "$OUT/fig-per-torus.svg"       "-5,4,3,-3,5,-1,5,-3,3,4"           torus   9      12   q
  "$R" "$OUT/fig-per-torus18.svg"     "-4,-2,2,4,-2,2,4,-2,-2,4,0,2,2,4"  torus   16     5
fi

if [[ $WHAT == coronas || $WHAT == all ]]; then
  echo "== maximal coronas (shipped CLI, no custom code; minutes per tile)"
  ( cd "$REPO" && cargo build --release --bin classify_tiles --features cli )
  C="$REPO/target/release/classify_tiles"
  corona() { "$C" --render="$1" --render-mode corona --render-out "$OUT/$2"; }
  corona "-4,2,4,-4,2,4,-2,4,-2,2,2,4" fig-heesch3    # the only Heesch-3 tile at perimeter <= 14
  corona "-1,2,3,1,2,1,4"              fig-heesch2a
  corona "-3,2,-3,4,3,0,3,-2,3,5"      fig-heesch2b
fi
echo "done"
