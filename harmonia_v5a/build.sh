#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  Harmonia — Build Script
#  Usage: ./build.sh [--clean] [--run]
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

CLEAN=0
RUN=0
for arg in "$@"; do
  [[ "$arg" == "--clean" ]] && CLEAN=1
  [[ "$arg" == "--run"   ]] && RUN=1
done

echo "═══════════════════════════════════════════════"
echo "  HARMONIA Build"
echo "═══════════════════════════════════════════════"

# ── 1. Install system dependencies ─────────────────────────────────────────
# On Debian/Ubuntu:
# apt install libfltk1.3-dev libgl1-mesa-dev libglu1-mesa-dev libasound2-dev cmake build-essential pkg-config python3-numpy python3-scipy

# ── 2. Install Python theory deps ──────────────────────────────────────────
echo ""
echo "▶ Checking Python dependencies..."
python3 -c "import numpy, json, math, itertools" && echo "  ✓ Python deps OK"

# ── 3. Configure ───────────────────────────────────────────────────────────
echo ""
echo "▶ Configuring with CMake..."
BUILD_DIR="build"
[[ $CLEAN -eq 1 ]] && rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake .. -DCMAKE_BUILD_TYPE=Release

# ── 4. Build ───────────────────────────────────────────────────────────────
echo ""
echo "▶ Building..."
make -j$(nproc)

echo ""
echo "════════════════════════════════════════════════"
echo "  ✓ Build complete: build/harmonia"
echo ""
echo "  Usage:"
echo "    ./build/harmonia"
echo ""
echo "  Shortcuts:"
echo "    Tonal Space: Left-click  node = add/toggle voice"
echo "                 Middle-drag      = rotate"
echo "                 Scroll           = zoom"
echo "    Theory:  'Suggest from current chord' = Markov suggestions"
echo "             'Suggest completion note'    = roughness-optimal note"
echo "             'Analyse current EDO'        = tuning analysis"
echo "════════════════════════════════════════════════"

cd ..

[[ $RUN -eq 1 ]] && exec ./build/harmonia
