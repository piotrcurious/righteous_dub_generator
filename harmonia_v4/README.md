# Harmonia — Psychoacoustic Multi-Voice Music Coder (V3 Psychoacoustic Edition)

A Linux FLTK/OpenGL application that demystifies tuning, harmony and voice-leading through the unified lens of **algebraic geometry**, **psychoacoustics** and **neural auditory cortex models**.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  FLTK/OpenGL Application (C++)                                  │
│                                                                 │
│  ┌──────────────┐  ┌─────────────────────┐  ┌───────────────┐ │
│  │ Voice Panel  │  │  Tonnetz Widget (GL) │  │ Theory Panel  │ │
│  │  • ADSR      │  │  • 5-limit lattice   │  │  • Next chord │ │
│  │  • Timbre    │  │  • Triad triangles   │  │  • Completion │ │
│  │  • Detune    │  │  • Abstract objects  │  │  • EDO errors │ │
│  │  • Roughness │  │  • Progression path  │  │  • Neural Psy  │ │
│  └──────┬───────┘  └─────────┬───────────┘  └───────┬───────┘ │
│         │                    │                       │          │
│  ┌──────▼────────────────────▼───────────────────────▼───────┐ │
│  │  Audio Engine (ALSA + additive synthesis thread)           │ │
│  │  • Plomp-Levelt roughness  • Virtual pitch autocorrelation │ │
│  │  • ADSR envelopes          • Chord identification          │ │
│  └─────────────────────────────────────────────────────────── ┘ │
└─────────────────────────────────────────────────────────────────┘
                        │  JSON / stdin-stdout
┌───────────────────────▼─────────────────────────────────────────┐
│  Python Theory Server (subprocess - server3.py)                 │
│  • Neural Auditory Cortex Model (Peripheral/Brainstem/Cortical) │
│  • 5/7-limit JI lattice navigation                              │
│  • Tymoczko orbifold voice-leading distances                    │
│  • Terhardt Subharmonic Summation (Virtual Pitch)               │
│  • Krumhansl-Kessler Tonal Hierarchy Analysis                   │
└─────────────────────────────────────────────────────────────────┘
```

---

## Prerequisites

```bash
# Ubuntu/Debian
sudo apt install libfltk1.3-dev libgl1-mesa-dev libglu1-mesa-dev \
                 libasound2-dev cmake build-essential pkg-config \
                 python3-numpy python3-scipy
```

---

## Build

```bash
chmod +x build.sh
./build.sh          # install deps, configure, compile
./build.sh --run    # build and launch immediately
```

---

## Concepts Implemented

### 1. Algebraic Geometry

**Pitch class space as 𝕊¹ (circle quotient)**
- Frequencies are mapped to pitch classes via `PC = ℝ⁺ / (f ~ 2f) ≅ ℝ/ℤ`
- 12-TET = ℤ₁₂, microtonal EDOs = ℤₙ

**Tonnetz (5-limit lattice)**
- Each Tonnetz node at `(a, b)` represents pitch class `(7a + 4b) mod 12`
- Horizontal axis = perfect fifths (×3/2)
- Diagonal axis = major thirds (×5/4)
- Triangles = major/minor triads (geometric encoding of chord structure)

**Orbifold chord space (Tymoczko 2006)**
- n-voice chords live in `Tⁿ/Sₙ` (torus quotient by permutation group)
- Voice-leading distance = geodesic length in this orbifold
- The Python server computes optimal voice-leading via Hungarian matching

**Maximally even scales**
- The diatonic scale is the unique 7-note subset of ℤ₁₂ that minimizes
  the evenness functional `E(S) = Σᵢⱼ |sᵢ−sⱼ − 12(i−j)/7|`

### 2. Psychoacoustics

**Sethares roughness model**
For each pair of harmonic partials `(fₐ, aₐ)` and `(f_b, a_b)`:
```
R_pair = aₐ · a_b · (Δf/CBW) · exp(1 − Δf/CBW)   for Δf ≤ CBW
```
where `CBW(f) = 24.7 · (4.37f/1000 + 1)` (Moore-Glasberg ERB).

This is computed in real time across all voice pairs and displayed as:
- Per-voice roughness (colour-coded in voice strip)
- Per-pair roughness (in Tonnetz as edge weights)
- Total aggregate roughness (spectrum bar)

**Virtual pitch via autocorrelation**
The combined spectrum (sum of all active harmonic series) is processed
on a log-frequency grid. Autocorrelation peaks identify the virtual
fundamental — the "missing bass note" implied by the chord gestalt.

**Abstract objects**
A chord is an "abstract object" when its virtual pitch and pitch-class
set are recognisable across voices even if no single voice carries the
complete information. Displayed as a pulsing gold halo in the Tonnetz.

**Timbre-dependent consonance**
Harmonic amplitudes per voice follow standard timbre models:
- `STRINGS:  amp[k] = exp(-0.15k) / k`  (rich but smooth)
- `SAWTOOTH: amp[k] = 1/k`              (bright, audible roughness)
- `BRASS:    amp[k] = k/6 · exp(-0.3(k-6))` (rising-then-falling)

### 3. Markov Chains & Combinatorics

**Progression Markov chain**
Transition probabilities derived from Bach 371 chorales:
```
I   → {IV:22%, V:30%, vi:18%, ii:10%, V7:8%, ...}
V7  → {I:75%, vi:14%, IV:6%, ...}
vii° → {I:65%, V:16%, vi:10%, ...}
```

Query: given current chord `(root, type, key)`, return ranked successors
with probability, Roman numeral, and Tonnetz delta coordinates.

**Completion suggestions**
Given active pitch classes `{p₁, ..., pₙ}`, rank all 12 pitch classes
by consonance score:
```
score(p) = 1 - ΔRoughness(p) / (n+1)
```
where `ΔRoughness` uses the Plomp-Levelt interval roughness table.

**EDO error analysis**
For EDO-n, compute the error vector for primes {3, 5, 7, 11, 13}:
```
error(p, n) = round(log₂(p) · n) · (1200/n) - log₂(p)·1200  [cents]
```

Key results:
- 12-TET: 5th error −2¢, M3 error +14¢ (usable but impure thirds)
- 31-TET: 5th error −5¢, M3 error +1¢, 7th error +1¢ (near-JI 7-limit)
- 53-TET: 5th error −0.07¢ (Pythagorean near-perfection)

---

## Key Interactions

| Action | Result |
|--------|--------|
| Click Tonnetz node | Add voice (or remove if already active) |
| Middle-drag Tonnetz | Pan the lattice |
| Scroll on Tonnetz | Zoom in/out |
| Timbre dropdown | Change harmonic series profile |
| Detune slider | Microtonal pitch adjustment (±50¢) |
| "Suggest from current chord" | Markov next-chord suggestions |
| "Suggest completion note" | Roughness-optimal voice addition |
| "Analyse current EDO" | JI error vector for current EDO setting |
| Click chord in browser | Highlight Tonnetz path + set progression |

---

## Signal Flow

```
Voice.frequency + Voice.harmonic_amp[]
    │
    ▼  (per sample, audio thread)
Σₖ amp[k] · sin(2π·k·f·t)  × ADSR envelope
    │
    ▼
Soft clip: y = x / (1+|x|)
    │
    ▼
ALSA PCM (S16LE, 44100 Hz, stereo)

Every 4096 samples (≈93ms):
    │
    ▼
Pairwise Sethares roughness across all voice pairs
    +
Virtual pitch (autocorrelation of summed log-spectrum)
    +
Chord identification (pitch-class → template matching)
    │
    ▼
AbstractObject → Tonnetz halo + UI label
```

---

## File Structure

```
harmonia/
├── CMakeLists.txt
├── build.sh
├── README.md
├── src/
│   ├── main.cpp              FLTK window, panels, idle loop
│   ├── voice.h               Voice data structure + timbre presets
│   ├── audio_engine.h/.cpp   ALSA synthesis, roughness, virtual pitch
│   ├── tonnetz_widget.h/.cpp OpenGL Tonnetz lattice display
│   └── theory_bridge.h/.cpp  Python subprocess IPC (JSON)
└── theory/
    └── server.py             Theory server:
                                - Markov chain progressions
                                - JI lattice (ji_ratio)
                                - Orbifold distance (Tymoczko)
                                - EDO error analysis
                                - Completion suggestions
```

---

## Extending

**Add a new timbre preset:**
In `voice.h`, add case to `setTimbre()` switch with custom `harmonic_amp[]` profile.

**Add a new theory command:**
In `theory/server.py`, add `elif cmd == "my_cmd":` handler in `handle()`.
In `theory_bridge.h/.cpp`, add corresponding `queryMyCmd()` method.

**Microtonal tuning:**
Use the Detune slider (±50¢) per voice, or set `Voice.frequency` directly to
any JI ratio by computing `f = root_hz * num/den`.

**Change EDO:**
Set the EDO spinner. Click "Analyse current EDO" to see the JI error vector.
The Tonnetz will re-project voice positions at the new resolution.

---

## References

- Tymoczko, D. (2006). The geometry of musical chords. *Science*, 313(5783), 72-74.
- Sethares, W. A. (1993). Local consonance and the relationship between timbre and scale. *JASA*, 94(3), 1218-1228.
- Plomp, R. & Levelt, W. J. M. (1965). Tonal consonance and critical bandwidth. *JASA*, 38(4), 548-560.
- Clough, J. & Douthett, J. (1991). Maximally even sets. *Journal of Music Theory*, 35(1/2), 93-173.
- Euler, L. (1739). *Tentamen novae theoriae musicae*. (Tonnetz precursor)
- Bellmund, J. L. S. et al. (2018). Navigating cognition: Spatial codes for human thinking. *Science*, 362.
