# Harmonia V5 — Psychoacoustic Multi-Voice Music Coder (Circular Tonal Space Edition)

A Linux FLTK/OpenGL application that demystifies tuning, harmony and voice-leading through the lens of **circular pitch class space**, **psychoacoustics** and **neural auditory cortex models**.

Version 5 introduces a **Circular Tonal Space** visualization that replaces the traditional Tonnetz to better represent the full tonal range of any Equal Division of the Octave (EDO).

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  FLTK/OpenGL Application (C++)                                  │
│                                                                 │
│  ┌──────────────┐  ┌─────────────────────┐  ┌───────────────┐ │
│  │ Voice Panel  │  │  Tonal Space (GL)   │  │ Theory Panel  │ │
│  │  • ADSR      │  │  • Circular Layout  │  │  • Analysis   │ │
│  │  • Timbre    │  │  • Octave Rings     │  │  • Neural Psy  │ │
│  │  • Detune    │  │  • Roughness lines  │  │               │ │
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
│  Python Theory Server (subprocess - server.py)                  │
│  • Neural Auditory Cortex Model (Peripheral/Brainstem/Cortical) │
│  • Circular Tonal Space mapping                                 │
│  • Tymoczko orbifold voice-leading distances                    │
│  • Terhardt Subharmonic Summation (Virtual Pitch)               │
│  • Krumhansl-Kessler Tonal Hierarchy Analysis                   │
└─────────────────────────────────────────────────────────────────┘
```

---

## Key Features in V5

### Circular Tonal Space
The new visualization arranges all available pitch classes of an EDO around a circle.
- **Concentric Rings:** Represent octaves (from octave 2 to 6).
- **Radial Slices:** Represent pitch classes.
- **Harmonic Connections:** Lines are drawn between nodes representing the nearest fifth and major third intervals, making the harmonic structure of any EDO immediately visible.
- **Roughness Heat:** Real-time roughness between active voices is displayed as glowing lines directly connecting the nodes in the circular space.

### EDO Agnostic
Unlike the 5-limit Tonnetz which is biased towards 12-EDO, the Circular Tonal Space natively supports any EDO from 5 to 72 divisions, ensuring all notes are equally accessible and visible.

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

## Key Interactions

| Action | Result |
|--------|--------|
| Click node | Toggle voice at specific Pitch Class and Octave |
| Middle-drag | Rotate the circular space |
| Scroll | Zoom in/out |
| Timbre dropdown | Change harmonic series profile |
| Detune slider | Microtonal pitch adjustment (±50¢) |
| Change EDO | Update the entire tonal space for the new resolution |

---

## File Structure

```
harmonia_v5/
├── CMakeLists.txt
├── build.sh
├── README.md
├── src/
│   ├── main.cpp                 FLTK window, panels, idle loop
│   ├── voice.h                  Voice data structure + timbre presets
│   ├── audio_engine.h/.cpp      ALSA synthesis, roughness, virtual pitch
│   ├── tonal_space_widget.h/.cpp Circular tonal space display
│   └── theory_bridge.h/.cpp     Python subprocess IPC (JSON)
└── theory/
    └── server.py                Theory server
```
