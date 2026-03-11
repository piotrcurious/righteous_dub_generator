#!/usr/bin/env python3
"""
Harmonia Theory Server  ·  Psychoacoustic Edition
══════════════════════════════════════════════════
Algebraic geometry + psychoacoustics + neural auditory cortex model.

Architecture
────────────
MODULE 0  PSYCHOACOUSTIC ENGINE
   0.1  Frequency / perceptual-scale conversions  (Bark, ERB, Greenwood)
   0.2  Cochlear mechanics  (critical bandwidth, basilar membrane)
   0.3  Inner hair cell model  (compression, half-wave rectification)
   0.4  Plomp–Levelt roughness  (exact CB-normalized smooth curve)
   0.5  Sethares spectral consonance  (full harmonic-series interaction)
   0.6  Auditory-nerve firing-rate model  (Sachs–Abbas sigmoid)
   0.7  Virtual pitch / Terhardt SHS  (brainstem periodicity extraction)
   0.8  Spectral masking  (Moore & Glasberg asymmetric excitation spread)
   0.9  Krumhansl–Kessler tonal hierarchy  (cortical expectation model)
   0.10 Auditory cortex integration  (peripheral + brainstem + cortical)

MODULE 1  SET-CLASS THEORY          (Allen Forte)
MODULE 2  NEO-RIEMANNIAN PLR        (Tonnetz transformations)
MODULE 3  FUNCTIONAL ANALYSIS       (algebraic containment)
MODULE 4  VOICE-LEADING ORBIFOLD    (Hungarian optimal assignment)
MODULE 5  TONNETZ TENSION           (geometric + psychoacoustic hybrid)
MODULE 6  PATTERN RECOGNITION       (sequence / PLR / circle-of-fifths)
MODULE 7  EDO / JI LATTICE          (prime-limit ratio analysis)
MODULE 8  VOICE COMPLETION          (psychoacoustic scoring)
MODULE 9  PIVOT CHORD SEARCH        (modulation theory)
   9.1  Roman numeral system            (diatonic + chromatic degrees)
   9.2  Diatonic chord inventory        (triads + seventh chords)
   9.3  Common-chord pivots             (shared diatonic membership)
   9.4  Secondary dominant pivots       (applied dominants V/x)
   9.5  Enharmonic pivots               (dim7 × 4, GerAug6 ↔ Dom7, aug+ × 3)
   9.6  Borrowed / modal-mixture pivots (parallel minor ♭III ♭VI ♭VII)
   9.7  Chromatic mediant pivots        (PLR-L / PLR-R cross-key)
   9.8  Pivot scoring                   (common-tone × KK × VL × tension)
   9.9  Modulation path search          (Dijkstra over dual-key chord graph)
"""

import sys
import json
import math
import itertools
from typing import List, Dict, Tuple, Optional, Set, FrozenSet, Any
from fractions import Fraction

NOTE_NAMES   = ['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
NOTE_NAMES_b = ['C','Db','D','Eb','E','F','Gb','G','Ab','A','Bb','B']


def nn(pc: int, prefer_flat: bool = False, edo: int = 12) -> str:
    """
    Return a human-readable name for pitch-class *pc* in the given EDO.

    EDO-12  : standard letter names (C, C#/Db, …)
    Other   : nearest 12-tone letter name with a ±N¢ deviation tag when
              the EDO step does not land exactly on a 12-tone semitone grid
              point (i.e. deviation > 2 cents).  This keeps labels compact
              (e.g. "D+39¢" for EDO-31 step 5) while remaining unique and
              musically meaningful for any EDO up to at least 74.
    """
    pc = pc % edo
    if edo == 12:
        return (NOTE_NAMES_b if prefer_flat else NOTE_NAMES)[pc]
    # Map to nearest 12-tone semitone (0–11)
    exact_semitones = pc * 12.0 / edo          # fractional semitones above C
    idx12 = round(exact_semitones) % 12
    base  = (NOTE_NAMES_b if prefer_flat else NOTE_NAMES)[idx12]
    # Deviation in cents from the nearest 12-tone grid point
    dev_cents = round((exact_semitones - round(exact_semitones)) * 100)
    if dev_cents == 0:
        return base
    sign = '+' if dev_cents > 0 else ''
    return f"{base}{sign}{dev_cents}¢"


# ──────────────────────────────────────────────────────────────────────────────
# Utilities: modular distances
# ──────────────────────────────────────────────────────────────────────────────
def mod_signed_dist(a: int, b: int, modulus: int = 12) -> int:
    d = (b - a) % modulus
    if d > modulus // 2:
        d -= modulus
    return d

def mod_abs_dist(a: int, b: int, modulus: int = 12) -> int:
    return abs(mod_signed_dist(a, b, modulus))


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 0  —  PSYCHOACOUSTIC ENGINE
# ──────────────────────────────────────────────────────────────────────────────

# ── 0.1  Frequency / perceptual-scale conversions ────────────────────────────

def pc_to_hz(pc: int, octave: int = 4, c4_hz: float = 261.625565, edo: int = 12) -> float:
    """Pitch-class + octave → frequency in Hz.  PC 0 = C, octave 4 → 261.625565 Hz."""
    return c4_hz * 2.0 ** (octave - 4 + pc / float(edo))

def midi_to_hz(midi: float) -> float:
    return 440.0 * 2.0 ** ((midi - 69.0) / 12.0)

def hz_to_midi(hz: float) -> float:
    return 69.0 + 12.0 * math.log2(max(hz, 1e-9) / 440.0)

def hz_to_bark(f: float) -> float:
    """
    Traunmüller (1990) Bark scale.
    Accurate to ±0.05 Bark over 20–20000 Hz.
    """
    if f <= 0.0:
        return 0.0
    return (26.81 * f / (1960.0 + f)) - 0.53

def bark_to_hz(z: float) -> float:
    """Inverse Bark scale (Traunmüller 1990)."""
    return 1960.0 * (z + 0.53) / (26.28 - (z + 0.53))

def hz_to_erb_number(f: float) -> float:
    """ERB-rate scale — number of ERBs below frequency f (Moore & Glasberg 1983)."""
    return 21.4 * math.log10(max(1.0, 1.0 + 0.00437 * f))

def erb_bandwidth(f: float) -> float:
    """
    Equivalent Rectangular Bandwidth at centre frequency f Hz.
    Moore & Glasberg (1983): ERB(f) = 24.7*(4.37*f/1000 + 1).
    """
    return 24.7 * (4.37 * f / 1000.0 + 1.0)

def critical_bandwidth_zwicker(f: float) -> float:
    """
    Critical bandwidth (Hz) — Zwicker & Terhardt (1980) approximation.
    More accurate than the older Plomp–Levelt table near 1 kHz.
    """
    return 25.0 + 75.0 * (1.0 + 1.4 * (f / 1000.0) ** 2) ** 0.69


# ── 0.2  Basilar membrane / cochlear place map ────────────────────────────────

def greenwood_place(f: float) -> float:
    """
    Greenwood (1990) cochlear frequency-place map.
    Returns basilar membrane position in mm from apex.
    Human cochlea ≈ 35 mm total length.
    Formula: x = (1/a) * log10(f/f0 + k),  a=0.06, f0=165.4 Hz, k=0.88
    """
    if f <= 0.0:
        return 0.0
    return (1.0 / 0.06) * math.log10(f / 165.4 + 0.88)

def greenwood_freq(x_mm: float) -> float:
    """Inverse Greenwood map: basilar membrane position → frequency."""
    return 165.4 * (10.0 ** (0.06 * x_mm) - 0.88)

def tonotopic_distance_mm(f1: float, f2: float) -> float:
    """
    Physical distance on the basilar membrane between two frequencies.
    This is the most fundamental psychoacoustic metric — closer on the BM
    means more interaction (roughness, masking, fusion).
    """
    return abs(greenwood_place(f1) - greenwood_place(f2))


# ── 0.3  Inner hair cell (IHC) model ─────────────────────────────────────────

def ihc_compression(amplitude: float,
                    level_db: float = 70.0,
                    knee_spl: float = 40.0,
                    compression_ratio: float = 0.25) -> float:
    """
    IHC input-output compressive nonlinearity (Meddis 1986 / Glasberg & Moore 1990).

    Implements the ~4:1 compression (ratio 0.25) observed on the basilar
    membrane above the compressive knee (~40 dB SPL at mid-frequency).

    Parameters
    ----------
    amplitude         : relative linear amplitude in [0, 1] (1 = loudest partial)
    level_db          : presentation level of the loudest partial in dB SPL
    knee_spl          : compression onset in dB SPL (Meddis 1986: ~40 dB)
    compression_ratio : output/input dB slope above knee (0.25 -> 4:1 ratio)

    Returns compressed linear amplitude (scale-compatible with the input).
    """
    if amplitude <= 0.0:
        return 0.0
    # Convert relative amplitude to dB SPL
    spl = level_db + 20.0 * math.log10(amplitude + 1e-12)
    if spl <= knee_spl:
        return amplitude                         # linear region
    # Compressed SPL above knee
    spl_out = knee_spl + (spl - knee_spl) * compression_ratio
    # Convert back to linear amplitude relative to reference (level_db)
    return 10.0 ** ((spl_out - level_db) / 20.0)

def ihc_halfwave(x: float) -> float:
    """Half-wave rectification (mechanical → neural transduction at IHC stereocilia)."""
    return max(0.0, x)

def ihc_lowpass(x: float, cutoff: float = 3500.0, fs: float = 22050.0) -> float:
    """
    First-order IHC low-pass filter (3.5 kHz cutoff, Meddis 1986).
    Approximated here as a scalar gain for spectral-domain use.
    Returns the fraction of energy passed for a pure tone of given cutoff.
    g = 1 / sqrt(1 + (f/fc)^2)
    """
    return 1.0 / math.sqrt(1.0 + (cutoff / max(fs, 1.0)) ** 2)


# ── 0.4  Plomp–Levelt roughness ───────────────────────────────────────────────

def plomp_levelt_roughness_pair(f1: float, f2: float,
                                a1: float = 1.0, a2: float = 1.0) -> float:
    """
    Sensory roughness between two pure sinusoids.

    Uses the Vassilakis (2001) re-parameterisation of the Sethares / Plomp-Levelt
    roughness model.  Key improvements over the plain Sethares formula:

      1. Critical bandwidth via ERB (Moore & Glasberg 1983) rather than the
         cruder Plomp (1964) linear approximation, giving better accuracy
         across the full audible range.

      2. Vassilakis amplitude-weighting factor:
             A = (a1·a2)^0.1 · [2·a_min/(a1+a2)]^3.11
         This correctly captures:
           · Joint loudness contribution  — (a1·a2)^0.1
           · Amplitude asymmetry penalty — [2·a_min/(a1+a2)]^3.11
             (maximum at equal amplitudes, zero when one tone vanishes)

      3. Shape constants b1=3.5, b2=5.75 (Sethares 1997 / Vassilakis 2001),
         giving peak roughness at s ≈ 0.25 CB with correct asymmetric fall-off.

    References
    ----------
    Vassilakis (2001) Perceptual and Physical Properties of Amplitude Fluctuation.
    Sethares (1997)   Tuning Timbre Spectrum Scale.
    Plomp & Levelt (1965) JASA.
    """
    if f1 <= 0.0 or f2 <= 0.0:
        return 0.0
    f_lo = min(f1, f2)
    # ERB-based critical bandwidth (more accurate than Plomp 1964 linear fit)
    cb = 24.7 * (4.37 * f_lo / 1000.0 + 1.0)
    s  = abs(f2 - f1) / cb
    if s == 0.0:
        return 0.0
    # Vassilakis amplitude weighting
    a_min = min(a1, a2)
    amp_weight = ((a1 * a2) ** 0.1) * ((2.0 * a_min / (a1 + a2 + 1e-12)) ** 3.11)
    roughness = amp_weight * (math.exp(-3.5 * s) - math.exp(-5.75 * s))
    return max(0.0, roughness)

def plomp_levelt_curve_at_f(f_center: float, n_points: int = 200
                            ) -> List[Tuple[float, float]]:
    """
    Returns (delta_hz, roughness) pairs for a P-L roughness curve centred at
    f_center, sweeping from 0 to 2·CB.  Useful for visualisation / debugging.
    """
    cb   = 0.021 * f_center + 19.0
    step = 2.0 * cb / n_points
    return [(i * step,
             plomp_levelt_roughness_pair(f_center, f_center + i * step))
            for i in range(n_points + 1)]


# ── 0.5  Sethares spectral consonance ────────────────────────────────────────

def _build_harmonic_series(f0: float, n: int,
                            rolloff: float = 0.88) -> Tuple[List[float], List[float]]:
    """
    Build harmonic series for a single pitch.
    Returns (frequencies, amplitudes) where amp_k = rolloff^(k-1).
    rolloff ≈ 0.88 matches a typical bowed string / piano tone;
    use 0.70 for brass, 0.95 for flute/sine.
    """
    freqs = [f0 * k for k in range(1, n + 1)]
    amps  = [rolloff ** (k - 1) for k in range(1, n + 1)]
    return freqs, amps

def sethares_roughness(freqs: List[float], amps: List[float]) -> float:
    """
    Total roughness of a multi-partial spectrum (Sethares 1993, 1997).
    Sums P-L roughness over all i < j pairs.
    freqs, amps must be co-indexed.
    """
    total = 0.0
    n = len(freqs)
    for i in range(n):
        for j in range(i + 1, n):
            total += plomp_levelt_roughness_pair(
                freqs[i], freqs[j], amps[i], amps[j])
    return total

def chord_roughness_psychoacoustic(
        pitch_classes : List[int],
        c4_hz         : float = 261.625565,
        octave        : int   = 4,
        n_harmonics   : int   = 8,
        rolloff       : float = 0.88,
        close_position: bool  = True,
        edo           : int   = 12,
        level_db      : float = 70.0,
) -> Dict[str, Any]:
    """
    Full spectral roughness of a chord via Vassilakis/Sethares model + IHC compression.

    Notes are voiced in close position (default) so that no note is more than
    an octave above the bass — this models normal keyboard voicing and avoids
    artificial roughness from wide spacing.

    Normalization: roughness is divided by the roughness of a reference minor-2nd
    dyad (one EDO step apart) using the same harmonic series and level.  This
    gives norm=1.0 for the roughest possible interval and correct relative ordering
    (e.g. m2 > M2 > M3 > P5 > P8=0).

    Returns
    -------
    roughness_raw          : raw Vassilakis/Sethares sum
    roughness_normalized   : 0–1 (m2-reference-scaled)
    consonance_score       : 1 – roughness_normalized
    roughness_per_pair     : list of (pcA, pcB, roughness) for transparency
    note_frequencies_hz    : fundamental Hz for each voiced pitch class
    """
    if not pitch_classes:
        return {"roughness_raw": 0.0, "roughness_normalized": 0.0,
                "consonance_score": 1.0, "roughness_per_pair": [],
                "note_frequencies_hz": []}

    pcs = sorted(set(p % edo for p in pitch_classes))

    # ── Voice in close position ────────────────────────────────────────────
    if close_position:
        note_freqs: List[float] = []
        oct_cur = octave
        prev_pc = -1
        for pc in pcs:
            if prev_pc >= 0 and pc <= prev_pc:
                oct_cur += 1
            note_freqs.append(c4_hz * 2.0 ** (oct_cur - 4 + pc / float(edo)))
            prev_pc = pc
    else:
        note_freqs = [c4_hz * 2.0 ** (pc / float(edo)) for pc in pcs]

    # ── Build & compress partial spectra ─────────────────────────────────
    def _build_compressed(f0: float) -> Tuple[List[float], List[float]]:
        fs, as_ = _build_harmonic_series(f0, n_harmonics, rolloff)
        # IHC compression: operate in true dB-SPL space
        as_ = [ihc_compression(a, level_db=level_db) for a in as_]
        return fs, as_

    all_freqs: List[float] = []
    all_amps:  List[float] = []
    note_partial_idx: List[Tuple[int, int]] = []
    for f0 in note_freqs:
        start = len(all_freqs)
        fs, as_ = _build_compressed(f0)
        all_freqs.extend(fs)
        all_amps.extend(as_)
        note_partial_idx.append((start, len(all_freqs)))

    # Normalise amplitudes globally so max = 1
    max_a = max(all_amps) if all_amps else 1.0
    all_amps = [a / max_a for a in all_amps]

    # ── Total roughness (all pairs) ───────────────────────────────────────
    raw = sethares_roughness(all_freqs, all_amps)

    # ── Per-note-pair roughness breakdown ─────────────────────────────────
    pair_roughness: List[Dict] = []
    for ni in range(len(pcs)):
        for nj in range(ni + 1, len(pcs)):
            si, ei = note_partial_idx[ni]
            sj, ej = note_partial_idx[nj]
            r = 0.0
            for ii in range(si, ei):
                for jj in range(sj, ej):
                    r += plomp_levelt_roughness_pair(
                        all_freqs[ii], all_freqs[jj],
                        all_amps[ii],  all_amps[jj])
            pair_roughness.append({
                "pc_a": pcs[ni], "pc_b": pcs[nj],
                "name_a": nn(pcs[ni], edo=edo), "name_b": nn(pcs[nj], edo=edo),
                "roughness": round(r, 5),
                "interval_semitones": mod_abs_dist(pcs[ni], pcs[nj], modulus=edo)
            })

    # ── Reference normalization: minor-2nd dyad (1 EDO step) ─────────────
    # Build the same-style spectrum for a m2 dyad at the same register.
    # This gives a proper "maximum roughness" reference so that a m2 ≈ 1.0
    # and larger intervals scale correctly.
    f_ref = note_freqs[0] if note_freqs else c4_hz
    f_ref2 = f_ref * 2.0 ** (1.0 / float(edo))
    ref_f, ref_a1 = _build_compressed(f_ref)
    ref_f2, ref_a2 = _build_compressed(f_ref2)
    # Merge and normalise
    rf_all = ref_f + ref_f2
    ra_all = ref_a1 + ref_a2
    ra_max = max(ra_all) if ra_all else 1.0
    ra_all = [a / ra_max for a in ra_all]
    r_ref   = sethares_roughness(rf_all, ra_all)
    n_pairs = max(1, len(pcs) * (len(pcs) - 1) // 2)

    norm = (raw / (n_pairs * r_ref)) if r_ref > 1e-12 else 0.0
    norm = min(1.0, norm)
    consonance = max(0.0, 1.0 - norm)

    return {
        "roughness_raw"        : round(raw, 5),
        "roughness_normalized" : round(norm, 4),
        "consonance_score"     : round(consonance, 4),
        "roughness_per_pair"   : sorted(pair_roughness,
                                        key=lambda x: -x["roughness"]),
        "note_frequencies_hz"  : [round(f, 2) for f in note_freqs],
    }


# ── 0.6  Auditory-nerve firing-rate model ────────────────────────────────────

def auditory_nerve_rate(stimulus_db : float,
                        threshold_db: float = 0.0,
                        r_spont     : float = 0.05,
                        r_max       : float = 1.0,
                        slope_db    : float = 0.1) -> float:
    """
    Sachs & Abbas (1974) sigmoidal rate-level function for an AN fibre.

    r(L) = r_spont + (r_max − r_spont) / [1 + exp(−slope · (L − L_sat))]

    where L_sat = threshold + 20 dB  (saturation knee ≈ 20 dB above threshold).
    Returns normalised firing rate in [r_spont, r_max].
    """
    excess = stimulus_db - threshold_db
    if excess < 0.0:
        return r_spont
    l_sat = 20.0  # dB above threshold
    return r_spont + (r_max - r_spont) / (
        1.0 + math.exp(-slope_db * (excess - l_sat)))

def an_rate_place_profile(freqs: List[float], amps: List[float],
                          level_db: float = 70.0,
                          n_cfs: int = 40) -> List[Dict[str, Any]]:
    """
    Auditory-nerve rate-place profile: neural excitation as a function of
    characteristic frequency (CF) along the tonotopic axis.

    Simulates n_cfs AN fibres with CFs logarithmically spaced from
    100 Hz to 8000 Hz.  Each fibre integrates input from all partials
    weighted by an asymmetric roex(p,r) auditory filter (Patterson et al. 1987)
    approximated on the ERB-number scale:

        W(g) = (1 + p_low*g)  * exp(-p_low*g)   for  g <= 0   (lower side)
        W(g) = (1 + p_high*g) * exp(-p_high*g)  for  g > 0    (upper side)

    where g = |f - CF| / CF  and  p_low=31, p_high=15  (typical values from
    Patterson 1987 for a 1 kHz CF; scaled here by ERB ratio).

    The asymmetry (lower skirt steeper than upper) correctly captures
    the well-known upward-spread of masking and the sharper low-frequency
    rejection of auditory filters.

    Returns a list of {"cf_hz", "erbnum", "firing_rate"} dicts.
    """
    cf_lo, cf_hi = 100.0, 8000.0
    cfs = [cf_lo * (cf_hi / cf_lo) ** (i / (n_cfs - 1)) for i in range(n_cfs)]
    profile = []
    for cf in cfs:
        erb_cf = hz_to_erb_number(cf)
        erb_bw = erb_bandwidth(cf)      # ERB at this CF

        # Sum input from all partials weighted by roex(p,r) filter
        total_power = 1e-30             # power floor
        for f_part, a_part in zip(freqs, amps):
            # Normalised deviation on ERB scale
            g = abs(f_part - cf) / max(erb_bw, 1.0)
            # Asymmetric roex: p_low=31 (steep lower skirt), p_high=15 (shallower upper)
            p = 31.0 if f_part <= cf else 15.0
            filter_gain = (1.0 + p * g) * math.exp(-p * g)
            # dB of this partial at this fibre
            part_db = level_db + 20.0 * math.log10(max(a_part, 1e-9))
            part_db_filtered = part_db + 20.0 * math.log10(max(filter_gain, 1e-9))
            total_power += 10.0 ** (part_db_filtered / 10.0)

        total_input_db = 10.0 * math.log10(total_power)

        # Absolute threshold varies with CF: ISO 226 simplified
        thresh = threshold_in_quiet(cf)
        rate = auditory_nerve_rate(total_input_db, threshold_db=thresh)
        profile.append({"cf_hz": round(cf, 1), "erbnum": round(erb_cf, 2),
                        "firing_rate": round(rate, 4)})
    return profile


# ── 0.7  Virtual pitch — Terhardt subharmonic summation (SHS) ────────────────

def virtual_pitch_strength(
        pitch_classes: List[int],
        c4_hz        : float = 261.625565,
        octave       : int   = 4,
        n_harmonics  : int   = 8,
        edo          : int   = 12,
) -> Dict[str, Any]:
    """
    Virtual (residue) pitch detection via Terhardt's Subharmonic Summation (1982).

    Algorithm
    ---------
    For each input tone f_n and each subharmonic rank k = 1..N_sub:
        F0_candidate = f_n / k

    For each candidate F0, the pitch salience is:
        score(F0) = sum_n  sum_k  w_harm(k) * w_dom(k*F0) * gaussian_match(f_n, k*F0)

    where:
        w_harm(k)    = 1/k                (harmonic weighting; lower harmonics dominate)
        w_dom(f)     = spectral dominance region weight (Terhardt 1972):
                       Gaussian centred at 700 Hz on a log-frequency axis, with
                       sigma_log = 1.5 octaves; ensures partials 500–2000 Hz
                       carry the most virtual-pitch information.
        gaussian_match = exp(-0.5*(f_n - k*F0)^2 / sigma_f^2)
                       with sigma_f = 0.012 * k * F0  (~20 cents tolerance)

    The normalised harmonicity (0-1) indicates how closely the chord resembles a
    single harmonic series — high values predict strong fusion and rootedness.

    Key improvements over the naive SHS implementation
    ---------------------------------------------------
    * Fine search grid: 1/8-step resolution (vs. full EDO step)
    * Candidates seeded from actual subharmonics of input tones — prevents
      trivial single-note matching from dominating
    * Spectral dominance weighting — mid-frequency partials count more
    * Tighter gaussian tolerance (20 cents instead of 35 cents)
    """
    if not pitch_classes:
        return {"virtual_pitch_hz": 0.0, "virtual_pitch_pc": 0,
                "virtual_pitch_name": "?", "harmonicity": 0.0,
                "harmonicity_label": "none", "top_candidates": []}

    pcs = sorted(set(p % edo for p in pitch_classes))
    note_freqs = [c4_hz * 2.0 ** (pc / float(edo)) for pc in pcs]

    n_sub = n_harmonics      # subharmonic ranks to consider per input tone
    n_eval = n_harmonics     # harmonic ranks to check per candidate F0

    # ── Build F0 candidate set from actual subharmonics ────────────────────────
    # For each input tone f_n, F0 = f_n/k for k = 1..n_sub
    # Cluster nearby candidates by snapping to a fine log-frequency grid
    # (1/8 EDO step resolution).
    f_min_candidate = min(note_freqs) / (n_sub + 0.5)
    f_max_candidate = max(note_freqs) * 1.02

    candidate_set: Set[float] = set()
    for f_n in note_freqs:
        for k in range(1, n_sub + 1):
            f0_cand = f_n / k
            if f_min_candidate < f0_cand < f_max_candidate:
                # Snap to 1/8-EDO-step grid to cluster near-identical candidates
                grid_idx = round(8 * edo * math.log2(f0_cand / f_min_candidate))
                snapped  = f_min_candidate * 2.0 ** (grid_idx / (8.0 * edo))
                candidate_set.add(snapped)

    if not candidate_set:
        candidate_set.add(note_freqs[0])

    # ── Spectral dominance weighting (Terhardt 1972) ───────────────────────────
    # Peak at ~700 Hz on log scale, sigma = 1.5 octaves
    _dom_f0    = 700.0
    _dom_sigma = 1.5   # octaves

    def spectral_dominance(f: float) -> float:
        if f <= 0.0:
            return 0.0
        log_ratio = math.log2(f / _dom_f0)
        return math.exp(-0.5 * (log_ratio / _dom_sigma) ** 2)

    # ── Score every candidate ─────────────────────────────────────────────────
    # Per-note best-k scoring with coherence weighting.
    # For each candidate F0, find the best harmonic rank k for each input tone,
    # then multiply by a coherence factor = (fraction of tones with meaningful
    # support)^2.  This prevents a single-note trivial match (F0 = one input note,
    # k=1) from dominating over a true common fundamental.
    _contribution_threshold = 0.04   # minimum per-note contribution to count

    candidates: List[Tuple[float, float]] = []
    for f0 in candidate_set:
        per_note_best: List[float] = []
        for f_note in note_freqs:
            best = 0.0
            for k in range(1, n_eval + 1):
                f_harm = k * f0
                if f_harm < f_note * 0.5 or f_harm > f_note * 2.0:
                    continue
                sigma_f = 0.012 * f_harm
                match   = math.exp(-0.5 * ((f_note - f_harm) / (sigma_f + 1e-9)) ** 2)
                w_harm  = 1.0 / k
                w_dom   = spectral_dominance(f_harm)
                best    = max(best, w_harm * w_dom * match)
            per_note_best.append(best)

        n_contributing = sum(1 for s in per_note_best if s > _contribution_threshold)
        coherence      = n_contributing / max(len(note_freqs), 1)
        total_score    = sum(per_note_best) * (coherence ** 2)
        candidates.append((f0, total_score))

    best_f0, best_score = max(candidates, key=lambda x: x[1])

    # ── Convert best F0 to pitch class ────────────────────────────────────────
    best_pc = round(edo * math.log2(max(best_f0, 1e-9) / c4_hz)) % edo

    # ── Harmonicity: calibrated normalization ────────────────────────────────────
    # The upper-bound reference is: the score achievable if each note were its
    # OWN fundamental (k=1, perfect match, coherence=1).  This is independent
    # of the search result and gives a meaningful ceiling:
    #
    #   ideal_upper = coherence_weight * sum_n { w_dom(f_n) }
    #
    # where coherence_weight = (n_notes / n_notes)^2 = 1.0 for a perfect chord.
    # Dividing by this yields harmonicity_raw in (0,1], where:
    #   ~1.0 → chord notes are exact harmonics of a common F0 (e.g. C-maj 4:5:6)
    #   ~0.3 → chord has some harmonic relationship (minor, dom7)
    #   ~0.1 → chord has no common harmonic series (cluster, aug, dim7)
    #
    # A calibration scale of 3.5 stretches the range so that the labels
    # ("high", "medium", "low") cover musically meaningful bands:
    #   high   (> 0.55) : pure intervals, major triads
    #   medium (> 0.25) : minor triads, dominant 7ths
    #   low    (≤ 0.25) : diminished, augmented, cluster
    _HARM_SCALE = 3.5

    ideal_upper = sum(spectral_dominance(f) for f in note_freqs)
    harmonicity_raw = best_score / max(ideal_upper, 1e-9)
    harmonicity     = min(1.0, harmonicity_raw * _HARM_SCALE)

    top5 = sorted(candidates, key=lambda x: -x[1])[:5]
    top5_out = []
    for f, s in top5:
        if f <= 0.0:
            continue
        pc_cand = round(edo * math.log2(f / c4_hz)) % edo
        top5_out.append({"f0_hz": round(f, 2), "pc": pc_cand,
                         "name": nn(pc_cand, edo=edo), "score": round(s, 4)})

    return {
        "virtual_pitch_hz"    : round(best_f0, 2),
        "virtual_pitch_pc"    : best_pc,
        "virtual_pitch_name"  : nn(best_pc, edo=edo),
        "harmonicity"         : round(harmonicity, 4),
        "harmonicity_label"   : ("high"   if harmonicity > 0.55 else
                                 "medium" if harmonicity > 0.25 else "low"),
        "top_candidates"      : top5_out,
    }


# ── 0.8  Spectral masking (Moore & Glasberg 1987) ────────────────────────────

def masking_excitation_db(masker_hz: float, masker_db: float,
                           probe_hz : float) -> float:
    """
    Excitation level (dB) at probe frequency due to a simultaneous masker tone.

    Uses the Moore & Glasberg (1987) asymmetric spreading function on the Bark
    scale with level-dependent upward-spread slope:

        Upward   slope  = -( 27 - 0.37 * max(0, masker_db - 40) )  dB/Bark
                          (slope shallows as level rises: more upward spread)
        Downward slope  = -40 dB/Bark  (downward spread, level-independent)

    The result is clamped to the absolute threshold of hearing so that masked
    tones below the quiet threshold are always flagged as inaudible.
    """
    if masker_hz <= 0.0 or probe_hz <= 0.0:
        return -200.0
    delta_bark = hz_to_bark(probe_hz) - hz_to_bark(masker_hz)
    if delta_bark >= 0.0:
        # Upward spread (probe above masker): level-dependent shallow slope
        slope = -(27.0 - 0.37 * max(0.0, masker_db - 40.0))
    else:
        # Downward spread: steeper, level-independent
        slope = -40.0
    excitation = masker_db + slope * abs(delta_bark)
    # Clamp to the absolute threshold of hearing (floor)
    return max(excitation, threshold_in_quiet(probe_hz))

def threshold_in_quiet(f: float) -> float:
    """
    Absolute threshold of hearing (ISO 226:2003 simplified fit).
    Returns approximate threshold in dB SPL.
    Valid ≈ 20–16000 Hz.
    """
    if f <= 0.0: return 120.0
    f_khz = f / 1000.0
    return (3.64 * f_khz ** -0.8
            - 6.5  * math.exp(-0.6 * (f_khz - 3.3) ** 2)
            + 1e-3 * f_khz ** 4)

def chord_masking_analysis(pitch_classes: List[int],
                           c4_hz       : float = 261.625565,
                           level_db    : float = 70.0,
                           edo         : int   = 12) -> Dict[str, Any]:
    """
    Psychoacoustic masking analysis for a chord.

    Each pitch class is treated as a pure tone at the given level.
    For each tone, the masking threshold is computed as:
        T_mask = max(T_quiet(f), max_over_j { excitation_j(f) })

    A tone is flagged 'masked' if its level is below T_mask.
    Sensation level = level_db − T_mask.

    Also reports the inter-tone CB overlap (tonotopic proximity) — the primary
    driver of roughness and fusion in the cochlea.
    """
    if not pitch_classes:
        return {"masked_tones": [], "tone_audibility": [], "cb_overlaps": []}

    pcs  = sorted(set(p % edo for p in pitch_classes))
    freqs = [c4_hz * 2.0 ** (pc / float(edo)) for pc in pcs]

    audibility = []
    for i, (pc, f_probe) in enumerate(zip(pcs, freqs)):
        t_quiet = threshold_in_quiet(f_probe)
        mask_th = t_quiet
        for j, f_mask in enumerate(freqs):
            if i != j:
                mask_th = max(mask_th,
                              masking_excitation_db(f_mask, level_db, f_probe))
        sl = level_db - mask_th
        audibility.append({
            "pc": pc, "name": nn(pc, edo=edo),
            "freq_hz": round(f_probe, 2),
            "masking_threshold_db": round(mask_th, 1),
            "sensation_level_db"  : round(sl, 1),
            "audible": sl > 0.0,
        })

    # Critical-band overlap between each pair
    cb_overlaps = []
    for i in range(len(pcs)):
        for j in range(i + 1, len(pcs)):
            cb_lo  = critical_bandwidth_zwicker(freqs[i])
            cb_hi  = critical_bandwidth_zwicker(freqs[j])
            delta  = abs(freqs[j] - freqs[i])
            avg_cb = 0.5 * (cb_lo + cb_hi)
            overlap = max(0.0, 1.0 - delta / avg_cb)   # 1 = full overlap, 0 = no
            bm_dist = tonotopic_distance_mm(freqs[i], freqs[j])
            cb_overlaps.append({
                "pc_a": pcs[i], "pc_b": pcs[j],
                "interval_semitones": mod_abs_dist(pcs[i], pcs[j], modulus=edo),
                "delta_hz": round(delta, 2),
                "cb_overlap_fraction": round(overlap, 3),
                "bm_distance_mm": round(bm_dist, 3),
            })

    masked = [t["pc"] for t in audibility if not t["audible"]]
    return {
        "masked_tones"     : masked,
        "masked_tone_names": [nn(p, edo=edo) for p in masked],
        "tone_audibility"  : audibility,
        "cb_overlaps"      : cb_overlaps,
    }


# ── 0.9  Krumhansl–Kessler tonal hierarchy ───────────────────────────────────

# Empirically measured probe-tone salience ratings (Krumhansl & Kessler 1982).
# Index i = degree i semitones above key root.
KK_MAJOR = [6.35, 2.23, 3.48, 2.33, 4.38, 4.09,
            2.52, 5.19, 2.39, 3.66, 2.29, 2.88]
KK_MINOR = [6.33, 2.68, 3.52, 5.38, 2.60, 3.53,
            2.54, 4.75, 3.98, 2.69, 3.34, 3.17]

def tonal_hierarchy_score(pitch_classes: List[int], key: int,
                          mode: str = "major", edo: int = 12) -> Dict[str, Any]:
    """
    Compute tonal stability of a chord using the Krumhansl-Kessler (1982)
    probe-tone profiles -- one of the most replicated findings in music cognition.

    The KK profiles reflect the *long-term memory* component of tonal tension,
    modelling how strongly each scale degree activates key-specific neural
    templates in auditory cortex.

    Normalization fix: uses the full range across both major and minor profiles
    (global min/max) so that I major chord correctly maps to ~1.0 ('very stable')
    rather than landing just below the threshold.

    For non-12 EDOs: linearly interpolates between adjacent 12-tone profile
    entries instead of rounding, giving smooth scores for microtonal intervals.

    Returns:
      kk_raw         : mean KK rating for the chord's PCs
      kk_normalized  : 0-1 (global min/max of combined M/m profiles)
      tonal_stability: same as kk_normalized, descriptive label
    """
    profile = KK_MAJOR if mode.lower() in ("major", "maj", "ionian") else KK_MINOR
    if not pitch_classes:
        return {"kk_raw": 0.0, "kk_normalized": 0.0,
                "tonal_stability": 0.0, "tonal_stability_label": "none"}
    pcs = list(set(p % edo for p in pitch_classes))

    def _kk_lookup(deg_edo: int) -> float:
        """Interpolated KK profile lookup for any EDO degree."""
        frac = (deg_edo % edo) * 12.0 / edo   # fractional 12-tone index
        lo   = int(frac) % 12
        hi   = (lo + 1) % 12
        t    = frac - int(frac)
        return profile[lo] * (1.0 - t) + profile[hi] * t

    raw = sum(_kk_lookup((pc - key) % edo) for pc in pcs) / len(pcs)

    # Normalization: use global min/max across both profiles for consistent scaling
    _kk_global_min = min(min(KK_MAJOR), min(KK_MINOR))   # 2.23
    _kk_global_max = max(max(KK_MAJOR), max(KK_MINOR))   # 6.35
    norm = (raw - _kk_global_min) / (_kk_global_max - _kk_global_min)
    norm = max(0.0, min(1.0, norm))

    label = ("very stable"   if norm > 0.70 else
             "stable"        if norm > 0.52 else
             "mildly stable" if norm > 0.28 else
             "unstable"      if norm > 0.10 else "very unstable")
    return {"kk_raw": round(raw, 3), "kk_normalized": round(norm, 4),
            "tonal_stability": round(norm, 4), "tonal_stability_label": label}


# ── 0.10  Auditory cortex integration model ───────────────────────────────────

def auditory_cortex_model(
        pitch_classes: List[int],
        key          : int   = 0,
        c4_hz        : float = 261.625565,
        octave       : int   = 4,
        level_db     : float = 70.0,
        n_harmonics  : int   = 8,
        edo          : int   = 12,
) -> Dict[str, Any]:
    """
    Three-level neural model of chord perception.

    ┌─────────────────────────────────────────────────────────────────────┐
    │  Level 1 — Peripheral (cochlea + IHC)                              │
    │   • Spectral roughness via Sethares + IHC compression               │
    │   • Critical-band masking (Moore & Glasberg)                        │
    │   • Basilar-membrane tonotopic distances                            │
    ├─────────────────────────────────────────────────────────────────────┤
    │  Level 2 — Brainstem (inferior colliculus / cochlear nucleus)       │
    │   • Virtual pitch via Terhardt SHS                                  │
    │   • Harmonicity — how closely chord resembles a harmonic series     │
    │   • Auditory-nerve rate-place profile                               │
    ├─────────────────────────────────────────────────────────────────────┤
    │  Level 3 — Cortical (auditory cortex / tonality networks)          │
    │   • Krumhansl–Kessler tonal hierarchy                               │
    │   • Spectral centroid (brightness)                                  │
    │   • Combined perceptual tension (weighted fusion of all levels)     │
    └─────────────────────────────────────────────────────────────────────┘

    Tension weighting (calibrated on common-practice chords):
        T = 0.40 · roughness  +  0.30 · (1−harmonicity)  +  0.30 · (1−tonalStability)
    """
    if not pitch_classes:
        return {}

    pcs = sorted(set(p % edo for p in pitch_classes))

    # ── Level 1: Peripheral ─────────────────────────────────────────────
    roughness_data = chord_roughness_psychoacoustic(
        pcs, c4_hz, octave, n_harmonics, edo=edo, level_db=level_db)
    masking_data   = chord_masking_analysis(pcs, c4_hz, level_db, edo=edo)

    # ── Level 2: Brainstem ──────────────────────────────────────────────
    vp_data = virtual_pitch_strength(pcs, c4_hz, octave, n_harmonics, edo=edo)

    # Build partial list for AN profile
    note_freqs = roughness_data["note_frequencies_hz"]
    all_f, all_a = [], []
    for f0 in note_freqs:
        fs, as_ = _build_harmonic_series(f0, n_harmonics)
        all_f.extend(fs); all_a.extend(as_)

    an_profile = an_rate_place_profile(all_f, all_a, level_db)

    # ── Level 3: Cortical ───────────────────────────────────────────────
    kk_data = tonal_hierarchy_score(pcs, key, edo=edo)

    # Spectral centroid (brightness indicator)
    centroid_hz = sum(note_freqs) / len(note_freqs) if note_freqs else c4_hz
    centroid_pc = round(edo * math.log2(centroid_hz / c4_hz)) % edo

    # ── Combined perceptual tension ─────────────────────────────────────
    # harmonicity from virtual_pitch_strength is already calibrated (0-1 scale,
    # where ~1.0 = perfect unison, ~0.6 = pure fifth, ~0.45 = major triad).
    t_rough  = roughness_data["roughness_normalized"]
    t_harm   = 1.0 - vp_data["harmonicity"]
    t_tonal  = 1.0 - kk_data["kk_normalized"]
    p_tension = 0.40 * t_rough + 0.30 * t_harm + 0.30 * t_tonal
    p_tension = round(min(1.0, max(0.0, p_tension)), 4)

    tension_label = ("at rest" if p_tension < 0.15 else
                     "mild"    if p_tension < 0.32 else
                     "tension" if p_tension < 0.52 else
                     "strong"  if p_tension < 0.70 else "maximal")

    roughness_label = ("smooth"   if t_rough < 0.15 else
                       "low"      if t_rough < 0.35 else
                       "moderate" if t_rough < 0.55 else
                       "high"     if t_rough < 0.75 else "harsh")

    # AN profile summary (peak activation CF and rate)
    peak_an = max(an_profile, key=lambda x: x["firing_rate"], default=None)

    return {
        "level_1_peripheral": {
            "roughness_normalized" : roughness_data["roughness_normalized"],
            "roughness_label"      : roughness_label,
            "consonance_score"     : roughness_data["consonance_score"],
            "roughness_per_pair"   : roughness_data["roughness_per_pair"],
            "masked_tones"         : masking_data["masked_tones"],
            "masked_tone_names"    : masking_data["masked_tone_names"],
            "cb_overlaps"          : masking_data["cb_overlaps"],
            "note_frequencies_hz"  : note_freqs,
        },
        "level_2_brainstem": {
            "virtual_pitch_hz"     : vp_data["virtual_pitch_hz"],
            "virtual_pitch_name"   : vp_data["virtual_pitch_name"],
            "harmonicity"          : vp_data["harmonicity"],
            "harmonicity_label"    : vp_data["harmonicity_label"],
            "vp_top_candidates"    : vp_data["top_candidates"],
            "an_peak_cf_hz"        : peak_an["cf_hz"]       if peak_an else None,
            "an_peak_firing_rate"  : peak_an["firing_rate"] if peak_an else None,
        },
        "level_3_cortical": {
            "kk_tonal_stability"   : kk_data["tonal_stability"],
            "kk_stability_label"   : kk_data["tonal_stability_label"],
            "kk_raw_score"         : kk_data["kk_raw"],
            "spectral_centroid_hz" : round(centroid_hz, 2),
            "spectral_centroid_pc" : centroid_pc,
            "spectral_centroid_name": nn(centroid_pc, edo=edo),
        },
        "perceptual_tension"      : p_tension,
        "perceptual_tension_label": tension_label,
        "tension_breakdown"       : {
            "roughness_component" : round(t_rough,  4),
            "harmonicity_component": round(t_harm,   4),
            "tonal_component"     : round(t_tonal,  4),
            "weights"             : {"roughness": 0.40, "harmonicity": 0.30, "tonal": 0.30},
        },
    }


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 1  —  SET-CLASS THEORY
# ──────────────────────────────────────────────────────────────────────────────

def interval_class(a: int, b: int, edo: int = 12) -> int:
    d = abs(a - b) % edo
    return min(d, edo - d)

def interval_vector(pcs: List[int], edo: int = 12) -> List[int]:
    if edo != 12: return [0]*6 # ICV is 12-tone specific
    icv = [0] * 6
    for i in range(len(pcs)):
        for j in range(i + 1, len(pcs)):
            ic = interval_class(pcs[i], pcs[j], edo)
            if 1 <= ic <= 6:
                icv[ic - 1] += 1
    return icv

def prime_form(pcs: List[int], edo: int = 12) -> Tuple:
    s = sorted(set(p % edo for p in pcs))
    if not s: return ()
    n = len(s)

    def normal_form(rotation):
        base = rotation[0]
        return tuple((p - base) % edo for p in rotation)

    candidates = []
    for i in range(n):
        rot = [s[(i + j) % n] for j in range(n)]
        candidates.append(normal_form(rot))
    inv = sorted([(edo - p) % edo for p in s])
    for i in range(n):
        rot = [inv[(i + j) % n] for j in range(n)]
        candidates.append(normal_form(rot))
    return min(candidates, key=lambda c: (c[-1], c))


FORTE_TABLE: Dict[Tuple, Tuple[str, str]] = {
    (0,1,2):       ("3-1",   "chromatic cluster"),
    (0,1,3):       ("3-2",   "minor + M2"),
    (0,1,4):       ("3-3",   "minor + m3"),
    (0,1,5):       ("3-4",   "minor + P4"),
    (0,1,6):       ("3-5",   "minor + TT"),
    (0,2,4):       ("3-6",   "whole-tone fragment"),
    (0,2,5):       ("3-7",   "quartal"),
    (0,2,6):       ("3-8",   "augmented 2nd"),
    (0,2,7):       ("3-9",   "suspended (sus2/sus4)"),
    (0,3,6):       ("3-10",  "diminished triad"),
    (0,3,7):       ("3-11",  "minor / major triad"),
    (0,4,7):       ("3-11",  "major / minor triad"),
    (0,4,8):       ("3-12",  "augmented triad"),
    (0,2,3,5):     ("4-10",  "minor pentatonic fragment"),
    (0,1,2,4):     ("4-2",   "minor 7 fragment"),
    (0,3,6,9):     ("4-28",  "diminished 7th"),
    (0,2,5,8):     ("4-27",  "half-diminished 7th"),
    (0,3,5,8):     ("4-26",  "minor 7th"),
    (0,1,4,8):     ("4-19",  "major 7th (no 5)"),
    (0,2,4,7):     ("4-22",  "dom 7 / min 7 fragment"),
    (0,2,4,8):     ("4-24",  "augmented + tone"),
    (0,2,6,8):     ("4-25",  "French aug 6th"),
    (0,1,5,8):     ("4-20",  "major 7th"),
    (0,3,4,7):     ("4-17",  "dim + maj chord"),
    (0,1,4,7):     ("4-18",  "minor-major 7th"),
    (0,2,4,6):     ("4-21",  "whole-tone tetrad"),
    (0,2,5,7):     ("4-23",  "quartal tetrad"),
    (0,1,6,7):     ("4-9",   "tritone cluster"),
    (0,1,3,7):     ("4-Z29", "Z-related: [0,1,3,7]"),
    (0,1,4,6):     ("4-Z15", "Z-related: [0,1,4,6]"),
    (0,2,4,7,9):   ("5-35",  "pentatonic scale"),
    (0,1,3,5,7):   ("5-23",  "minor pentatonic inverse"),
    (0,2,4,5,7):   ("5-34",  "diatonic core"),
    (0,1,3,5,8):   ("5-25",  "min6/9 no root"),
    (0,2,3,5,7):   ("5-24",  "diatonic pentatonic variant"),
    (0,2,4,6,8):   ("5-33",  "whole-tone penta"),
    (0,1,2,4,7):   ("5-20",  "hexatonic penta"),
}


def set_class_info(pcs: List[int], edo: int = 12) -> Dict[str, Any]:
    pcs = sorted(set(p % edo for p in pcs))
    if not pcs:
        return {"prime_form": [], "interval_vector": [], "forte": "?",
                "common_name": "empty"}
    if edo != 12:
        return {"prime_form": list(prime_form(pcs, edo)), "interval_vector": [],
                "forte": "?", "common_name": "non-12 EDO", "cardinality": len(pcs)}

    pf = prime_form(pcs, 12)
    icv = interval_vector(pcs, 12)
    forte_name, common = FORTE_TABLE.get(pf, ("?", "unknown"))
    z_related = ("4-Z29" if "Z15" in forte_name else
                 "4-Z15" if "Z29" in forte_name else None)
    inv_pf   = prime_form([(12 - p) % 12 for p in pcs], 12)
    t_sym    = [n for n in range(1, 12)
                if sorted(set((p + n) % 12 for p in pcs)) == pcs]
    IC_NAMES = ["m2/M7","M2/m7","m3/M6","M3/m6","P4/P5","TT"]
    ic_desc  = ", ".join(f"{cnt}x{IC_NAMES[i]}"
                         for i, cnt in enumerate(icv) if cnt > 0)
    return {
        "prime_form": list(pf), "interval_vector": icv,
        "forte": forte_name, "common_name": common,
        "z_related": z_related, "inv_symmetric": (pf == inv_pf),
        "transpositional_symmetry": t_sym,
        "cardinality": len(pcs),
        "ic_description": ic_desc or "no intervals",
    }


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 2  —  NEO-RIEMANNIAN TRANSFORMATIONS
# ──────────────────────────────────────────────────────────────────────────────

def triad_pcs(root: int, quality: str, edo: int = 12) -> List[int]:
    q = quality.lower()
    def d(s12): return round(s12 * edo / 12.0)
    if q in ('maj', 'major', ''):  return [root%edo, (root+d(4))%edo, (root+d(7))%edo]
    if q in ('min', 'minor', 'm'): return [root%edo, (root+d(3))%edo, (root+d(7))%edo]
    if q in ('dim', 'diminished'): return [root%edo, (root+d(3))%edo, (root+d(6))%edo]
    if q in ('aug', 'augmented'):  return [root%edo, (root+d(4))%edo, (root+d(8))%edo]
    return [root%edo, (root+d(4))%edo, (root+d(7))%edo]

def identify_triad(pcs: List[int], edo: int = 12) -> Optional[Tuple[int, str]]:
    s = set(p % edo for p in pcs)
    def d(s12): return round(s12 * edo / 12.0)
    for root in range(edo):
        for q, ivs in [('maj',[4,7]),('min',[3,7]),('dim',[3,6]),('aug',[4,8])]:
            if s == {root%edo, (root+d(ivs[0]))%edo, (root+d(ivs[1]))%edo}:
                return (root, q)
    return None

def plr_P(root, quality, edo=12):
    if quality=='maj': return (root,'min')
    if quality=='min': return (root,'maj')
    return (root, quality)

def plr_L(root, quality, edo=12):
    d4 = round(4 * edo / 12.0)
    if quality=='maj': return ((root+d4)%edo,'min')
    if quality=='min': return ((root-d4)%edo,'maj')
    return (root, quality)

def plr_R(root, quality, edo=12):
    d9 = round(9 * edo / 12.0)
    d3 = round(3 * edo / 12.0)
    if quality=='maj': return ((root+d9)%edo,'min')
    if quality=='min': return ((root+d3)%edo,'maj')
    return (root, quality)

def plr_N(root, quality, edo=12):
    r,q = plr_L(root, quality, edo); r,q = plr_P(r,q, edo); return plr_L(r,q, edo)

def plr_S(root, quality, edo=12):
    r,q = plr_L(root, quality, edo); r,q = plr_P(r,q, edo); return plr_R(r,q, edo)

def plr_H(root, quality, edo=12):
    r,q = plr_L(root, quality, edo); r,q = plr_P(r,q, edo); return plr_L(r,q, edo)

PLR_OPS = {
    'P': (plr_P, "Parallel",           "Change 3rd ±1st. Keeps root+fifth."),
    'L': (plr_L, "Leittonwechsel",     "Leading-tone exchange. Root ±1 semitone."),
    'R': (plr_R, "Relative",           "Relative maj/min. Fifth ±whole-tone."),
    'N': (plr_N, "Nebenverwandt",      "LPL: fifth of major <-> root of minor."),
    'S': (plr_S, "Slide",              "LPR: shares third only; root+fifth both move."),
    'H': (plr_H, "Hexatonic Pole",     "LPL: maximum Tonnetz distance (tritone of roots)."),
}


def plr_transform(root: int, quality: str, op: str, edo: int = 12) -> Dict[str, Any]:
    fn, name, desc = PLR_OPS.get(op, (None, op, "unknown"))
    if fn is None: return {"error": f"unknown operation {op}"}
    new_root, new_quality = fn(root, quality, edo)
    pcs_from = triad_pcs(root, quality, edo)
    pcs_to   = triad_pcs(new_root, new_quality, edo)
    common   = sorted(set(pcs_from) & set(pcs_to))
    changed_from = sorted(set(pcs_from) - set(pcs_to))
    changed_to   = sorted(set(pcs_to)   - set(pcs_from))
    motions = []
    for cf in changed_from:
        if changed_to:
            ct   = min(changed_to, key=lambda x: mod_abs_dist(cf, x, modulus=edo))
            diff = mod_signed_dist(cf, ct, modulus=edo)
            motions.append({"from_pc": cf, "to_pc": ct, "semitones": diff,
                            "from_name": nn(cf, edo=edo), "to_name": nn(ct, edo=edo)})
    return {
        "op": op, "op_name": name, "op_description": desc,
        "from": {"root": root, "quality": quality,
                 "label": nn(root, edo=edo) + ("" if quality=="maj" else "m"),
                 "pcs": pcs_from},
        "to":   {"root": new_root, "quality": new_quality,
                 "label": nn(new_root, edo=edo) + ("" if new_quality=="maj" else "m"),
                 "pcs": pcs_to},
        "common_tones": common,
        "common_tone_names": [nn(p, edo=edo) for p in common],
        "voice_motions": motions,
        "vl_cost": sum(abs(m["semitones"]) for m in motions),
        "tonnetz_delta": list(tonnetz_projection_interval(new_root, edo)),
    }

def all_plr_neighbors(root, quality, edo=12):
    return [plr_transform(root, quality, op, edo) for op in ['P','L','R']]

def plr_path(root, quality, target_root, target_quality, max_depth=6, edo=12):
    from collections import deque
    start, goal = (root, quality), (target_root, target_quality)
    if start == goal: return []
    queue, visited = deque([(start, [])]), {start}
    while queue:
        state, path = queue.popleft()
        if len(path) >= max_depth: continue
        r, q = state
        for op in ['P','L','R']:
            nr, nq = PLR_OPS[op][0](r, q, edo)
            ns = (nr, nq)
            if ns == goal: return path + [op]
            if ns not in visited:
                visited.add(ns)
                queue.append((ns, path + [op]))
    return None


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 3  —  FUNCTIONAL ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────

MODES = {
    'Ionian':     [0,2,4,5,7,9,11],
    'Dorian':     [0,2,3,5,7,9,10],
    'Phrygian':   [0,1,3,5,7,8,10],
    'Lydian':     [0,2,4,6,7,9,11],
    'Mixolydian': [0,2,4,5,7,9,10],
    'Aeolian':    [0,2,3,5,7,8,10],
    'Locrian':    [0,1,3,5,6,8,10],
}

def diatonic_set(key: int, mode: str = 'Ionian', edo: int = 12) -> List[int]:
    ivs = MODES.get(mode, MODES['Ionian'])
    return [(key + round(iv * edo / 12.0)) % edo for iv in ivs]

def tritone_of_key(key: int, edo: int = 12) -> Tuple[int, int]:
    lt = round(11 * edo / 12.0)
    sd = round(5 * edo / 12.0)
    return ((key + lt) % edo, (key + sd) % edo)


def functional_analysis(pcs: List[int], key: int, edo: int = 12) -> Dict[str, Any]:
    pcs_set = set(p % edo for p in pcs)
    diatonic = set((key + round(s * edo / 12.0)) % edo
                   for s in [0,2,4,5,7,9,11])
    lt, sd = tritone_of_key(key, edo)
    tonic    = key % edo

    def d(s12): return round(s12 * edo / 12.0)
    dom_root = (key + d(7)) % edo

    in_diatonic    = pcs_set.issubset(diatonic)
    chromatic_tones = sorted(pcs_set - diatonic)
    has_tritone    = (lt in pcs_set and sd in pcs_set)

    degree_labels_12 = {
        0:'1', 1:'b2', 2:'2', 3:'b3', 4:'3', 5:'4',
        6:'#4/b5', 7:'5', 8:'b6', 9:'6', 10:'b7', 11:'7'
    }
    scale_degrees = []
    for pc in sorted(pcs_set):
        deg_edo = (pc - key) % edo
        deg_12 = round(deg_edo * 12.0 / edo)
        scale_degrees.append({
            "pc": pc, "name": nn(pc, edo=edo),
            "degree": deg_edo,
            "degree_label": degree_labels_12.get(deg_12, str(deg_edo))
        })

    if has_tritone:
        function        = "DOMINANT"
        function_reason = (
            f"Contains the diatonic tritone [{nn(lt, edo=edo)},{nn(sd, edo=edo)}] (near-7 and near-4 above {nn(key, edo=edo)}). "
            f"Contrary-motion resolution: {nn(lt, edo=edo)}→{nn(tonic, edo=edo)} "
            f"and {nn(sd, edo=edo)}→{nn((key+d(4))%edo, edo=edo)}.")
        tension_level   = 4 if len(pcs_set) >= 4 else 3
    elif tonic in pcs_set and not has_tritone:
        if (key+d(7))%edo in pcs_set or (key+d(4))%edo in pcs_set:
            function = "TONIC"
            function_reason = (
                f"Contains {nn(tonic, edo=edo)} (1) with consonant support and no tritone. "
                f"Maximally stable: low roughness, tonic groundedness.")
            tension_level = 0
        else:
            function        = "TONIC-WEAK"
            function_reason = f"Contains {nn(tonic, edo=edo)} without tritone, sparse consonant support."
            tension_level   = 1
    elif sd in pcs_set and lt not in pcs_set:
        function        = "SUBDOMINANT"
        function_reason = (
            f"Contains 4th degree {nn(sd, edo=edo)} without leading tone {nn(lt, edo=edo)}. "
            f"Pre-dominant or plagal function.")
        tension_level   = 2
    elif (key+d(9))%edo in pcs_set and not has_tritone:
        function        = "TONIC-SUBST"
        function_reason = f"vi-type: tonic substitute by common-tone relation."
        tension_level   = 1
    elif (key+d(2))%edo in pcs_set and not has_tritone:
        function        = "PREDOMINANT"
        function_reason = f"ii-type: supertonic creates forward motion toward V."
        tension_level   = 2
    else:
        function        = "MEDIANT"
        function_reason = "Neither tonic stability nor dominant tension — modal colour."
        tension_level   = 1

    tendency_tones = []
    if lt in pcs_set:
        tendency_tones.append({"pc": lt, "name": nn(lt, edo=edo), "role": "leading tone (7)",
            "tendency": f"resolves UP to {nn(tonic, edo=edo)} (1)", "force": "strong"})
    if sd in pcs_set:
        tendency_tones.append({"pc": sd, "name": nn(sd, edo=edo), "role": "subdominant (4)",
            "tendency": f"resolves DOWN to {nn((key+d(4))%edo, edo=edo)} (3)",
            "force": "strong" if has_tritone else "moderate"})

    containing_modes = []
    for mode_name, mode_ivs in MODES.items():
        for k2 in range(edo):
            mode_set = set((k2 + round(iv * edo / 12.0)) % edo for iv in mode_ivs)
            if pcs_set and pcs_set.issubset(mode_set):
                containing_modes.append({"key": nn(k2, edo=edo), "mode": mode_name})
                break
        if len(containing_modes) >= 3: break

    return {
        "function": function, "function_reason": function_reason,
        "tension_level": tension_level, "has_tritone": has_tritone,
        "tritone": [lt, sd] if has_tritone else [],
        "tritone_names": [nn(lt), nn(sd)] if has_tritone else [],
        "in_diatonic": in_diatonic, "chromatic_tones": chromatic_tones,
        "scale_degrees": scale_degrees, "tendency_tones": tendency_tones,
        "containing_modes": containing_modes,
    }


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 4  —  VOICE-LEADING ORBIFOLD  (Hungarian / Munkres)
# ──────────────────────────────────────────────────────────────────────────────

def linear_sum_assignment(cost_matrix: List[List[float]]) -> Tuple[List[int], List[int]]:
    """Pure-Python O(n³) Hungarian algorithm."""
    n = len(cost_matrix)
    if n == 0: return [], []
    cost = [list(map(float, row)) for row in cost_matrix]
    u, v, p, way = ([0.0]*(n+1), [0.0]*(n+1), [0]*(n+1), [0]*(n+1))
    for i in range(1, n+1):
        p[0] = i
        minv  = [float('inf')] * (n+1)
        used  = [False] * (n+1)
        j0    = 0
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], float('inf'), 0
            for j in range(1, n+1):
                if not used[j]:
                    cur = cost[i0-1][j-1] - u[i0] - v[j]
                    if cur < minv[j]: minv[j] = cur; way[j] = j0
                    if minv[j] < delta: delta = minv[j]; j1 = j
            for j in range(n+1):
                if used[j]: u[p[j]] += delta; v[j] -= delta
                else:        minv[j] -= delta
            j0 = j1
            if p[j0] == 0: break
        while True:
            j1 = way[j0]; p[j0] = p[j1]; j0 = j1
            if j0 == 0: break
    row_ind = [-1] * n
    for j in range(1, n+1):
        if p[j] != 0: row_ind[p[j]-1] = j-1
    return list(range(n)), row_ind


def orbifold_voice_leading(chord_a: List[int], chord_b: List[int], edo: int = 12) -> Dict[str, Any]:
    """Optimal voice-leading distance in T^n/S_n (Hungarian exact assignment)."""
    a = [p % edo for p in chord_a] if chord_a else []
    b = [p % edo for p in chord_b] if chord_b else []
    n = max(len(a), len(b))
    if not a or not b:
        return {"distance": 0, "motions": [], "smoothness": "rest",
                "contrary_motions": 0, "parallel_motions": 0, "oblique_motions": 0,
                "description": "no motion"}
    if len(a) < n: a = list(a) + [a[-1]] * (n - len(a))
    if len(b) < n: b = list(b) + [b[-1]] * (n - len(b))
    cost = [[mod_abs_dist(ai, bj, modulus=edo) for bj in b] for ai in a]
    try:
        rows, cols = linear_sum_assignment(cost)
        pairs = [(i, cols[i]) for i in range(n)]
    except Exception:
        assigned_b = set()
        pairs = []
        for i in range(n):
            candidates = sorted(((cost[i][j],j) for j in range(n)
                                 if j not in assigned_b))
            if candidates: _, j = candidates[0]; assigned_b.add(j); pairs.append((i,j))
    match = [(a[i], b[j], mod_signed_dist(a[i], b[j], modulus=edo)) for i,j in pairs]
    motions = [{"from_pc": m[0], "from_name": nn(m[0], edo=edo),
                "to_pc":   m[1], "to_name":   nn(m[1], edo=edo),
                "semitones": m[2],
                "motion_type": ("oblique" if m[2]==0 else
                                "step"    if abs(m[2])<=2 else
                                "leap"    if abs(m[2])<=5 else "distant")}
               for m in match]
    n_contrary = sum(1 for i in range(len(motions))
                     for j in range(i+1,len(motions))
                     if motions[i]["semitones"] * motions[j]["semitones"] < 0)
    n_parallel = sum(1 for i in range(len(motions))
                     for j in range(i+1,len(motions))
                     if motions[i]["semitones"] == motions[j]["semitones"] != 0)
    n_oblique  = sum(1 for m in motions if m["semitones"] == 0)
    total_dist = sum(abs(m["semitones"]) for m in motions)
    smoothness = "smooth" if total_dist<=3 else "moderate" if total_dist<=6 else "disjunct"
    return {"distance": total_dist, "motions": motions,
            "contrary_motions": n_contrary, "parallel_motions": n_parallel,
            "oblique_motions": n_oblique, "smoothness": smoothness,
            "description": f"Min. {total_dist} semitone{'s' if total_dist!=1 else ''} — {smoothness}"}


def resolution_paths(root: int, quality: str, key: int, edo: int=12) -> List[Dict[str,Any]]:
    pcs  = triad_pcs(root, quality, edo=edo)
    func = functional_analysis(pcs, key, edo)
    results: List[Dict[str,Any]] = []

    if func["has_tritone"]:
        lt, sd = func["tritone"]
        d4 = round(4 * edo / 12.0)
        for tgt_root, tgt_q, mode_label in [
                (key,'maj',"tonic major"), (key,'min',"tonic minor")]:
            tgt_pcs = triad_pcs(tgt_root, tgt_q, edo=edo)
            vl = orbifold_voice_leading(pcs, tgt_pcs, edo=edo)
            results.append({
                "target_pcs": tgt_pcs,
                "target_label": nn(tgt_root, edo=edo)+("" if tgt_q=="maj" else "m"),
                "target_root": tgt_root, "target_quality": tgt_q,
                "rule": "TRITONE_RESOLUTION", "rule_class": "necessity",
                "explanation": (
                    f"Tritone [{nn(lt, edo=edo)},{nn(sd, edo=edo)}] resolves inward: "
                    f"{nn(lt, edo=edo)}→{nn(key, edo=edo)} (+1) and {nn(sd, edo=edo)}→{nn((key+d4)%edo, edo=edo)} (−1)."),
                "voice_leading": vl, "priority": 1})
        vi_root  = (key + round(9 * edo / 12.0)) % edo
        vi_pcs   = triad_pcs(vi_root, 'min', edo=edo)
        results.append({
            "target_pcs": vi_pcs, "target_label": nn(vi_root, edo=edo)+"m",
            "target_root": vi_root, "target_quality": "min",
            "rule": "DECEPTIVE_CADENCE", "rule_class": "necessity",
            "explanation": f"Deceptive cadence → vi: {nn(lt, edo=edo)}→{nn(key, edo=edo)} but bass steps up.",
            "voice_leading": orbifold_voice_leading(pcs, vi_pcs, edo=edo), "priority": 2})

    for op in ['P','L','R']:
        t = plr_transform(root, quality, op, edo)
        if "error" in t: continue
        tgt_root, tgt_q = t["to"]["root"], t["to"]["quality"]
        tgt_pcs = triad_pcs(tgt_root, tgt_q, edo=edo)
        common  = sorted(set(pcs) & set(tgt_pcs))
        results.append({
            "target_pcs": tgt_pcs, "target_label": nn(tgt_root, edo=edo) + ("" if tgt_q=="maj" else "m"),
            "target_root": tgt_root, "target_quality": tgt_q,
            "rule": f"PLR_{op}", "rule_class": "minimal_motion",
            "explanation": (
                f"{op} ({t['op_name']}): {t['op_description']} "
                f"Common tones retained: {[nn(p, edo=edo) for p in common]}."),
            "voice_leading": orbifold_voice_leading(pcs, tgt_pcs, edo=edo), "priority": 3})

    diat = diatonic_set(key, edo=edo)
    root_pc = root % edo
    if root_pc in diat:
        root_idx = diat.index(root_pc)
        for step_dir, dirn in [(-1,"descending"),(+1,"ascending")]:
            ni = (root_idx + step_dir) % 7
            nr = diat[ni]
            third, fifth = diat[(ni+2)%7], diat[(ni+4)%7]
            t_iv = (third-nr)%edo; f_iv = (fifth-nr)%edo

            # Map intervals back to 12-tone qualities
            t_12 = round(t_iv * 12.0 / edo)
            f_12 = round(f_iv * 12.0 / edo)

            nq = ('maj' if t_12==4 and f_12==7 else
                  'min' if t_12==3 and f_12==7 else
                  'dim' if t_12==3 and f_12==6 else 'maj')
            tgt_pcs = triad_pcs(nr, nq, edo=edo)
            results.append({
                "target_pcs": tgt_pcs,
                "target_label": nn(nr, edo=edo)+("" if nq=="maj" else "m"),
                "target_root": nr, "target_quality": nq,
                "rule": f"DIATONIC_STEP_{dirn.upper()}", "rule_class": "scalar",
                "explanation": f"Scale-step {dirn}: {nn(root_pc, edo=edo)}→{nn(nr, edo=edo)} in {nn(key, edo=edo)} major.",
                "voice_leading": orbifold_voice_leading(pcs, tgt_pcs, edo=edo), "priority": 4})

    seen, unique = set(), []
    for r in sorted(results, key=lambda x: (x["priority"], x["voice_leading"]["distance"])):
        k = tuple(sorted(r["target_pcs"]))
        if k not in seen: seen.add(k); unique.append(r)
    return unique[:8]


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 5  —  TONNETZ TENSION  (geometric + psychoacoustic hybrid)
# ──────────────────────────────────────────────────────────────────────────────

def tonnetz_projection_interval(pc: int, edo: int = 12) -> Tuple[int, int]:
    """
    Map a pitch class to (fifths, thirds) Tonnetz lattice coordinates.
    Searches for (a,b) minimising |a·steps5 + b·steps3 − pc| (mod edo)
    with a Manhattan-distance tiebreaker.
    """
    steps_fifth = round(edo * math.log2(3 / 2))
    steps_third = round(edo * math.log2(5 / 4))
    target = pc % edo
    best, best_score = (0, 0), float('inf')
    for a in range(-edo, edo + 1):
        for b in range(-edo//2, edo//2 + 1):
            residue  = (a * steps_fifth + b * steps_third) % edo
            mismatch = min((residue - target) % edo, (target - residue) % edo)
            score    = (0 if mismatch == 0 else 1000) + abs(a) + abs(b)*1.4
            if score < best_score: best_score = score; best = (a, b)
    return best

def tonnetz_projection(iv: int, edo: int = 12) -> Tuple[int, int]:
    return tonnetz_projection_interval(iv, edo)

def tonic_region_centroid(key: int, edo: int = 12) -> Tuple[float, float]:
    vi_steps  = round(9  * edo / 12)
    iii_steps = round(4  * edo / 12)
    region    = [key%edo, (key+vi_steps)%edo, (key+iii_steps)%edo]
    xs = [tonnetz_projection_interval(r, edo)[0] for r in region]
    ys = [tonnetz_projection_interval(r, edo)[1] for r in region]
    return sum(xs)/3.0, sum(ys)/3.0


def tonnetz_tension(pcs: List[int], key: int, edo: int = 12,
                    include_psychoacoustic: bool = True) -> Dict[str, Any]:
    """
    Hybrid tension: geometric Tonnetz distance (algebraic) blended with
    psychoacoustic roughness and tonal hierarchy (perceptual).

    Geometric component  (weight 0.5): centroid distance from tonic region
    Roughness component  (weight 0.3): Sethares spectral consonance
    KK tonal component   (weight 0.2): Krumhansl–Kessler stability

    When include_psychoacoustic=False, falls back to pure geometric tension
    (legacy behaviour).
    """
    if not pcs:
        return {"tension": 0.0, "tension_label": "rest", "distance": 0.0,
                "chord_centroid": [0.0, 0.0], "tonic_centroid": [0.0, 0.0]}

    # ── Geometric (Tonnetz) ─────────────────────────────────────────────
    cx = cy = 0.0
    for pc in pcs:
        x, y = tonnetz_projection_interval(pc, edo)
        cx += x; cy += y
    cx /= len(pcs); cy /= len(pcs)
    tcx, tcy = tonic_region_centroid(key, edo)
    dist = math.sqrt((cx - tcx)**2 + (cy - tcy)**2)
    geom_norm = min(1.0, dist / 4.0)

    if not include_psychoacoustic:
        norm = geom_norm
        label = ("at rest" if norm<0.15 else "mild" if norm<0.35 else
                 "tension" if norm<0.55 else "strong" if norm<0.75 else "maximal")
        return {"tension": round(norm,3), "tension_label": label,
                "distance": round(dist,3),
                "chord_centroid": [round(cx,2), round(cy,2)],
                "tonic_centroid": [round(tcx,2), round(tcy,2)]}

    # ── Psychoacoustic components ────────────────────────────────────────
    roughness_data = chord_roughness_psychoacoustic(pcs, edo=edo)
    kk_data        = tonal_hierarchy_score(pcs, key, edo=edo)

    roughness_t = roughness_data["roughness_normalized"]
    kk_t        = 1.0 - kk_data["kk_normalized"]

    # Hybrid tension
    tension = 0.50 * geom_norm + 0.30 * roughness_t + 0.20 * kk_t
    tension = round(min(1.0, max(0.0, tension)), 4)

    label = ("at rest" if tension<0.15 else "mild" if tension<0.30 else
             "tension" if tension<0.50 else "strong" if tension<0.70 else "maximal")

    return {
        "tension": tension, "tension_label": label,
        "distance": round(dist, 3),
        "chord_centroid": [round(cx,2), round(cy,2)],
        "tonic_centroid": [round(tcx,2), round(tcy,2)],
        "geometric_component" : round(geom_norm, 4),
        "roughness_component" : round(roughness_t, 4),
        "kk_tonal_component"  : round(kk_t, 4),
        "consonance_score"    : roughness_data["consonance_score"],
        "kk_stability"        : kk_data["tonal_stability"],
    }


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 6  —  PATTERN RECOGNITION
# ──────────────────────────────────────────────────────────────────────────────

def detect_sequence(chord_sequence: List[Dict], edo: int = 12) -> Dict[str, Any]:
    if len(chord_sequence) < 3:
        return {"sequence_type": "none", "reason": "Too short"}
    roots     = [c.get("root", 0) for c in chord_sequence]
    intervals = [(roots[i+1] - roots[i]) % edo for i in range(len(roots)-1)]
    if len(set(intervals)) == 1:
        iv = intervals[0]
        iv12 = round(iv * 12.0 / edo)
        iv_name = {1:"chromatic ascent", 2:"whole-tone ascent",
                   3:"minor-3rd cycle", 4:"major-3rd cycle (hexatonic)",
                   5:"ascending fourths", 7:"descending fourths (circle of 5ths)",
                   9:"minor-3rd descent", 10:"whole-tone descent",
                   11:"chromatic descent"}.get(iv12 if edo==12 else -1, f"constant +{iv} steps")
        period = edo // math.gcd(edo, iv)
        return {"sequence_type": "transposition", "interval": iv,
                "interval_name": iv_name, "period": period,
                "description": f"T_{iv} transposition sequence. Orbit period = {period} steps.",
                "algebraic_explanation": (
                    f"T_{iv} ∈ Z_{edo} has order {period}. "
                    f"The sequence is an orbit of this cyclic group action on pitch-class space.")}
    ops_used = []
    for i in range(len(chord_sequence)-1):
        ra,qa = chord_sequence[i].get("root",0), chord_sequence[i].get("quality","maj")
        rb,qb = chord_sequence[i+1].get("root",0), chord_sequence[i+1].get("quality","maj")
        found = False
        for op in ['P','L','R']:
            nr,nq = PLR_OPS[op][0](ra,qa)
            if nr==rb and nq==qb: ops_used.append(op); found=True; break
        if not found: ops_used.append('?')
    if '?' not in ops_used and len(set(ops_used))==1:
        op = ops_used[0]
        return {"sequence_type": f"PLR_{op}_sequence", "operation": op, "period": 2,
                "description": (f"Iterated {op}-transform ({PLR_OPS[op][1]}). "
                                f"{op}² = identity, so this oscillates between two chords.")}
    if sum(1 for iv in intervals if iv==5) >= len(intervals)-1:
        return {"sequence_type": "circle_of_fifths",
                "description": "Descending-fifth root motion — traverses the diatonic circle."}
    return {"sequence_type": "irregular", "operations": ops_used}


def modal_mixture_analysis(chord_pcs: List[int], key: int, edo: int = 12) -> Dict[str, Any]:
    major_pcs  = set(diatonic_set(key, edo=edo))
    chord_set  = set(p % edo for p in chord_pcs)
    chromatic  = chord_set - major_pcs
    borrowed   = [m for m in MODES
                  if m != 'Ionian' and
                  chord_set.issubset(set(diatonic_set(key, m, edo=edo)))]
    pivot_keys = []
    for k2 in range(edo):
        if k2 != key and chord_set.issubset(set(diatonic_set(k2, edo=edo))):
            fa2 = functional_analysis(list(chord_set), k2, edo=edo)
            pivot_keys.append({"key": nn(k2, edo=edo), "function": fa2["function"]})
    return {
        "in_major": not chromatic,
        "chromatic_tones": sorted(chromatic),
        "chromatic_names": [nn(p, edo=edo) for p in sorted(chromatic)],
        "borrowed_from_modes": borrowed, "pivot_keys": pivot_keys[:4],
        "is_modal_mixture": bool(chromatic and borrowed),
    }


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 7  —  EDO / JI LATTICE
# ──────────────────────────────────────────────────────────────────────────────

def ji_ratio(fifths: int, thirds: int, sevenths: int = 0) -> Tuple[int, int, float]:
    log2_val = (fifths * math.log2(3) + thirds * math.log2(5) +
                sevenths * math.log2(7))
    cents = (log2_val % 1.0) * 1200.0
    num = 3**max(0,fifths)  * 5**max(0,thirds)  * 7**max(0,sevenths)
    den = 3**max(0,-fifths) * 5**max(0,-thirds) * 7**max(0,-sevenths)
    while num >= 2*den: den *= 2
    while den > num:   num *= 2
    g = math.gcd(num, den); num //= g; den //= g
    return num, den, cents


def edo_analysis(edo: int) -> Dict[str, Any]:
    step = 1200.0 / edo
    result = {"edo": edo, "step_cents": round(step, 4), "primes": {}, "intervals": []}
    for p in [3,5,7,11,13]:
        ji_c  = math.log2(p) * 1200.0 % 1200.0
        steps = round(ji_c / step)
        err   = steps * step - ji_c
        result["primes"][str(p)] = {
            "ji_cents": round(ji_c,3), "edo_steps": steps,
            "edo_cents": round(steps*step,3), "error_cents": round(err,3)}
    for name, num, den in [
            ("octave",2,1), ("fifth",3,2), ("fourth",4,3),
            ("major_third",5,4), ("minor_third",6,5), ("major_sixth",5,3),
            ("harmonic_7th",7,4), ("11th_harm",11,8)]:
        ji_c  = math.log2(num/den) * 1200.0
        steps = round(ji_c / step)
        err   = steps * step - ji_c
        result["intervals"].append({
            "name": name, "ratio": f"{num}/{den}",
            "ji_cents": round(ji_c,2), "edo_steps": steps,
            "edo_cents": round(steps*step,2), "error_cents": round(err,2)})
    e5 = [abs(result["primes"][str(p)]["error_cents"]) for p in [3,5]]
    e7 = [abs(result["primes"][str(p)]["error_cents"]) for p in [3,5,7]]
    result["consonance_score"]      = round(100-min(sum(e5),100),1)
    result["seven_limit_score"]     = round(100-min(sum(e7),100),1)
    return result


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 8  —  VOICE COMPLETION  (psychoacoustic scoring)
# ──────────────────────────────────────────────────────────────────────────────

def suggest_completion(pitch_classes: List[int], key: int,
                       edo: int = 12,
                       c4_hz: float = 261.625565) -> List[Dict[str, Any]]:
    """
    Suggest pitch classes to add to the existing chord, scored using the full
    psychoacoustic model (Sethares roughness + Terhardt virtual pitch +
    KK tonal stability) instead of the primitive interval-roughness table.

    Scoring (additive, all normalized 0–1 higher = better):
        40% — spectral consonance (1 − roughness_normalized)
        25% — tonal stability (KK profile for key)
        20% — harmonicity / virtual pitch strength
        15% — voice-leading smoothness to nearest resolution target

    Structural reasons are preserved for transparency.
    """
    existing = set(p % edo for p in pitch_classes)
    diat = set(diatonic_set(key, edo=edo))

    # JI reference cents: for each degree above the key, report the cents
    # of the nearest standard JI interval so the caller can gauge the EDO's
    # fidelity.  Works for any EDO by mapping the step to cents directly.
    _JI_TARGETS = [
        (0.0,), (111.7, 100.0), (203.9, 200.0), (315.6, 300.0),
        (386.3, 400.0), (498.0, 500.0), (590.2, 600.0), (702.0, 700.0),
        (813.7, 800.0), (884.4, 900.0), (969.0, 1000.0), (1088.3, 1100.0),
    ]
    def _ji_cents_for_step(step: int) -> float:
        """Return nearest JI reference cents for an EDO step above the key."""
        edo_cents = (step % edo) * 1200.0 / edo
        idx12 = round(edo_cents / 100) % 12
        targets = _JI_TARGETS[idx12]
        return min(targets, key=lambda c: abs(c - edo_cents))

    suggestions = []
    for pc in range(edo):
        if pc in existing:
            continue
        test_pcs = sorted(existing | {pc})

        # ── Psychoacoustic scoring ────────────────────────────────────
        roughness_data = chord_roughness_psychoacoustic(test_pcs, c4_hz, edo=edo)
        kk_data        = tonal_hierarchy_score(test_pcs, key, edo=edo)
        vp_data        = virtual_pitch_strength(test_pcs, c4_hz, edo=edo)

        consonance   = roughness_data["consonance_score"]
        kk_stability = kk_data["kk_normalized"]
        harmonicity  = vp_data["harmonicity"]

        # Voice-leading smoothness component (0–1): prefer step motion
        vl_cost    = sum(mod_abs_dist(pc, e, modulus=edo) for e in existing)
        vl_smooth  = max(0.0, 1.0 - vl_cost / (6.0 * max(1, len(existing))))

        score = (0.40 * consonance +
                 0.25 * kk_stability +
                 0.20 * harmonicity  +
                 0.15 * vl_smooth)

        # ── Structural reasons ────────────────────────────────────────
        triad_id = identify_triad(test_pcs, edo=edo)
        fa       = functional_analysis(test_pcs, key, edo)
        tension  = tonnetz_tension(test_pcs, key, edo)

        reasons = []
        if triad_id:
            r, q = triad_id
            reasons.append(f"Completes {nn(r, edo=edo)} {'maj' if q=='maj' else q} triad")
        lt, sd = tritone_of_key(key, edo)
        if pc == lt:
            reasons.append(f"Adds leading tone {nn(lt, edo=edo)} → dominant function")
        if ((pc==lt and sd in existing) or (pc==sd and lt in existing)):
            reasons.append("Creates diatonic tritone → dominant tension")
        if not existing:
            reasons.append("Provides root")
        if pc not in diat:
            # modal mixture needs edo support if used here
            pass
        if edo == 12:
            sc = set_class_info(test_pcs, 12)
            if sc.get("forte","?") != "?":
                reasons.append(f"Forms {sc['forte']} ({sc['common_name']})")
        if vp_data["virtual_pitch_name"] != "?":
            reasons.append(
                f"Virtual root: {vp_data['virtual_pitch_name']} "
                f"(harmonicity {vp_data['harmonicity_label']})")

        suggestions.append({
            "pc": pc, "name": nn(pc, edo=edo),
            "score": round(score, 4),
            # Psychoacoustic breakdown
            "consonance_score"   : round(consonance, 4),
            "roughness_normalized": roughness_data["roughness_normalized"],
            "kk_stability"       : round(kk_stability, 4),
            "harmonicity"        : round(harmonicity, 4),
            "vl_smoothness"      : round(vl_smooth, 4),
            # Context
            "cents_ji"           : round(_ji_cents_for_step((pc - key) % edo), 1),
            "in_key"             : pc in diat,
            "function_if_added"  : fa["function"],
            "tension_if_added"   : tension["tension"],
            "virtual_pitch_if_added": vp_data["virtual_pitch_name"],
            "structural_reasons" : reasons,
            "completes_triad"    : triad_id is not None,
        })

    suggestions.sort(key=lambda x: -x["score"])
    return suggestions[:8]


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 9  —  PIVOT CHORD SEARCH  (modulation theory)
# ──────────────────────────────────────────────────────────────────────────────

# ── 9.1  Roman numeral system ─────────────────────────────────────────────────

_DIAT_SEMITONES = [0, 2, 4, 5, 7, 9, 11]
_DIAT_DEGREE    = {s: i for i, s in enumerate(_DIAT_SEMITONES)}
_ROMAN_UP       = ['I',  'II',  'III',  'IV',  'V',  'VI',  'VII']
_ROMAN_LO       = ['i',  'ii',  'iii',  'iv',  'v',  'vi',  'vii']

_CHROM_MAP: Dict[int, Tuple[str, int, bool]] = {
    1:  ('b', 1, True),
    3:  ('b', 2, True),
    6:  ('#', 3, False),
    8:  ('b', 5, True),
    10: ('b', 6, True),
}


def roman_numeral(semitones_from_key: int, quality: str, edo: int = 12) -> str:
    """
    Convert (interval from tonic in EDO steps, chord quality) to Roman numeral.
    For non-12 EDO the step count is first mapped to the nearest 12-tone
    semitone so that the standard diatonic/chromatic degree table still applies.
    """
    # Map EDO interval to nearest 12-tone semitone for the lookup
    s = round(semitones_from_key * 12.0 / edo) % 12
    q = quality.lower()
    suffix = {'dim': 'o', 'aug': '+', 'maj7': 'M7', 'min7': 'm7',
              'dom7': '7', 'hdim7': 'o7', 'dim7': 'o7'}.get(q, '')

    if s in _DIAT_DEGREE:
        idx  = _DIAT_DEGREE[s]
        base = _ROMAN_UP[idx] if q in ('maj','aug','maj7','dom7','') else _ROMAN_LO[idx]
        return base + suffix
    elif s in _CHROM_MAP:
        prefix, idx, is_flat = _CHROM_MAP[s]
        if s == 6 and q in ('dim','dim7'):
            prefix = '#'; idx = 3
        elif s == 6:
            prefix = 'b'; idx = 4
        base = _ROMAN_UP[idx] if q in ('maj','aug','maj7','dom7','') else _ROMAN_LO[idx]
        return prefix + base + suffix
    else:
        return f"({s}){suffix}"


def roman_label(root_pc: int, quality: str, key: int, edo: int = 12) -> str:
    """Roman numeral of a chord (root_pc, quality) relative to key."""
    return roman_numeral((root_pc - key) % edo, quality, edo)


# ── 9.2  Diatonic chord inventory ─────────────────────────────────────────────

_SEVENTH_TYPES: Dict[Tuple[int,int,int], str] = {
    (4, 7, 11): 'maj7',
    (4, 7, 10): 'dom7',
    (3, 7, 10): 'min7',
    (3, 6, 10): 'hdim7',
    (3, 6,  9): 'dim7',
    (4, 8, 11): 'augmaj7',
    (4, 8, 10): 'aug7',
}


def _chord_quality_from_intervals(third_iv: int, fifth_iv: int) -> str:
    if third_iv == 4 and fifth_iv == 7: return 'maj'
    if third_iv == 3 and fifth_iv == 7: return 'min'
    if third_iv == 3 and fifth_iv == 6: return 'dim'
    if third_iv == 4 and fifth_iv == 8: return 'aug'
    return 'other'


def diatonic_chord_inventory(
        key: int,
        mode: str = 'Ionian',
        include_sevenths: bool = True,
        include_borrowed: bool = False,
        edo: int = 12,
) -> List[Dict[str, Any]]:
    """
    Build the full diatonic chord inventory for a key.
    """
    scale = diatonic_set(key, mode, edo)
    inventory: List[Dict[str, Any]] = []

    def d(s12): return round(s12 * edo / 12.0)

    for i in range(7):
        root    = scale[i]
        third   = scale[(i + 2) % 7]
        fifth   = scale[(i + 4) % 7]
        seventh = scale[(i + 6) % 7] if include_sevenths else None

        third_iv = (third - root) % edo
        fifth_iv = (fifth - root) % edo

        # Map back to 12 for quality guess
        t12 = round(third_iv * 12.0 / edo)
        f12 = round(fifth_iv * 12.0 / edo)
        q   = _chord_quality_from_intervals(t12, f12)
        pcs = sorted({root, third, fifth})

        seventh_type = None
        seventh_pcs  = pcs
        if include_sevenths and seventh is not None:
            sev_iv       = (seventh - root) % edo
            s12 = round(sev_iv * 12.0 / edo)
            seventh_type = _SEVENTH_TYPES.get((t12, f12, s12))
            seventh_pcs  = sorted({root, third, fifth, seventh})

        deg_semi = (root - key) % edo
        rom      = roman_numeral(deg_semi, q, edo)
        func     = functional_analysis(pcs, key, edo)["function"]
        kk       = tonal_hierarchy_score(pcs, key, edo=edo)["tonal_stability"]
        ten      = tonnetz_tension(pcs, key, edo, include_psychoacoustic=False)["tension"]

        entry: Dict[str, Any] = {
            "root"            : root,
            "quality"         : q,
            "pcs"             : pcs,
            "label"           : nn(root, edo=edo) + ('' if q=='maj' else
                                            'm'  if q=='min' else
                                            'o'  if q=='dim' else '+'),
            "roman"           : rom,
            "degree_semitones": deg_semi,
            "scale_degree"    : i + 1,
            "function"        : func,
            "kk_stability"    : round(kk, 4),
            "tension"         : round(ten, 4),
            "source"          : "diatonic",
        }
        if include_sevenths and seventh_type:
            entry["seventh_type"] = seventh_type
            entry["pcs_with_7"]   = seventh_pcs
            entry["roman_7"]      = roman_numeral(deg_semi, seventh_type, edo)

        inventory.append(entry)

    if include_borrowed:
        borrowed_defs = [
            (3,  'maj', '',  'bIII  - borrowed from Aeolian / parallel minor'),
            (8,  'maj', '',  'bVI   - borrowed from Aeolian / parallel minor'),
            (10, 'maj', '',  'bVII  - borrowed from Mixolydian / parallel minor'),
            (5,  'min', 'm', 'iv    - minor subdominant (borrowed)'),
            (2,  'dim', 'o', 'iio   - borrowed from Aeolian'),
            (1,  'maj', '',  'bII   - Neapolitan (borrowed)'),
        ]
        existing_sets = {frozenset(c['pcs']) for c in inventory}
        for deg12, q, suffix, src in borrowed_defs:
            root   = (key + d(deg12)) % edo
            pcs    = triad_pcs(root, q, edo)
            pcs_fs = frozenset(pcs)
            if pcs_fs in existing_sets:
                continue
            deg_semi = (root - key) % edo
            rom      = roman_numeral(deg_semi, q, edo)
            func     = functional_analysis(pcs, key, edo)["function"]
            kk       = tonal_hierarchy_score(pcs, key, edo=edo)["tonal_stability"]
            ten      = tonnetz_tension(pcs, key, edo, include_psychoacoustic=False)["tension"]
            inventory.append({
                "root"            : root,
                "quality"         : q,
                "pcs"             : pcs,
                "label"           : nn(root, edo=edo) + suffix,
                "roman"           : rom,
                "degree_semitones": deg_semi,
                "scale_degree"    : None,
                "function"        : func,
                "kk_stability"    : round(kk, 4),
                "tension"         : round(ten, 4),
                "source"          : "borrowed",
                "source_desc"     : src,
            })
    return inventory


def _inventory_map(key: int, include_borrowed: bool = False, edo: int = 12) -> Dict[FrozenSet, Dict]:
    """FrozenSet[pcs] -> chord-info dict for fast pivot lookup."""
    return {frozenset(c['pcs']): c
            for c in diatonic_chord_inventory(key, include_borrowed=include_borrowed, edo=edo)}


# ── 9.3  Pivot scoring ────────────────────────────────────────────────────────

def _pivot_score(common_tones: int, kk_from: float, kk_to: float,
                 vl_to_v: int, tension_from: float, tension_to: float) -> float:
    """
    Composite pivot quality score (0-1).

    Weights:
        40% - common-tone count (normalised to [0,1])
        20% - KK stability in key_from (smooth departure)
        15% - KK stability in key_to   (establishes new key)
        15% - VL distance to V of key_to  (lower = better)
        10% - tension balance (similar tension in both contexts = smoother pivot)
    """
    ct_norm   = min(1.0, common_tones / 3.0)
    vl_norm   = max(0.0, 1.0 - vl_to_v / 8.0)
    t_balance = 1.0 - abs(tension_from - tension_to)
    return round(
        0.40 * ct_norm   +
        0.20 * kk_from   +
        0.15 * kk_to     +
        0.15 * vl_norm   +
        0.10 * t_balance,
        4)


# ── 9.4  Common-chord pivots ──────────────────────────────────────────────────

def common_chord_pivots(
        key_from        : int,
        key_to          : int,
        include_borrowed: bool = True,
        edo: int = 12,
) -> List[Dict[str, Any]]:
    """
    Find all chords whose pitch-class set is diatonic in *both* key_from
    and key_to (optionally including borrowed chords from parallel minor).

    Each result:
        pcs, label, roman_from, roman_to, function_from, function_to,
        kk_stability_from/to, tension_from/to,
        vl_to_dominant_of_target, pivot_score, interpretation
    """
    inv_from    = _inventory_map(key_from, include_borrowed, edo)
    inv_to      = _inventory_map(key_to,   include_borrowed, edo)
    d7 = round(7 * edo / 12.0)
    v_of_to_pcs = triad_pcs((key_to + d7) % edo, 'maj', edo)
    results     = []

    for pcs_fs in inv_from:
        if pcs_fs not in inv_to:
            continue
        cf       = inv_from[pcs_fs]
        ct       = inv_to  [pcs_fs]
        pcs_list = sorted(pcs_fs)
        vl_to_v  = orbifold_voice_leading(pcs_list, v_of_to_pcs, edo)["distance"]
        score    = _pivot_score(len(pcs_list), cf["kk_stability"], ct["kk_stability"],
                                vl_to_v, cf["tension"], ct["tension"])
        results.append({
            "type"                     : "common_chord",
            "pcs"                      : pcs_list,
            "label"                    : cf["label"],
            "root"                     : cf["root"],
            "quality"                  : cf["quality"],
            "roman_from"               : cf["roman"],
            "roman_to"                 : ct["roman"],
            "function_from"            : cf["function"],
            "function_to"              : ct["function"],
            "source_from"              : cf.get("source","diatonic"),
            "source_to"                : ct.get("source","diatonic"),
            "kk_stability_from"        : cf["kk_stability"],
            "kk_stability_to"          : ct["kk_stability"],
            "tension_from"             : cf["tension"],
            "tension_to"               : ct["tension"],
            "vl_to_dominant_of_target" : vl_to_v,
            "pivot_score"              : score,
            "interpretation"           : (
                f"{cf['roman']} in {nn(key_from, edo=edo)} ({cf['function'].lower()}) "
                f"-> {ct['roman']} in {nn(key_to, edo=edo)} ({ct['function'].lower()})"),
        })

    results.sort(key=lambda x: -x["pivot_score"])
    return results


# ── 9.5  Secondary dominant pivots ────────────────────────────────────────────

def secondary_dominant_pivots(key_from: int, key_to: int, edo: int = 12) -> List[Dict[str, Any]]:
    """
    Find V7/x applied-dominant chords that are functional in key_from and
    resolve naturally to a diatonic chord in key_to.

    The strongest case: V7/x in key_from IS the V7 of key_to -> direct
    dominant approach to the new tonic.

    Works for any EDO by scaling all intervals proportionally from 12-EDO.
    """
    def d(s12): return round(s12 * edo / 12.0)

    scale_from    = diatonic_set(key_from, edo=edo)
    v7_of_to_root = (key_to + d(7)) % edo
    inv_to        = _inventory_map(key_to, include_borrowed=True, edo=edo)
    results       = []

    for i, deg_root in enumerate(scale_from):
        secdom_root = (deg_root + d(7)) % edo
        # Secondary dominant 7th chord intervals (0, M3, P5, m7)
        secdom_pcs  = [(secdom_root + d(iv)) % edo for iv in [0, 4, 7, 10]]
        roman_from  = f"V7/{roman_numeral((deg_root - key_from) % edo, 'maj', edo)}"

        if secdom_root == v7_of_to_root:
            pcs_s    = sorted(set(secdom_pcs))
            vl_to_i  = orbifold_voice_leading(pcs_s, triad_pcs(key_to, 'maj', edo), edo)["distance"]
            kk_from  = tonal_hierarchy_score(pcs_s, key_from, edo=edo)["tonal_stability"]
            kk_to    = tonal_hierarchy_score(pcs_s, key_to,   edo=edo)["tonal_stability"]
            ten_from = tonnetz_tension(pcs_s, key_from, edo, include_psychoacoustic=False)["tension"]
            ten_to   = tonnetz_tension(pcs_s, key_to,   edo, include_psychoacoustic=False)["tension"]
            score    = _pivot_score(2, kk_from, kk_to, vl_to_i, ten_from, ten_to)
            results.append({
                "type"                    : "secondary_dominant",
                "pcs"                     : pcs_s,
                "label"                   : nn(secdom_root, edo=edo) + "7",
                "root"                    : secdom_root,
                "quality"                 : "dom7",
                "roman_from"              : roman_from,
                "roman_to"                : "V7",
                "function_from"           : "DOMINANT (applied)",
                "function_to"             : "DOMINANT",
                "applied_to_scale_degree" : i + 1,
                "applied_to_name"         : nn(deg_root, edo=edo),
                "kk_stability_from"       : round(kk_from, 4),
                "kk_stability_to"         : round(kk_to,   4),
                "vl_to_tonic_of_target"   : vl_to_i,
                "pivot_score"             : score,
                "interpretation"          : (
                    f"{roman_from} in {nn(key_from, edo=edo)} IS V7 of {nn(key_to, edo=edo)}. "
                    f"Strongest applied-dominant pivot -> resolves directly to new tonic."),
            })
        else:
            # Subsets: try full V7 then three omit-one subsets
            secdom_full = sorted(set(secdom_pcs))
            for sub_pcs in [secdom_full,
                            sorted({secdom_pcs[0], secdom_pcs[1], secdom_pcs[2]}),
                            sorted({secdom_pcs[0], secdom_pcs[1], secdom_pcs[3]}),
                            sorted({secdom_pcs[0], secdom_pcs[2], secdom_pcs[3]})]:
                sub_pcs = sorted(set(sub_pcs))
                if frozenset(sub_pcs) in inv_to:
                    ct_info  = inv_to[frozenset(sub_pcs)]
                    kk_from  = tonal_hierarchy_score(secdom_full, key_from, edo=edo)["tonal_stability"]
                    vl_to_v  = orbifold_voice_leading(
                        sub_pcs, triad_pcs((key_to + d(7)) % edo, 'maj', edo), edo)["distance"]
                    ten_from = tonnetz_tension(secdom_full, key_from, edo, include_psychoacoustic=False)["tension"]
                    score    = _pivot_score(len(sub_pcs), kk_from, ct_info["kk_stability"],
                                            vl_to_v, ten_from, ct_info["tension"])
                    results.append({
                        "type"                    : "secondary_dominant",
                        "pcs"                     : sub_pcs,
                        "label"                   : nn(secdom_root, edo=edo) + "7",
                        "root"                    : secdom_root,
                        "quality"                 : "dom7",
                        "roman_from"              : roman_from,
                        "roman_to"                : ct_info["roman"],
                        "function_from"           : "DOMINANT (applied)",
                        "function_to"             : ct_info["function"],
                        "applied_to_scale_degree" : i + 1,
                        "applied_to_name"         : nn(deg_root, edo=edo),
                        "kk_stability_from"       : round(kk_from, 4),
                        "kk_stability_to"         : ct_info["kk_stability"],
                        "vl_to_dominant_of_target": vl_to_v,
                        "pivot_score"             : score,
                        "interpretation"          : (
                            f"{roman_from} in {nn(key_from, edo=edo)} = "
                            f"{ct_info['roman']} in {nn(key_to, edo=edo)}."),
                    })
                    break

    results.sort(key=lambda x: -x["pivot_score"])
    return results


# ── 9.6  Enharmonic pivots ────────────────────────────────────────────────────

def enharmonic_pivots(key_from: int, key_to: int, edo: int = 12) -> List[Dict[str, Any]]:
    """
    Three classic enharmonic pivot families, generalized for any EDO.

    (A) Symmetric 4-note chord (dim7 generalization)
        Requires edo divisible by 4.  Step = edo//4.
        Each note can serve as the leading tone of a different key.

    (B) German Augmented Sixth <-> Dominant Seventh  (EDO-12 only)
        Enharmonic identity between Ger+6 and a dom7 a semitone above.

    (C) Symmetric 3-note chord (augmented-triad generalization)
        Requires edo divisible by 3.  Step = edo//3.
        Each note can be III+ or V+ in three different keys.
    """
    results = []

    def d(s12): return round(s12 * edo / 12.0)
    lt_step = d(1)   # semitone analog in this EDO (leading-tone motion)

    # ── (A) Symmetric 4-note chords ───────────────────────────────────────
    if edo % 4 == 0:
        sym4 = edo // 4
        seen4: Set[FrozenSet] = set()
        for base_root in range(sym4):
            chord_fs = frozenset((base_root + k * sym4) % edo for k in range(4))
            if chord_fs in seen4:
                continue
            seen4.add(chord_fs)
            for lt in sorted(chord_fs):
                tonic_implied = (lt + lt_step) % edo
                is_from = (lt == (key_from + d(11)) % edo)
                is_to   = (tonic_implied == key_to)
                if not (is_from or is_to):
                    continue
                pcs_list = sorted(chord_fs)
                kk_from  = tonal_hierarchy_score(pcs_list, key_from, edo=edo)["kk_normalized"]
                kk_to    = tonal_hierarchy_score(pcs_list, key_to,   edo=edo)["kk_normalized"]
                vl_to_v  = orbifold_voice_leading(
                    pcs_list, triad_pcs((key_to + d(7)) % edo, 'maj', edo), edo)["distance"]
                ten_from = tonnetz_tension(pcs_list, key_from, edo, include_psychoacoustic=False)["tension"]
                ten_to   = tonnetz_tension(pcs_list, key_to,   edo, include_psychoacoustic=False)["tension"]
                score    = _pivot_score(3, kk_from, kk_to, vl_to_v, ten_from, ten_to)
                results.append({
                    "type"                    : "enharmonic_dim7",
                    "pcs"                     : pcs_list,
                    "label"                   : nn(lt, edo=edo) + "o7",
                    "root_as_leading_tone"    : lt,
                    "implied_tonic"           : tonic_implied,
                    "implied_tonic_name"      : nn(tonic_implied, edo=edo),
                    "roman_from"              : f"viio7 in {nn(key_from, edo=edo)} (LT={nn(lt, edo=edo)})",
                    "roman_to"                : f"viio7 in {nn(tonic_implied, edo=edo)} (LT={nn(lt, edo=edo)})",
                    "kk_stability_from"       : round(kk_from, 4),
                    "kk_stability_to"         : round(kk_to,   4),
                    "vl_to_dominant_of_target": vl_to_v,
                    "pivot_score"             : score,
                    "interpretation"          : (
                        f"Sym-4 chord {[nn(p, edo=edo) for p in pcs_list]}: "
                        f"LT {nn(lt, edo=edo)} -> tonic {nn(tonic_implied, edo=edo)}. "
                        f"Any of the {4} symmetry-equivalent notes can be reinterpreted as LT."),
                })

    # ── (B) German Aug6 <-> Dom7 (EDO-12 only) ────────────────────────────
    if edo == 12:
        ger_root  = (key_from + 8) % 12
        ger_pcs   = sorted({ger_root, key_from % 12,
                            (key_from + 3) % 12, (key_from + 6) % 12})
        dom7_key  = (key_from + 1) % 12
        kk_from_g = tonal_hierarchy_score(ger_pcs, key_from, edo=12)["kk_normalized"]
        kk_to_g   = tonal_hierarchy_score(ger_pcs, key_to,   edo=12)["kk_normalized"]
        vl_ger    = orbifold_voice_leading(ger_pcs, triad_pcs((key_to + 7) % 12, 'maj', 12), 12)["distance"]
        ten_from_g = tonnetz_tension(ger_pcs, key_from, 12, include_psychoacoustic=False)["tension"]
        ten_to_g   = tonnetz_tension(ger_pcs, key_to,   12, include_psychoacoustic=False)["tension"]
        score_g    = _pivot_score(2, kk_from_g, kk_to_g, vl_ger, ten_from_g, ten_to_g)
        results.append({
            "type"                    : "enharmonic_german_aug6",
            "pcs"                     : ger_pcs,
            "label"                   : f"Ger+6/{nn(key_from, edo=12)} = {nn(ger_root, edo=12)}7",
            "root"                    : ger_root,
            "roman_from"              : f"Ger+6 (bVI aug6th in {nn(key_from, edo=12)})",
            "roman_to"                : f"V7 in {nn(dom7_key, edo=12)} (enharmonic Dom7)",
            "dom7_resolves_to_key"    : nn(dom7_key, edo=12),
            "kk_stability_from"       : round(kk_from_g, 4),
            "kk_stability_to"         : round(kk_to_g,   4),
            "vl_to_dominant_of_target": vl_ger,
            "pivot_score"             : score_g,
            "interpretation"          : (
                f"Ger+6 in {nn(key_from, edo=12)} = pcs {[nn(p, edo=12) for p in ger_pcs]}. "
                f"Enharmonically = {nn(ger_root, edo=12)}7 (dom7). "
                f"As dom7 -> resolves to {nn(dom7_key, edo=12)} major. "
                f"Classic pivot: {nn(key_from, edo=12)} major -> {nn(dom7_key, edo=12)} major."),
        })

    # ── (C) Symmetric 3-note chords ───────────────────────────────────────
    if edo % 3 == 0:
        sym3 = edo // 3
        seen3: Set[FrozenSet] = set()
        for base_root in range(sym3):
            chord_fs = frozenset((base_root + k * sym3) % edo for k in range(3))
            if chord_fs in seen3:
                continue
            seen3.add(chord_fs)
            for root in sorted(chord_fs):
                for tonic_impl, role in [
                        ((root - d(4)) % edo, 'III+'),
                        ((root - d(7)) % edo, 'V+')]:
                    if tonic_impl not in {key_from, key_to}:
                        continue
                    pcs_list  = sorted(chord_fs)
                    kk_from_a = tonal_hierarchy_score(pcs_list, key_from, edo=edo)["tonal_stability"]
                    kk_to_a   = tonal_hierarchy_score(pcs_list, key_to,   edo=edo)["tonal_stability"]
                    vl_a      = orbifold_voice_leading(
                        pcs_list, triad_pcs((key_to + d(7)) % edo, 'maj', edo), edo)["distance"]
                    ten_from_a = tonnetz_tension(pcs_list, key_from, edo, include_psychoacoustic=False)["tension"]
                    ten_to_a   = tonnetz_tension(pcs_list, key_to,   edo, include_psychoacoustic=False)["tension"]
                    score_a    = _pivot_score(2, kk_from_a, kk_to_a, vl_a, ten_from_a, ten_to_a)
                    results.append({
                        "type"                    : "enharmonic_augmented_triad",
                        "pcs"                     : pcs_list,
                        "label"                   : nn(root, edo=edo) + "+",
                        "root"                    : root,
                        "roman_from"              : roman_label(root, 'aug', key_from, edo),
                        "roman_to"                : roman_label(root, 'aug', key_to,   edo),
                        "role"                    : role,
                        "kk_stability_from"       : round(kk_from_a, 4),
                        "kk_stability_to"         : round(kk_to_a,   4),
                        "vl_to_dominant_of_target": vl_a,
                        "pivot_score"             : score_a,
                        "interpretation"          : (
                            f"Sym-3 chord {[nn(p, edo=edo) for p in pcs_list]} as {role}: "
                            f"enharmonic pivot via equal-M3 division of octave."),
                    })

    results.sort(key=lambda x: -x["pivot_score"])
    return results


# ── 9.7  Borrowed / modal-mixture pivots ──────────────────────────────────────

def borrowed_pivots(key_from: int, key_to: int, edo: int = 12) -> List[Dict[str, Any]]:
    """
    Find chords that are *borrowed* in key_from (from any mode of that tonic)
    but *diatonic* in key_to, or vice versa.

    Searches all 7 modal rotations of key_from for chords not in Ionian(key_from)
    that ARE present in Ionian(key_to).
    """
    inv_to  = _inventory_map(key_to, include_borrowed=False, edo=edo)
    ionian_from = {frozenset(c['pcs']) for c in diatonic_chord_inventory(key_from, edo=edo)}
    results = []

    for mode_name in MODES:
        mode_scale = diatonic_set(key_from, mode_name, edo=edo)
        for i in range(7):
            root  = mode_scale[i]
            third = mode_scale[(i+2) % 7]
            fifth = mode_scale[(i+4) % 7]
            t_iv  = (third - root) % edo
            f_iv  = (fifth - root) % edo

            t12 = round(t_iv * 12.0 / edo)
            f12 = round(f_iv * 12.0 / edo)
            q     = _chord_quality_from_intervals(t12, f12)
            pcs   = sorted({root, third, fifth})
            pcs_fs= frozenset(pcs)

            if pcs_fs in ionian_from:
                continue
            if pcs_fs not in inv_to:
                continue

            ct_info  = inv_to[pcs_fs]
            kk_from  = tonal_hierarchy_score(pcs, key_from, edo=edo)["tonal_stability"]
            kk_to    = ct_info["kk_stability"]
            d7 = round(7 * edo / 12.0)
            vl_to_v  = orbifold_voice_leading(
                pcs, triad_pcs((key_to + d7) % edo, 'maj', edo), edo)["distance"]
            ten_from = tonnetz_tension(pcs, key_from, edo, include_psychoacoustic=False)["tension"]
            ten_to   = ct_info["tension"]
            score    = _pivot_score(len(pcs), kk_from, kk_to, vl_to_v, ten_from, ten_to)
            rom_from = roman_label(root, q, key_from, edo)

            results.append({
                "type"                    : "borrowed_pivot",
                "pcs"                     : pcs,
                "label"                   : nn(root, edo=edo) + ('' if q=='maj' else
                                                        'm' if q=='min' else 'o'),
                "root"                    : root,
                "quality"                 : q,
                "roman_from"              : rom_from,
                "roman_to"                : ct_info["roman"],
                "function_from"           : functional_analysis(pcs, key_from, edo)["function"],
                "function_to"             : ct_info["function"],
                "source_mode"             : mode_name,
                "kk_stability_from"       : round(kk_from, 4),
                "kk_stability_to"         : round(kk_to, 4),
                "vl_to_dominant_of_target": vl_to_v,
                "tension_from"            : round(ten_from, 4),
                "tension_to"              : round(ten_to,   4),
                "pivot_score"             : score,
                "interpretation"          : (
                f"{rom_from} borrowed from {nn(key_from, edo=edo)} {mode_name} "
                f"= {ct_info['roman']} diatonic in {nn(key_to, edo=edo)}."),
            })

    seen: Set[FrozenSet] = set()
    unique = []
    for r in sorted(results, key=lambda x: -x["pivot_score"]):
        k = frozenset(r["pcs"])
        if k not in seen:
            seen.add(k); unique.append(r)
    return unique[:12]


# ── 9.8  Chromatic mediant pivots ─────────────────────────────────────────────

def chromatic_mediant_pivots(key_from: int, key_to: int, edo: int = 12) -> List[Dict[str, Any]]:
    """
    BFS through the PLR Tonnetz graph from the tonic of key_from (depth <= 3).
    Reports any intermediate chord that has a diatonic role in key_to or IS
    the tonic/dominant of key_to.

    Characteristic of late-Romantic chromatic mediant modulations (Schubert,
    Brahms, Wagner) where keys related by M3 or m3 share at most 1 tone but
    connect through smooth PLR voice-leading.
    """
    from collections import deque
    results  = []
    inv_to   = _inventory_map(key_to, include_borrowed=True, edo=edo)
    d7 = round(7 * edo / 12.0)
    v_of_to  = triad_pcs((key_to + d7) % edo, 'maj', edo)

    queue   = deque([(key_from % edo, 'maj', [])])
    visited = {(key_from % edo, 'maj')}

    while queue:
        r, q, ops = queue.popleft()
        if len(ops) >= 3:
            continue
        for op in ['P', 'L', 'R']:
            nr, nq = PLR_OPS[op][0](r, q, edo)
            ns = (nr, nq)
            if ns in visited:
                continue
            visited.add(ns)
            new_ops  = ops + [op]
            pcs      = triad_pcs(nr, nq, edo)
            pcs_fs   = frozenset(pcs)
            role_to  = inv_to.get(pcs_fs)
            is_tonic = sorted(pcs) == sorted(triad_pcs(key_to, 'maj', edo))
            is_dom   = sorted(pcs) == sorted(v_of_to)

            if role_to or is_tonic or is_dom:
                kk_from  = tonal_hierarchy_score(pcs, key_from, edo=edo)["tonal_stability"]
                kk_to    = (role_to["kk_stability"] if role_to else
                            tonal_hierarchy_score(pcs, key_to, edo=edo)["tonal_stability"])
                vl_cost  = orbifold_voice_leading(pcs, v_of_to, edo)["distance"]
                ten_from = tonnetz_tension(pcs, key_from, edo, include_psychoacoustic=False)["tension"]
                ten_to   = tonnetz_tension(pcs, key_to, edo, include_psychoacoustic=False)["tension"]
                tn_kf    = tonnetz_projection_interval(key_from % edo, edo)
                tn_nr    = tonnetz_projection_interval(nr, edo)
                tn_dist  = math.sqrt((tn_kf[0]-tn_nr[0])**2 + (tn_kf[1]-tn_nr[1])**2)
                score    = _pivot_score(3 - len(new_ops), kk_from, kk_to,
                                        vl_cost, ten_from, ten_to)
                roman_to = (role_to["roman"] if role_to else
                            "I" if is_tonic else "V")
                results.append({
                    "type"                    : "chromatic_mediant",
                    "pcs"                     : sorted(pcs),
                    "label"                   : nn(nr, edo=edo) + ("" if nq=="maj" else "m"),
                    "root"                    : nr,
                    "quality"                 : nq,
                    "plr_path"                : new_ops,
                    "plr_path_str"            : "->".join(new_ops),
                    "roman_from"              : roman_label(nr, nq, key_from, edo),
                    "roman_to"                : roman_to,
                    "tonnetz_distance"        : round(tn_dist, 3),
                    "kk_stability_from"       : round(kk_from, 4),
                    "kk_stability_to"         : round(kk_to,   4),
                    "vl_to_dominant_of_target": vl_cost,
                    "tension_from"            : round(ten_from, 4),
                    "tension_to"              : round(ten_to,   4),
                    "pivot_score"             : score,
                    "interpretation"          : (
                        f"PLR {new_ops} from {nn(key_from, edo=edo)} tonic -> "
                        f"{nn(nr, edo=edo)}{'m' if nq=='min' else ''} = {roman_to} "
                        f"in {nn(key_to, edo=edo)} (Tonnetz dist {round(tn_dist,2)})."),
                })
            queue.append((nr, nq, new_ops))

    results.sort(key=lambda x: -x["pivot_score"])
    return results[:8]


# ── 9.9  Key-distance helpers ──────────────────────────────────────────────────

def _extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Extended Euclidean algorithm: returns (gcd, x, y) s.t. a·x + b·y = gcd."""
    if a == 0:
        return b, 0, 1
    g, x, y = _extended_gcd(b % a, a)
    return g, y - (b // a) * x, x

def cof_distance(key_a: int, key_b: int, edo: int = 12) -> int:
    """
    Steps on the circle of fifths between two keys, for any EDO.

    The CoF position of a pitch class p is  (p · inv5) mod edo,
    where inv5 is the modular inverse of steps_fifth = round(edo · log2(3/2)).
    When gcd(steps_fifth, edo) ≠ 1 the fifth doesn't span all edo steps, so
    we fall back to plain chromatic distance.
    """
    steps_fifth = round(edo * math.log2(3.0 / 2.0))
    g, x, _ = _extended_gcd(steps_fifth % edo, edo)
    if g != 1:
        # Fifth sub-group doesn't cover the full octave (e.g. EDO 6, 9, 10...)
        d = abs((key_b - key_a) % edo)
        return min(d, edo - d)
    inv5   = x % edo
    pos_a  = (key_a * inv5) % edo
    pos_b  = (key_b * inv5) % edo
    d      = abs(pos_b - pos_a)
    return min(d, edo - d)

def _cof_relationship(d: int, edo: int = 12) -> str:
    """Describe the circle-of-fifths distance between two keys."""
    if edo == 12:
        return {0: "same key",
                1: "closely related (1 step on CoF)",
                2: "moderately related (2 steps)",
                3: "distantly related (3 steps)",
                4: "remote (4 steps)",
                5: "very remote (5 steps)",
                6: "tritone / maximally remote"}.get(d, "remote")
    # General: express as fraction of the maximum possible CoF distance
    max_d = edo // 2
    if d == 0:
        return "same key"
    frac = d / max(max_d, 1)
    if frac <= 0.15: return f"closely related ({d} CoF steps)"
    if frac <= 0.35: return f"moderately related ({d} CoF steps)"
    if frac <= 0.55: return f"distantly related ({d} CoF steps)"
    if frac <= 0.80: return f"remote ({d} CoF steps)"
    return f"maximally remote ({d} CoF steps)"

def common_tones_between_keys(key_a: int, key_b: int, edo: int = 12) -> List[int]:
    """Pitch classes diatonic in BOTH key_a and key_b (major scales)."""
    return sorted(set(diatonic_set(key_a, edo=edo)) & set(diatonic_set(key_b, edo=edo)))


# ── 9.10  Modulation path search (Dijkstra) ────────────────────────────────────

def modulation_path_search(key_from: int, key_to: int,
                           max_depth: int = 6, edo: int = 12) -> List[Dict[str, Any]]:
    """
    Dijkstra search for the minimum voice-leading cost modulation path:
        I(key_from)  -->  [diatonic chords in key_from]
                     -->  [pivot chord: diatonic in BOTH keys]
                     -->  [diatonic chords in key_to]
                     -->  I(key_to)

    Graph nodes: diatonic triads of key_from + key_to (plus borrowed chords).
    Edge weight: orbifold voice-leading distance between adjacent chords.
    Pivot edge:  zero-cost reinterpretation when a chord is in both inventories.

    The search enforces that the pivot is not taken on the very first step
    (the music must first establish key_from) and the goal is to reach
    the major tonic of key_to after the pivot.

    State: (phase, root, quality)
        phase 0 = still in key_from
        phase 1 = crossed pivot, now in key_to
    """
    import heapq

    inv_from_list = diatonic_chord_inventory(key_from, include_borrowed=True, edo=edo)
    inv_to_list   = diatonic_chord_inventory(key_to,   include_borrowed=True, edo=edo)

    from_map: Dict[Tuple[int,str], Dict] = {
        (c["root"], c["quality"]): c for c in inv_from_list}
    to_map:   Dict[Tuple[int,str], Dict] = {
        (c["root"], c["quality"]): c for c in inv_to_list}

    pivot_nodes: Set[Tuple[int,str]] = set(from_map) & set(to_map)

    goal_ck  = (key_to % edo, 'maj')
    start_ck = (key_from % edo, 'maj')

    # heap: (total_vl_cost, phase, root, quality, path_list)
    start_step = {
        "root": key_from % edo, "quality": 'maj',
        "label": nn(key_from, edo=edo),
        "roman": roman_label(key_from, 'maj', key_from, edo), "key_context": "from",
        "function": "TONIC",
        "vl_from_prev": 0, "cumulative_cost": 0,
    }
    heap: List = [(0, 0, 0, key_from % edo, 'maj', [start_step])]
    _cnt = [1]
    visited: Set[Tuple[int, int, str]] = set()
    best_path: Optional[List[Dict]] = None
    best_cost = float('inf')

    while heap:
        cost, _seq, phase, cur_r, cur_q, path = heapq.heappop(heap)
        state = (phase, cur_r, cur_q)
        if state in visited or cost >= best_cost:
            continue
        visited.add(state)

        cur_pcs = triad_pcs(cur_r, cur_q, edo)
        ck      = (cur_r, cur_q)

        # Check goal
        if phase == 1 and ck == goal_ck:
            if cost < best_cost:
                best_cost = cost; best_path = path
            continue

        if len(path) >= max_depth:
            continue

        if phase == 0:
            # Move within key_from
            for (nr, nq), ninfo in from_map.items():
                if (0, nr, nq) in visited:
                    continue
                vl       = orbifold_voice_leading(cur_pcs, triad_pcs(nr, nq, edo), edo)["distance"]
                new_cost = cost + vl
                if new_cost >= best_cost:
                    continue
                step = {
                    "root": nr, "quality": nq,
                    "label": nn(nr, edo=edo) + ("" if nq=="maj" else "m"),
                    "roman": ninfo["roman"], "key_context": "from",
                    "function": ninfo["function"],
                    "vl_from_prev": vl, "cumulative_cost": new_cost,
                }
                heapq.heappush(heap, (new_cost, _cnt[0], 0, nr, nq, path + [step])); _cnt[0]+=1

            # Take pivot (must have played at least 2 chords in key_from)
            if ck in pivot_nodes and len(path) >= 2:
                pivot_step = dict(path[-1])
                pivot_step["key_context"]       = "pivot"
                pivot_step["roman_as_pivot_from"]= from_map[ck]["roman"]
                pivot_step["roman_as_pivot_to"]  = to_map[ck]["roman"]
                pivot_step["pivot_note"]         = (
                    f"{from_map[ck]['roman']} in {nn(key_from, edo=edo)} "
                    f"= {to_map[ck]['roman']} in {nn(key_to, edo=edo)}")
                heapq.heappush(heap,
                    (cost, _cnt[0], 1, cur_r, cur_q, path[:-1] + [pivot_step])); _cnt[0]+=1

        elif phase == 1:
            # Move within key_to
            for (nr, nq), ninfo in to_map.items():
                if (1, nr, nq) in visited:
                    continue
                vl       = orbifold_voice_leading(cur_pcs, triad_pcs(nr, nq, edo), edo)["distance"]
                new_cost = cost + vl
                if new_cost >= best_cost:
                    continue
                step = {
                    "root": nr, "quality": nq,
                    "label": nn(nr, edo=edo) + ("" if nq=="maj" else "m"),
                    "roman": ninfo["roman"], "key_context": "to",
                    "function": ninfo["function"],
                    "vl_from_prev": vl, "cumulative_cost": new_cost,
                }
                heapq.heappush(heap, (new_cost, _cnt[0], 1, nr, nq, path + [step])); _cnt[0]+=1

    return best_path or []


# ── 9.11  Master pivot search ─────────────────────────────────────────────────

def pivot_search(
        key_from           : int,
        key_to             : int,
        include_borrowed   : bool = True,
        include_secondaries: bool = True,
        include_enharmonic : bool = True,
        include_chromatic  : bool = True,
        max_per_type       : int  = 6,
        edo                : int  = 12,
) -> Dict[str, Any]:
    """
    Master pivot-chord search: runs all pivot-finding strategies and returns
    a unified, de-duplicated, ranked result set.

    Parameters
    ----------
    key_from / key_to  : integer pitch classes (0=C, 1=C#, ..., 11=B)
    include_*          : toggle each pivot family
    max_per_type       : cap results per family before merging

    Returns
    -------
    {
      key_from, key_to, key_from_name, key_to_name,
      cof_distance, cof_relationship,
      common_scale_tones, common_scale_names, n_common_tones,
      diatonic_inventory_from, diatonic_inventory_to,
      pivots_by_type: {
        common_chord: [...],
        secondary_dominant: [...],
        enharmonic: [...],
        borrowed: [...],
        chromatic_mediant: [...]
      },
      all_pivots: [...],   # merged, de-duplicated, ranked by pivot_score
      best_pivot: {...},   # highest scoring
      modulation_path: [...]  # Dijkstra-optimal chord sequence
    }
    """
    cof_d  = cof_distance(key_from, key_to, edo)
    common = common_tones_between_keys(key_from, key_to, edo)

    pivots_by_type: Dict[str, List] = {}

    pivots_by_type["common_chord"] = common_chord_pivots(
        key_from, key_to, include_borrowed=include_borrowed, edo=edo)[:max_per_type]

    if include_secondaries:
        pivots_by_type["secondary_dominant"] = secondary_dominant_pivots(
            key_from, key_to, edo=edo)[:max_per_type]

    if include_enharmonic:
        pivots_by_type["enharmonic"] = enharmonic_pivots(
            key_from, key_to, edo=edo)[:max_per_type]

    if include_borrowed:
        pivots_by_type["borrowed"] = borrowed_pivots(
            key_from, key_to, edo=edo)[:max_per_type]

    if include_chromatic:
        pivots_by_type["chromatic_mediant"] = chromatic_mediant_pivots(
            key_from, key_to, edo=edo)[:max_per_type]

    # Merge and de-duplicate
    seen_pcs: Set[FrozenSet] = set()
    all_pivots: List[Dict]   = []
    for type_pivots in pivots_by_type.values():
        for p in type_pivots:
            k = frozenset(p["pcs"])
            if k not in seen_pcs:
                seen_pcs.add(k); all_pivots.append(p)
    all_pivots.sort(key=lambda x: -x["pivot_score"])

    best = all_pivots[0] if all_pivots else None

    # Modulation path
    mod_path = modulation_path_search(key_from, key_to, edo=edo)

    return {
        "key_from"              : key_from,
        "key_to"                : key_to,
        "key_from_name"         : nn(key_from, edo=edo),
        "key_to_name"           : nn(key_to, edo=edo),
        "cof_distance"          : cof_d,
        "cof_relationship"      : _cof_relationship(cof_d, edo),
        "common_scale_tones"    : common,
        "common_scale_names"    : [nn(p, edo=edo) for p in common],
        "n_common_tones"        : len(common),
        "diatonic_inventory_from": [
            {"roman": c["roman"], "label": c["label"], "pcs": c["pcs"],
             "function": c["function"]}
            for c in diatonic_chord_inventory(key_from, edo=edo)],
        "diatonic_inventory_to": [
            {"roman": c["roman"], "label": c["label"], "pcs": c["pcs"],
             "function": c["function"]}
            for c in diatonic_chord_inventory(key_to, edo=edo)],
        "pivots_by_type"        : pivots_by_type,
        "all_pivots"            : all_pivots,
        "best_pivot"            : best,
        "modulation_path"       : mod_path,
    }



# ──────────────────────────────────────────────────────────────────────────────
# MODULE 11 — RANK-R LATTICE TEMPERAMENT
# ──────────────────────────────────────────────────────────────────────────────

def find_lattice_mapping(primes: List[int], rank: int) -> List[List[int]]:
    """
    Find an optimal rank-r mapping matrix.
    Currently uses known good mappings for standard prime sets.
    """
    if primes == [2, 3, 5] and rank == 2:
        return [[1, 1, 0], [0, 1, 4]]
    if primes == [2, 3, 5, 7] and rank == 2:
        return [[1, 1, 0, -2], [0, 1, 4, 10]]

    # Generic fallback: use mappings from EDO-12 and EDO-31
    m1 = [round(12 * math.log2(p)) for p in primes]
    if rank == 1: return [m1]

    m2 = [round(31 * math.log2(p)) for p in primes]
    # Simple linear independence check
    is_parallel = True
    if len(primes) > 0:
        ratio = m1[0]/m2[0] if m2[0] != 0 else None
        for i in range(1, len(primes)):
            cur_ratio = m1[i]/m2[i] if m2[i] != 0 else None
            if cur_ratio != ratio:
                is_parallel = False
                break
    if is_parallel:
        m2 = [round(19 * math.log2(p)) for p in primes]

    return [m1, m2]

def solve_lattice_generators(M: List[List[int]], primes: List[int]) -> List[float]:
    """Solve for optimal (Tenney-L2) generator sizes given mapping M."""
    nP = len(primes)
    r  = len(M)
    v  = [math.log2(p) for p in primes]
    w_weight = [1.0/max(1e-3, abs(vi)) for vi in v]

    # M_W = M * diag(W)
    MW = [[M[i][j] * w_weight[j] for j in range(nP)] for i in range(r)]
    # A = MW * M^T (r x r)
    A = [[sum(MW[i][k] * M[j][k] for k in range(nP)) for j in range(r)] for i in range(r)]
    # B = MW * V (r x 1)
    B = [sum(MW[i][k] * v[k] for k in range(nP)) for i in range(r)]

    if r == 1:
        return [B[0] / A[0][0]] if A[0][0] != 0 else [1.0/12.0]
    if r == 2:
        det = A[0][0]*A[1][1] - A[0][1]*A[1][0]
        if abs(det) > 1e-12:
            return [(B[0]*A[1][1] - B[1]*A[0][1]) / det,
                    (A[0][0]*B[1] - A[1][0]*B[0]) / det]
    return [1.0/12.0] * r

def get_lattice_tuning_info(primes: List[int], rank: int) -> Dict[str, Any]:
    M = find_lattice_mapping(primes, rank)
    w = solve_lattice_generators(M, primes)

    nP = len(primes)
    v_target = [math.log2(p) for p in primes]
    v_actual = [sum(M[i][j] * w[i] for i in range(rank)) for j in range(nP)]
    errors_cents = [(v_actual[j] - v_target[j]) * 1200.0 for j in range(nP)]

    return {
        "primes": primes,
        "rank": rank,
        "mapping": M,
        "generators_log2": w,
        "generators_cents": [wi * 1200.0 for wi in w],
        "prime_errors_cents": [round(e, 4) for e in errors_cents],
        "max_error_cents": round(max(abs(e) for e in errors_cents), 4)
    }

def get_lattice_coords(log2_f: float, generators_log2: List[float], bounds: int = 10) -> List[int]:
    """Find lattice coordinates (k1, k2, ...) that best approximate log2_f."""
    rank = len(generators_log2)
    best_k = [0] * rank
    min_err = float('inf')

    if rank == 1:
        w1 = generators_log2[0]
        k1 = round(log2_f / w1)
        return [k1]

    if rank == 2:
        w1, w2 = generators_log2[0], generators_log2[1]
        # Search in a range around the naive projection
        for k1 in range(-bounds, bounds + 1):
            for k2 in range(-bounds, bounds + 1):
                err = abs(k1*w1 + k2*w2 - log2_f)
                if err < min_err:
                    min_err = err
                    best_k = [k1, k2]
        return best_k

    return best_k

def get_lattice_notes(generators_log2: List[float], bounds: List[int]) -> List[Dict[str, Any]]:
    """
    Generate all notes within specified lattice bounds.
    bounds is a list of [min, max] for each generator.
    """
    rank = len(generators_log2)
    ranges = [range(bounds[i*2], bounds[i*2+1] + 1) for i in range(rank)]
    notes = []
    for coords in itertools.product(*ranges):
        log2_f = sum(c * w for c, w in zip(coords, generators_log2))
        notes.append({
            "coords": list(coords),
            "log2_f": log2_f,
            "freq_hz": 261.625565 * (2.0 ** log2_f)
        })
    return notes


# ──────────────────────────────────────────────────────────────────────────────
# DISPATCH
# ──────────────────────────────────────────────────────────────────────────────

def _detect_quality_from_pcs(pcs: List[int], edo: int = 12) -> str:
    """
    Quality detection from a pitch-class set.
    Tests each pc as a potential root (handling all inversions).
    Returns the quality of the best-matching common chord type.
    Falls back to 'maj' if nothing matches.
    """
    if not pcs:
        return "maj"
    pcs_set = sorted(set(p % edo for p in pcs))

    def _d(s12: float) -> int:
        return round(s12 * edo / 12.0)

    TEMPLATES = [
        ("maj",   [0, _d(4), _d(7)]),
        ("min",   [0, _d(3), _d(7)]),
        ("dom7",  [0, _d(4), _d(7), _d(10)]),
        ("maj7",  [0, _d(4), _d(7), _d(11)]),
        ("min7",  [0, _d(3), _d(7), _d(10)]),
        ("dim",   [0, _d(3), _d(6)]),
        ("aug",   [0, _d(4), _d(8)]),
        ("hdim7", [0, _d(3), _d(6), _d(10)]),
        ("dim7",  [0, _d(3), _d(6), _d(9)]),
    ]
    for root_cand in pcs_set:
        ivs = sorted((p - root_cand) % edo for p in pcs_set)
        for qual, template in TEMPLATES:
            if ivs == sorted(template):
                return qual
    # Fallback: check if it has a minor third above any candidate root
    for root_cand in pcs_set:
        ivs = sorted((p - root_cand) % edo for p in pcs_set)
        if _d(3) in ivs and _d(4) not in ivs:
            return "min"
    return "maj"


def _detect_root_from_pcs(pcs: List[int], edo: int = 12) -> int:
    """
    Root detection from a pitch-class set.
    Returns the pc that best fits as a chord root by matching standard templates.
    Falls back to the lowest pc.
    """
    if not pcs:
        return 0
    pcs_set = sorted(set(p % edo for p in pcs))

    def _d(s12: float) -> int:
        return round(s12 * edo / 12.0)

    TEMPLATES = [
        ("maj",   [0, _d(4), _d(7)]),
        ("min",   [0, _d(3), _d(7)]),
        ("dom7",  [0, _d(4), _d(7), _d(10)]),
        ("maj7",  [0, _d(4), _d(7), _d(11)]),
        ("min7",  [0, _d(3), _d(7), _d(10)]),
        ("dim",   [0, _d(3), _d(6)]),
        ("aug",   [0, _d(4), _d(8)]),
        ("hdim7", [0, _d(3), _d(6), _d(10)]),
        ("dim7",  [0, _d(3), _d(6), _d(9)]),
    ]
    for root_cand in pcs_set:
        ivs = sorted((p - root_cand) % edo for p in pcs_set)
        for _, template in TEMPLATES:
            if ivs == sorted(template):
                return root_cand
    return pcs_set[0]   # lowest pc fallback


def handle(cmd: Dict[str, Any]) -> Dict[str, Any]:
    c    = cmd.get("cmd", "")
    key  = cmd.get("key", 0)
    edo  = cmd.get("edo", 12)
    pcs  = cmd.get("pcs", [])
    root = cmd.get("root", _detect_root_from_pcs(pcs, edo) if pcs else 0)
    # Auto-detect quality from pcs when not explicitly provided
    qual = cmd.get("quality", _detect_quality_from_pcs(pcs, edo) if pcs else "maj")
    c4hz = cmd.get("c4_hz", 261.625565)

    # ── Command aliases for backwards-compatibility ────────────────────────────
    _ALIASES = {
        "functional_analysis"  : "analyze_chord",
        "set_class_analysis"   : "set_class",
        "interval_vector"      : "set_class",
        "voice_leading"        : "orbifold_distance",
        "recognize_pattern"    : "detect_sequence",
        "edo_ji_lattice"       : "ji_lattice",
        "analyze"              : "analyze_chord",
        "voice_completion"     : "suggest_completion",
        "chord_completion"     : "suggest_completion",
    }
    c = _ALIASES.get(c, c)

    if c == "analyze_chord":
        if not pcs: pcs = triad_pcs(root, qual, edo)
        res = {"result": "analyze_chord"}
        res.update(set_class_info(pcs, edo))
        res.update(functional_analysis(pcs, key, edo))
        res.update(tonnetz_tension(pcs, key, edo))
        return res

    elif c == "psychoacoustic_analysis":
        """
        Full three-level auditory cortex analysis of a chord.
        Accepts: pcs, key, c4_hz (default 261.625565), octave (default 4),
                 level_db (default 70), n_harmonics (default 8)
        """
        octave     = cmd.get("octave", 4)
        level_db   = cmd.get("level_db", 70.0)
        n_harmonics= cmd.get("n_harmonics", 8)
        if not pcs: pcs = triad_pcs(root, qual, edo)
        return {
            "result": "psychoacoustic_analysis",
            **auditory_cortex_model(pcs, key, c4hz, octave, level_db, n_harmonics, edo)
        }

    elif c == "spectral_roughness":
        """
        Sethares spectral roughness for a chord.
        Accepts: pcs, c4_hz, octave, n_harmonics, rolloff
        """
        if not pcs: pcs = triad_pcs(root, qual, edo)
        return {
            "result": "spectral_roughness",
            **chord_roughness_psychoacoustic(
                pcs, c4hz,
                cmd.get("octave", 4),
                cmd.get("n_harmonics", 8),
                cmd.get("rolloff", 0.88),
                edo=edo)
        }

    elif c == "virtual_pitch":
        """Virtual pitch detection via Terhardt SHS."""
        if not pcs: pcs = triad_pcs(root, qual, edo)
        return {
            "result": "virtual_pitch",
            **virtual_pitch_strength(pcs, c4hz, cmd.get("octave", 4), edo=edo)
        }

    elif c == "masking_analysis":
        """Spectral masking analysis."""
        if not pcs: pcs = triad_pcs(root, qual, edo)
        return {
            "result": "masking_analysis",
            **chord_masking_analysis(pcs, c4hz, cmd.get("level_db", 70.0), edo=edo)
        }

    elif c == "roughness_curve":
        """
        Plomp-Levelt roughness curve centred at a reference frequency.
        Accepts: ref_hz (default 261.625565), n_points (default 150)
        Returns list of {delta_hz, roughness} for plotting.
        """
        ref_hz  = cmd.get("ref_hz", 261.625565)
        n_pts   = cmd.get("n_points", 150)
        curve   = plomp_levelt_curve_at_f(ref_hz, n_pts)
        return {"result": "roughness_curve", "ref_hz": ref_hz,
                "points": [{"delta_hz": round(d,2), "roughness": round(r,5)}
                           for d,r in curve]}

    elif c == "tonal_hierarchy":
        """Krumhansl-Kessler tonal stability."""
        if not pcs: pcs = triad_pcs(root, qual, edo)
        return {
            "result": "tonal_hierarchy",
            **tonal_hierarchy_score(pcs, key, cmd.get("mode", "major"), edo=edo)
        }

    elif c == "plr_transform":
        # Accept both 'op' and 'transform' as the operation key
        op = cmd.get("op") or cmd.get("transform") or "P"
        return {"result": "plr_transform",
                **plr_transform(root, qual, op, edo)}

    elif c == "plr_neighbors":
        return {"result": "plr_neighbors",
                "neighbors": all_plr_neighbors(root, qual, edo)}

    elif c == "plr_path":
        path = plr_path(root, qual,
                        cmd.get("root_b", 0), cmd.get("quality_b","maj"), edo=edo)
        return {"result": "plr_path", "path": path or [],
                "length": len(path) if path else -1,
                "reachable": path is not None}

    elif c == "resolution_paths":
        return {"result": "resolution_paths",
                "paths": resolution_paths(root, qual, key, edo)}

    elif c == "suggest_completion":
        # Accept both 'pcs' and 'pitch_classes' as input key
        input_pcs = cmd.get("pcs", cmd.get("pitch_classes", []))
        return {"result": "suggest_completion",
                "suggestions": suggest_completion(input_pcs, key, edo, c4hz)}

    elif c == "orbifold_distance":
        vl = orbifold_voice_leading(cmd.get("chord_a",[]), cmd.get("chord_b",[]), edo)
        return {"result": "orbifold_distance", **vl}

    elif c == "detect_sequence":
        # Accept either List[Dict] with 'root'/'quality' keys,
        # or List[List[int]] (pcs lists) — convert the latter automatically
        raw_chords = cmd.get("chords", [])
        chord_dicts = []
        for ch in raw_chords:
            if isinstance(ch, dict):
                chord_dicts.append(ch)
            elif isinstance(ch, (list, tuple)):
                ch_pcs = list(ch)
                q = _detect_quality_from_pcs(ch_pcs, edo)
                rt = _detect_root_from_pcs(ch_pcs, edo)
                chord_dicts.append({"root": rt, "quality": q, "pcs": ch_pcs})
        return {"result": "detect_sequence",
                **detect_sequence(chord_dicts, edo)}

    elif c == "modal_mixture":
        return {"result": "modal_mixture",
                **modal_mixture_analysis(pcs, key, edo)}

    elif c == "tonnetz_tension":
        t = tonnetz_tension(pcs, key, edo)
        # Expose both key names: 'tension' (internal) and 'tonnetz_tension' (documented)
        return {"result": "tonnetz_tension",
                "tonnetz_tension": t.get("tension"),
                **t}

    elif c == "ji_lattice":
        num, den, cents = ji_ratio(cmd.get("fifths",0),
                                   cmd.get("thirds",0),
                                   cmd.get("sevenths",0))
        return {"result": "ji_lattice",
                "ratio": f"{num}/{den}", "cents": round(cents, 3)}

    elif c == "edo_analysis":
        return {"result": "edo_analysis", **edo_analysis(edo)}

    elif c == "get_lattice_tuning":
        p = cmd.get("primes", [2, 3, 5])
        r = cmd.get("rank", 2)
        return {"result": "get_lattice_tuning", **get_lattice_tuning_info(p, r)}

    elif c == "get_lattice_coords":
        f = cmd.get("log2_f", 0.0)
        w = cmd.get("generators", [1.0, 0.5849])
        return {"result": "get_lattice_coords", "coords": get_lattice_coords(f, w)}

    elif c == "get_lattice_notes":
        w = cmd.get("generators", [1.0, 0.5849])
        b = cmd.get("bounds", [-2, 2, -5, 6]) # e.g. octaves -2..2, fifths -5..6
        return {"result": "get_lattice_notes", "notes": get_lattice_notes(w, b)}

    elif c == "tonnetz_projection":
        df, dt = tonnetz_projection_interval(cmd.get("interval", 0), edo)
        return {"result": "tonnetz_projection", "df": df, "dt": dt}

    elif c == "set_class":
        if not pcs: pcs = triad_pcs(root, qual, edo)
        return {"result": "set_class", **set_class_info(pcs, edo)}

    elif c == "pivot_search":
        """
        Full pivot-chord search between two keys.
        Required: key_from (int pc), key_to (int pc)
        Optional: include_borrowed, include_secondaries,
                  include_enharmonic, include_chromatic (all bool, default True)
                  max_per_type (int, default 6)
        """
        return {
            "result": "pivot_search",
            **pivot_search(
                cmd.get("key_from", key),
                cmd.get("key_to", 0),
                include_borrowed    = cmd.get("include_borrowed",    True),
                include_secondaries = cmd.get("include_secondaries", True),
                include_enharmonic  = cmd.get("include_enharmonic",  True),
                include_chromatic   = cmd.get("include_chromatic",   True),
                max_per_type        = cmd.get("max_per_type", 6),
                edo                 = edo,
            )
        }

    elif c == "diatonic_inventory":
        """
        List all diatonic chords (with optional borrowed) for a key.
        Required: key (int pc)
        Optional: mode (str, default Ionian), include_borrowed (bool),
                  include_sevenths (bool)
        """
        return {
            "result": "diatonic_inventory",
            "key": key, "key_name": nn(key, edo=edo),
            "mode": cmd.get("mode", "Ionian"),
            "chords": diatonic_chord_inventory(
                key,
                mode             = cmd.get("mode", "Ionian"),
                include_sevenths = cmd.get("include_sevenths", True),
                include_borrowed = cmd.get("include_borrowed", False),
                edo              = edo,
            )
        }

    elif c == "common_chord_pivots":
        """
        Find chords diatonic in both key_from and key_to.
        Required: key_from, key_to
        """
        return {
            "result"  : "common_chord_pivots",
            "key_from": cmd.get("key_from", 0),
            "key_to"  : cmd.get("key_to", 0),
            "pivots"  : common_chord_pivots(
                cmd.get("key_from", 0),
                cmd.get("key_to",   0),
                include_borrowed = cmd.get("include_borrowed", True),
                edo              = edo,
            )
        }

    elif c == "modulation_path":
        """
        Dijkstra-optimal modulation path from I(key_from) to I(key_to).
        Required: key_from, key_to
        Optional: max_depth (int, default 6)
        """
        return {
            "result"  : "modulation_path",
            "key_from": cmd.get("key_from", 0),
            "key_to"  : cmd.get("key_to",   0),
            "path"    : modulation_path_search(
                cmd.get("key_from", 0),
                cmd.get("key_to",   0),
                max_depth = cmd.get("max_depth", 6),
                edo       = edo,
            )
        }

    elif c == "roman_numeral":
        """
        Compute Roman numeral for a chord relative to a key.
        Required: root (int pc), quality (str), key (int pc)
        """
        return {
            "result"       : "roman_numeral",
            "root"         : root,
            "quality"      : qual,
            "key"          : key,
            "roman"        : roman_label(root, qual, key, edo),
            "degree_semitones": (root - key) % edo,
        }

    elif c == "ping":
        return {
            "result": "pong", "status": "ok",
            "modules": [
                "psychoacoustic_engine", "set_class", "plr_neo_riemannian",
                "functional_analysis", "voice_leading_orbifold",
                "tonnetz_tension_hybrid", "pattern_recognition",
                "edo_ji_lattice", "voice_completion_psychoacoustic",
                "pivot_search",
            ],
            "commands": [
                # Psychoacoustic engine
                "psychoacoustic_analysis", "spectral_roughness",
                "virtual_pitch", "masking_analysis", "roughness_curve",
                "tonal_hierarchy",
                # Chord analysis
                "analyze_chord",        # alias: functional_analysis, analyze
                "set_class",            # alias: set_class_analysis, interval_vector
                "roman_numeral",
                # Neo-Riemannian / PLR
                "plr_transform",        # op or transform key; P, L, R, N, S, H
                "plr_neighbors",
                "plr_path",
                "resolution_paths",
                # Voice leading
                "orbifold_distance",    # alias: voice_leading; keys: chord_a, chord_b
                # Tonnetz
                "tonnetz_tension",
                "tonnetz_projection",
                # Harmony / completion
                "suggest_completion",   # alias: voice_completion, chord_completion; key: pcs
                "modal_mixture",
                # Pattern recognition
                "detect_sequence",      # alias: recognize_pattern; key: chords
                # EDO / JI
                "edo_analysis",
                "ji_lattice",           # alias: edo_ji_lattice; keys: fifths, thirds, sevenths
                # Pivot / modulation
                "pivot_search",
                "common_chord_pivots",
                "modulation_path",
                "diatonic_inventory",
            ],
        }

    return {"result": "error", "message": f"unknown command: {c}"}


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    sys.stderr.write(
        "Harmonia Theory Server (Psychoacoustic Edition) started\n"
        "  Peripheral : Sethares roughness + IHC compression + Moore masking\n"
        "  Brainstem  : Terhardt virtual pitch (SHS) + AN rate-place profile\n"
        "  Cortical   : Krumhansl–Kessler tonal hierarchy + Tonnetz hybrid\n"
    )
    sys.stderr.flush()
    for line in sys.stdin:
        line = line.strip()
        if not line: continue
        try:
            req    = json.loads(line)
            tag    = req.get("tag")
            result = handle(req)
            if tag: result["tag"] = tag
            print(json.dumps(result), flush=True)
        except Exception as e:
            import traceback
            sys.stderr.write(traceback.format_exc())
            err_resp = {"result": "error", "message": str(e)}
            try:
                if 'req' in locals() and req.get("tag"):
                    err_resp["tag"] = req["tag"]
            except Exception:
                pass
            print(json.dumps(err_resp), flush=True)


if __name__ == "__main__":
    main()
