import math
from typing import List, Dict, Tuple, Optional, Any
from .utils import nn

def hz_to_bark(f: float) -> float:
    if f <= 0.0: return 0.0
    return (26.81 * f / (1960.0 + f)) - 0.53

def critical_bandwidth_zwicker(f: float) -> float:
    return 25.0 + 75.0 * (1.0 + 1.4 * (f / 1000.0) ** 2) ** 0.69

def greenwood_place(f: float) -> float:
    if f <= 0.0: return 0.0
    return (1.0 / 0.06) * math.log10(f / 165.4 + 0.88)

def tonotopic_distance_mm(f1: float, f2: float) -> float:
    return abs(greenwood_place(f1) - greenwood_place(f2))

def plomp_levelt_roughness_pair(f1: float, f2: float, a1: float = 1.0, a2: float = 1.0) -> float:
    if f1 <= 0.0 or f2 <= 0.0: return 0.0
    f_lo = min(f1, f2)
    cb   = 0.021 * f_lo + 19.0
    s    = abs(f2 - f1) / cb
    roughness = math.exp(-3.5 * s) - math.exp(-5.75 * s)
    return max(0.0, a1 * a2 * roughness)

def sethares_roughness(freqs: List[float], amps: List[float]) -> float:
    total = 0.0
    n = len(freqs)
    for i in range(n):
        for j in range(i + 1, n):
            total += plomp_levelt_roughness_pair(freqs[i], freqs[j], amps[i], amps[j])
    return total

def chord_roughness_psychoacoustic(freqs: List[float], n_harmonics: int = 8, rolloff: float = 0.88, edo: int = 12) -> Dict[str, Any]:
    if not freqs: return {"roughness_raw": 0.0, "consonance_score": 1.0}
    note_freqs = sorted(freqs)
    pcs = [int(round(math.log2(f / 261.63) * edo)) % edo for f in note_freqs]
    all_freqs, all_amps, note_partial_idx = [], [], []
    for f0 in note_freqs:
        start = len(all_freqs)
        for k in range(1, n_harmonics + 1):
            all_freqs.append(f0 * k)
            all_amps.append(rolloff ** (k - 1))
        note_partial_idx.append((start, len(all_freqs)))
    raw = sethares_roughness(all_freqs, all_amps)
    f_ref = note_freqs[0]
    r_ref = plomp_levelt_roughness_pair(f_ref, f_ref * 2.0 ** (1.0 / edo))
    n_pairs = max(1, len(note_freqs) * (len(note_freqs) - 1) // 2)
    norm = raw / (n_pairs * n_harmonics ** 2 * r_ref) if r_ref > 0 else 0.0
    return {"roughness_raw": round(raw, 5), "roughness_normalized": round(min(1.0, norm), 4), "consonance_score": round(max(0.0, 1.0 - norm), 4)}

def virtual_pitch_strength(freqs: List[float], n_harmonics: int = 8, edo: int = 12) -> Dict[str, Any]:
    if not freqs: return {"harmonicity": 0.0}
    f_lo, f_hi = min(freqs) / 4.0, max(freqs) * 1.05
    n_steps = max(1, round(12 * math.log2(f_hi / f_lo)))
    best_f0, best_score = 0.0, -1.0
    for i in range(n_steps + 1):
        f0 = f_lo * 2.0 ** (i / 12.0)
        score = sum((1.0/k) * math.exp(-0.5 * ((fn - k*f0) / (0.02*k*f0))**2) for fn in freqs for k in range(1, n_harmonics+1))
        if score > best_score: best_score = score; best_f0 = f0
    harmonicity = min(1.0, best_score / (sum(1.0/k for k in range(1, n_harmonics+1)) * len(freqs)))
    pc = int(round(math.log2(best_f0 / 261.63) * edo)) % edo
    return {"virtual_pitch_hz": round(best_f0, 2), "virtual_pitch_name": nn(pc, edo=edo), "harmonicity": round(harmonicity, 4)}

def auditory_cortex_model(freqs: List[float], key: int = 0, edo: int = 12) -> Dict[str, Any]:
    rough = chord_roughness_psychoacoustic(freqs, edo=edo)
    vp = virtual_pitch_strength(freqs, edo=edo)
    return {"perceptual_tension": round(0.5 * rough["roughness_normalized"] + 0.5 * (1.0 - vp["harmonicity"]), 4), "roughness": rough, "virtual_pitch": vp}
