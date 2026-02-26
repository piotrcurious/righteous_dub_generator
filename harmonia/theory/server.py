#!/usr/bin/env python3
"""
Harmonia Theory Server  ·  Structural Edition
══════════════════════════════════════════════
Replaces Markov probabilities with explicit algebraic structure.

Theoretical foundations:
  1. Set-class theory (interval vectors, prime form, Forte numbers, Z-relation)
  2. Neo-Riemannian group theory (PLR transformations on the Tonnetz)
  3. Diatonic functional analysis (function from algebraic containment)
  4. Voice-leading necessity (resolution rules as algebraic constraints)
  5. Tonnetz geometry (tension as distance from tonic region)
  6. Harmonic patterns (sequences, axes of symmetry, modal mixture)
  7. EDO/JI lattice analysis
"""

import sys, json, math, itertools
from typing import List, Dict, Tuple, Optional, Set, FrozenSet

NOTE_NAMES = ['C','D♭','D','E♭','E','F','G♭','G','A♭','A','B♭','B']

def nn(pc:int, prefer_flat=True, edo:int=12) -> str:
    if edo == 12:
        return NOTE_NAMES[pc % 12]
    return f"[{pc}]"


# ══════════════════════════════════════════════════════════════════════════════
#  MODULE 1 — SET-CLASS THEORY
#  Every chord is characterized by its interval-class vector (ICV).
#  ICV is a structural invariant under transposition AND inversion.
#  It reveals exactly which intervals appear and how often.
# ══════════════════════════════════════════════════════════════════════════════

def interval_class(a: int, b: int) -> int:
    """Interval class between two pitch classes (1-6)."""
    d = abs(a - b) % 12
    return min(d, 12 - d)

def interval_vector(pcs: List[int]) -> List[int]:
    """
    Interval-class vector (ICV) — the structural fingerprint of a chord.
    ICV[i] = number of interval-class (i+1) dyads in the chord.
    ICV is invariant under transposition and inversion.

    Example: major triad {0,4,7}
      pairs: (0,4)=ic4, (0,7)=ic5, (4,7)=ic3
      ICV = [0,0,1,1,1,0]
    """
    icv = [0] * 6
    for i in range(len(pcs)):
        for j in range(i+1, len(pcs)):
            ic = interval_class(pcs[i], pcs[j])
            if 1 <= ic <= 6:
                icv[ic-1] += 1
    return icv

def prime_form(pcs: List[int]) -> Tuple:
    """
    Allen Forte prime form: most compact, left-packed normal form
    invariant under transposition and inversion.
    """
    s = sorted(set(p % 12 for p in pcs))
    if not s: return ()
    n = len(s)

    def normal_form(rotation):
        t = [(p - rotation[0]) % 12 for p in rotation]
        return tuple(t)

    candidates = []
    for i in range(n):
        rot = [s[(i+j) % n] for j in range(n)]
        candidates.append(normal_form(rot))
    inv = sorted([(12 - p) % 12 for p in s])
    for i in range(n):
        rot = [inv[(i+j) % n] for j in range(n)]
        candidates.append(normal_form(rot))

    def pack_score(c):
        return (c[-1], c)

    return min(candidates, key=pack_score)

# Forte set-class table
FORTE_TABLE: Dict[Tuple, Tuple[str,str]] = {
    (0,1,2): ("3-1",  "chromatic cluster"),
    (0,1,3): ("3-2",  "minor + M2"),
    (0,1,4): ("3-3",  "minor + m3"),
    (0,1,5): ("3-4",  "minor + P4"),
    (0,1,6): ("3-5",  "minor + TT"),
    (0,2,4): ("3-6",  "whole-tone fragment"),
    (0,2,5): ("3-7",  "quartal"),
    (0,2,6): ("3-8",  "augmented 2nd"),
    (0,2,7): ("3-9",  "suspended (sus2/sus4)"),
    (0,3,6): ("3-10", "diminished triad"),
    (0,3,7): ("3-11", "minor / major triad"),
    (0,4,7): ("3-11", "major / minor triad"),
    (0,4,8): ("3-12", "augmented triad"),
    (0,2,3,5): ("4-10", "minor pentatonic fragment"),
    (0,1,2,4): ("4-2",  "minor 7 fragment"),
    (0,3,6,9): ("4-28", "diminished 7th"),
    (0,2,5,8): ("4-27", "half-diminished 7th"),
    (0,3,5,8): ("4-26", "minor 7th"),
    (0,1,4,8): ("4-19", "major 7th (no 5)"),
    (0,2,4,7): ("4-22", "dom 7 / min 7 fragment"),
    (0,2,4,8): ("4-24", "augmented + tone"),
    (0,2,6,8): ("4-25", "French aug 6th"),
    (0,1,5,8): ("4-20", "major 7th"),
    (0,3,4,7): ("4-17", "dim + maj chord"),
    (0,1,4,7): ("4-18", "minor-major 7th"),
    (0,2,4,6): ("4-21", "whole-tone tetrad"),
    (0,2,5,7): ("4-23", "quartal tetrad"),
    (0,1,6,7): ("4-9",  "tritone cluster"),
    (0,1,3,7): ("4-Z29","Z-related: [0,1,3,7]"),
    (0,1,4,6): ("4-Z15","Z-related: [0,1,4,6]"),
    (0,2,4,7,9): ("5-35", "pentatonic scale"),
    (0,1,3,5,7): ("5-23", "minor pentatonic inverse"),
    (0,2,4,5,7): ("5-34", "diatonic core"),
    (0,1,3,5,8): ("5-25", "min6/9 no root"),
    (0,2,3,5,7): ("5-24", "diatonic pentatonic variant"),
    (0,2,4,6,8): ("5-33", "whole-tone penta"),
    (0,1,2,4,7): ("5-20", "hexatonic penta"),
}

def set_class_info(pcs: List[int]) -> Dict:
    """Full structural analysis of a pitch-class set."""
    pcs = sorted(set(p % 12 for p in pcs))
    if not pcs:
        return {"prime_form":[], "interval_vector":[], "forte":"?", "common_name":"empty"}
    pf  = prime_form(pcs)
    icv = interval_vector(pcs)
    forte_name, common = FORTE_TABLE.get(pf, ("?", "unknown"))

    z_related = None
    if "Z15" in forte_name: z_related = "4-Z29"
    if "Z29" in forte_name: z_related = "4-Z15"

    inv_pf = prime_form([(12-p)%12 for p in pcs])
    inv_symmetric = (pf == inv_pf)

    t_sym = []
    for n in range(1, 12):
        if sorted(set((p+n)%12 for p in pcs)) == pcs:
            t_sym.append(n)

    IC_NAMES = ["m2/M7","M2/m7","m3/M6","M3/m6","P4/P5","TT"]
    ic_desc = ", ".join(f"{cnt}x{IC_NAMES[i]}" for i, cnt in enumerate(icv) if cnt > 0)

    return {
        "prime_form": list(pf),
        "interval_vector": icv,
        "forte": forte_name,
        "common_name": common,
        "z_related": z_related,
        "inv_symmetric": inv_symmetric,
        "transpositional_symmetry": t_sym,
        "cardinality": len(pcs),
        "ic_description": ic_desc or "no intervals",
    }


# ══════════════════════════════════════════════════════════════════════════════
#  MODULE 2 — NEO-RIEMANNIAN TRANSFORMATIONS
#
#  P, L, R are involutions on the 24 consonant triads.
#  Each changes exactly ONE voice by a minimal interval.
#  Their group structure: the dihedral group D_6 acting on triads.
#  This gives ALL triadic connections with explicit voice-leading derivation.
# ══════════════════════════════════════════════════════════════════════════════

def triad_pcs(root: int, quality: str) -> List[int]:
    q = quality.lower()
    if q in ('maj','major',''):    return [root%12, (root+4)%12, (root+7)%12]
    if q in ('min','minor','m'):   return [root%12, (root+3)%12, (root+7)%12]
    if q in ('dim','diminished'):  return [root%12, (root+3)%12, (root+6)%12]
    if q in ('aug','augmented'):   return [root%12, (root+4)%12, (root+8)%12]
    return [root%12, (root+4)%12, (root+7)%12]

def identify_triad(pcs: List[int]) -> Optional[Tuple[int,str]]:
    s = set(p%12 for p in pcs)
    for root in range(12):
        for q, ivs in [('maj',[4,7]),('min',[3,7]),('dim',[3,6]),('aug',[4,8])]:
            if s == {root%12, (root+ivs[0])%12, (root+ivs[1])%12}:
                return (root, q)
    return None

def plr_P(root: int, quality: str) -> Tuple[int,str]:
    """
    P (Parallel): C maj <-> C min.
    Keeps root and fifth. Changes third by 1 semitone.
    This is the minimal perturbation within the same root dyad.
    """
    if quality == 'maj': return (root, 'min')
    if quality == 'min': return (root, 'maj')
    return (root, quality)

def plr_L(root: int, quality: str) -> Tuple[int,str]:
    """
    L (Leittonwechsel): C maj <-> E min.
    Keeps third and fifth. Root descends 1 semitone (leading-tone exchange).
    Derived from: the only 1-semitone move keeping 2 notes of a major triad.
    """
    if quality == 'maj':
        return ((root + 4) % 12, 'min')
    if quality == 'min':
        return ((root - 4) % 12, 'maj')
    return (root, quality)

def plr_R(root: int, quality: str) -> Tuple[int,str]:
    """
    R (Relative): C maj <-> A min.
    Keeps root and third. Fifth moves by whole-tone.
    Derived from: the unique 2-semitone move keeping 2 notes of a major triad.
    """
    if quality == 'maj':
        return ((root + 9) % 12, 'min')
    if quality == 'min':
        return ((root + 3) % 12, 'maj')
    return (root, quality)

def plr_N(root: int, quality: str) -> Tuple[int,str]:
    """N (Nebenverwandt) = L then P then L. Fifth of major <-> root of minor."""
    r,q = plr_L(root,quality); r,q = plr_P(r,q); return plr_L(r,q)

def plr_S(root: int, quality: str) -> Tuple[int,str]:
    """S (Slide) = L then P then R. Shares only the third."""
    r,q = plr_L(root,quality); r,q = plr_P(r,q); return plr_R(r,q)

def plr_H(root: int, quality: str) -> Tuple[int,str]:
    """H (Hexatonic Pole) = L then P then L. Antipodal in Tonnetz."""
    r,q = plr_L(root,quality); r,q = plr_P(r,q); return plr_L(r,q)

PLR_OPS = {
    'P': (plr_P, "Parallel",       "Change 3rd ±1st. Keeps root+fifth."),
    'L': (plr_L, "Leittonwechsel", "Leading-tone exchange. Root ±1 semitone."),
    'R': (plr_R, "Relative",       "Relative maj/min. Fifth ±whole-tone."),
    'N': (plr_N, "Nebenverwandt",  "LPL: fifth of major <-> root of minor."),
    'S': (plr_S, "Slide",          "LPR: shares third only; root+fifth both move."),
    'H': (plr_H, "Hexatonic Pole", "LPL: maximum Tonnetz distance (tritone of roots)."),
}

def plr_transform(root: int, quality: str, op: str) -> Dict:
    fn, name, desc = PLR_OPS.get(op, (None,op,"unknown"))
    if fn is None:
        return {"error": f"unknown operation {op}"}

    new_root, new_quality = fn(root, quality)
    pcs_from = triad_pcs(root, quality)
    pcs_to   = triad_pcs(new_root, new_quality)

    common       = sorted(set(pcs_from) & set(pcs_to))
    changed_from = sorted(set(pcs_from) - set(pcs_to))
    changed_to   = sorted(set(pcs_to)   - set(pcs_from))

    motions = []
    for cf in changed_from:
        if changed_to:
            ct   = min(changed_to, key=lambda x: abs((x-cf)%12 - (cf-x)%12))
            diff = (ct - cf) % 12
            if diff > 6: diff -= 12
            motions.append({"from_pc": cf, "to_pc": ct, "semitones": diff,
                            "from_name": nn(cf), "to_name": nn(ct)})

    return {
        "op": op, "op_name": name, "op_description": desc,
        "from": {"root": root,     "quality": quality,
                 "label": nn(root)     + ("" if quality     == "maj" else "m"),
                 "pcs": pcs_from},
        "to":   {"root": new_root, "quality": new_quality,
                 "label": nn(new_root) + ("" if new_quality == "maj" else "m"),
                 "pcs": pcs_to},
        "common_tones":      common,
        "common_tone_names": [nn(p) for p in common],
        "voice_motions":     motions,
        "vl_cost":           sum(abs(m["semitones"]) for m in motions),
        "tonnetz_delta":     list(tonnetz_projection_interval(new_root)),
    }

def all_plr_neighbors(root: int, quality: str) -> List[Dict]:
    return [plr_transform(root, quality, op) for op in ['P','L','R']]

def plr_path(root: int, quality: str, target_root: int, target_quality: str,
             max_depth: int = 6) -> Optional[List[str]]:
    """BFS on the PLR graph to find shortest path between two triads."""
    from collections import deque
    start = (root, quality); goal = (target_root, target_quality)
    if start == goal: return []
    queue = deque([(start, [])])
    visited = {start}
    while queue:
        state, path = queue.popleft()
        if len(path) >= max_depth: continue
        r, q = state
        for op in ['P','L','R']:
            nr, nq = PLR_OPS[op][0](r, q)
            ns = (nr, nq)
            if ns == goal: return path + [op]
            if ns not in visited:
                visited.add(ns); queue.append((ns, path + [op]))
    return None


# ══════════════════════════════════════════════════════════════════════════════
#  MODULE 3 — FUNCTIONAL ANALYSIS FROM ALGEBRAIC CONTAINMENT
#
#  Function is determined by structural relationships, not statistics:
#  DOMINANT  = contains the diatonic tritone (unique in the key)
#  TONIC     = contains scale degree 1 without the tritone
#  SUBDOM    = contains 4th scale degree, excludes leading tone
# ══════════════════════════════════════════════════════════════════════════════

MODES = {
    'Ionian':     [0,2,4,5,7,9,11],
    'Dorian':     [0,2,3,5,7,9,10],
    'Phrygian':   [0,1,3,5,7,8,10],
    'Lydian':     [0,2,4,6,7,9,11],
    'Mixolydian': [0,2,4,5,7,9,10],
    'Aeolian':    [0,2,3,5,7,8,10],
    'Locrian':    [0,1,3,5,6,8,10],
}

def diatonic_set(key: int, mode: str = 'Ionian') -> List[int]:
    ivs = MODES.get(mode, MODES['Ionian'])
    return [(key + iv) % 12 for iv in ivs]

def tritone_of_key(key: int) -> Tuple[int,int]:
    """
    The diatonic tritone: scale degrees 7 and 4 (leading tone and subdominant).
    Appears EXACTLY ONCE in the major scale — its resolution is unique.
    """
    return ((key + 11) % 12, (key + 5) % 12)

def functional_analysis(pcs: List[int], key: int, edo: int = 12) -> Dict:
    """
    Determine harmonic function from algebraic structure alone.
    Every statement is derivable, not statistical.
    """
    pcs_set  = set(p % edo for p in pcs)
    # Scale diatonic set to EDO if not 12
    if edo == 12:
        diatonic = set(diatonic_set(key))
    else:
        # Simple heuristic for other EDOs: just use scaled major scale steps
        major_steps = [0, 2, 4, 5, 7, 9, 11]
        scaled_steps = [round(s * edo / 12) for s in major_steps]
        diatonic = set((key + s) % edo for s in scaled_steps)
    lt, sd   = tritone_of_key(key)      # leading-tone (7th), subdominant (4th)
    tonic    = key % 12
    dom_root = (key + 7) % 12

    in_diatonic    = pcs_set.issubset(diatonic)
    chromatic_tones = sorted(pcs_set - diatonic)
    has_tritone    = (lt in pcs_set and sd in pcs_set)

    degree_labels = {
        0:'1', 1:'b2', 2:'2', 3:'b3', 4:'3', 5:'4',
        6:'#4/b5', 7:'5', 8:'b6', 9:'6', 10:'b7', 11:'7'
    }
    scale_degrees = [
        {"pc": pc, "name": nn(pc),
         "degree": (pc-key)%12,
         "degree_label": degree_labels.get((pc-key)%12, str((pc-key)%12))}
        for pc in sorted(pcs_set)
    ]

    # ── Functional determination by algebraic rule
    if has_tritone and edo == 12:
        function = "DOMINANT"
        function_reason = (
            f"Contains the diatonic tritone [{nn(lt)},{nn(sd)}] (7 and 4 above {nn(key)}). "
            f"The major scale contains exactly one tritone — this one. "
            f"Its unique contrary-motion resolution: {nn(lt)}->{nn(tonic)} (up 1) "
            f"and {nn(sd)}->{nn((key+4)%12)} (down 1). "
            f"No other diatonic resolution of this tritone exists."
        )
        tension_level = 4 if len(pcs_set) >= 4 else 3
    elif tonic in pcs_set and not has_tritone:
        if (key+7)%12 in pcs_set or (key+4)%12 in pcs_set:
            function       = "TONIC"
            function_reason = (
                f"Contains {nn(tonic)} (1) with consonant support "
                f"({nn((key+4)%12) if (key+4)%12 in pcs_set else nn((key+7)%12)} = "
                f"{'3' if (key+4)%12 in pcs_set else '5'}) and no tritone. "
                f"Maximally stable: low roughness, tonic groundedness, no tendency tones."
            )
            tension_level = 0
        else:
            function       = "TONIC-WEAK"
            function_reason = f"Contains {nn(tonic)} without tritone, but sparse consonant support."
            tension_level = 1
    elif sd in pcs_set and lt not in pcs_set:
        function       = "SUBDOMINANT"
        function_reason = (
            f"Contains 4th degree {nn(sd)} without leading tone {nn(lt)}. "
            f"Moves away from tonic by subdominant motion (third-relation or P4 above tonic). "
            f"Pre-dominant or plagal function: can resolve to V (through V) or I (plagally)."
        )
        tension_level = 2
    elif (key+9)%12 in pcs_set and not has_tritone:
        function       = "TONIC-SUBST"
        function_reason = (
            f"vi-type chord: shares {nn(tonic)} and {nn((key+4)%12)} with I. "
            f"Tonic substitute by common-tone relation (deceptive resolution target)."
        )
        tension_level = 1
    elif (key+2)%12 in pcs_set and not has_tritone:
        function       = "PREDOMINANT"
        function_reason = (
            f"ii-type: supertonic {nn((key+2)%12)} (2) creates forward motion toward V. "
            f"Approaches dominant through upper-fifth relation. "
            f"Stronger than IV because it adds 2 — the note dominant seventh also contains."
        )
        tension_level = 2
    else:
        function       = "MEDIANT"
        function_reason = "Neither tonic stability nor dominant tension — modal color."
        tension_level = 1

    # ── Tendency tones
    tendency_tones = []
    if lt in pcs_set:
        tendency_tones.append({
            "pc": lt, "name": nn(lt), "role": "leading tone (7)",
            "tendency":   f"resolves UP 1 semitone to {nn(tonic)} (1)",
            "force":      "strong",
            "derivation": "Semitone below tonic = maximal upward gravitational pull. "
                          "Defined structurally: it IS the note 1 semitone below 1."
        })
    if sd in pcs_set:
        tendency_tones.append({
            "pc": sd, "name": nn(sd), "role": "subdominant (4)",
            "tendency":   f"resolves DOWN 1 semitone to {nn((key+4)%12)} (3)",
            "force":      "strong" if has_tritone else "moderate",
            "derivation": ("Forms tritone with leading tone — both resolve inward by semitone "
                           "(contrary motion). Algebraically unique." if has_tritone
                           else "Fourth degree in diatonic context tends stepwise to mediant.")
        })

    # ── Which modes contain this chord?
    containing_modes = []
    for mode_name, mode_ivs in MODES.items():
        for k2 in range(12):
            mode_set = set((k2+iv)%12 for iv in mode_ivs)
            if pcs_set and pcs_set.issubset(mode_set):
                containing_modes.append({"key": nn(k2), "mode": mode_name})
                break
        if len(containing_modes) >= 3: break

    return {
        "function":         function,
        "function_reason":  function_reason,
        "tension_level":    tension_level,
        "has_tritone":      has_tritone,
        "tritone":          [lt, sd] if has_tritone else [],
        "tritone_names":    [nn(lt), nn(sd)] if has_tritone else [],
        "in_diatonic":      in_diatonic,
        "chromatic_tones":  chromatic_tones,
        "scale_degrees":    scale_degrees,
        "tendency_tones":   tendency_tones,
        "containing_modes": containing_modes,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  MODULE 4 — VOICE-LEADING NECESSITY
#
#  Resolution rules are algebraic constraints, not correlational weights.
#  Priority order:
#    1. Tritone resolution (algebraically unique in key)
#    2. PLR minimal-motion neighbors
#    3. Diatonic bass-step motion
#    4. Common-tone maximization
# ══════════════════════════════════════════════════════════════════════════════

def orbifold_voice_leading(chord_a: List[int], chord_b: List[int]) -> Dict:
    """Optimal voice-leading in T^n/S_n (Tymoczko orbifold geodesic)."""
    n = max(len(chord_a), len(chord_b))
    a = list(chord_a) + [chord_a[-1]] * (n - len(chord_a)) if chord_a else []
    b = list(chord_b) + [chord_b[-1]] * (n - len(chord_b)) if chord_b else []

    if not a or not b:
        return {"distance": 0, "motions": [], "smoothness": "rest",
                "contrary_motions": 0, "parallel_motions": 0, "oblique_motions": 0,
                "description": "no motion"}

    best_d = float('inf'); best_match = []
    for perm in itertools.permutations(range(n)):
        d = 0; match = []
        for i, j in enumerate(perm):
            diff = (b[j] - a[i]) % 12
            if diff > 6: diff -= 12
            d += abs(diff); match.append((a[i], b[j], diff))
        if d < best_d: best_d = d; best_match = match

    motions = [{"from_pc": m[0], "from_name": nn(m[0]),
                "to_pc": m[1],   "to_name":   nn(m[1]),
                "semitones": m[2],
                "motion_type": ("oblique" if m[2]==0 else
                                "step" if abs(m[2])<=2 else
                                "leap" if abs(m[2])<=5 else "distant")}
               for m in best_match]

    n_contrary = sum(1 for i in range(len(motions))
                     for j in range(i+1, len(motions))
                     if motions[i]["semitones"] * motions[j]["semitones"] < 0)
    n_parallel = sum(1 for i in range(len(motions))
                     for j in range(i+1, len(motions))
                     if motions[i]["semitones"] == motions[j]["semitones"] != 0)
    n_oblique  = sum(1 for m in motions if m["semitones"] == 0)

    smoothness = "smooth" if best_d<=3 else "moderate" if best_d<=6 else "disjunct"

    return {
        "distance": best_d, "motions": motions,
        "contrary_motions": n_contrary, "parallel_motions": n_parallel,
        "oblique_motions":  n_oblique,  "smoothness":       smoothness,
        "description": f"Min. {best_d} semitone{'s' if best_d!=1 else ''} — {smoothness}",
    }

def resolution_paths(root: int, quality: str, key: int, edo: int = 12) -> List[Dict]:
    """
    Enumerate ALL structurally valid resolutions, each with its algebraic justification.
    These are not probabilities — they are derived from containment, interval structure,
    and minimal voice-leading geometry.
    """
    pcs  = triad_pcs(root, quality)
    func = functional_analysis(pcs, key, edo)
    results = []

    # ── Rule 1: Tritone resolution (if dominant function)
    if func["has_tritone"]:
        lt, sd = func["tritone"]
        for tgt_root, tgt_q, mode_label in [
            (key, 'maj', "tonic major"),
            (key, 'min', "tonic minor")
        ]:
            tgt_pcs = triad_pcs(tgt_root, tgt_q)
            vl = orbifold_voice_leading(pcs, tgt_pcs)
            results.append({
                "target_pcs": tgt_pcs,
                "target_label": nn(tgt_root) + ("" if tgt_q=="maj" else "m"),
                "target_root": tgt_root, "target_quality": tgt_q,
                "rule": "TRITONE_RESOLUTION",
                "rule_class": "necessity",
                "explanation": (
                    f"Tritone [{nn(lt)},{nn(sd)}] resolves inward by contrary semitone motion: "
                    f"{nn(lt)}->{nn(key)} (+1) and {nn(sd)}->{nn((key+4)%12)} (-1). "
                    f"This is the unique diatonic resolution — no other resolution "
                    f"exists within the scale that moves both tones by <=1 semitone."
                ),
                "voice_leading": vl, "priority": 1,
            })
        # Deceptive: V -> vi (leading tone still resolves, but 5th moves up instead)
        vi_root = (key+9)%12
        vi_pcs  = triad_pcs(vi_root, 'min')
        vl = orbifold_voice_leading(pcs, vi_pcs)
        results.append({
            "target_pcs": vi_pcs,
            "target_label": nn(vi_root)+"m",
            "target_root": vi_root, "target_quality": "min",
            "rule": "DECEPTIVE_CADENCE",
            "rule_class": "necessity",
            "explanation": (
                f"Deceptive cadence: {nn(lt)}->{nn(key)} (leading tone resolves normally), "
                f"but {nn((key+7)%12)}->{nn((key+9)%12)} (fifth steps up whole-tone). "
                f"The root {nn(key)} is present, denying full tonic arrival."
            ),
            "voice_leading": vl, "priority": 2,
        })

    # ── Rule 2: PLR minimal-motion neighbors (always structurally valid)
    for op in ['P','L','R']:
        t = plr_transform(root, quality, op)
        tgt_root, tgt_q = t["to"]["root"], t["to"]["quality"]
        tgt_pcs  = triad_pcs(tgt_root, tgt_q)
        vl       = orbifold_voice_leading(pcs, tgt_pcs)
        common   = sorted(set(pcs) & set(tgt_pcs))
        results.append({
            "target_pcs": tgt_pcs,
            "target_label": t["to"]["label"],
            "target_root": tgt_root, "target_quality": tgt_q,
            "rule": f"PLR_{op}",
            "rule_class": "minimal_motion",
            "explanation": (
                f"{op} ({t['op_name']}): {t['op_description']} "
                f"Common tones retained: {[nn(p) for p in common]}. "
                f"Single-voice motion: {t['voice_motions'][0]['from_name']}->"
                f"{t['voice_motions'][0]['to_name']} if available."
                if t.get('voice_motions') else t['op_description']
            ),
            "voice_leading": vl, "priority": 3,
        })

    # ── Rule 3: Diatonic bass-step motion
    diat = diatonic_set(key)
    if root%12 in diat:
        root_idx = diat.index(root%12)
        for step, dirn in [(-1, "descending"), (+1, "ascending")]:
            ni      = (root_idx + step) % 7
            nr      = diat[ni]
            # Quality from diatonic thirds
            third   = diat[(ni+2)%7]; fifth = diat[(ni+4)%7]
            t_iv    = (third-nr)%12;  f_iv  = (fifth-nr)%12
            nq      = 'maj' if t_iv==4 and f_iv==7 else \
                      'min' if t_iv==3 and f_iv==7 else \
                      'dim' if t_iv==3 and f_iv==6 else 'maj'
            tgt_pcs = triad_pcs(nr, nq)
            vl      = orbifold_voice_leading(pcs, tgt_pcs)
            results.append({
                "target_pcs": tgt_pcs,
                "target_label": nn(nr) + ("" if nq=="maj" else "m"),
                "target_root": nr, "target_quality": nq,
                "rule": f"DIATONIC_STEP_{dirn.upper()}",
                "rule_class": "scalar",
                "explanation": (
                    f"Scale-step {dirn}: root {nn(root%12)} ({root_idx+1}) -> "
                    f"{nn(nr)} ({ni+1}) in {nn(key)} major. "
                    f"Smooth bass, no voice-leading obligation."
                ),
                "voice_leading": vl, "priority": 4,
            })

    # Deduplicate and sort
    seen = set(); unique = []
    for r in sorted(results, key=lambda x: (x["priority"], x["voice_leading"]["distance"])):
        k = tuple(sorted(r["target_pcs"]))
        if k not in seen:
            seen.add(k); unique.append(r)
    return unique[:8]


# ══════════════════════════════════════════════════════════════════════════════
#  MODULE 5 — TONNETZ TENSION AS GEOMETRIC DISTANCE
# ══════════════════════════════════════════════════════════════════════════════

def tonnetz_projection_interval(pc: int, edo: int = 12) -> Tuple[int,int]:
    """Map pitch class to (fifths, thirds) lattice coordinates."""
    fifth_steps = round(math.log2(1.5) * edo)
    third_steps = round(math.log2(1.25) * edo)
    best = (0, 0); best_d = 999.0
    # Search a reasonable range around the origin
    for a in range(-8, 9):
        for b in range(-6, 7):
            if (a*fifth_steps + b*third_steps) % edo == pc % edo:
                d = abs(a) + abs(b)*1.4
                if d < best_d: best_d = d; best = (a, b)
    return best

# Keep old name for compat
def tonnetz_projection(iv: int, edo: int = 12) -> Tuple[int,int]:
    return tonnetz_projection_interval(iv, edo)

def tonic_region_centroid(key: int, edo: int = 12) -> Tuple[float, float]:
    """Centroid of the tonic region: I, vi, iii in the Tonnetz."""
    # Scale intervals to EDO
    vi_steps = round(9 * edo / 12)
    iii_steps = round(4 * edo / 12)
    region_roots = [key % edo, (key + vi_steps) % edo, (key + iii_steps) % edo]
    xs = [tonnetz_projection_interval(r, edo)[0] for r in region_roots]
    ys = [tonnetz_projection_interval(r, edo)[1] for r in region_roots]
    return sum(xs)/3.0, sum(ys)/3.0

def tonnetz_tension(pcs: List[int], key: int, edo: int = 12) -> Dict:
    if not pcs:
        return {"tension": 0.0, "tension_label": "rest", "distance": 0.0,
                "chord_centroid": [0.0,0.0], "tonic_centroid": [0.0,0.0]}
    cx = cy = 0.0
    for pc in pcs:
        x, y = tonnetz_projection_interval(pc, edo)
        cx += x; cy += y
    cx /= len(pcs); cy /= len(pcs)
    tcx, tcy = tonic_region_centroid(key, edo)
    dist = math.sqrt((cx-tcx)**2 + (cy-tcy)**2)
    norm = min(1.0, dist / 4.0)
    label = ("at rest" if norm < 0.15 else "mild" if norm < 0.35 else
             "tension" if norm < 0.55 else "strong" if norm < 0.75 else "maximal")
    return {
        "tension": round(norm,3), "tension_label": label,
        "distance": round(dist,3),
        "chord_centroid": [round(cx,2), round(cy,2)],
        "tonic_centroid": [round(tcx,2), round(tcy,2)],
    }


# ══════════════════════════════════════════════════════════════════════════════
#  MODULE 6 — PATTERN RECOGNITION (STRUCTURAL)
# ══════════════════════════════════════════════════════════════════════════════

def detect_sequence(chord_sequence: List[Dict]) -> Dict:
    """Detect algebraic pattern in a chord sequence."""
    if len(chord_sequence) < 3:
        return {"sequence_type": "none", "reason": "Too short"}

    roots = [c.get("root", 0) for c in chord_sequence]
    intervals = [(roots[i+1]-roots[i])%12 for i in range(len(roots)-1)]

    if len(set(intervals)) == 1:
        iv = intervals[0]
        interval_name = {
            1:"chromatic ascent", 2:"whole-tone ascent", 3:"minor-3rd cycle",
            4:"major-3rd cycle (hexatonic)", 5:"ascending fourths",
            7:"descending fourths (circle of 5ths)",
            9:"minor-3rd descent", 10:"whole-tone descent", 11:"chromatic descent"
        }.get(iv, f"constant +{iv} semitones")
        period = 12 // math.gcd(12, iv)
        return {
            "sequence_type": "transposition",
            "interval": iv, "interval_name": interval_name,
            "period": period,
            "description": f"T_{iv} transposition sequence. Completes after {period} steps.",
            "algebraic_explanation": (
                f"The operator T_{iv} acts on pitch-class space as a cyclic group element "
                f"of order {period}. The sequence is an orbit of this action."
            ),
        }

    # PLR sequence?
    ops_used = []
    for i in range(len(chord_sequence)-1):
        ra, qa = chord_sequence[i].get("root",0), chord_sequence[i].get("quality","maj")
        rb, qb = chord_sequence[i+1].get("root",0), chord_sequence[i+1].get("quality","maj")
        found = False
        for op in ['P','L','R']:
            nr, nq = PLR_OPS[op][0](ra, qa)
            if nr==rb and nq==qb:
                ops_used.append(op); found=True; break
        if not found: ops_used.append('?')

    if '?' not in ops_used and len(set(ops_used)) == 1:
        op = ops_used[0]
        return {
            "sequence_type": f"PLR_{op}_sequence",
            "operation": op, "period": 2,
            "description": f"Iterated {op}-transform ({PLR_OPS[op][1]}). "
                           f"Since {op}^2 = identity, this oscillates between two chords.",
        }

    fifth_descents = sum(1 for iv in intervals if iv==5)
    if fifth_descents >= len(intervals)-1:
        return {"sequence_type": "circle_of_fifths",
                "description": "Descending-fifth root motion — traverses the diatonic circle."}

    return {"sequence_type": "irregular", "operations": ops_used}

def modal_mixture_analysis(chord_pcs: List[int], key: int) -> Dict:
    """Structural analysis of modal borrowing."""
    major_pcs = set(diatonic_set(key))
    chord_set = set(p%12 for p in chord_pcs)
    chromatic = chord_set - major_pcs

    borrowed_from = []
    for mode_name in MODES:
        mode_set = set(diatonic_set(key, mode_name))
        if chord_set.issubset(mode_set) and mode_name != 'Ionian':
            borrowed_from.append(mode_name)

    pivot_keys = []
    for k2 in range(12):
        if k2 == key: continue
        if chord_set.issubset(set(diatonic_set(k2))):
            fa2 = functional_analysis(list(chord_set), k2)
            pivot_keys.append({"key": nn(k2), "function": fa2["function"]})

    return {
        "in_major": not chromatic,
        "chromatic_tones": sorted(chromatic),
        "chromatic_names": [nn(p) for p in sorted(chromatic)],
        "borrowed_from_modes": borrowed_from,
        "pivot_keys": pivot_keys[:4],
        "is_modal_mixture": bool(chromatic and borrowed_from),
    }


# ══════════════════════════════════════════════════════════════════════════════
#  MODULE 7 — EDO / JI LATTICE
# ══════════════════════════════════════════════════════════════════════════════

def ji_ratio(fifths: int, thirds: int, sevenths: int = 0):
    log2_val = (fifths*math.log2(3) + thirds*math.log2(5) +
                sevenths*math.log2(7)) % 1.0
    cents = log2_val * 1200.0
    num = 3**max(0,fifths)*5**max(0,thirds)*7**max(0,sevenths)
    den = 3**max(0,-fifths)*5**max(0,-thirds)*7**max(0,-sevenths)
    while num >= 2*den: den *= 2
    while den > num:    num *= 2
    def gcd(a,b): return a if b==0 else gcd(b,a%b)
    g = gcd(num,den); return num//g, den//g, cents

def edo_analysis(edo: int) -> Dict:
    step = 1200.0/edo
    result = {"edo":edo, "step_cents":round(step,4), "primes":{}, "intervals":[]}
    for p in [3,5,7,11,13]:
        ji_c = math.log2(p)*1200.0 % 1200.0
        steps = round(ji_c/step)
        err = steps*step - ji_c
        result["primes"][str(p)] = {
            "ji_cents":round(ji_c,3),"edo_steps":steps,
            "edo_cents":round(steps*step,3),"error_cents":round(err,3)}
    for name,num,den in [("octave",2,1),("fifth",3,2),("fourth",4,3),
                          ("major_third",5,4),("minor_third",6,5),("major_sixth",5,3),
                          ("harmonic_7th",7,4),("11th_harm",11,8)]:
        ji_c = math.log2(num/den)*1200.0
        steps = round(ji_c/step); err = steps*step - ji_c
        result["intervals"].append({"name":name,"ratio":f"{num}/{den}",
            "ji_cents":round(ji_c,2),"edo_steps":steps,
            "edo_cents":round(steps*step,2),"error_cents":round(err,2)})
    e5 = [abs(result["primes"][str(p)]["error_cents"]) for p in [3,5]]
    e7 = [abs(result["primes"][str(p)]["error_cents"]) for p in [3,5,7]]
    result["consonance_score"]  = round(100-min(sum(e5),100),1)
    result["seven_limit_score"] = round(100-min(sum(e7),100),1)
    return result


# ══════════════════════════════════════════════════════════════════════════════
#  MODULE 8 — VOICE COMPLETION WITH STRUCTURAL REASONING
# ══════════════════════════════════════════════════════════════════════════════

INTERVAL_ROUGHNESS = {
    0:0.00,1:0.90,2:0.70,3:0.20,4:0.12,5:0.08,
    6:0.55,7:0.05,8:0.15,9:0.13,10:0.25,11:0.30,
}

def suggest_completion(pitch_classes: List[int], key: int, edo: int = 12) -> List[Dict]:
    existing = set(p % edo for p in pitch_classes)
    if edo == 12:
        diat = set(diatonic_set(key))
    else:
        major_steps = [0, 2, 4, 5, 7, 9, 11]
        diat = set((key + round(s * edo / 12)) % edo for s in major_steps)
    ji_map   = {0:0.0,2:203.9,4:386.3,5:498.0,7:702.0,9:884.4,11:1088.3,
                3:315.6,6:590.2,8:813.7,10:969.0,1:111.7}

    suggestions = []
    for pc in range(12):
        if pc in existing: continue
        roughness_delta = sum(
            INTERVAL_ROUGHNESS.get(min(abs(pc-e)%12, 12-abs(pc-e)%12), 0.3)
            for e in existing)

        test_list = sorted(existing | {pc})
        triad_id  = identify_triad(test_list)
        fa        = functional_analysis(test_list, key, edo)
        tension   = tonnetz_tension(test_list, key, edo)

        reasons = []
        if triad_id:
            r,q = triad_id
            reasons.append(f"Completes {nn(r)}{'m' if q=='min' else ''} {q} triad")
        lt = (key+11)%12; sd = (key+5)%12
        if pc == lt: reasons.append(f"Adds leading tone {nn(lt)} (7) -> dominant function")
        if (pc==lt and sd in existing) or (pc==sd and lt in existing):
            reasons.append("Creates diatonic tritone -> dominant tension")
        if not existing: reasons.append("Provides root")
        if pc not in diat:
            mm = modal_mixture_analysis(test_list, key)
            if mm["borrowed_from_modes"]:
                reasons.append(f"Modal mixture from {mm['borrowed_from_modes'][0]}")
        sc = set_class_info(test_list)
        if sc.get("forte","?") != "?":
            reasons.append(f"Forms {sc['forte']} ({sc['common_name']})")

        score = max(0.0, 1.0 - roughness_delta/(len(existing)+1))
        suggestions.append({
            "pc": pc, "name": nn(pc),
            "score": round(score,4),
            "roughness_delta": round(roughness_delta,4),
            "cents_ji": round(ji_map.get((pc-key)%12, 0.0),1),
            "in_key": pc in diat,
            "function_if_added": fa["function"],
            "tension_if_added":  tension["tension"],
            "structural_reasons": reasons,
            "completes_triad": triad_id is not None,
        })

    suggestions.sort(key=lambda x: (-x["score"], x["roughness_delta"]))
    return suggestions[:8]


# ══════════════════════════════════════════════════════════════════════════════
#  DISPATCH
# ══════════════════════════════════════════════════════════════════════════════

def handle(cmd: Dict) -> Dict:
    c = cmd.get("cmd","")

    if c == "analyze_chord":
        pcs  = cmd.get("pcs",[])
        root = cmd.get("root", pcs[0] if pcs else 0)
        qual = cmd.get("quality","maj")
        key  = cmd.get("key",0)
        edo  = cmd.get("edo", 12)
        if not pcs: pcs = triad_pcs(root, qual)
        res = {"result":"analyze_chord"}
        res.update(set_class_info(pcs))
        res.update(functional_analysis(pcs, key, edo))
        res.update(tonnetz_tension(pcs, key, edo))
        return res

    elif c == "plr_transform":
        return {"result":"plr_transform",
                **plr_transform(cmd.get("root",0), cmd.get("quality","maj"),
                                cmd.get("op","P"))}

    elif c == "plr_neighbors":
        nbrs = all_plr_neighbors(cmd.get("root",0), cmd.get("quality","maj"))
        return {"result":"plr_neighbors", "neighbors": nbrs}

    elif c == "plr_path":
        path = plr_path(cmd.get("root_a",0), cmd.get("quality_a","maj"),
                        cmd.get("root_b",0), cmd.get("quality_b","maj"))
        return {"result":"plr_path", "path": path or [],
                "length": len(path) if path else -1,
                "reachable": path is not None}

    elif c == "resolution_paths":
        edo = cmd.get("edo", 12)
        return {"result":"resolution_paths",
                "paths": resolution_paths(cmd.get("root",0), cmd.get("quality","maj"),
                                          cmd.get("key",0), edo)}

    elif c == "suggest_completion":
        edo = cmd.get("edo", 12)
        return {"result":"suggest_completion",
                "suggestions": suggest_completion(cmd.get("pitch_classes",[]),
                                                  cmd.get("key",0), edo)}

    elif c == "orbifold_distance":
        vl = orbifold_voice_leading(cmd.get("chord_a",[]), cmd.get("chord_b",[]))
        return {"result":"orbifold_distance", **vl}

    elif c == "detect_sequence":
        return {"result":"detect_sequence",
                **detect_sequence(cmd.get("chords",[]))}

    elif c == "modal_mixture":
        return {"result":"modal_mixture",
                **modal_mixture_analysis(cmd.get("pcs",[]), cmd.get("key",0))}

    elif c == "tonnetz_tension":
        return {"result":"tonnetz_tension",
                **tonnetz_tension(cmd.get("pcs",[]), cmd.get("key",0), cmd.get("edo", 12))}

    elif c == "ji_lattice":
        num,den,cents = ji_ratio(cmd.get("fifths",0),cmd.get("thirds",0),
                                  cmd.get("sevenths",0))
        return {"result":"ji_lattice","ratio":f"{num}/{den}","cents":round(cents,3)}

    elif c == "edo_analysis":
        return {"result":"edo_analysis", **edo_analysis(cmd.get("edo",12))}

    elif c == "tonnetz_projection":
        df,dt = tonnetz_projection_interval(cmd.get("interval",0))
        return {"result":"tonnetz_projection","df":df,"dt":dt}

    elif c == "ping":
        return {"result":"pong","status":"ok",
                "modules":["set_class","plr","functional","voice_leading",
                           "tension","patterns","edo","completion"]}

    return {"result":"error","message":f"unknown command: {c}"}


def main():
    sys.stderr.write("Harmonia Theory Server (Structural Edition) started\n")
    sys.stderr.flush()
    for line in sys.stdin:
        line = line.strip()
        if not line: continue
        try:
            req = json.loads(line)
            tag = req.get("tag")
            result = handle(req)
            if tag:
                result["tag"] = tag
            print(json.dumps(result), flush=True)
        except Exception as e:
            import traceback
            sys.stderr.write(traceback.format_exc())
            err_resp = {"result":"error","message":str(e)}
            if 'req' in locals() and req.get("tag"):
                err_resp["tag"] = req.get("tag")
            print(json.dumps(err_resp), flush=True)

if __name__ == "__main__":
    main()
