import math
from typing import List, Dict, Tuple, Optional, Any

NOTE_NAMES   = ['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
NOTE_NAMES_b = ['C','Db','D','Eb','E','F','Gb','G','Ab','A','Bb','B']

def nn(pc: int, prefer_flat: bool = False, edo: int = 12) -> str:
    if edo == 12:
        return (NOTE_NAMES_b if prefer_flat else NOTE_NAMES)[pc % 12]
    return f"{pc}"

def mod_signed_dist(a: int, b: int, modulus: int = 12) -> int:
    d = (b - a) % modulus
    if d > modulus // 2:
        d -= modulus
    return d

def mod_abs_dist(a: int, b: int, modulus: int = 12) -> int:
    return abs(mod_signed_dist(a, b, modulus))
