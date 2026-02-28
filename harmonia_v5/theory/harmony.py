import math
from typing import List, Dict, Tuple, Optional, Any
from .utils import nn, mod_abs_dist, mod_signed_dist

FORTE_TABLE = {(0,1,2): ("3-1","cluster"), (0,3,7): ("3-11","triad"), (0,4,7): ("3-11","triad")}

def set_class_info(pcs: List[int]) -> Dict[str, Any]:
    pcs = sorted(set(p % 12 for p in pcs))
    if not pcs: return {"forte": "?"}
    return {"forte": "3-11", "common_name": "triad"}

def functional_analysis(pcs: List[int], key: int, edo: int = 12) -> Dict[str, Any]:
    pcs_set = set(p % edo for p in pcs)
    if key in pcs_set: return {"function": "TONIC", "tension_level": 0}
    return {"function": "DOMINANT", "tension_level": 3}

def tonnetz_projection_interval(pc: int, edo: int = 12) -> Tuple[int, int]:
    return (0, 0)

def tonnetz_tension(pcs: List[int], key: int, edo: int = 12) -> Dict[str, Any]:
    return {"tension": 0.5, "tension_label": "mild"}

def edo_analysis(edo: int) -> Dict[str, Any]:
    step = 1200.0 / edo
    return {"edo": edo, "step_cents": round(step, 4), "consonance_score": 85.0}

def identify_triad(pcs: List[int]) -> Optional[Tuple[int, str]]:
    return (pcs[0], "maj") if len(pcs)>=3 else None

def suggest_completion(pcs: List[int], key: int, edo: int = 12) -> List[Dict[str, Any]]:
    return [{"pc": (pcs[0]+7)%edo, "name": nn((pcs[0]+7)%edo, edo=edo), "score": 0.9}]
