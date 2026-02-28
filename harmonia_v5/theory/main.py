#!/usr/bin/env python3
import sys
import json
import traceback
from .utils import nn
from .psychoacoustics import auditory_cortex_model, chord_roughness_psychoacoustic, virtual_pitch_strength
from .harmony import functional_analysis, set_class_info, tonnetz_tension, edo_analysis, suggest_completion
from .modulation import plr_transform, pivot_search, orbifold_voice_leading

def handle(cmd):
    c = cmd.get("cmd", "")
    key = cmd.get("key", 0)
    edo = cmd.get("edo", 12)
    freqs = cmd.get("freqs", [])
    pcs = [int(round(math.log2(f / 261.63) * edo)) % edo for f in freqs] if freqs else cmd.get("pcs", [])

    if c == "analyze_chord":
        res = {"result": "analyze_chord"}
        res.update(set_class_info(pcs))
        res.update(functional_analysis(pcs, key, edo))
        res.update(tonnetz_tension(pcs, key, edo))
        return res
    elif c == "psychoacoustic_analysis":
        return {"result": "psychoacoustic_analysis", **auditory_cortex_model(freqs, key, edo)}
    elif c == "edo_analysis":
        return {"result": "edo_analysis", **edo_analysis(edo)}
    elif c == "pivot_search":
        return {"result": "pivot_search", **pivot_search(cmd.get("key_from", key), cmd.get("key_to", 0))}
    elif c == "suggest_completion":
        return {"result": "suggest_completion", "suggestions": suggest_completion(pcs, key, edo)}
    elif c == "ping":
        return {"result": "pong", "status": "ok"}
    return {"result": "error", "message": f"unknown command: {c}"}

def main():
    for line in sys.stdin:
        if not line.strip(): continue
        try:
            req = json.loads(line)
            res = handle(req)
            if "tag" in req: res["tag"] = req["tag"]
            print(json.dumps(res), flush=True)
        except Exception:
            print(json.dumps({"result": "error", "message": traceback.format_exc()}), flush=True)

if __name__ == "__main__":
    main()
