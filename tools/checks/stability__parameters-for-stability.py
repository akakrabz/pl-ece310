"""Independent check of questions/stability/parameters-for-stability.

For many values of K the checker evaluates b(K), a(K) from params["model"], finds poles and zeros
with numpy.roots, cancels coincident pole/zero pairs, and calls the system stable iff every
surviving pole has magnitude < 1. The answer (interval + endpoint membership, single value, or set
option) must agree with that verdict on the whole grid, at the endpoints and just inside / outside them;
every distractor of a multiple choice must disagree with it somewhere."""

import html
import re
from fractions import Fraction

import numpy as np


def _coefs(lst, K):
    return [float(Fraction(c0)) + float(Fraction(c1)) * K for c0, c1 in lst]


def _roots(c):
    """Roots in z of sum_k c[k] z^-k (multiply by z^(len-1))."""
    c = list(c)
    while len(c) > 1 and abs(c[-1]) < 1e-14:
        c.pop()                                    # trailing zero coefficients = roots at z = 0
    return list(np.roots(c)) if len(c) > 1 else []


def stable(params, K):
    b = _coefs(params["model"]["b"], K)
    a = _coefs(params["model"]["a"], K)
    poles, zeros = _roots(a), _roots(b)
    for zr in zeros:
        for i, p in enumerate(poles):
            if abs(p - zr) < 1e-7:
                poles.pop(i)
                break
    return all(abs(p) < 1 - 1e-9 for p in poles)


def _grid(lo, hi, extra=()):
    pts = list(np.linspace(lo - 4, hi + 4, 321))
    for e in (lo, hi, *extra):
        pts += [float(e), e - 1e-3, e + 1e-3]
    return pts


def _norm(s):
    return re.sub(r"\s+", "", html.unescape(s or ""))


def _member(d, K):
    if d.get("outside"):
        return abs(K) > d["c"]
    inside = False
    if d["lo"] is not None:
        inside = (d["lo"] < K < d["hi"]) or (d["lo_in"] and K == d["lo"]) or (d["hi_in"] and K == d["hi"])
    return inside or (d["point"] is not None and abs(K - float(Fraction(d["point"]))) < 1e-12)


def check(params, correct):
    probs = []
    fam = params["family"]
    if fam == "interval":
        lo, hi = correct["K_lo"], correct["K_hi"]
        if (lo, hi) != (params["lo"], params["hi"]):
            probs.append("endpoint params disagree with correct answers")
        opt = [c for c in params["ends_choices"] if _norm(c["text"]) == _norm(correct["K_ends"]["html"])]
        if len(opt) != 1:
            return probs + ["cannot identify the correct endpoint option"]
        lo_in, hi_in = opt[0]["lo_in"], opt[0]["hi_in"]
        for K in _grid(lo, hi):
            want = (lo < K < hi) or (lo_in and K == lo) or (hi_in and K == hi)
            if stable(params, K) != want:
                probs.append(f"K = {K:.4f}: numpy says stable={stable(params, K)}, answer says {want}")
                break
    elif fam == "single":
        Kv = correct["K_val"]
        if not stable(params, Kv):
            probs.append(f"K = {Kv} does not make the system stable")
        for K in _grid(Kv - 3, Kv + 3):
            if abs(K - Kv) > 1e-6 and stable(params, K):
                probs.append(f"K = {K:.4f} is also stable (answer is not unique)")
                break
        for w in params["Kwrong"]:
            if stable(params, float(Fraction(w))):
                probs.append(f"distractor value {w} is stable")
    else:
        opts = [c for c in params["K_choices"] if _norm(c["text"]) == _norm(correct["K_set"]["html"])]
        if len(opts) != 1:
            return ["cannot identify the correct option"]
        s = opts[0]["set"]
        Kc = float(Fraction(params["Kc"]))
        c = max(abs(o["set"].get("c") or o["set"].get("hi") or 0) for o in params["K_choices"])
        grid = _grid(-c, c, extra=(Kc,))
        for K in grid:
            if stable(params, K) != _member(s, K):
                probs.append(f"K = {K:.4f}: numpy says stable={stable(params, K)}, option {opts[0]['text']} says {_member(s, K)}")
                break
        for o in params["K_choices"]:
            if o is opts[0]:
                continue
            if all(_member(o["set"], K) == stable(params, K) for K in grid):
                probs.append(f"distractor {o['text']} is also correct")
        if len({_norm(o["text"]) for o in params["K_choices"]}) != len(params["K_choices"]):
            probs.append("duplicate option texts")
    return probs


def submissions(params, correct):
    fam = params["family"]
    if fam == "interval":
        lo, hi = int(correct["K_lo"]), int(correct["K_hi"])
        key = {_norm(o["html"]): o["key"] for o in params["K_ends"]}
        texts = {(c["lo_in"], c["hi_in"]): c["text"] for c in params["ends_choices"]}
        right = (params["lo_in"], params["hi_in"])
        cases = [
            ({"K_lo": f"{2 * lo}/2", "K_hi": f"{hi}.0"}, {"K_lo": 1, "K_hi": 1}),           # equivalent forms
            ({"K_lo": f" {lo} ", "K_hi": f"{3 * hi}/3"}, {"K_lo": 1, "K_hi": 1}),
            ({"K_lo": "0"}, {"K_lo": 0 if lo != 0 else 1}),                                 # "K must be positive"
            ({"K_lo": str(lo - 1), "K_hi": str(hi + 1)}, {"K_lo": 0, "K_hi": 0}),           # off by one
            ({"K_hi": "abc"}, {"K_hi": "invalid"}),
            ({"K_ends": key[_norm(texts[right])]}, {"K_ends": 1}),
        ]
        if hi - lo != 2:                                                                  # |p| < 1 read as |K| < 1
            cases.append(({"K_lo": "-1", "K_hi": "1"}, {"K_lo": 0 if lo != -1 else 1, "K_hi": 0 if hi != 1 else 1}))
        # endpoint slips: "<=" everywhere (marginal = stable), or ignoring the cancellation at an endpoint
        for wrong in {(True, True), (False, False)} - {right}:
            cases.append(({"K_ends": key[_norm(texts[wrong])]}, {"K_ends": 0}))
        return cases
    if fam == "single":
        Kv = Fraction(params["Kval"])
        cases = [
            ({"K_val": f"{Kv.numerator}/{Kv.denominator}" if Kv.denominator != 1 else f"{Kv.numerator}/1"}, {"K_val": 1}),
            ({"K_val": f"{float(Kv):.8f}"}, {"K_val": 1}),
        ]
        if Kv.denominator != 1:                                                           # 4 significant digits suffice
            cases.append(({"K_val": f"{float(Kv):.4g}"}, {"K_val": 1}))
        for w in params["Kwrong"]:            # sign error, z instead of z^-1, cancelling the stable pole
            cases.append(({"K_val": str(Fraction(w))}, {"K_val": 0}))
        cases.append(({"K_val": "K"}, {"K_val": "invalid"}))
        return cases
    key = {_norm(o["html"]): o["key"] for o in params["K_set"]}
    cases = []
    for c in params["K_choices"]:
        cases.append(({"K_set": key[_norm(c["text"])]}, {"K_set": 1 if c["correct"] else 0}))
    return cases
