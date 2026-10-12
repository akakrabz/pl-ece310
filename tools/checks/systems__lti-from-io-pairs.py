"""Independent check of questions/systems/lti-from-io-pairs.

The data are re-classified from scratch: linearity is disproved iff rank[X | Y] > rank[X] (SymPy,
exact; rows = pairs, columns = a common index window), time-invariance iff some input is a pure
shift of another while its output is not the same shift of the other output, LTI iff some pair
violates y = x * h with h read off pair 1 (a unit impulse at n = m). Exactly one of the four
options must describe the result, and it must be the one marked correct. Part (b) is recomputed
with np.convolve."""

import numpy as np
import sympy

TEXT = {
    "lti": "consistent with an LTI system",
    "notL": "The system is not linear",
    "notTI": "The system is not time-invariant",
    "joint": "cannot be LTI, but these pairs show neither",
}


def _arr(vals, start, lo, hi):
    out = np.zeros(hi - lo, dtype=int)
    for i, v in enumerate(vals):
        out[start - lo + i] = v
    return out


def _trim(vals, start):
    a = np.trim_zeros(np.array(vals, dtype=int), "f")
    start += len(vals) - len(a)
    a = np.trim_zeros(a, "b")
    return tuple(int(v) for v in a), start


def _classify(params):
    P = params["pairs"]
    lo = min(min(p["nx"], p["ny"]) for p in P) - 12
    hi = max(max(p["nx"] + len(p["x"]), p["ny"] + len(p["y"])) for p in P) + 12
    X = sympy.Matrix([list(_arr(p["x"], p["nx"], lo, hi)) for p in P])
    Y = sympy.Matrix([list(_arr(p["y"], p["ny"], lo, hi)) for p in P])
    l_dis = X.row_join(Y).rank() > X.rank()
    ti_dis = False
    for i, pi in enumerate(P):
        for j, pj in enumerate(P):
            if i == j:
                continue
            xi, si = _trim(pi["x"], pi["nx"])
            xj, sj = _trim(pj["x"], pj["nx"])
            if xi == xj:
                k = sj - si
                yi, ti = _trim(pi["y"], pi["ny"])
                yj, tj = _trim(pj["y"], pj["ny"])
                if not (yi == yj and tj - ti == k):
                    ti_dis = True
    x1, s1 = _trim(P[0]["x"], P[0]["nx"])
    assert x1 == (1,), "pair 1 is not a unit impulse"
    h, hs = np.array(P[0]["y"]), P[0]["ny"] - s1
    lti_dis = False
    for p in P:
        y = np.convolve(p["x"], h)
        got = _trim(list(y), p["nx"] + hs)
        if got != _trim(p["y"], p["ny"]):
            lti_dis = True
    if not lti_dis:
        cls = "lti"
    elif l_dis and ti_dis:
        cls = "both"
    elif l_dis:
        cls = "notL"
    elif ti_dis:
        cls = "notTI"
    else:
        cls = "joint"
    return cls, h, hs, X.rank()


def check(params, correct):
    probs = []
    cls, h, hs, rx = _classify(params)
    opts = params["verdict"]
    marked = correct["verdict"]["html"]
    matching = [o for o in opts if cls in TEXT and TEXT[cls] in o["html"]]
    if cls == "both":
        probs.append("data disprove both properties: no option is right")
    if len(matching) != 1:
        probs.append(f"{len(matching)} options describe class {cls}")
    elif matching[0]["html"] != marked:
        probs.append(f"class {cls} but marked option is {marked!r}")
    if len(opts) != 4:
        probs.append(f"{len(opts)} options")
    if params["family"] != cls:
        probs.append(f"family {params['family']} but data classify as {cls}")
    if params["rel"] in ("none", "scaledshift", "sumshift", "shift") and rx != 3:
        probs.append("inputs are linearly dependent but the solution says they are independent")
    w = np.convolve(params["v"], h)
    ws = params["nv"] + hs
    got = np.array(correct["wrow"]["_value"][0])
    if correct["wfirst"] != ws or got.shape != w.shape or not np.allclose(got, w):
        probs.append(f"part (b): ({correct['wfirst']}, {got}) != ({ws}, {w})")
    if w[0] == 0 or w[-1] == 0:
        probs.append("part (b) output has a zero end sample")
    if params["m"] == 0:
        probs.append("m = 0: the shift trap of part (b) is absent")
    return probs


def _row(v, sep=", "):
    return "[" + sep.join(str(int(round(float(t)))) for t in v) + "]"


def submissions(params, correct):
    cls, h, hs, _ = _classify(params)
    w = np.convolve(params["v"], h)
    ws = params["nv"] + hs
    y1s = params["pairs"][0]["ny"]
    cases = [
        ({"verdict": correct["verdict"]["key"]}, {"verdict": 1}),
        ({"wrow": _row(w, " ")}, {"wrow": 1}),                                     # MATLAB style
        ({"wfirst": f"+{ws}" if ws >= 0 else f" {ws} "}, {"wfirst": 1}),
        ({"wfirst": str(params["nv"] + y1s)}, {"wfirst": 0}),                      # used y1 as h (shift not undone)
        ({"wrow": _row(list(w) + [0])}, {"wrow": 0}),
        ({"wfirst": "0.5"}, {"wfirst": "invalid"}),
    ]
    corr = np.convolve(params["v"], h[::-1])                                       # h not flipped (correlation)
    if corr.shape != w.shape or not np.allclose(corr, w):
        cases.append(({"wrow": _row(corr)}, {"wrow": 0}))
    else:
        off = list(w)
        off[0] += 1
        cases.append(({"wrow": _row(off)}, {"wrow": 0}))
    for o in params["verdict"]:
        if o["key"] != correct["verdict"]["key"]:
            cases.append(({"verdict": o["key"]}, {"verdict": 0}))
    return cases
