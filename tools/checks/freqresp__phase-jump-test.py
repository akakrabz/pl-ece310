"""Independent check of questions/freqresp/phase-jump-test.

Each system's H_k(w) is rebuilt from the DISPLAYED text (impulse response or difference equation ->
taps -> DTFT sum; closed forms -> a numpy expression translated from the LaTeX) and compared with
params' spec. The +-pi jumps are counted from the phase itself: on a dense grid that never hits a
zero, a step of the principal phase whose wrapped size is about pi is a +-pi jump (a 2 pi wrap wraps
to about 0, a linear phase moves by a tiny amount). Checked: the checkbox answer is exactly the set of
systems with jumps; 1-3 such systems; a touched zero (min |H| = 0, no jump) among the others whenever
one exists, and a never-zero system whenever two exist; the designated count; the smallest positive
jump (or an empty answer when there is none), where |H| must vanish exactly. Probes: net-correct
partial scores (one touch system ticked too, only the touch system, everything); zeros of |H|
counted as jumps; only the jumps in (0, pi); 2 pi wraps counted; for (c) the empty box, a touched
zero, the second jump, 0."""

import math
import re
from fractions import Fraction

import numpy as np

GRID = np.pi * (-1 + (2 * np.arange(20000) + 1) / 20000)


def _terms(tex, sym):
    pat = re.compile(r"([+-]?)\s*(\d*)\s*" + re.escape(sym) + r"\[n(?:([+-])(\d+))?\]")
    out = {}
    for sg, c, s2, k in pat.findall(tex):
        coef = (-1 if sg == "-" else 1) * (int(c) if c else 1)
        shift = 0 if not k else (int(k) if s2 == "+" else -int(k))
        out[-shift] = coef
    return out


def _closed_fn(tex):
    s = tex.split("=", 1)[1]
    s = s.split(r",\quad")[0]
    s = s.replace(r"\left(", "(").replace(r"\right)", ")").replace(r"\,", " ")
    s = re.sub(r"\\tfrac\{(\d+)\}\{(\d+)\}", r"(\1/\2)", s)
    s = s.replace(r"e^{j|\omega|}", "exp(1j*abs(w))").replace(r"e^{-j|\omega|}", "exp(-1j*abs(w))")
    s = re.sub(r"e\^\{-j(\d*)\\omega\}", lambda m: f"exp(-1j*{m.group(1) or 1}*w)", s)
    s = re.sub(r"\\(cos|sin)\((\d*)\\omega\)", lambda m: f"{m.group(1)}({m.group(2) or 1}*w)", s)
    s = re.sub(r"^\s*j\s*", "1j*", s)
    s = re.sub(r"(\d|\))\s*(exp|cos|sin|\()", r"\1*\2", s)
    s = re.sub(r"\)\s*(exp|cos|sin|\()", r")*\1", s)
    assert re.fullmatch(r"[\d\s()+\-*/.jwexpcosinab]+", s), (tex, s)
    return lambda w: eval(s, {"__builtins__": {}}, {"w": w, "exp": np.exp, "cos": np.cos, "sin": np.sin, "abs": np.abs})


def H_fn(sysp):
    tex = sysp["def_tex"]
    if tex.startswith("y[n] ="):
        taps = _terms(tex[6:], "x")
    elif re.match(r"h_\d\[n\] =", tex):
        taps = _terms(tex.split("=", 1)[1], r"\delta")
    else:
        f = _closed_fn(tex)
        return lambda w: np.asarray(f(np.asarray(w, float)), complex) * np.ones(np.shape(w))
    n = np.array(sorted(taps), float)
    v = np.array([taps[k] for k in sorted(taps)], float)
    return lambda w: np.exp(-1j * np.outer(np.atleast_1d(w), n)) @ v


def H_spec(spec, w):
    w = np.asarray(w, float)
    k = spec["kind"]
    if "taps" in spec:
        n = np.array([t[0] for t in spec["taps"]], float)
        v = np.array([t[1] for t in spec["taps"]], float)
        return np.exp(-1j * np.outer(np.atleast_1d(w), n)) @ v
    if k == "quad":
        a0, a1, a2 = spec["a"]
        return np.exp(-1j * spec["M"] * w) * (a0 + a1 * np.cos(w) + a2 * np.cos(2 * w))
    if k == "cosK":
        return np.exp(-1j * spec["M"] * w) * np.cos(spec["K"] * w)
    if k == "sinK":
        return 1j * np.exp(-1j * spec["M"] * w) * np.sin(spec["K"] * w)
    b = float(Fraction(spec["b"]))
    return (b + np.cos(w)) * np.exp(1j * spec["s"] * np.abs(w))


def flips(H):
    ang = np.angle(H(GRID))
    d = np.angle(np.exp(1j * np.diff(ang)))
    idx = np.nonzero(np.abs(d) > np.pi / 2)[0]
    locs = [(GRID[i] + GRID[i + 1]) / 2 / np.pi for i in idx]
    wraps = int(np.sum((np.abs(np.diff(ang)) > 1.5 * np.pi) & (np.abs(d) < np.pi / 2)))
    return locs, wraps


def mag_zeros(H):
    """(locations of the zeros of |H| on (-pi, pi) as multiples of pi, smallest exact |H| found):
    grid minima, snapped to the nearest fraction with denominator <= 24 and evaluated there."""
    mag = np.abs(H(GRID))
    z = [GRID[i] / np.pi for i in range(1, len(GRID) - 1) if mag[i] <= mag[i - 1] and mag[i] <= mag[i + 1] and mag[i] < 1e-2]
    snapped = [float(Fraction(x).limit_denominator(24)) for x in z]
    vals = [abs(H(np.array([np.pi * x]))[0]) for x in snapped] + [abs(H(np.array([np.pi]))[0]), float(mag.min())]
    return z, float(min(vals))


def check(params, correct):
    probs = []
    systems = params["systems"]
    jumpers, nonj = [], []
    for s in systems:
        H = H_fn(s)
        if np.max(np.abs(H(GRID[::37]) - H_spec(s["spec"], GRID[::37]))) > 1e-9:
            probs.append(f"system {s['idx']}: displayed {s['def_tex']!r} differs from spec {s['spec']}")
        locs, _ = flips(H)
        if len(locs) != s["njumps"] or (len(locs) > 0) != s["jumps"]:
            probs.append(f"system {s['idx']}: {len(locs)} jumps on the grid, params say {s['njumps']} / {s['jumps']}")
        (jumpers if locs else nonj).append((s, H))
    right = {o["html"].strip() for o in correct["jumps"]}
    want = {f"System {s['idx']}" for s, _ in jumpers}
    if right != want:
        probs.append(f"checkbox answer {right}, jumping systems {want}")
    if not 1 <= len(jumpers) <= 3:
        probs.append(f"{len(jumpers)} systems with jumps")
    mins = [mag_zeros(H)[1] for _, H in nonj]
    if nonj and not any(m < 1e-9 for m in mins):
        probs.append("no touched-zero system among the non-jumping ones")
    if len(nonj) >= 2 and not any(m > 0.1 for m in mins):
        probs.append("no never-zero system among the non-jumping ones")
    des = next(s for s in systems if s["idx"] == params["k"])
    H = H_fn(des)
    locs, _ = flips(H)
    if correct["njumps"] != len(locs):
        probs.append(f"designated: {len(locs)} jumps, answer {correct['njumps']}")
    pos = [x for x in locs if x > 1e-6]
    if pos:
        if correct["wjump"] == "" or abs(pos[0] - correct["wjump"]) > 2e-3:
            probs.append(f"first positive jump at {pos[0]}, answer {correct['wjump']!r}")
        elif abs(H(np.array([np.pi * correct["wjump"]]))[0]) > 1e-9:
            probs.append("|H| is not zero at the answer of (c)")
    elif correct["wjump"] != "":
        probs.append(f"no positive jump, but (c) answer {correct['wjump']}")
    for s in systems:
        for t in (s["def_tex"], s["A_tex"], s.get("derive_tex", "")):
            for pat in (r"(?<![\d.])1\\(cos|sin|delta)", r"\+ -", r"- -", r"e\^\{-j1\\", r"(?<![\d.{])1x\["):
                if re.search(pat, t):
                    probs.append(f"display slip /{pat}/ in {t!r}")
    return probs


def _num(x):
    q = Fraction(x).limit_denominator(1000)
    if abs(float(q) - x) < 1e-12:
        return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return f"{x:.6g}"


def submissions(params, correct):
    cases = []
    keys = {o["html"].strip(): o["key"] for o in params["jumps"]}
    right = [o["key"] for o in correct["jumps"]]
    c = len(right)
    cases.append(({"jumps": right}, {"jumps": 1}))
    touch = [keys[f"System {s['idx']}"] for s in params["systems"] if s["cat"] == "touch"]
    if touch:
        cases.append(({"jumps": right + touch[:1]}, {"jumps": (c - 1) / c}))      # a touched zero taken for a jump
        cases.append(({"jumps": touch[:1]}, {"jumps": 0}))
    cases.append(({"jumps": list(keys.values())}, {"jumps": max(0.0, (c - (4 - c)) / c)}))
    des = next(s for s in params["systems"] if s["idx"] == params["k"])
    H = H_fn(des)
    n = correct["njumps"]
    cases.append(({"njumps": str(n)}, {"njumps": 1}))
    zs, _ = mag_zeros(H)
    if len(zs) != n:
        cases.append(({"njumps": str(len(zs))}, {"njumps": 0}))                    # every zero of |H| counted
    locs, wraps = flips(H)
    half = len([x for x in locs if x > 1e-6])
    if half != n:
        cases.append(({"njumps": str(half)}, {"njumps": 0}))                       # only (0, pi)
    if wraps:
        cases.append(({"njumps": str(n + wraps)}, {"njumps": 0}))                  # 2 pi wraps counted too
    cases.append(({"njumps": "2.5"}, {"njumps": "invalid"}))
    w = correct["wjump"]
    tz = [z for z in zs if z > 1e-3]
    if w == "":
        cases.append(({"wjump": ""}, {"wjump": 1}))
        cases.append(({"wjump": "  "}, {"wjump": 1}))
        if tz:
            cases.append(({"wjump": _num(float(Fraction(tz[0]).limit_denominator(24)))}, {"wjump": 0}))   # a touched zero
        cases.append(({"wjump": "0"}, {"wjump": 0}))
    else:
        cases.append(({"wjump": _num(w)}, {"wjump": 1}))
        cases.append(({"wjump": f"{w:.4f}"}, {"wjump": 1}))
        cases.append(({"wjump": ""}, {"wjump": 0}))
        pos = [x for x in locs if x > 1e-6]
        if len(pos) > 1:
            cases.append(({"wjump": _num(float(Fraction(pos[1]).limit_denominator(24)))}, {"wjump": 0}))  # not the smallest
        cases.append(({"wjump": "0"}, {"wjump": 0}))
    return cases
