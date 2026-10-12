"""Independent check of questions/sampling/sinusoid-aliasing.

(a) the analog sinusoid is sampled with numpy and compared with A cos(omega n + phi) from the
answers; (b) f_min by brute force over every integer frequency 0 <= f <= f_s/2, fitting
a cos + b sin at that frequency to the samples (exact fit with amplitude A); (c) every listed option
is sampled and compared with the stated x[n]."""

import math
import re
from fractions import Fraction as F

import numpy as np

N = np.arange(64)


def _samples(kind, A, f, fs, th):
    arg = 2 * np.pi * f * N / fs + np.pi * th
    return A * (np.sin(arg) if kind == "sin" else np.cos(arg))


def _xa_samples(params):
    a = params["a"]
    return _samples(a["kind"], a["A"], float(F(a["f0"])), float(F(a["fs"])), float(F(a["th"])))


def fmin_brute(params):
    a = params["a"]
    fs = float(F(a["fs"]))
    x = _xa_samples(params)
    for f in range(0, int(fs // 2) + 1):
        M = np.column_stack([np.cos(2 * np.pi * f * N / fs), np.sin(2 * np.pi * f * N / fs)])
        coef, *_ = np.linalg.lstsq(M, x, rcond=None)
        if np.max(np.abs(M @ coef - x)) < 1e-7 and abs(np.hypot(*coef) - a["A"]) < 1e-7:
            return f
    return None


def _target_b(params):
    b = params["b"]
    fn = "sin" if b["kind"] == "sin" else "cos"
    return fn, b["A"] * (np.sin if fn == "sin" else np.cos)(np.pi * float(F(b["w0"])) * N + np.pi * float(F(b["ph"])))


def _opt_value(html):
    m = re.search(r"([fT]) = (-?)(\\tfrac\{(\d+)\}\{(\d+)\}|-?[\d.]+)", html)
    sign = -1 if m.group(2) else 1
    if m.group(4):
        return m.group(1), sign * F(int(m.group(4)), int(m.group(5)))
    return m.group(1), sign * F(m.group(3))


def option_ok(params, html):
    b = params["b"]
    fn, target = _target_b(params)
    var, v = _opt_value(html)
    ph = float(F(b["ph"]))
    if var == "f":
        x = _samples(fn, b["A"], float(v), float(F(b["fs"])), ph)
    else:
        x = _samples(fn, b["A"], float(F(b["f0"])), 1 / float(v), ph)
    return bool(np.max(np.abs(x - target)) < 1e-7)


def check(params, correct):
    probs = []
    a = params["a"]
    x = _xa_samples(params)
    w, phi = correct["w"], correct["phi"]
    y = a["A"] * np.cos(np.pi * w * N + np.pi * phi)
    if np.max(np.abs(x - y)) > 1e-7:
        probs.append(f"samples of x_a differ from {a['A']} cos({w} pi n + {phi} pi)")
    if not (0 < w < 1):
        probs.append(f"omega = {w} pi not in (0, pi)")
    if not (-1 < phi < 1) or abs(phi) < 1e-12:
        probs.append(f"phi = {phi} pi not in (-pi, pi) or zero")
    fm = fmin_brute(params)
    if fm is None or abs(fm - correct["fmin"]) > 1e-9:
        probs.append(f"f_min = {correct['fmin']} but brute force gives {fm}")
    marked = {o["html"] for o in correct["cands"]}
    for o in params["cands"]:
        if option_ok(params, o["html"]) != (o["html"] in marked):
            probs.append(f"option {o['html']} marked {o['html'] in marked}, sampling says {not o['html'] in marked}")
    if not 1 <= len(marked) <= 3:
        probs.append(f"{len(marked)} correct options (expected 1-3)")
    if params["b"]["kind"] != "sin":
        fs_or = [o for o in params["cands"] if o["html"] not in marked]
        if not fs_or:
            probs.append("no distractor")
    return probs


def _s(q):
    q = F(q).limit_denominator(10000)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def _differs(u, v):
    return abs(u - v) > 2 * (1e-6 + 1e-3 * abs(v))


def submissions(params, correct):
    a = params["a"]
    fs, f0, th = F(a["fs"]), F(a["f0"]), F(a["th"])
    r = 2 * f0 / fs
    w, phi = F(correct["w"]).limit_denominator(1000), F(correct["phi"]).limit_denominator(1000)
    thc = th - F(1, 2) if a["kind"] == "sin" else th
    cases = [
        ({"w": _s(w), "phi": _s(phi)}, {"w": 1, "phi": 1}),
        ({"w": f"{float(w):.4g}", "phi": f"{float(phi):.4g}"}, {"w": 1, "phi": 1}),
        ({"w": repr(float(w) * math.pi)}, {"w": 0}),                                   # radians, not multiples of pi
        ({"phi": _s(-phi)}, {"phi": 0}),                                               # sign of the phase
        ({"fmin": f"{correct['fmin']:.6g}"}, {"fmin": 1}),
        ({"fmin": f"{2 * int(correct['fmin'])}/2"}, {"fmin": 1}),
    ]
    if r > 1 and _differs(float(r), float(w)):
        cases.append(({"w": _s(r)}, {"w": 0}))                                         # not reduced to (0, pi)
    if r < 1:
        cases.append(({"w": _s(2 - r)}, {"w": 0}))                                     # "aliased" to f_s - f0 anyway
    m = r % 2
    if m > 1:                                                                          # folded: forgot the phase flip
        p0 = thc
        while p0 <= -1:
            p0 += 2
        while p0 > 1:
            p0 -= 2
        if _differs(float(p0), float(phi)):
            cases.append(({"phi": _s(p0)}, {"phi": 0}))
    if a["kind"] == "sin" and _differs(float(th), float(phi)) and _differs(float(-th), float(phi)):
        cases.append(({"phi": _s(th)}, {"phi": 0}))                                    # forgot sin = cos(. - pi/2)
    if _differs(float(f0), correct["fmin"]):
        cases.append(({"fmin": _s(f0)}, {"fmin": 0}))                                  # the original f0
    if _differs(float(fs) - correct["fmin"], correct["fmin"]):
        cases.append(({"fmin": _s(fs - F(int(correct['fmin'])))}, {"fmin": 0}))        # f_s - f
    opts = params["cands"]
    good = [o["key"] for o in opts if option_ok(params, o["html"])]
    bad = [o["key"] for o in opts if not option_ok(params, o["html"])]
    c = len(good)
    cases += [
        ({"cands": good}, {"cands": 1}),
        ({"cands": list(reversed(good))}, {"cands": 1}),
        ({"cands": good + bad[:1]}, {"cands": (c - 1) / c}),
        ({"cands": good[:1]}, {"cands": 1 / c}),
        ({"cands": bad}, {"cands": 0}),
    ]
    return cases
