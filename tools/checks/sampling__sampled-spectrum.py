"""Independent check of questions/sampling/sampled-spectrum.

X_a(Omega) is rebuilt in numpy from the stored shape parameters (in rad/s, not the u = omega/pi
units of the library), and X_d(omega) = (1/T) sum_k X_a((omega - 2 pi k)/T) is summed over
|k| <= 60 at omega -+ 1e-9 (both one-sided limits must agree: no query sits on a jump)."""

import math
from fractions import Fraction

import numpy as np

EPS = 1e-9


def _xa(spec):
    A, B = float(Fraction(spec["A"])), float(Fraction(spec["B"]))
    WB = 2 * math.pi * B
    shape = spec["shape"]

    def f(W):
        W = np.abs(np.asarray(W, float))
        inside = W < WB
        if shape == "tri":
            return np.where(inside, A * (1 - W / WB), 0.0)
        if shape == "rect":
            return np.where(inside, A, 0.0)
        if shape == "vee":
            return np.where(inside, A * W / WB, 0.0)
        W1 = float(Fraction(spec["r"])) * WB
        return np.where(W <= W1, A, np.where(inside, A * (WB - W) / (WB - W1), 0.0))
    return f


def xd(spec, w):
    T = 1 / float(Fraction(spec["fs"]))
    k = np.arange(-60, 61)
    return float(np.sum(_xa(spec)((w - 2 * np.pi * k) / T)) / T)


def xd_value(spec, u):
    w = math.pi * u
    lo, hi = xd(spec, w - EPS), xd(spec, w + EPS)
    return lo, hi


def central(spec, u):
    T = 1 / float(Fraction(spec["fs"]))
    return float(_xa(spec)(np.pi * u / T)) / T


def max_copies(spec, u):
    T = 1 / float(Fraction(spec["fs"]))
    return max(float(_xa(spec)((np.pi * u - 2 * np.pi * k) / T)) / T for k in range(-8, 9))


def _close(a, b, scale=1.0):
    return abs(a - b) <= 1e-6 * max(1.0, abs(b), scale)


def _scale(spec):
    return float(Fraction(spec["A"]) * Fraction(spec["fs"]))


def check(params, correct):
    probs = []
    spec = params["spec"]
    B, fs = Fraction(spec["B"]), Fraction(spec["fs"])
    b = 2 * B / fs
    if not _close(correct["wB"], float(b)):
        probs.append(f"omega_B {correct['wB']} != Omega_B T / pi = {b}")
    for name, u in (("X0", 0.0), ("Xpi", 1.0), ("X1", float(Fraction(params["u1"])))):
        lo, hi = xd_value(spec, u)
        if not _close(lo, hi, _scale(spec)):
            probs.append(f"{name}: X_d jumps at u = {u} ({lo} vs {hi}): ambiguous query")
        if not _close(correct[name], (lo + hi) / 2, _scale(spec)):
            probs.append(f"{name} = {correct[name]} but the copy sum gives {lo}, {hi} at u = {u}")
    reg = "above" if fs > 2 * B else ("at" if fs == 2 * B else "below")
    marked = correct["regime"]["html"]
    want = {"above": "Above", "at": "Exactly at", "below": "Below"}[reg]
    if not marked.strip().startswith(want):
        probs.append(f"regime marked {marked!r}, expected {reg} (f_s = {fs}, 2B = {2 * B})")
    if len(params["regime"]) != 3:
        probs.append("regime needs 3 options")
    u1 = Fraction(params["u1"])
    if not (0 < u1 < 1):
        probs.append(f"u1 = {u1} outside (0, 1)")
    # the figure must label the band edge and must not show the answers' digital frequencies
    if f"{2 * B}π".replace("1π", "π") not in params["fig_html"] and f"{2 * B}π" not in params["fig_html"]:
        probs.append("band edge tick missing in the figure")
    return probs


def _differs(wrong, right):
    return abs(wrong - right) > 2 * (1e-6 + 1e-3 * abs(right))


def _num(x):
    return repr(float(x)) if not float(x).is_integer() else str(int(x))


def _frac_str(x):
    q = Fraction(x).limit_denominator(1000)
    return f"{q.numerator}/{q.denominator}" if q.denominator != 1 else str(q.numerator)


def _key(params, start):
    return next(o["key"] for o in params["regime"] if o["html"].strip().startswith(start))


def submissions(params, correct):
    spec = params["spec"]
    B, fs = Fraction(spec["B"]), Fraction(spec["fs"])
    T = 1 / fs
    b = float(2 * B / fs)
    u1 = float(Fraction(params["u1"]))
    cases = [
        ({"wB": _frac_str(b)}, {"wB": 1}),                                   # exact fraction
        ({"wB": f"{b:.4g}"}, {"wB": 1}),                                     # 4 significant digits
        ({"wB": _num(float(2 * B * fs))}, {"wB": 0}),                        # Omega_B / T instead of Omega_B T
        ({"wB": _num(b / 2)}, {"wB": 0}),                                    # B T (dropped the 2)
        ({"X0": f"{correct['X0']:.5g}", "Xpi": _frac_str(correct["Xpi"])}, {"X0": 1, "Xpi": 1}),
        ({"X1": f"{2 * Fraction(correct['X1']).limit_denominator(1000).numerator}/"
                f"{2 * Fraction(correct['X1']).limit_denominator(1000).denominator}"}, {"X1": 1}),
    ]
    for name, u in (("X0", 0.0), ("Xpi", 1.0), ("X1", u1)):
        right = correct[name]
        if right != 0:                                                       # forgot the 1/T
            cases.append(({name: _num(right * float(T))}, {name: 0}))
        c0 = central(spec, u)
        lo, hi = xd(spec, math.pi * u - EPS), xd(spec, math.pi * u + EPS)
        sc = _scale(spec)
        if _close(lo, hi, sc) and _differs(c0, right) and _close(central(spec, u - 1e-7), central(spec, u + 1e-7), 1e3 * sc):
            cases.append(({name: _num(c0)}, {name: 0}))                      # only the copy at 0, no overlap added
        m = max_copies(spec, u)
        if _differs(m, right) and _differs(m, c0) and _close(max_copies(spec, u - 1e-7), max_copies(spec, u + 1e-7), 1e3 * sc):
            cases.append(({name: _num(m)}, {name: 0}))                       # copies drawn, not added
    reg = "above" if fs > 2 * B else ("at" if fs == 2 * B else "below")
    starts = {"above": "Above", "at": "Exactly at", "below": "Below"}
    for other in starts:
        cases.append(({"regime": _key(params, starts[other])}, {"regime": 1 if other == reg else 0}))
    cases.append(({"X0": "1/0"}, {"X0": "invalid"}))
    return cases
