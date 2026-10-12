"""Independent check of questions/dtft/modulation.

Y_d is re-evaluated with numpy from the modulation property, using a periodic X_d built with np.mod
(no copy bookkeeping as in server.py); for triangle spectra (absolutely summable x[n]) Y_d is also
computed from the DEFINITION, sum_n x[n] m[n] e^{-j w n} with x[n] = (A W/2) sinc^2(W pi n/2) and the
modulating sequence m[n], truncated at |n| <= 6000. Points must not sit on a jump of a rectangle
(left and right limits agree), and the band edges are the centre of the copy +- W."""

import collections
import html
import math
import re
from fractions import Fraction

import numpy as np
from ece310 import svg


def fr(s):
    return float(Fraction(s))


def X_per(shape, A, W, w):
    """2pi-periodic base spectrum; w in radians (array)."""
    u = np.mod(np.asarray(w, float) + math.pi, 2 * math.pi) - math.pi
    inside = np.abs(u) <= W * math.pi + 1e-12
    if shape == "rect":
        return np.where(inside, float(A), 0.0)
    return np.where(inside, A * (1 - np.abs(u) / (W * math.pi)), 0.0)


def Y_np(params, w, variant="true"):
    shape, kind, A = params["shape"], params["kind"], params["A"]
    W, w0 = fr(params["W"]), fr(params["w0"]) * math.pi
    X = lambda v: X_per(shape, A, W, v)                         # noqa: E731
    Xn = lambda v: np.where(np.abs(v) <= W * math.pi + 1e-12,   # noqa: E731  non-periodic (no wrap-around)
                            X_per(shape, A, W, v), 0.0)
    if kind == "cos":
        if variant == "half":
            return X(w - w0) + X(w + w0)
        if variant == "nowrap":
            return 0.5 * Xn(w - w0) + 0.5 * Xn(w + w0)
        if variant == "onecopy":
            return 0.5 * X(w - w0)
        return 0.5 * X(w - w0) + 0.5 * X(w + w0)
    if kind == "cos2":
        if variant == "half":
            return 0.5 * X(w) + 0.5 * X(w - 2 * w0) + 0.5 * X(w + 2 * w0)
        if variant == "nodc":
            return 0.25 * X(w - 2 * w0) + 0.25 * X(w + 2 * w0)
        if variant == "nowrap":
            return 0.5 * Xn(w) + 0.25 * Xn(w - 2 * w0) + 0.25 * Xn(w + 2 * w0)
        return 0.5 * X(w) + 0.25 * X(w - 2 * w0) + 0.25 * X(w + 2 * w0)
    if variant == "half":
        return 0.5 * X(w - math.pi)
    if variant == "noshift":
        return X(w)
    return X(w - math.pi)


def dtft_direct(params, w):
    """Triangle only: Y_d(w) from the definition."""
    A, W = params["A"], fr(params["W"])
    n = np.arange(-6000, 6001)
    half = W * math.pi * n / 2
    safe = np.where(n == 0, 1.0, half)
    x = (A * W / 2) * np.where(n == 0, 1.0, (np.sin(safe) / safe) ** 2)
    w0 = fr(params["w0"]) * math.pi
    m = {"cos": np.cos(w0 * n), "cos2": np.cos(w0 * n) ** 2, "alt": np.cos(math.pi * n)}[params["kind"]]
    return float(np.real(np.sum(x * m * np.exp(-1j * w * n))))


def _width(label):
    """Estimated width in px at font-size 13 (sans-serif): digits 7.3, '/' 3.8, minus 7.6, pi and others 7.5."""
    return sum(7.3 if ch.isdigit() else 3.8 if ch == "/" else 7.6 if ch == "\u2212" else 7.5 for ch in label)


def label_overlaps(fig_html):
    """Text labels in one row (same y) whose estimated extents (plus a 3 px gap) collide."""
    rows = collections.defaultdict(list)
    for m in re.finditer(r'<text x="([\d.]+)" y="([\d.]+)" text-anchor="middle">([^<]*)</text>', fig_html):
        rows[m.group(2)].append((float(m.group(1)), html.unescape(m.group(3))))
    out = []
    for items in rows.values():
        items.sort()
        for (x1, l1), (x2, l2) in zip(items, items[1:]):
            need = (_width(l1) + _width(l2)) / 2 + 3
            if x2 - x1 < need:
                out.append(f"figure labels {l1!r} and {l2!r} are {x2 - x1:.1f} px apart (need {need:.1f})")
    return out


def check(params, correct):
    probs = []
    pts = [fr(s) for s in params["pts"]]
    if len(set(pts)) != 3 or not all(0 <= v <= 1 for v in pts):
        probs.append(f"bad points {pts}")
    for i, v in enumerate(pts, start=1):
        w = v * math.pi
        y = float(Y_np(params, w))
        if abs(y - correct[f"y{i}"]) > 1e-9:
            probs.append(f"Y_d({v} pi) = {y}, answer {correct[f'y{i}']}")
        lo, hi = float(Y_np(params, w - 1e-7)), float(Y_np(params, w + 1e-7))
        if params["shape"] == "rect" and abs(lo - hi) > 1e-9:
            probs.append(f"point {v} pi sits on a jump ({lo} vs {hi})")
        variants = {"cos": ["half", "nowrap", "onecopy"], "cos2": ["half", "nodc", "nowrap"],
                    "alt": ["half", "noshift"]}[params["kind"]]
        if not any(_differs(float(Y_np(params, w, var)), y) for var in variants):
            probs.append(f"no classic mistake differs from Y_d({v} pi) = {y}")
        if params["shape"] == "tri":
            d = dtft_direct(params, w)
            if abs(d - y) > 2e-3 * max(1.0, abs(y)):
                probs.append(f"direct DTFT sum {d} != {y} at {v} pi")
    W, w0 = fr(params["W"]), fr(params["w0"])
    centre = {"cos": w0, "cos2": 2 * w0, "alt": 1.0}[params["kind"]]
    if not 0 < centre <= 1:
        probs.append(f"copy centre {centre} not in (0, 1]")
    if abs(correct["wl"] - (centre - W)) > 1e-12 or abs(correct["wh"] - (centre + W)) > 1e-12:
        probs.append(f"edges ({correct['wl']}, {correct['wh']}) != ({centre - W}, {centre + W})")
    for lab in (svg.fmt_pi(Fraction(params["W"])), svg.fmt_pi(-Fraction(params["W"]))):
        if f">{lab}<" not in params["fig_html"]:
            probs.append(f"figure lacks the tick {lab}")
    if any(abs(correct[k] - 13 / 7) < 1e-3 or abs(correct[k] - 2 / 7) < 1e-3 for k in ("y1", "y2", "y3", "wl", "wh")):
        probs.append("an answer equals a format example (13/7 or 2/7)")
    probs += label_overlaps(params["fig_html"])
    if len({round(correct[f"y{i}"], 9) for i in (1, 2, 3)}) < 2:
        probs.append("all three requested values of Y_d are equal")
    if params["kind"] == "cos" and abs(w0 - W) < 1e-12:
        probs.append("w0 == W")
    return probs


def _differs(a, b):
    return abs(a - b) > 2 * (1e-6 + 1e-3 * abs(b))


def submissions(params, correct):
    cases = []
    kind = params["kind"]
    variants = {"cos": ["half", "nowrap", "onecopy"], "cos2": ["half", "nodc", "nowrap"], "alt": ["half", "noshift"]}[kind]
    for i, s in enumerate(params["pts"], start=1):
        y = correct[f"y{i}"]
        q = Fraction(y).limit_denominator(1000)
        cases.append(({f"y{i}": f"{y:.4g}"}, {f"y{i}": 1}))
        cases.append(({f"y{i}": f"{3 * q.numerator}/{3 * q.denominator}"}, {f"y{i}": 1}))
        for var in variants:
            wrong = float(Y_np(params, fr(s) * math.pi, var))
            if _differs(wrong, y):
                cases.append(({f"y{i}": f"{wrong:.6g}"}, {f"y{i}": 0}))
    wl, wh = correct["wl"], correct["wh"]
    W, w0 = fr(params["W"]), fr(params["w0"])
    cases.append(({"wl": f"{wl:.4g}", "wh": f"{wh:.4g}"}, {"wl": 1, "wh": 1}))
    c = (wl + wh) / 2
    probes = [(c - W / 2, c + W / 2)]                          # cutoff read as the full width
    if kind == "cos2":
        probes.append((w0 - W, w0 + W))                         # copy at w0 instead of 2 w0
    if wh > 1:
        probes.append((wl, wh - 2))                             # folded the upper edge
    for a, b in probes:
        exp = {"wl": 0 if _differs(a, wl) else 1, "wh": 0 if _differs(b, wh) else 1}
        if 0 in exp.values():
            cases.append(({"wl": f"{a:.6g}", "wh": f"{b:.6g}"}, exp))
    return cases
