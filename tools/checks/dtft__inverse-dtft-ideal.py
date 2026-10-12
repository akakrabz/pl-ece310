"""Independent check of questions/dtft/inverse-dtft-ideal.

x[n] is recomputed from the DEFINITION of the inverse DTFT, integrating G(w) e^{-j w n0} e^{j w n} piece by
piece over the plotted levels (exact antiderivatives of the exponential, numpy), not from the sinc
formula; the generated answer strings are parsed with SymPy and compared with it at n0 +- 1..40; the
rectangle decomposition is compared with the levels on a frequency grid; the figure must carry a
tick label at every cutoff."""

import collections
import html
import math
import re
from fractions import Fraction

import numpy as np
from checklib import sympy
from ece310 import svg

N = sympy.Symbol("n")


def _levels(params):
    return [(float(Fraction(a)), float(Fraction(b)), float(Fraction(h))) for a, b, h in params["levels"]]


def inv_dtft(params, n):
    """(1/2pi) int_{-pi}^{pi} G(w) e^{-j w n0} e^{j w n} dw, G piecewise constant (levels in units of pi)."""
    k = n - params["n0"]
    tot = 0j
    for a, b, h in _levels(params):
        if k == 0:
            tot += h * (b - a) * math.pi
        else:
            tot += h * (np.exp(1j * b * math.pi * k) - np.exp(1j * a * math.pi * k)) / (1j * k)
    return tot / (2 * math.pi)


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


def _eval(expr_str, n):
    loc = {"n": N, "pi": sympy.pi, "sinc": sympy.sinc}
    e = sympy.sympify(expr_str.replace("^", "**"), locals=loc)
    return complex(e.subs(N, n).evalf(30))


def _differs(a, b):
    return abs(float(a) - float(b)) > 2 * (1e-6 + 1e-3 * abs(float(b)))


def check(params, correct):
    probs = []
    n0, m = params["n0"], params["m"]
    ks = [k for k in range(-40, 41) if k != 0]
    for k in ks:
        t = inv_dtft(params, n0 + k)
        if abs(t.imag) > 1e-9:
            probs.append(f"x[{n0 + k}] = {t} is not real")
            break
        got = _eval(params["x_sym"], n0 + k)
        if abs(got - t.real) > 1e-9:
            probs.append(f"x_sym gives {got} at n = {n0 + k}, inverse DTFT {t.real}")
            break
        if params.get("mod_sym") and abs(_eval(params["mod_sym"], n0 + k) - t.real) > 1e-9:
            probs.append(f"modulated form wrong at n = {n0 + k}")
            break
    x0 = inv_dtft(params, n0).real
    if abs(x0 - correct["x0"]) > 1e-12:
        probs.append(f"x[n0] = {x0}, answer {correct['x0']}")
    xm = inv_dtft(params, n0 + m).real
    if abs(xm - correct["xm"]) > 1e-12:
        probs.append(f"x[n0+m] = {xm}, answer {correct['xm']}")
    if abs(abs(xm) - 1.414) < 0.01 or abs(x0 - 13 / 7) < 1e-3:
        probs.append("an answer equals a format example (1.414 or 13/7)")
    if abs(xm) < 1e-3:
        probs.append("x[n0+m] is (nearly) zero")
    if params["n1"] != n0 + m:
        probs.append("n1 != n0 + m")
    if Fraction(correct["x0"]).limit_denominator(100) != Fraction(sum(
            Fraction(h) * (Fraction(b) - Fraction(a)) for a, b, h in params["levels"]) / 2):
        probs.append("x[n0] is not the exact average height")
    # decomposition into centred rectangles agrees with the plotted levels
    w = np.linspace(-0.999, 0.999, 4001)
    lev = np.zeros_like(w)
    for a, b, h in _levels(params):
        lev += np.where((w > a) & (w < b), h, 0.0)
    rec = np.zeros_like(w)
    for h, W in params["rects"]:
        rec += float(Fraction(h)) * (np.abs(w) < float(Fraction(W)))
    edges = sorted({abs(a) for a, _, _ in _levels(params)} | {abs(b) for _, b, _ in _levels(params)})
    mask = np.all([np.abs(np.abs(w) - e) > 2e-3 for e in edges], axis=0)
    if np.max(np.abs(lev - rec)[mask]) > 1e-12:
        probs.append("rectangle decomposition does not match the plotted levels")
    probs += label_overlaps(params["fig_html"])
    for e in edges:
        if 0 < e < 1:
            for lab in (svg.fmt_pi(Fraction(e).limit_denominator(12)), svg.fmt_pi(-Fraction(e).limit_denominator(12))):
                if f">{lab}<" not in params["fig_html"]:
                    probs.append(f"figure lacks the tick label {lab}")
    return probs


# ----------------------------------------------------------------------------- student-style strings
def _arg(n0):
    return "n" if n0 == 0 else f"(n{'-' if n0 > 0 else '+'}{abs(n0)})"


def _rects(params):
    return [(Fraction(h), Fraction(W)) for h, W in params["rects"]]


def _sinc_sum(rects, n0, mode="std", den="pi"):
    a = _arg(n0)
    parts = []
    for h, W in rects:
        if W == 1:
            continue
        if mode == "sinc":            # h W sinc(W pi n)
            parts.append(f"({h * W})*sinc({W}*pi*{a})")
        elif mode == "divx":          # h W sin(W pi n)/(W pi n)
            parts.append(f"({h * W})*sin({W}*pi*{a})/({W}*pi*{a})")
        else:
            parts.append(f"({h})*sin(pi*{a}*{W})/({den}*{a})")
    return " + ".join(parts)


def submissions(params, correct):
    n0, fam = params["n0"], params["family"]
    rects = _rects(params)
    cases = [
        ({"x": _sinc_sum(rects, n0, "sinc")}, {"x": 1}),            # course sinc notation
        ({"x": _sinc_sum(rects, n0, "divx")}, {"x": 1}),            # (W/pi)-weighted sin(x)/x form
        ({"x": _sinc_sum(rects, n0, den="1")}, {"x": 0}),           # forgot the pi in pi n
        ({"x": "sin(pi*n/4)/(pi*m)"}, {"x": "invalid"}),
    ]
    if fam == "bp":
        cases.append(({"x": params["mod_sym"]}, {"x": 1}))           # 2 cos(wm n) x low-pass
        cases.append(({"x": params["mod_sym"].replace("2*", "", 1)}, {"x": 0}))   # forgot the 2
    if fam == "hp":
        A, W = rects[0][0], 1 - rects[1][1]
        a = _arg(n0)
        cases.append(({"x": f"(-1)^({a})*({A})*sin({W}*pi*{a})/(pi*{a})"}, {"x": 1}))   # band centred at pi
        cases.append(({"x": f"({A})*sin({rects[1][1]}*pi*{a})/(pi*{a})"}, {"x": 0}))    # no minus sign
    if fam == "stair":
        (d, wa), (A2, wb) = rects
        A1 = d + A2
        a = _arg(n0)
        cases.append(({"x": f"({A1})*sin({wa}*pi*{a})/(pi*{a}) + ({A2})*sin({wb}*pi*{a})/(pi*{a})"}, {"x": 0}))
    if any(abs(h) != 1 for h, W in rects if W != 1):                 # heights dropped
        cases.append(({"x": _sinc_sum([(1 if h > 0 else -1, W) for h, W in rects], n0)}, {"x": 0}))
    if n0 != 0:                                                      # delay with the wrong sign
        cases.append(({"x": _sinc_sum(rects, -n0)}, {"x": 0}))
        cases.append(({"x": _sinc_sum(rects, 0)}, {"x": 0}))         # ignored the phase factor
    # x[n0]
    x0 = Fraction(correct["x0"]).limit_denominator(100)
    cases.append(({"x0": f"{float(x0):.4g}"}, {"x0": 1}))
    cases.append(({"x0": f"{2 * x0.numerator}/{2 * x0.denominator}"}, {"x0": 1}))
    for wrong in {float(x0) * math.pi, float(max(abs(h) for h, _ in rects)),
                  float(abs(sum(h * W for h, W in rects if W != 1)))}:          # area/2 not /2pi; peak; all-pass forgotten
        if _differs(wrong, x0):
            cases.append(({"x0": f"{wrong:.6g}"}, {"x0": 0}))
    # x[n0 + m]
    xm = correct["xm"]
    cases += [({"xm": f"{xm:.4g}"}, {"xm": 1}), ({"xm": f"{xm:.7f}"}, {"xm": 1}),
              ({"xm": f"{-xm:.6g}"}, {"xm": 0}), ({"xm": f"{xm * math.pi:.6g}"}, {"xm": 0})]
    if n0 != 0:
        g_undelayed = inv_dtft(dict(params, n0=0), n0 + params["m"]).real
        if _differs(g_undelayed, xm):
            cases.append(({"xm": f"{g_undelayed:.6g}"}, {"xm": 0}))
    return cases
