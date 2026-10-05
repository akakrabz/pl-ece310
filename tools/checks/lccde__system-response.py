"""Independent check of questions/lccde/system-response: scipy.signal.lfilter on the input sequence
x[n] = a^n (n = 0..39) versus the closed form y[n] (SymPy) and the numbers y[0], y[1]."""

from fractions import Fraction

import numpy as np
from checklib import N as NSYM
from checklib import S
from scipy import signal


def _setup(params, L=40):
    b = [float(Fraction(v)) for v in params["b"]]
    a = [float(Fraction(v)) for v in params["a"]]
    ain = float(Fraction(params["a_in"]))
    x = ain ** np.arange(L)
    return b, a, ain, signal.lfilter(b, a, x)


def _closed(expr, L=40):
    f = S(expr)
    return np.array([complex(f.subs(NSYM, k).evalf(30)).real for k in range(L)])


def _exact(params, L=40):
    """Exact recursion from rest with Fractions (lfilter in floating point cannot follow a cancelled
    unstable mode for long: round-off excites it)."""
    b = [Fraction(v) for v in params["b"]]
    a = [Fraction(v) for v in params["a"]]
    ain = Fraction(params["a_in"])
    x = [ain ** k for k in range(L)]
    y = []
    for k in range(L):
        acc = sum(b[i] * x[k - i] for i in range(len(b)) if k - i >= 0)
        acc -= sum(a[i] * y[k - i] for i in range(1, len(a)) if k - i >= 0)
        y.append(acc / a[0])
    return y


def check(params, correct):
    probs = []
    b, a, ain, y = _setup(params, 12)
    yc = _closed(correct["y"])
    if not np.all(np.abs(yc[:12] - y) <= 1e-9 * np.maximum(1.0, np.abs(y))):
        probs.append(f"closed form {yc[:5]} != lfilter {y[:5]}")
    ye = _exact(params)
    if any(abs(float(Fraction(v)) - c) > 1e-9 * max(1.0, abs(float(v))) for v, c in zip(ye, yc)):
        probs.append("closed form differs from the exact recursion for some n < 40")
    if not np.isclose(correct["y0"], y[0], rtol=1e-12, atol=1e-12) or not np.isclose(correct["y1"], y[1], rtol=1e-12, atol=1e-12):
        probs.append(f"y0, y1 = {correct['y0']}, {correct['y1']} but lfilter gives {y[0]}, {y[1]}")
    # the advertised cancellation really happens
    zb = np.roots(b) if len(b) > 1 else np.array([])
    pa = np.roots(a)
    kind = params["kind"]
    if kind == "xcancel" and not np.any(np.abs(zb - ain) < 1e-9):
        probs.append("xcancel: H has no zero at the input pole")
    if kind == "hcancel" and not any(np.min(np.abs(pa - z)) < 1e-9 for z in zb):
        probs.append("hcancel: B and A share no root")
    if kind in ("plain1", "plain2") and (any(np.min(np.abs(pa - z)) < 1e-9 for z in zb) or np.any(np.abs(zb - ain) < 1e-9)):
        probs.append("plain variant with a cancellation")
    return probs


def _term(A, p, how):
    A, p = Fraction(A), Fraction(p)
    if how == "shift":      # A p^n = (A p) p^(n-1)
        return f"({A * p})*({p})^(n-1)"
    if how == "inverse":    # p^n = (1/p)^(-n)
        return f"({A})*({1 / p})^(-n)"
    return f"({A})*({p})^n"


def submissions(params, correct):
    poles = [Fraction(v) for v in params["poles_y"]]
    res = [Fraction(v) for v in params["res_y"]]
    ain = Fraction(params["a_in"])
    cases = []
    eq1 = " + ".join(_term(A, p, "shift") for A, p in zip(res, poles))
    eq2 = " + ".join(_term(A, p, "inverse" if p != 1 else "plain") for A, p in zip(res, poles))
    cases.append(({"y": eq1}, {"y": 1}))
    cases.append(({"y": eq2}, {"y": 1}))
    # mistakes: residues swapped between the poles; every residue with the wrong sign; input mode dropped
    if len(set(res)) > 1:
        sw = res[1:] + res[:1]
        cases.append(({"y": " + ".join(_term(A, p, "plain") for A, p in zip(sw, poles))}, {"y": 0}))
    cases.append(({"y": " + ".join(_term(-A, p, "plain") for A, p in zip(res, poles))}, {"y": 0}))
    if ain in poles and len(poles) > 1:
        keep = [(A, p) for A, p in zip(res, poles) if p != ain]
        cases.append(({"y": " + ".join(_term(A, p, "plain") for A, p in keep)}, {"y": 0}))
    cases.append(({"y": "0.5^n"}, {"y": "invalid"}))
    # expressions that evaluate to NaN / infinity must never pass the integer-sample grader
    for bad in ("0/0", "zoo", "1/(n-n)", "oo-oo"):
        cases.append(({"y": bad}, {"y": 0}))
    # numbers: fractions and decimals; sign error on the feedback terms for y[1]
    y0, y1 = Fraction(correct["y0"]).limit_denominator(10**6), Fraction(correct["y1"]).limit_denominator(10**6)
    cases.append(({"y0": f"{y0.numerator}/{y0.denominator}", "y1": repr(float(y1))}, {"y0": 1, "y1": 1}))
    # 4-significant-digit decimals score 1 (course rule: rtol 1e-3)
    cases.append(({"y0": f"{float(y0):.4g}", "y1": f"{float(y1):.4g}"}, {"y0": 1, "y1": 1}))
    b = [float(Fraction(v)) for v in params["b"]]
    a = [float(Fraction(v)) for v in params["a"]]
    a_bad = [a[0]] + [-v for v in a[1:]]
    yb = signal.lfilter(b, a_bad, float(ain) ** np.arange(3))
    if not np.isclose(yb[1], float(y1)):
        cases.append(({"y1": repr(float(yb[1]))}, {"y1": 0}))
    # the impulse response value instead of the output value
    h = signal.lfilter(b, a, np.r_[1.0, np.zeros(2)])
    if not np.isclose(h[1], float(y1)):
        cases.append(({"y1": repr(float(h[1]))}, {"y1": 0}))
    return cases
