"""Independent check of questions/ztransform/improper-long-division.

Direct terms and residues from scipy.signal.residuez (not poly.pdivmod), h[0..3] from
scipy.signal.lfilter (the causal recursion), and the decomposition re-evaluated at test points."""

from fractions import Fraction

import numpy as np
from checklib import impulse_response, signal


def _ba(params):
    b = [float(Fraction(x)) for x in params["num"]]
    poles = [float(Fraction(q)) for q in params["poles"]]
    a = np.poly(poles)
    return b, a, poles


def check(params, correct):
    probs = []
    b, a, poles = _ba(params)
    r, p, k = signal.residuez(b, a)
    L = params["L"]
    want_C = [correct["C0"]] + ([correct["C1"]] if L == 2 else [])
    if len(k) != L or not np.allclose(k, want_C, atol=1e-8):
        probs.append(f"direct terms {k} != {want_C}")
    names = ["A1", "A2"][: len(poles)]
    for name, q in zip(names, poles):
        i = int(np.argmin(np.abs(p - q)))
        if abs(p[i] - q) > 1e-7 or abs(r[i] - correct[name]) > 1e-7:
            probs.append(f"{name} = {correct[name]} but residuez gives {r[i]} at {p[i]}")
    h = impulse_response(b, a, 4)
    for n in range(4):
        if abs(h[n] - correct[f"h{n}"]) > 1e-9 * max(1, abs(h[n])):
            probs.append(f"h[{n}] = {correct[f'h{n}']} but lfilter gives {h[n]}")
    if len(b) < len(a):
        probs.append("H is proper")
    for z0 in (1.7 + 2.2j, -3.1 + 0.4j):     # decomposition reproduces H(z)
        H = np.polyval(b[::-1], 1 / z0) / np.polyval(a[::-1], 1 / z0)
        dec = sum(c * z0 ** (-i) for i, c in enumerate(want_C))
        dec += sum(correct[nm] / (1 - q / z0) for nm, q in zip(names, poles))
        if abs(H - dec) > 1e-9 * max(1, abs(H)):
            probs.append(f"decomposition differs from H at {z0}")
    return probs


def plain(x):
    q = Fraction(x).limit_denominator(10000)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def alt(x):
    """decimal if it terminates, else a non-reduced fraction."""
    q = Fraction(x).limit_denominator(10000)
    d = q.denominator
    while d % 2 == 0:
        d //= 2
    while d % 5 == 0:
        d //= 5
    if d == 1 and q.denominator != 1:
        return repr(float(q))
    return f"{2 * q.numerator}/{2 * q.denominator}"


def differs(wrong, right):
    """True if pl-number-input (relabs, rtol 1e-3, atol 1e-6) would mark `wrong` incorrect, with margin."""
    return abs(float(wrong) - float(right)) > 2 * (1e-6 + 1e-3 * abs(float(right)))


def sig4(x):
    return f"{float(x):.4g}"


def submissions(params, correct):
    b, a, poles = _ba(params)
    r, p, k = signal.residuez(b, a)
    names = [nm for nm in ("C0", "C1", "A1", "A2", "h0", "h1", "h2", "h3") if nm in correct]
    cases = [
        ({nm: alt(correct[nm]) for nm in names}, {nm: 1 for nm in names}),     # 0.75, 6/8 ...
        ({nm: plain(correct[nm]) for nm in names}, {nm: 1 for nm in names}),   # reduced fractions
        ({"h1": "1.2.3"}, {"h1": "invalid"}),
    ]
    dec = {nm: sig4(correct[nm]) for nm in names if Fraction(correct[nm]).limit_denominator(10000).denominator != 1}
    if dec:                                       # 4 significant digits
        cases.append((dec, {nm: 1 for nm in dec}))
    sumA = float(np.sum(r).real)
    if differs(b[0], correct["C0"]):              # "divided from the front": C0 = b0 = h[0]
        cases.append(({"C0": plain(b[0])}, {"C0": 0}))
    if differs(sumA, correct["h0"]):              # forgot the delta terms in h[0]
        cases.append(({"h0": plain(sumA)}, {"h0": 0}))
    if "C1" in correct:                           # forgot C1 delta[n-1] in h[1]
        exp1 = sum((ri * pi).real for ri, pi in zip(r, p))
        if differs(exp1, correct["h1"]):
            cases.append(({"h1": plain(exp1)}, {"h1": 0}))
    if "A2" in correct and differs(correct["A2"], correct["A1"]) and differs(correct["A1"], correct["A2"]):
        cases.append(({"A1": plain(correct["A2"]), "A2": plain(correct["A1"])}, {"A1": 0, "A2": 0}))
    return cases
