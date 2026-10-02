"""Independent check of questions/ztransform/complex-poles.

Poles from numpy.roots, samples from scipy.signal.lfilter (causal recursion), and the closed form
(the correct pl-symbolic-input answer) evaluated numerically against lfilter for n = 0..39."""

import cmath
import math
from fractions import Fraction

import numpy as np
import sympy
from checklib import N, impulse_response

TH = {"pi/2": math.pi / 2, "pi/3": math.pi / 3, "2*pi/3": 2 * math.pi / 3}
THS = {"pi/2": "pi*n/2", "pi/3": "pi*n/3", "2*pi/3": "2*pi*n/3"}
ALIAS = {"pi/2": "3*pi*n/2", "pi/3": "5*pi*n/3", "2*pi/3": "4*pi*n/3"}     # (2*pi - theta) n


def differs(wrong, right):
    """True if pl-number-input (relabs, rtol 1e-3, atol 1e-6) would mark `wrong` incorrect, with margin."""
    return abs(float(wrong) - float(right)) > 2 * (1e-6 + 1e-3 * abs(float(right)))


def sig4(x):
    return f"{float(x):.4g}"


def _ba(params):
    return [float(Fraction(x)) for x in params["b"]], [float(Fraction(x)) for x in params["a"]]


def closed_fn(expr_json):
    e = sympy.sympify(expr_json["_value"], locals={"n": N})
    return sympy.lambdify(N, e, "math")


def check(params, correct):
    probs = []
    b, a = _ba(params)
    r, th = float(Fraction(params["r"])), TH[params["theta"]]
    roots = np.roots(a)
    want = [r * cmath.exp(1j * th), r * cmath.exp(-1j * th)]
    if not all(min(abs(x - w) for x in roots) < 1e-9 for w in want):
        probs.append(f"poles {roots} are not {r} e^(+-j{th})")
    h = impulse_response(b, a, 40)
    for n in range(3):
        if abs(h[n] - correct[f"h{n}"]) > 1e-9 * max(1, abs(h[n])):
            probs.append(f"h[{n}] = {correct[f'h{n}']} but lfilter gives {h[n]}")
    f = closed_fn(correct["hn"])
    scale = abs(float(Fraction(params["A"]))) + 2 * abs(float(Fraction(params["k"]))) + 1
    for n in range(40):
        v = float(f(n))
        if abs(v - h[n]) > 1e-9 * scale * max(1.0, r**n):
            probs.append(f"closed form at n={n}: {v} vs lfilter {h[n]}")
            break
    return probs


def plain(x):
    q = Fraction(x).limit_denominator(10000)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def alt(x):
    q = Fraction(x).limit_denominator(10000)
    d = q.denominator
    while d % 2 == 0:
        d //= 2
    while d % 5 == 0:
        d //= 5
    if d == 1 and q.denominator != 1:
        return repr(float(q))
    return f"{2 * q.numerator}/{2 * q.denominator}"


def rpow(r):
    r = Fraction(r)
    if r == 1:
        return "1"
    return f"({plain(r)})^n" if r.denominator != 1 else f"{r.numerator}^n"


def submissions(params, correct):
    b, a = _ba(params)
    r = Fraction(params["r"])
    A, k, sq = Fraction(params["A"]), Fraction(params["k"]), params["sqrt3"]
    ths = THS[params["theta"]]
    Bs = f"({plain(k)})*sqrt(3)" if sq else f"({plain(k)})"          # B as typed
    B2 = 3 * k * k if sq else k * k
    M2 = A * A + B2
    rp = rpow(r)
    cases = [
        ({"h0": alt(correct["h0"]), "h1": alt(correct["h1"]), "h2": alt(correct["h2"])}, {"h0": 1, "h1": 1, "h2": 1}),
        ({"h0": plain(correct["h0"]), "h1": plain(correct["h1"]), "h2": plain(correct["h2"])}, {"h0": 1, "h1": 1, "h2": 1}),
        ({"hn": f"0.5^n*cos({ths})"}, {"hn": "invalid"}),
        ({"hn": "0/0"}, {"hn": 0}),                                        # NaN must not earn credit
        ({"hn": "zoo"}, {"hn": 0}),                                        # complex infinity neither
        # aliased angle (2*pi - theta): cos(.) unchanged, sin(.) flips sign, equal for every integer n
        ({"hn": f"{rp}*(({plain(A)})*cos({ALIAS[params['theta']]}) - {Bs}*sin({ALIAS[params['theta']]}))"}, {"hn": 1}),
    ]
    dec = {nm: sig4(correct[nm]) for nm in ("h0", "h1", "h2")
           if Fraction(correct[nm]).limit_denominator(10000).denominator != 1}
    if dec:                                                                # 4 significant digits
        cases.append((dec, {nm: 1 for nm in dec}))
    # amplitude-phase form  M r^n cos(theta n - phi)
    if A == 0:
        ap = f"{Bs}*{rp}*cos({ths} - pi/2)"
    else:
        ratio = f"({plain(k / A)})*sqrt(3)" if sq else f"({plain(k / A)})"
        ap = f"{'-' if A < 0 else ''}sqrt({plain(M2)})*{rp}*cos({ths} - atan({ratio}))"
    cases.append(({"hn": ap}, {"hn": 1}))
    # r^n split into num^n/den^n, sqrt(3) factors rewritten, terms reversed
    num_pow = f"{r.numerator}^n" if r.numerator != 1 else "1"
    den_pow = f"/{r.denominator}^n" if r.denominator != 1 else ""
    Bs2 = f"({plain(3 * k)})/sqrt(3)" if sq else f"({plain(k)})"
    cases.append(({"hn": f"({Bs2}*sin({ths}) + ({plain(A)})*cos({ths}))*{num_pow}{den_pow}"}, {"hn": 1}))
    # mistakes
    Bval = float(k) * (math.sqrt(3) if sq else 1.0)
    if abs(Bval - float(A)) > 1e-9:                                       # cos and sin parts swapped
        cases.append(({"hn": f"{rp}*({Bs}*cos({ths}) + ({plain(A)})*sin({ths}))"}, {"hn": 0}))
    if r != 1:                                                            # forgot r^n
        cases.append(({"hn": f"({plain(A)})*cos({ths}) + {Bs}*sin({ths})"}, {"hn": 0}))
    if differs(b[1], correct["h1"]):                                      # h[1] = b1 (ignored the feedback)
        cases.append(({"h1": plain(b[1])}, {"h1": 0}))
    wrong2 = a[1] * correct["h1"] - a[2] * correct["h0"]                   # recursion with the wrong sign of a1
    if differs(wrong2, correct["h2"]):
        cases.append(({"h2": plain(wrong2)}, {"h2": 0}))
    return cases
