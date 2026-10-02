"""Independent check of questions/ztransform/all-possible-rocs.

SymPy cancels the displayed H(z) (so a cancelled pole is detected independently), the poles come
from the reduced denominator, residues from scipy.signal.residuez, and h[n] for the stated ROC is
rebuilt from them and verified against H(z) by a truncated z-transform sum inside that ROC."""

import math
import re
from fractions import Fraction

import numpy as np
import sympy
from checklib import signal, zsum

INF = math.inf
W = sympy.Symbol("w")


def parse_roc(html):
    s = html.replace("$", "").strip()
    s = re.sub(r"\\tfrac\{(\d+)\}\{(\d+)\}", r"\1/\2", s)
    z = r"\lvert z\rvert"
    if s.startswith(z + r" \gt "):
        return float(Fraction(s[len(z + r" \gt "):])), INF
    if s.startswith(z + r" \lt "):
        return 0.0, float(Fraction(s[len(z + r" \lt "):]))
    lo, mid, hi = s.split(r" \lt ")
    assert mid == z, s
    return float(Fraction(lo)), float(Fraction(hi))


def reduced(params):
    """(b, a, poles) of the cancelled H(z) in powers of z^-1 (a[0] = 1)."""
    num = sum(sympy.Rational(c) * W**k for k, c in enumerate(params["shown_num"]))
    den = sympy.Mul(*[1 - sympy.Rational(q) * W for q in params["den_poles"]])
    n_red, d_red = sympy.fraction(sympy.cancel(num / den))
    a = [sympy.Rational(x) for x in sympy.Poly(d_red, W).all_coeffs()[::-1]]
    b = [sympy.Rational(x) for x in sympy.Poly(n_red, W).all_coeffs()[::-1]]
    b = [x / a[0] for x in b]
    a = [x / a[0] for x in a]
    poles = [1 / r for r in sympy.Poly(d_red, W).all_roots()]
    return [float(x) for x in b], [float(x) for x in a], [float(q) for q in poles]


def true_rocs(poles):
    mags = sorted({float(Fraction(abs(q)).limit_denominator(1000)) for q in poles})
    return [(0.0, mags[0])] + list(zip(mags, mags[1:])) + [(mags[-1], INF)]


def h_values(params):
    """h[n] for the stated ROC, from scipy residues and our own side rule."""
    b, a, poles = reduced(params)
    r, p, k = signal.residuez(b, a)
    lo = float(Fraction(params["spec_inner"]))
    hi = INF if params["spec_outer"] is None else float(Fraction(params["spec_outer"]))

    def h(n):
        v = 0.0
        for ri, pi in zip(r, p):
            if abs(pi) <= lo + 1e-9 and n >= 0:
                v += (ri * pi**n).real
            elif abs(pi) >= hi - 1e-9 and n <= -1:
                v -= (ri * pi**n).real
        return v
    return h, b, a, poles, (lo, hi), (r, p, k)


def differs(wrong, right):
    """True if pl-number-input (relabs, rtol 1e-3, atol 1e-6) would mark `wrong` incorrect, with margin."""
    return abs(float(wrong) - float(right)) > 2 * (1e-6 + 1e-3 * abs(float(right)))


def sig4(x):
    return f"{float(x):.4g}"


def cancelled(params, poles):
    """Displayed denominator poles that are not poles of the reduced H (found by SymPy)."""
    shown = [float(Fraction(q)) for q in params["den_poles"]]
    return [q for q in shown if min(abs(q - t) for t in poles) > 1e-9]


def check(params, correct):
    probs = []
    h, b, a, poles, (lo, hi), (r, p, k) = h_values(params)
    for c in cancelled(params, poles):
        if lo < abs(c) < hi:
            probs.append(f"part (c) ROC {(lo, hi)} crosses the cancelled pole |{c}|")
    nopt = len(params["stab"])
    if nopt <= len(true_rocs(poles)):
        probs.append(f"{nopt} options for {len(true_rocs(poles))} possible ROCs: counting options answers (a)")
    if all(parse_roc(o["html"]) in true_rocs(poles) for o in params["stab"]):
        probs.append("no invalid region among the options")
    if len(k) and np.any(np.abs(k) > 1e-9):
        probs.append(f"reduced H is not proper: direct terms {k}")
    if any(abs(abs(q) - 1) < 1e-9 for q in poles):
        probs.append("pole on the unit circle")
    rocs = true_rocs(poles)
    if correct["nroc"] != len(rocs):
        probs.append(f"nroc {correct['nroc']} != {len(rocs)} (poles {poles})")
    if len({round(abs(q), 9) for q in poles}) != len(poles):
        probs.append("pole magnitudes not distinct")
    stable = [rr for rr in rocs if rr[0] < 1 < rr[1]]
    if parse_roc(correct["stab"]["html"]) != stable[0]:
        probs.append(f"marked stable ROC {correct['stab']['html']} != {stable[0]}")
    shown = [parse_roc(o["html"]) for o in params["stab"]]
    if not all(rr in shown for rr in rocs):
        probs.append("a possible ROC is missing from the options")
    if len(set(shown)) != len(shown):
        probs.append("duplicate options")
    if (lo, hi) not in rocs or (lo, hi) == stable[0]:
        probs.append(f"part (c) ROC {(lo, hi)} invalid or equal to the stable one")
    for name, n in (("hm2", -2), ("h0", 0), ("h2", 2)):
        if abs(h(n) - correct[name]) > 1e-8 * max(1, abs(h(n))):
            probs.append(f"{name} = {correct[name]} but independent h[{n}] = {h(n)}")
    # our own h must really be the inverse of H in that ROC
    rad = (1.5 * lo + 0.3) if hi == INF else (0.6 * hi if lo == 0 else math.sqrt(lo * hi))
    for z0 in (rad, rad * complex(math.cos(1.1), math.sin(1.1))):
        Hz = np.polyval(b[::-1], 1 / z0) / np.polyval(a[::-1], 1 / z0)
        got = zsum(h, z0, -400, 400)
        if abs(got - Hz) > 1e-8 * max(1, abs(Hz)):
            probs.append(f"independent h does not invert H at z = {z0}: {got} vs {Hz}")
    return probs


def frac_alt(x):
    q = Fraction(x).limit_denominator(10000)
    d = q.denominator
    while d % 2 == 0:
        d //= 2
    while d % 5 == 0:
        d //= 5
    if d == 1 and q.denominator != 1:
        return repr(float(q))
    return f"{3 * q.numerator}/{3 * q.denominator}"


def plain(x):
    q = Fraction(x).limit_denominator(10000)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def submissions(params, correct):
    h, b, a, poles, (lo, hi), (r, p, k) = h_values(params)
    m = len(poles)
    cases = [
        ({"nroc": str(2 ** m)}, {"nroc": 0}),                       # 2^(number of poles)
        ({"nroc": str(m)}, {"nroc": 0}),                            # number of poles
        ({"nroc": "3.5"}, {"nroc": "invalid"}),
        ({"hm2": frac_alt(correct["hm2"]), "h0": frac_alt(correct["h0"]), "h2": frac_alt(correct["h2"])},
         {"hm2": 1, "h0": 1, "h2": 1}),                             # decimal / non-reduced fractions
        ({"hm2": plain(correct["hm2"]), "h0": plain(correct["h0"]), "h2": plain(correct["h2"])},
         {"hm2": 1, "h0": 1, "h2": 1}),                             # reduced fractions
    ]
    dec = {nm: sig4(correct[nm]) for nm in ("hm2", "h0", "h2")
           if Fraction(correct[nm]).limit_denominator(10000).denominator != 1}
    if dec:                                                         # 4 significant digits
        cases.append((dec, {nm: 1 for nm in dec}))
    if len(params["den_poles"]) > m:                                # counted the cancelled pole too
        cases.append(({"nroc": str(m + 2)}, {"nroc": 0}))
    # stability: picking the causal ROC (or a bogus ring through the cancelled pole)
    stable = [rr for rr in true_rocs(poles) if rr[0] < 1 < rr[1]][0]
    for o in params["stab"]:
        rr = parse_roc(o["html"])
        if rr != stable and (rr[1] == INF or rr[0] < 1 < rr[1]):
            cases.append(({"stab": o["key"]}, {"stab": 0}))
    # h mistakes
    left = [(ri, pi) for ri, pi in zip(r, p) if abs(pi) >= hi - 1e-9]
    if left:
        wrong_sign = sum((ri * pi**-2).real for ri, pi in left)          # forgot -A p^n u[-n-1]'s minus
        if differs(wrong_sign, correct["hm2"]):
            cases.append(({"hm2": plain(wrong_sign)}, {"hm2": 0}))
        all_sum = sum(ri.real for ri in r)                              # h[0] of the causal system
        if differs(all_sum, correct["h0"]):
            cases.append(({"h0": plain(all_sum)}, {"h0": 0}))
        causal_h2 = sum((ri * pi**2).real for ri, pi in zip(r, p))
        if differs(causal_h2, correct["h2"]):
            cases.append(({"h2": plain(causal_h2)}, {"h2": 0}))
    else:                                                               # causal ROC asked: anti-causal value
        anti_m2 = -sum((ri * pi**-2).real for ri, pi in zip(r, p))
        if differs(anti_m2, correct["hm2"]):
            cases.append(({"hm2": plain(anti_m2)}, {"hm2": 0}))
    return cases
