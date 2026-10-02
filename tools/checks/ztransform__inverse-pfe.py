"""Independent check of questions/ztransform/inverse-pfe.

A1, A2 from scipy.signal.residuez; the ROC from the stated property with our own rule; h[n] from
the correct pl-symbolic-input answers, verified by summing sum_n h[n] z^-n at points inside the
expected ROC and comparing with H(z) there (so sidedness AND signs are checked)."""

import math
import re
from fractions import Fraction

import numpy as np
import sympy
from checklib import N, signal, zsum

INF = math.inf


def _poles(params):
    return [Fraction(p) for p in params["poles"]]


def parse_roc(html, mags=None):
    """'$\\tfrac{1}{2} \\lt \\lvert z\\rvert \\lt 2$' -> (0.5, 2.0); exterior -> (r, inf).
    mags = (|p1|, |p2|) resolves the symbolic options of the multiplied-out variants."""
    s = html.replace("$", "").strip()
    if mags is not None:
        s = s.replace(r"\lvert p_1\rvert", str(mags[0])).replace(r"\lvert p_2\rvert", str(mags[1]))
    s = re.sub(r"\\tfrac\{(\d+)\}\{(\d+)\}", r"\1/\2", s)
    z = r"\lvert z\rvert"
    if s.startswith(z + r" \gt "):
        return float(Fraction(s[len(z + r" \gt "):])), INF
    if s.startswith(z + r" \lt "):
        return 0.0, float(Fraction(s[len(z + r" \lt "):]))
    lo, mid, hi = s.split(r" \lt ")
    assert mid == z, s
    return float(Fraction(lo)), float(Fraction(hi))


def expected_roc(poles, prop):
    m1, m2 = sorted(abs(float(q)) for q in poles)
    if prop == "causal":
        return m2, INF
    if prop == "left":
        return 0.0, m1
    for lo, hi in ((0.0, m1), (m1, m2), (m2, INF)):
        if lo < 1 < hi:
            return lo, hi
    raise AssertionError("no stable ROC")


def test_radius(lo, hi):
    if hi == INF:
        return 1.5 * lo + 0.3
    if lo == 0:
        return 0.6 * hi
    return math.sqrt(lo * hi)


def seq_fn(expr_json):
    e = sympy.sympify(expr_json["_value"], locals={"n": N})
    return sympy.lambdify(N, e, "math")


def mags_of(poles):
    """(|p1|, |p2|) as exact fractions, smaller first (independent of server's labelling)."""
    m = sorted(abs(Fraction(q)) for q in poles)
    return m[0], m[1]


def differs(wrong, right):
    """True if pl-number-input (relabs, rtol 1e-3, atol 1e-6) would mark `wrong` incorrect, with margin."""
    return abs(float(wrong) - float(right)) > 2 * (1e-6 + 1e-3 * abs(float(right)))


def sig4(x):
    """4-significant-digit decimal, as a student would type it."""
    return f"{float(x):.4g}"


def residues(params):
    """{pole: residue} from scipy (independent of the cover-up code in server.py)."""
    b = [float(Fraction(c)) for c in params["num"]]
    a = np.poly([float(q) for q in _poles(params)])
    r, p, k = signal.residuez(b, a)
    return r, p, k


def check(params, correct):
    probs = []
    poles = _poles(params)
    r, p, k = residues(params)
    if len(k) and np.any(np.abs(k) > 1e-9):
        probs.append(f"unexpected direct terms {k}")
    for name, q in (("A1", poles[0]), ("A2", poles[1])):
        i = int(np.argmin(np.abs(p - float(q))))
        if abs(r[i] - correct[name]) > 1e-7 or abs(p[i] - float(q)) > 1e-7:
            probs.append(f"{name} = {correct[name]} but residuez gives {r[i]} at pole {p[i]}")
    if not abs(poles[0]) < abs(poles[1]):
        probs.append("p1 must have the smaller magnitude")
    if any(abs(q) == 1 for q in poles):
        probs.append("pole on the unit circle")
    if correct["A1"] == correct["A2"]:
        probs.append("A1 == A2 (swap probe meaningless)")
    if params["multiplied"]:
        roots = sorted(np.roots(np.poly([float(q) for q in poles])), key=abs)
        if not (np.allclose([correct["p1"], correct["p2"]], np.real(roots)) and np.allclose(np.imag(roots), 0)):
            probs.append(f"poles {correct.get('p1')}, {correct.get('p2')} vs roots {roots}")
    lo, hi = expected_roc(poles, params["prop"])
    opts = params["roc"]
    mags = mags_of(poles)
    if len(opts) != 3:
        probs.append(f"{len(opts)} ROC options, expected 3")
    if params["multiplied"]:
        for o in opts:   # the options must not give the pole magnitudes away
            if re.search(r"\d", o["html"].replace("p_1", "").replace("p_2", "")):
                probs.append(f"multiplied-out variant shows a number in ROC option {o['html']}")
    if sorted(parse_roc(o["html"], mags) for o in opts) != sorted(
            [(0.0, float(min(abs(q) for q in poles))), tuple(sorted(float(abs(q)) for q in poles)),
             (float(max(abs(q) for q in poles)), INF)]):
        probs.append("ROC options are not the three possible ROCs")
    if parse_roc(correct["roc"]["html"], mags) != (lo, hi):
        probs.append(f"marked ROC {correct['roc']['html']} but expected {(lo, hi)} for {params['prop']}")
    # h[n]: z-transform of the answer must equal H(z) inside the expected ROC
    hp, hn = seq_fn(correct["hpos"]), seq_fn(correct["hneg"])

    def h(n):
        return float(hp(n)) if n >= 0 else float(hn(n))

    b = [float(Fraction(c)) for c in params["num"]]
    rad = test_radius(lo, hi)
    for z0 in (rad, rad * complex(math.cos(0.7), math.sin(0.7))):
        Hz = (b[0] + b[1] / z0) / ((1 - float(poles[0]) / z0) * (1 - float(poles[1]) / z0))
        got = zsum(h, z0, -400, 400)
        if abs(got - Hz) > 1e-8 * max(1, abs(Hz)):
            probs.append(f"sum h[n] z^-n = {got} != H({z0}) = {Hz}")
    return probs


# ---------------------------------------------------------------- student-style strings
def plain(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def pow_str(p):
    p = Fraction(p)
    s = plain(p)
    return f"{s}^n" if (p > 0 and p.denominator == 1) else f"({s})^n"


def terms_str(terms):
    """'2*(1/2)^n - (3/4)*(-2)^n'; '0' if empty."""
    out = ""
    for c, p in terms:
        c = Fraction(c)
        body = pow_str(p) if abs(c) == 1 else f"({plain(abs(c))})*{pow_str(p)}"
        out += ("-" if c < 0 else "") + body if not out else (" - " if c < 0 else " + ") + body
    return out or "0"


def terms_alt(terms):
    """Same sum, written differently: reversed order, p^n split into (-1)^n a^n / b^n,
    coefficient as integer numerator over denominator, and b^(-n) for unit numerators."""
    out = []
    for c, p in reversed(terms):
        c, p = Fraction(c), Fraction(p)
        parts = [str(abs(c.numerator))]
        if p < 0:
            parts.append("(-1)^n")
        if abs(p.numerator) != 1:
            parts.append(f"{abs(p.numerator)}^n")
        t = "*".join(parts)
        if p.denominator != 1:
            t += f"*{p.denominator}^(-n)" if abs(p.numerator) == 1 else f"/{p.denominator}^n"
        if c.denominator != 1:
            t += f"/{c.denominator}"
        out.append(("-" if c < 0 else "+") + "(" + t + ")")
    s = "".join(out)
    return (s[1:] if s.startswith("+") else s) or "0"


def frac_alt(x):
    """An equivalent way to type the number x: decimal if it terminates, else a non-reduced fraction."""
    q = Fraction(x).limit_denominator(1000)
    d = q.denominator
    while d % 2 == 0:
        d //= 2
    while d % 5 == 0:
        d //= 5
    if d == 1 and q.denominator != 1:
        return repr(float(q))
    return f"{2 * q.numerator}/{2 * q.denominator}"


def submissions(params, correct):
    poles = _poles(params)
    r, p, k = residues(params)
    A = []
    for q in poles:
        i = int(np.argmin(np.abs(p - float(q))))
        A.append(Fraction(float(r[i].real)).limit_denominator(1000))
    lo, hi = expected_roc(poles, params["prop"])
    right = [(a, q) for a, q in zip(A, poles) if abs(float(q)) <= lo + 1e-12]
    left = [(-a, q) for a, q in zip(A, poles) if abs(float(q)) >= hi - 1e-12]
    key_of = {parse_roc(o["html"], mags_of(poles)): o["key"] for o in params["roc"]}
    causal_key = next(v for (l, h), v in key_of.items() if h == INF)

    cases = [
        ({"A1": frac_alt(A[0]), "A2": frac_alt(A[1])}, {"A1": 1, "A2": 1}),      # 0.75 / 6/8 forms
        ({"hpos": terms_alt(right), "hneg": terms_alt(left)}, {"hpos": 1, "hneg": 1}),   # rearranged
        ({"hpos": terms_str(right), "hneg": terms_str(left)}, {"hpos": 1, "hneg": 1}),   # plain form
        ({"hpos": "0.5^n"}, {"hpos": "invalid"}),                                 # decimals refused
    ]
    dec = {nm: sig4(v) for nm, v in (("A1", A[0]), ("A2", A[1])) if Fraction(v).denominator != 1}
    if params["multiplied"]:
        dec.update({nm: sig4(q) for nm, q in (("p1", poles[0]), ("p2", poles[1])) if q.denominator != 1})
    if dec:                                                                       # 4 significant digits
        cases.append((dec, {nm: 1 for nm in dec}))
    if differs(A[1], A[0]) and differs(A[0], A[1]):                              # residues swapped
        cases.append(({"A1": plain(A[1]), "A2": plain(A[0])}, {"A1": 0, "A2": 0}))
    if right:   # (a/b)^n typed as b^(-n) where possible, terms in reverse order
        cases.append(({"hpos": " + ".join(f"({plain(a)})*{pow_str(q)}" for a, q in reversed(right))}, {"hpos": 1}))
    if left:    # forgot the minus sign of the left-sided pair
        cases.append(({"hneg": terms_str([(-a, q) for a, q in left])}, {"hneg": 0}))
    if (lo, hi)[1] != INF:   # answered the causal system instead
        causal = [(a, q) for a, q in zip(A, poles)]
        cases.append(({"roc": causal_key, "hpos": terms_str(causal), "hneg": "0"},
                      {"roc": 0, "hpos": 0, "hneg": 0}))
    else:                    # causal asked: the both-left-sided (anti-causal) answer
        anti_key = next(v for (l, h), v in key_of.items() if l == 0)
        cases.append(({"roc": anti_key, "hpos": "0", "hneg": terms_str([(-a, q) for a, q in zip(A, poles)])},
                      {"roc": 0, "hpos": 0, "hneg": 0}))
    if params["multiplied"]:
        cases.append(({"p1": plain(poles[1]), "p2": plain(poles[0])}, {"p1": 0, "p2": 0}))          # swapped
        if differs(1 / poles[0], poles[0]) and differs(1 / poles[1], poles[1]):
            cases.append(({"p1": plain(1 / poles[0]), "p2": plain(1 / poles[1])}, {"p1": 0, "p2": 0}))  # roots in z^-1
    return cases
