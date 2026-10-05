"""Independent check of questions/convolution/infinite-convolution.

Rebuilds x and h as truncated float arrays, convolves them with np.convolve (exact on the
compared range because both are right-sided) and compares with the closed form (parsed with
SymPy) at many n, plus zeros before n0 and the requested sample."""

from fractions import Fraction as F

import numpy as np
from checklib import S, N

NMIN, LEN = -6, 240          # arrays cover n = NMIN .. NMIN + LEN - 1


def _signals(p):
    n = np.arange(NMIN, NMIN + LEN)
    if p["template"] in ("ab", "aa"):
        a, b = float(F(p["a"])), float(F(p["b"]))
        x = p["A"] * a ** n.astype(float) * (n >= p["k1"])
        h = p["B"] * b ** n.astype(float) * (n >= p["k2"])
    else:
        x = np.zeros(LEN)
        for i, ci in enumerate(p["c"]):
            x[p["nx"] + i - NMIN] = ci
        b = float(F(p["b"]))
        h = p["B"] * b ** n.astype(float) * (n >= p["k"])
    return x, h


def _numeric_y(p):
    x, h = _signals(p)
    y = np.convolve(x, h)            # y[i] is at n = 2*NMIN + i, exact for n < NMIN + LEN + NMIN
    return lambda m: float(y[m - 2 * NMIN])


def _eval(expr_str, m):
    return float(S(expr_str).subs(N, m))


def check(params, correct):
    probs = []
    p = params
    yn = _numeric_y(p)
    n0 = correct["n0"]
    form = correct["yform"]["_value"]
    if "." in form:
        probs.append(f"closed form contains a float: {form}")
    # first nonzero index
    first = next(m for m in range(2 * NMIN, 40) if abs(yn(m)) > 1e-12)
    if first != n0:
        probs.append(f"n0 = {n0} but numeric y first nonzero at {first}")
    valid_from = n0 + (p["L"] - 1 if p["template"] == "finite" else 0)
    expr = S(form)
    for m in list(range(valid_from, valid_from + 12)) + [valid_from + 25]:
        want = yn(m)
        got = float(expr.subs(N, m))
        if abs(got - want) > 1e-9 * max(1.0, abs(want)):
            probs.append(f"closed form at n={m}: {got} vs numeric {want}")
            break
    m = n0 + p["c_off"]
    if abs(correct["yc"] - yn(m)) > 1e-9 * max(1.0, abs(yn(m))):
        probs.append(f"y[n0+{p['c_off']}] = {correct['yc']} vs numeric {yn(m)}")
    if p["template"] == "finite":
        m1 = n0 + p["c_off"]
        if m1 >= valid_from or abs(float(expr.subs(N, m1)) - yn(m1)) < 1e-9:
            probs.append("finite template: the sample in (c) is not a transient sample")
    if p["template"] == "ab" and abs(F(p["a"])) == abs(F(p["b"])):
        probs.append("|a| = |b| in the a != b template")
    if p["template"] == "aa" and n0 == 0:
        probs.append("a = b template without a shift (no-shift probe would coincide)")
    return probs


# ----------------------------------------------------------------------------- submissions
def _pb(q):
    """Base as a student types it: 2, (1/2), (-1/3), (-2)."""
    q = F(q)
    s = str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return s if (q > 0 and q.denominator == 1) else f"({s})"


def _pn(q):
    q = F(q)
    s = str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return f"({s})"


def _differs(s1, s2, ns):
    e1, e2 = S(s1.replace("^", "**")), S(s2.replace("^", "**"))
    return any(abs(complex(e1.subs(N, m)) - complex(e2.subs(N, m))) > 1e-9 * max(1.0, abs(complex(e2.subs(N, m))))
               for m in ns)


def submissions(params, correct):
    p = params
    n0 = correct["n0"]
    form = correct["yform"]["_value"]
    t = p["template"]
    eq, bad = [], []
    if t == "ab":
        a, b, A, B, k1, k2 = F(p["a"]), F(p["b"]), F(p["A"]), F(p["B"]), p["k1"], p["k2"]
        Ca = A * B * a ** (1 - k2) * b ** k2 / (a - b)
        Cb = A * B * a ** k1 * b ** (1 - k1) / (b - a)
        # key's shifted form  C a^{n0} a^{n-n0}
        eq.append(f"{_pn(Ca * a ** n0)}*{_pb(a)}^(n-({n0})) + {_pn(Cb * b ** n0)}*{_pb(b)}^(n-({n0}))")
        # reciprocal bases: (1/2)^n written as 2^(-n)
        eq.append(f"{_pn(Ca)}*{_pb(1 / a)}^(-n) + {_pn(Cb)}*{_pb(1 / b)}^(-n)")
        # common factor AB/(b-a)
        eq.append(f"{_pn(A * B / (b - a))}*({_pn(a ** k1 * b ** (1 - k1))}*{_pb(b)}^n - {_pn(a ** (1 - k2) * b ** k2)}*{_pb(a)}^n)")
        # mistakes: identity without the shifts; a and b swapped (shift attached to the wrong base); sign
        bad.append(f"{_pn(A * B)}*({_pb(a)}^(n+1) - {_pb(b)}^(n+1))/{_pn(a - b)}")
        Sb = A * B * b ** (1 - k2) * a ** k2 / (b - a)        # coefficient of b^n after swapping
        Sa = A * B * b ** k1 * a ** (1 - k1) / (a - b)        # coefficient of a^n after swapping
        bad.append(f"{_pn(Sb)}*{_pb(b)}^n + {_pn(Sa)}*{_pb(a)}^n")
        bad.append(f"-({form})")
    elif t == "aa":
        a, AB = F(p["a"]), F(p["A"] * p["B"])
        eq.append(f"{_pn(AB)}*n*{_pb(a)}^n + {_pn(AB * (1 - n0))}*{_pb(a)}^n")
        eq.append(f"{_pn(AB)}*(n-{_pn(n0 - 1)})*{_pb(1 / a)}^(-n)")
        bad.append(f"{_pn(AB)}*(n+1)*{_pb(a)}^n")              # identity without the shift
        bad.append(f"{_pn(AB)}*(n-{_pn(n0)})*{_pb(a)}^n")       # term count off by one
        bad.append(f"-({form})")
    else:
        b, B, c, nx = F(p["b"]), F(p["B"]), p["c"], p["nx"]
        C = B * sum(ci * b ** (-(nx + i)) for i, ci in enumerate(c))
        eq.append(" + ".join(f"{_pn(ci * B)}*{_pb(b)}^(n-({nx + i}))" for i, ci in enumerate(c)))
        eq.append(f"{_pn(C)}*{_pb(1 / b)}^(-n)")
        bad.append(f"{_pn(B * sum(c))}*{_pb(b)}^n")                                    # copies not shifted
        bad.append(" + ".join(f"{_pn(ci * B)}*{_pb(b)}^(n+({nx + i}))" for i, ci in enumerate(c)))  # wrong direction
        bad.append(f"-({form})")
    ns = list(range(n0 + 3, n0 + 9))
    cases = [({"yform": s}, {"yform": 1}) for s in eq]
    cases += [({"yform": s}, {"yform": 0}) for s in bad if _differs(s, form, ns)]
    cases += [({"n0": str(n0 + 1)}, {"n0": 0}), ({"n0": str(n0 - 1)}, {"n0": 0}),
              ({"yc": p["yc_frac"]}, {"yc": 1}),
              ({"n0": "1.5"}, {"n0": "invalid"})]
    # (c): a 4-significant-digit decimal is accepted (rtol 1e-3); the sample one index early is not
    yc = correct["yc"]
    if yc != int(yc):
        cases.append(({"yc": f"{yc:.4g}"}, {"yc": 1}))
    yn = _numeric_y(p)
    wrong = yn(n0 + p["c_off"] - 1)
    if abs(wrong - yc) > 1e-2 * max(abs(yc), 1e-4):
        cases.append(({"yc": repr(wrong)}, {"yc": 0}))
    return cases
