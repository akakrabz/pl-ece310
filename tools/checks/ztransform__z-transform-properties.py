"""Independent check of questions/ztransform/z-transform-properties.

The given pair is x[n] = A a^n u[n] (ROC |z| > |a|) or x[n] = -A a^n u[-n-1] (ROC |z| < |a|). This file
builds y[n] sample by sample from x[n] and the operation (its own code; the self-convolution closed form
is first confirmed by direct convolution sums), compares Y(z) with truncated sums sum_n y[n] z^-n
(checklib.zsum) inside the ROC, derives the ROC of y from its samples (side, first/last index) and the
boundary radius, confirms that boundary numerically (the truncated sums settle just on the ROC side of
it and blow up just outside), and parses every displayed option back into a region."""

import cmath
import html
import re
from fractions import Fraction as F

import sympy
from checklib import S, Z, zsum


def x_sample(p, n):
    A, a = F(p["A"]), F(p["a"])
    if p["side"] == "R":
        return A * a ** n if n >= 0 else F(0)
    return -A * a ** n if n <= -1 else F(0)


def y_sample(p, n):
    op, A, a = p["op"], F(p["A"]), F(p["a"])
    if op == "shift":
        return x_sample(p, n - p["k"])
    if op == "mult_n":
        return n * x_sample(p, n)
    if op == "scale":
        return F(p["c"]) ** n * x_sample(p, n)
    if op == "reverse":
        return x_sample(p, -n)
    if op == "conv":    # closed form, confirmed against direct sums in check()
        if p["side"] == "R":
            return A * A * (n + 1) * a ** n if n >= 0 else F(0)
        return -A * A * (n + 1) * a ** n if n <= -2 else F(0)
    if op == "diff":
        return x_sample(p, n) - x_sample(p, n - 1)
    raise ValueError(op)


def panel_y_sample(terms, n):
    """The answer panel's time-domain y[n] (terms c * P(n) * b^(n-s) * step), evaluated at n."""
    tot = F(0)
    for t in terms:
        if (t["side"] == "R" and n < t["e"]) or (t["side"] == "L" and n > t["e"]):
            continue
        tot += F(t["c"]) * {"1": 1, "n": n, "n+1": n + 1}[t["poly"]] * F(t["b"]) ** (n - t["s"])
    return tot


def boundary(p):
    a = abs(F(p["a"]))
    if p["op"] == "scale":
        return a * abs(F(p["c"]))
    if p["op"] == "reverse":
        return 1 / a
    return a


def region(p):
    nz = [n for n in range(-70, 71) if y_sample(p, n) != 0]
    first, last = min(nz), max(nz)
    to_right = any(y_sample(p, n) != 0 for n in range(55, 71))
    to_left = any(y_sample(p, n) != 0 for n in range(-70, -54))
    assert to_right != to_left
    r = boundary(p)
    if to_right:
        return (r, None, False, first >= 0)
    return (F(0), r, last <= 0, False)


def _num(s):
    m = re.fullmatch(r"\\[td]?frac\{(\d+)\}\{(\d+)\}", s)
    if m:
        return F(int(m[1]), int(m[2]))
    if re.fullmatch(r"\d+", s):
        return F(int(s))
    if s == r"\infty":
        return None
    raise ValueError(f"cannot parse radius {s!r}")


def parse_roc(text):
    s = html.unescape(text).replace("$", "")
    s = re.sub(r"\s+", "", s).replace(r"\lvertz\rvert", "Z")
    m = re.fullmatch(r"Z\\gt(.+)", s)
    if m:
        return (_num(m[1]), None, False, True)
    m = re.fullmatch(r"Z\\lt(.+)", s)
    if m:
        return (F(0), _num(m[1]), True, False)
    m = re.fullmatch(r"(.+?)\\ltZ\\lt(.+)", s)
    if m:
        return (_num(m[1]), _num(m[2]), False, False)
    raise ValueError(f"cannot parse ROC text {text!r}")


def converges(p, rho):
    zz = rho * cmath.exp(0.41j)
    try:
        s1 = zsum(lambda n: y_sample(p, n), zz, -300, 300)     # probes sit at 0.8 r / 1.25 r: tail ~ 0.8^300
        s2 = zsum(lambda n: y_sample(p, n), zz, -600, 600)
    except OverflowError:
        return False
    return abs(s2 - s1) <= 1e-7 * max(1.0, abs(s2))


def test_point(reg):
    lo, hi = reg[0], reg[1]
    rho = float(lo) * 1.6 + 0.5 if hi is None else float(hi) * 0.6
    return rho


def evalz(expr, zz):
    return complex(S(expr.replace("^", "**")).subs(Z, zz).evalf(30))


def check(params, correct):
    probs = []
    p = params
    if p["op"] == "conv":   # closed form of x*x against direct convolution sums
        for n in range(-12, 13):
            direct = sum((x_sample(p, k) * x_sample(p, n - k) for k in range(-40, 41)), F(0))
            if direct != y_sample(p, n):
                probs.append(f"self-convolution closed form wrong at n = {n}")
                break
    for n in range(-15, 16):
        if panel_y_sample(p["y_terms"], n) != y_sample(p, n):
            probs.append(f"answer panel's time-domain y[n] is wrong at n = {n}")
            break
    reg = region(p)
    r = float(boundary(p))
    # numerical confirmation of the boundary circle
    for rho, want in ((0.8 * r, reg[1] is not None), (1.25 * r, reg[1] is None)):
        if converges(p, rho) != want:
            probs.append(f"|z| = {rho:.4g}: convergence disagrees with the ROC {reg}")
    for ang in (0.3, 2.4):
        zz = test_point(reg) * cmath.exp(1j * ang)
        want = zsum(lambda n: y_sample(p, n), zz, -400, 400)  # test point ratio <= 0.63
        got = evalz(correct["Y"], zz)
        if abs(got - want) > 1e-7 * max(1.0, abs(want)):
            probs.append(f"Y({zz:.3f}) = {got} but the truncated sum gives {want}")
    opts = p["roc_choices"]
    if len(opts) < 4:
        probs.append(f"only {len(opts)} ROC options")
    try:
        regions = [parse_roc(o["text"]) for o in opts]
    except ValueError as exc:
        return probs + [str(exc)]
    if len(set(regions)) != len(regions):
        probs.append("two options describe the same region")
    good = [g for g, o in zip(regions, opts) if o["correct"]]
    if len(good) != 1 or good[0] != reg:
        probs.append(f"marked {good} but the ROC of y is {reg}")
    if parse_roc(p["roc_tex"]) != reg:
        probs.append("answer panel ROC differs")
    return probs


def _key_for(params, text):
    norm = lambda s: re.sub(r"\s+", "", html.unescape(s))
    for o in params["roc"]:
        if norm(o["html"]) == norm(text):
            return o["key"]
    raise KeyError(text)


def submissions(params, correct):
    p = params
    op, A, a, k, c = p["op"], F(p["A"]), F(p["a"]), p["k"], F(p["c"])
    reg = region(p)
    zz = test_point(reg) * cmath.exp(0.77j)
    ref = evalz(correct["Y"], zz)
    X = f"({A})/(1-({a})*z^(-1))"
    if op == "shift":
        pos = f"({A})*z^({1 - k})/(z-({a}))"
        bad = [f"({A})*z^({k})/(1-({a})*z^(-1))", X]                                   # z^{+k}; shift forgotten
    elif op == "mult_n":
        pos = f"({A * a})*z/(z-({a}))^2"
        bad = [f"-({A * a})*z^(-1)/(1-({a})*z^(-1))^2",                               # +z dX/dz (sign error)
               f"({A * a})*z^(-2)/(1-({a})*z^(-1))^2"]                                # -dX/dz (factor z missing)
    elif op == "scale":
        pos = f"({A})*z/(z-({a * c}))"
        bad = [f"({A})/(1-({a / c})*z^(-1))", X]                                       # X(cz); unscaled
    elif op == "reverse":
        pos = f"({-A / a})/(z-({1 / a}))"
        bad = [f"({A})/(1-({1 / a})*z^(-1))", X]                                       # pole inverted without the rest; unchanged
    elif op == "conv":
        pos = f"({A * A})*z^2/(z-({a}))^2"
        bad = [f"2*{X}", f"({A})/(1-({a})*z^(-2))"]                                    # X + X; X(z^2)
    else:
        pos = f"({A})*(z-1)/(z-({a}))"
        bad = [f"({A})*(1-z)/(1-({a})*z^(-1))", f"({A})*(1+z^(-1))/(1-({a})*z^(-1))"]  # (1 - z) X; sign of the delay
    cases = [({"Y": pos}, {"Y": 1}),                                                    # positive powers of z
             ({"Y": str(sympy.factor(sympy.cancel(S(correct["Y"]))))}, {"Y": 1})]       # one factored fraction
    for b in bad:
        if abs(evalz(b, zz) - ref) > 1e-6 * max(1.0, abs(ref)):
            cases.append(({"Y": b}, {"Y": 0}))
    for o in p["roc_choices"]:
        cases.append(({"roc": _key_for(p, o["text"])}, {"roc": 1 if o["correct"] else 0}))
    return cases
