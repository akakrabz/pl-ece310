"""Independent check of questions/signals/geometric-sums.

(a) the limits and the base are parsed from the displayed LaTeX and the sum is computed term by
    term (exact Fractions for finite sums, 3000-term float sums for infinite ones, numpy for the
    roots of unity); (b) the defining sum is evaluated by brute force at several n and compared
    with the stored SymPy answer; (c) convergence is tested numerically at sample values of a in
    each of the regions |a| < 1/c, 1/c < |a| < c, |a| > c and matched with the marked option."""

import math
import re
from fractions import Fraction

import numpy as np
from checklib import N, S as Sym

INF = None


# ----------------------------------------------------------------------------- LaTeX parsing
def parse_base(t):
    """'\\left(\\tfrac{2}{3}\\right)' / '\\left(-\\tfrac{1}{2}\\right)' / '2' / '\\left(-2\\right)' -> Fraction."""
    t = t.strip()
    m = re.fullmatch(r"\\left\((-?)\\tfrac\{(\d+)\}\{(\d+)\}\\right\)", t)
    if m:
        return Fraction(int(m.group(2)), int(m.group(3))) * (-1 if m.group(1) else 1)
    m = re.fullmatch(r"\\left\((-?\d+)\\right\)", t)
    if m:
        return Fraction(int(m.group(1)))
    m = re.fullmatch(r"(\d+)", t)
    if m:
        return Fraction(int(m.group(1)))
    raise ValueError(f"cannot parse base {t!r}")


def parse_lim(t):
    t = t.strip()
    if t == r"\infty":
        return math.inf
    if t == r"-\infty":
        return -math.inf
    return int(t)


SUM = r"\\sum_\{n=([^}]+)\}\^\{([^}]+)\} ?"


def term_sum(f, lo, hi, cap=3000):
    """sum_{n=lo}^{hi} f(n), infinite limits truncated at +-cap."""
    lo = -cap if lo == -math.inf else lo
    hi = cap if hi == math.inf else hi
    return sum(f(n) for n in range(lo, hi + 1))


def value_a(tex):
    m = re.fullmatch(SUM + r"(.+)\^\{n\}", tex)
    if m:
        lo, hi, a = parse_lim(m.group(1)), parse_lim(m.group(2)), parse_base(m.group(3))
        if lo != -math.inf and hi != math.inf:
            return term_sum(lambda n: a ** n, lo, hi), ("fin", a, lo, hi)
        return term_sum(lambda n: float(a) ** n, lo, hi), ("inf" if hi == math.inf else "left", a, lo, hi)
    m = re.fullmatch(SUM + r"(.+)\^\{\\lvert n\\rvert\}", tex)
    if m:
        a = parse_base(m.group(3))
        return term_sum(lambda n: float(a) ** abs(n), -math.inf, math.inf), ("two_abs", a, None, None)
    m = re.fullmatch(SUM + r"\\left\((.+)\^\{n\}\\,u\[n\] \+ (.+)\^\{n\}\\,u\[(-n|-n - 1)\]\\right\)", tex)
    if m:
        a, b, arg = parse_base(m.group(3)), parse_base(m.group(4)), m.group(5)
        top = 0 if arg == "-n" else -1
        val = term_sum(lambda n: float(a) ** n, 0, math.inf) + term_sum(lambda n: float(b) ** n, -math.inf, top)
        return val, ("two_pair", a, b, top)
    m = re.fullmatch(r"\\sum_\{n=0\}\^\{(\d+)\} e\^\{(-?)j(\d*)\\pi n(?:/(\d+))?\}", tex)
    if m:
        Nn = int(m.group(1)) + 1
        r = Fraction(int(m.group(3) or 1), int(m.group(4) or 1)) * (-1 if m.group(2) else 1)
        n = np.arange(Nn)
        val = complex(np.sum(np.exp(1j * math.pi * float(r) * n)))
        if abs(val.imag) > 1e-9:
            raise ValueError("complex roots-of-unity sum")
        return val.real, ("roots", r, Nn, None)
    raise ValueError(f"cannot parse (a) {tex!r}")


def value_b(tex, valid, n):
    """Brute-force s[n] from the displayed sum (truncated at 3000 terms)."""
    m = re.fullmatch(r"\\sum_\{k=([^}]+)\}\^\{([^}]+)\} (.+)\^\{k\}\\,(.+)\^\{n-k\}", tex)
    if m:
        a, b = parse_base(m.group(3)), parse_base(m.group(4))
        lo, hi = parse_lim(m.group(1)), m.group(2)
        assert lo == 0 and hi == "n"
        return sum(a ** k * b ** (n - k) for k in range(0, n + 1))
    m = re.fullmatch(r"\\sum_\{k=([^}]+)\}\^\{([^}]+)\} (.+)\^\{k\}", tex)
    if not m:
        raise ValueError(f"cannot parse (b) {tex!r}")
    a = parse_base(m.group(3))
    lo = n if m.group(1) == "n" else parse_lim(m.group(1))
    hi = n if m.group(2) == "n" else parse_lim(m.group(2))
    if lo != -math.inf and hi != math.inf:
        return sum(a ** k for k in range(lo, hi + 1))
    return term_sum(lambda k: float(a) ** k, lo, hi)


def b_points(valid):
    if valid.startswith(r"n \ge"):
        m0 = int(valid.split(r"\ge")[1])
        return list(range(m0, m0 + 7))
    return list(range(-3, 6))


# ----------------------------------------------------------------------------- (c)
def converges(spec, a):
    """Numerical convergence test of the series in conv_spec at the value a."""
    c = spec["c"]
    base = a / c if spec["base"] == "a/c" else c * a
    ex = spec["exp"]
    e = {"n": lambda n: n, "-n": lambda n: -n, "abs": lambda n: abs(n), "-abs": lambda n: -abs(n)}[ex]
    lo, hi = spec["lo"], spec["hi"]
    ends = []
    if hi is None:
        ends.append(4000)
    if lo is None:
        ends.append(-4000)
    lb = math.log(abs(base))
    return all(e(n) * lb < -50 for n in ends)      # |term| at the far ends below e^-50


def option_regions(html, c):
    """Option text -> set of regions {1: |a| < 1/c, 2: 1/c < |a| < c, 3: |a| > c} it describes."""
    if "every" in html:
        return {1, 2, 3}
    if "no value" in html:
        return set()
    m = re.search(r"\\lvert a\\rvert \\(lt|gt) (\\tfrac\{1\}\{(\d+)\}|(\d+))", html)
    small = m.group(3) is not None
    k = int(m.group(3) or m.group(4))
    assert k == c, (html, c)
    if m.group(1) == "lt":
        return {1} if small else {1, 2}
    return {2, 3} if small else {3}


def regions_converging(spec):
    c = spec["c"]
    out = set()
    for reg, mag in ((1, 1 / (2 * c)), (2, 1.0), (3, 2.0 * c)):
        res = {converges(spec, s * mag) for s in (1, -1)}
        assert len(res) == 1
        if res.pop():
            out.add(reg)
    return out


def check(params, correct):
    probs = []
    val, info = value_a(params["sumA_tex"])
    if not math.isclose(float(val), correct["S"], rel_tol=1e-9, abs_tol=1e-12):
        probs.append(f"(a) S = {correct['S']} but the terms sum to {float(val)} ({info})")
    q = Fraction(correct["S"]).limit_denominator(1000)
    if abs(float(q) - correct["S"]) > 1e-12 or q.denominator > 64:
        probs.append(f"(a) S = {correct['S']} is not a simple fraction")
    if info[0] != "roots" and correct["S"] == 0:
        probs.append("(a) S = 0")
    if len(_a_mistakes(params, correct)) < 2:
        probs.append("(a) fewer than two distinct classic mistakes")
    expr = Sym(correct["closed"]["_value"])
    for n in b_points(params["sumB_valid_tex"]):
        want = float(value_b(params["sumB_tex"], params["sumB_valid_tex"], n))
        got = complex(expr.subs(N, n).evalf(30))
        if abs(got - want) > 1e-9 * max(1.0, abs(want)):
            probs.append(f"(b) s[{n}] = {got} from the closed form but {want} from the sum")
            break
    if "." in correct["closed"]["_value"]:
        probs.append("(b) float in the symbolic answer")
    regs = regions_converging(params["conv_spec"])
    marked = option_regions(correct["conv"]["html"], params["conv_spec"]["c"])
    if regs != marked:
        probs.append(f"(c) converges in regions {regs} but the marked option {correct['conv']['html']} says {marked}")
    texts = [option_regions(o["html"], params["conv_spec"]["c"]) for o in params["conv"]]
    if sum(t == regs for t in texts) != 1 or len(texts) != 6:
        probs.append("(c) option list does not contain exactly one correct description")
    return probs


# ----------------------------------------------------------------------------- probes
def _num(x):
    return f"{float(x):.6g}"


def _a_mistakes(params, correct):
    _, info = value_a(params["sumA_tex"])
    S = correct["S"]
    fam = info[0]
    c = []
    if fam == "fin":
        _, a, lo, hi = info
        Nt = hi - lo + 1
        c += [(1 - a ** Nt) / (1 - a), (a ** lo - a ** hi) / (1 - a), (a ** lo - a ** (hi + 1)) / (1 + a)]
    elif fam == "inf":
        _, a, lo, _ = info
        c += [1 / (1 - a), a ** (lo + 1) / (1 - a), a ** (-lo) / (1 - a)]   # start at 0 / one late / a^{-n1}
    elif fam == "left":
        _, a, _, hi = info
        c += [a ** hi / (1 - a), a ** hi / (a - 1), 1 / (1 - 1 / a)]
    elif fam == "two_abs":
        a = info[1]
        c += [2 / (1 - a), 1 / (1 - a)]
    elif fam == "two_pair":
        _, a, b, top = info
        right = 1 / (1 - a)
        c += [right + (1 / (b - 1) if top == 0 else b / (b - 1)), right + 1 / (1 - b), right]
    else:
        Nn = info[2]
        c += [Nn, 1] if abs(S) < 1e-12 else [0, 1]
    out = []
    for v in c:
        v = float(v)
        if abs(v - S) > 2e-2 * max(1e-2, abs(S)) and all(abs(v - o) > 1e-9 for o in out):
            out.append(v)
    return out


def _p(x):
    x = Fraction(x)
    return f"({x.numerator})" if x.denominator == 1 else f"({x.numerator}/{x.denominator})"


def _pq(x, e):
    """x^e written as p^e/q^e."""
    x = Fraction(x)
    return f"(({x.numerator})^({e})/{x.denominator}^({e}))"


def _differs(s1, s2, ns):
    e1, e2 = Sym(s1.replace("^", "**")), Sym(s2.replace("^", "**"))
    return any(abs(complex(e1.subs(N, m).evalf(30)) - complex(e2.subs(N, m).evalf(30)))
               > 1e-6 * max(1.0, abs(complex(e2.subs(N, m).evalf(30)))) for m in ns)


def submissions(params, correct):
    cases = []
    S = correct["S"]
    q = Fraction(S).limit_denominator(1000)
    cases.append(({"S": f"{q.numerator}/{q.denominator}"}, {"S": 1}))
    cases.append(({"S": f"{3 * q.numerator}/{3 * q.denominator}"}, {"S": 1}))
    if q.denominator != 1:
        cases.append(({"S": f"{S:.4g}"}, {"S": 1}))
    for v in _a_mistakes(params, correct):
        cases.append(({"S": _num(v)}, {"S": 0}))

    spec = params["sumB_spec"]
    fam, a = spec["fam"], Fraction(spec["a"])
    A = _p(a)
    if fam == "G1":
        eq = [f"(1 - {A}^(n+1))/(1 - {A})", f"({A}^(n+1) - 1)/({A} - 1)", f"(1 - {_pq(a, 'n+1')})/(1 - {A})"]
        bad = [f"(1 - {A}^n)/(1 - {A})", f"1/(1 - {A})"]
    elif fam == "G2":
        eq = [f"{A}^n/(1 - {A})", f"{_pq(a, 'n')}/(1 - {A})", f"{A}^(n+1)/({A} - {A}^2)"]
        bad = [f"1/(1 - {A})", f"{A}^(n+1)/(1 - {A})"]
    elif fam == "G3":
        b = Fraction(spec["b"])
        B = _p(b)
        eq = [f"({A}^(n+1) - {B}^(n+1))/({A} - {B})", f"{B}^n*(1 - ({A}/{B})^(n+1))/(1 - {A}/{B})",
              f"({_pq(a, 'n+1')} - {_pq(b, 'n+1')})/({A} - {B})"]
        bad = [f"({A}^(n+1) - {B}^(n+1))/({B} - {A})", f"({A}^n - {B}^n)/({A} - {B})"]
    elif fam == "G4":
        eq = [f"{A}^(n+1)/({A} - 1)", f"{A}^n/(1 - 1/{A})", f"{_pq(a, 'n+1')}/({A} - 1)"]
        bad = [f"{A}^n/(1 - {A})", f"{A}^n/({A} - 1)"]
    else:
        m = spec["m"]
        eq = [f"({A}^({m}) - {A}^(n+1))/(1 - {A})", f"{A}^({m})*(1 - {A}^(n-({m})+1))/(1 - {A})",
              f"({A}^(n+1) - {A}^({m}))/({A} - 1)"]
        bad = [f"(1 - {A}^(n+1))/(1 - {A})", f"({A}^({m}) - {A}^n)/(1 - {A})"]
    ns = b_points(params["sumB_valid_tex"])
    form = correct["closed"]["_value"]
    cases += [({"closed": s}, {"closed": 1}) for s in eq]
    cases += [({"closed": s}, {"closed": 0}) for s in bad if _differs(s, form, ns)]
    cases.append(({"closed": "0.5^n"}, {"closed": "invalid"}))

    right = correct["conv"]["key"]
    for o in params["conv"]:
        cases.append(({"conv": o["key"]}, {"conv": 1 if o["key"] == right else 0}))
    return cases
