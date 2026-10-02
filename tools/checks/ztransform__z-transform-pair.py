"""Independent check of questions/ztransform/z-transform-pair.

x[n] is rebuilt sample by sample from the term specs with this file's own code (not ece310.zt).
X(z) is compared with truncated sums sum_n x[n] z^-n (checklib.zsum) at points inside the ROC.
The ROC is derived here from the samples: the radii from the pieces that extend to +infinity /
-infinity and their exponential bases, z = 0 / z = infinity membership from the SUMMED signal's
last / first nonzero index. Every multiple-choice option's displayed LaTeX is parsed back into a
region and compared with it."""

import cmath
import html
import math
import re
from fractions import Fraction as F

import sympy
from checklib import S, Z, zsum

COSV = {"pi/2": [F(1), F(0), F(-1), F(0)],
        "pi/3": [F(1), F(1, 2), F(-1, 2), F(-1), F(-1, 2), F(1, 2)],
        "2*pi/3": [F(1), F(-1, 2), F(-1, 2)],
        "pi": [F(1), F(-1)]}


# ----------------------------------------------------------------------------- signal
def sample(spec, n):
    t = spec["type"]
    if t == "right_exp":
        return F(spec["A"]) * F(spec["a"]) ** n if n >= spec["k"] else F(0)
    if t == "left_exp":
        return F(spec["B"]) * F(spec["b"]) ** n if n <= spec["m"] else F(0)
    if t == "n_exp":
        return F(spec["A"]) * n * F(spec["a"]) ** n if n >= 0 else F(0)
    if t == "finite":
        i = n - spec["start"]
        return F(spec["values"][i]) if 0 <= i < len(spec["values"]) else F(0)
    if t == "cos":
        cv = COSV[spec["w"]]
        return F(spec["A"]) * F(spec["r"]) ** n * cv[n % len(cv)] if n >= 0 else F(0)
    raise ValueError(t)


def signal(specs):
    return lambda n: sum((sample(s, n) for s in specs), F(0))


def base_radius(spec):
    t = spec["type"]
    return abs(F(spec["a"])) if t in ("right_exp", "n_exp") else abs(F(spec["b"])) if t == "left_exp" else abs(F(spec["r"]))


# ----------------------------------------------------------------------------- regions (inner, outer|None, has0, hasinf) or "EMPTY"
def term_region(spec):
    nz = [n for n in range(-70, 71) if sample(spec, n) != 0]
    first, last = min(nz), max(nz)
    to_right = any(sample(spec, n) != 0 for n in range(55, 71))
    to_left = any(sample(spec, n) != 0 for n in range(-70, -54))
    assert not (to_right and to_left)
    if to_right:
        return (base_radius(spec), None, False, first >= 0)
    if to_left:
        return (F(0), base_radius(spec), last <= 0, False)
    return (F(0), None, last <= 0, first >= 0)


def closed_form(specs):
    """ROC of the sum: (inner, outer|None, has0, hasinf) or 'EMPTY'."""
    x = signal(specs)
    regs = [term_region(sp) for sp in specs]
    rights = [r[0] for r in regs if r[1] is None and r[0] > 0]
    lefts = [r[1] for r in regs if r[1] is not None]
    if any(x(n) != 0 for n in range(55, 71)) != bool(rights) or any(x(n) != 0 for n in range(-70, -54)) != bool(lefts):
        raise AssertionError("infinite tails of the pieces cancel: closed form does not apply")
    inner = max(rights, default=F(0))
    outer = min(lefts) if lefts else None
    if outer is not None and inner >= outer:
        return "EMPTY"
    return (inner, outer, not rights and all(x(n) == 0 for n in range(1, 71)),
            not lefts and all(x(n) == 0 for n in range(-70, 0)))


def intersect(regions):
    inner = max(r[0] for r in regions)
    outs = [r[1] for r in regions if r[1] is not None]
    outer = min(outs) if outs else None
    if outer is not None and inner >= outer:
        return "EMPTY"
    return (inner, outer, all(r[2] for r in regions), all(r[3] for r in regions))


def test_radius(reg):
    lo, hi = float(reg[0]), reg[1]
    if hi is None:
        return 1.5 * lo + 0.6 if lo > 0 else 1.3
    hi = float(hi)
    return 0.6 * hi if lo == 0 else math.sqrt(lo * hi)


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
    if "notexist" in s or "empty" in s:
        return "EMPTY"
    if s in (r"\text{all}z", r"\text{all}z\text{(including}0\text{and}\infty)"):
        return (F(0), None, True, True)
    m = re.fullmatch(r"Z\\gt(.+)", s)
    if m:
        return (_num(m[1]), None, False, True)
    m = re.fullmatch(r"Z\\lt(.+)", s)
    if m:
        R = _num(m[1])
        return (F(0), R, True, False)
    m = re.fullmatch(r"(.+?)\\ltZ\\lt(.+)", s)
    if m:
        return (_num(m[1]), _num(m[2]), False, False)
    raise ValueError(f"cannot parse ROC text {text!r}")


# ----------------------------------------------------------------------------- transforms written here
def zinv_form(spec):
    """This term's X(z) in powers of z^-1 (own derivation of the table pairs)."""
    t = spec["type"]
    if t == "right_exp":
        A, a, k = F(spec["A"]), F(spec["a"]), spec["k"]
        return f"({A * a ** k})*z^({-k})/(1-({a})*z^(-1))"
    if t == "left_exp":
        B, b, m = F(spec["B"]), F(spec["b"]), spec["m"]
        return f"({-B * b ** (m + 1)})*z^({-(m + 1)})/(1-({b})*z^(-1))"
    if t == "n_exp":
        A, a = F(spec["A"]), F(spec["a"])
        return f"({A * a})*z^(-1)/(1-({a})*z^(-1))^2"
    if t == "finite":
        return "(" + "+".join(f"({v})*z^({-(spec['start'] + i)})" for i, v in enumerate(spec["values"])) + ")"
    if t == "cos":
        A, r, c = F(spec["A"]), F(spec["r"]), COSV[spec["w"]][1]
        return f"(({A})-({A * r * c})*z^(-1))/(1-({2 * r * c})*z^(-1)+({r * r})*z^(-2))"
    raise ValueError(t)


def pos_form(spec):
    """The same term written in positive powers of z."""
    t = spec["type"]
    if t == "right_exp":
        A, a, k = F(spec["A"]), F(spec["a"]), spec["k"]
        return f"({A * a ** k})*z^({1 - k})/(z-({a}))"
    if t == "left_exp":
        B, b, m = F(spec["B"]), F(spec["b"]), spec["m"]
        return f"({-B * b ** (m + 1)})*z^({-m})/(z-({b}))"
    if t == "n_exp":
        A, a = F(spec["A"]), F(spec["a"])
        return f"({A * a})*z/(z-({a}))^2"
    if t == "finite":
        ns = [spec["start"] + i for i in range(len(spec["values"]))]
        N = max(0, max(ns))
        return "((" + "+".join(f"({v})*z^({N - n})" for v, n in zip(spec["values"], ns)) + f")/z^({N}))"
    if t == "cos":
        A, r, c = F(spec["A"]), F(spec["r"]), COSV[spec["w"]][1]
        return f"((({A})*z^2-({A * r * c})*z)/(z^2-({2 * r * c})*z+({r * r})))"
    raise ValueError(t)


def mistakes_for(spec):
    """Classic wrong transforms of one term: (label, expression)."""
    t = spec["type"]
    out = [("sign", f"(-({zinv_form(spec)}))")]
    if t == "right_exp":
        A, a, k = F(spec["A"]), F(spec["a"]), spec["k"]
        c = A * a ** k
        if k != 0:
            out.append(("no z^-k", f"({c})/(1-({a})*z^(-1))"))
            out.append(("no a^k", f"({A})*z^({-k})/(1-({a})*z^(-1))"))
            out.append(("z^+k", f"({c})*z^({k})/(1-({a})*z^(-1))"))
        if abs(a) != 1:
            out.append(("1/a pole", f"({c})*z^({-k})/(1-({1 / a})*z^(-1))"))
    elif t == "left_exp":
        B, b, m = F(spec["B"]), F(spec["b"]), spec["m"]
        c = -B * b ** (m + 1)
        if m + 1 != 0:
            out.append(("no z^-k", f"({c})/(1-({b})*z^(-1))"))
            out.append(("no b^k", f"({-B})*z^({-(m + 1)})/(1-({b})*z^(-1))"))
        out.append(("1/b pole", f"({c})*z^({-(m + 1)})/(1-({1 / b})*z^(-1))"))
    elif t == "n_exp":
        A, a = F(spec["A"]), F(spec["a"])
        out.append(("not squared", f"({A * a})*z^(-1)/(1-({a})*z^(-1))"))
        out.append(("no a", f"({A})*z^(-1)/(1-({a})*z^(-1))^2"))
    elif t == "finite":
        out.append(("z^+n", "(" + "+".join(f"({v})*z^({spec['start'] + i})" for i, v in enumerate(spec["values"])) + ")"))
    elif t == "cos":
        A, r, c = F(spec["A"]), F(spec["r"]), COSV[spec["w"]][1]
        out.append(("num sign", f"(({A})+({A * r * c})*z^(-1))/(1-({2 * r * c})*z^(-1)+({r * r})*z^(-2))"))
    return out


def evalz(expr, zz):
    return complex(S(expr.replace("^", "**")).subs(Z, zz).evalf(30))


# ----------------------------------------------------------------------------- check
def check(params, correct):
    probs = []
    specs = params["specs"]
    reg = closed_form(specs)
    if reg == "EMPTY":
        return ["the ROC is empty (this question promises a z-transform)"]
    if intersect([term_region(s) for s in specs]) != reg:
        probs.append(f"pieces cancel at an end of the signal: the intersection of the pieces' ROCs is not the ROC {reg}")
    inf_radii = [base_radius(s) for s in specs if s["type"] != "finite"]
    if len(set(inf_radii)) != len(inf_radii):
        probs.append(f"repeated pole magnitude {inf_radii}: the ROC would be ambiguous")
    # X(z) against truncated sums inside the ROC
    x = signal(specs)
    rho = test_radius(reg)
    for ang in (0.37, 2.2):
        zz = rho * cmath.exp(1j * ang)
        want = zsum(x, zz)
        got = evalz(correct["X"], zz)
        mine = sum(evalz(zinv_form(s), zz) for s in specs)
        if abs(got - want) > 1e-7 * max(1.0, abs(want)):
            probs.append(f"X({zz:.3f}) = {got} but the truncated sum gives {want}")
        if abs(mine - want) > 1e-7 * max(1.0, abs(want)):
            probs.append(f"checker's own transform disagrees with the truncated sum at {zz:.3f}")
    # the ROC options
    opts = params["roc_choices"]
    if len(opts) < 4:
        probs.append(f"only {len(opts)} ROC options")
    regions = []
    for o in opts:
        try:
            regions.append(parse_roc(o["text"]))
        except ValueError as exc:
            probs.append(str(exc))
            return probs
    if len(set(regions)) != len(regions):
        probs.append("two options describe the same region")
    good = [r for r, o in zip(regions, opts) if o["correct"]]
    if len(good) != 1:
        probs.append(f"{len(good)} options marked correct")
    elif good[0] != reg:
        probs.append(f"marked ROC {good[0]} but the signal's ROC is {reg}")
    if reg not in regions:
        probs.append("the true ROC is not among the options")
    if parse_roc(params["roc_tex"]) != reg:
        probs.append("answer panel ROC differs from the true ROC")
    return probs


def _key_for(params, text):
    norm = lambda s: re.sub(r"\s+", "", html.unescape(s))
    for o in params["roc"]:
        if norm(o["html"]) == norm(text):
            return o["key"]
    raise KeyError(text)


def submissions(params, correct):
    specs = params["specs"]
    reg = closed_form(specs)
    zz = test_radius(reg) * cmath.exp(0.9j)
    ref = evalz(correct["X"], zz)

    def differs(expr):
        return abs(evalz(expr, zz) - ref) > 1e-6 * max(1.0, abs(ref))

    cases = []
    # equivalent forms (score 1)
    cases.append(({"X": " + ".join(pos_form(s) for s in specs)}, {"X": 1}))                     # positive powers of z
    one = sympy.factor(sympy.cancel(S(correct["X"])))
    cases.append(({"X": str(one)}, {"X": 1}))                                                  # one combined fraction
    cases.append(({"X": " + ".join(zinv_form(s) for s in reversed(specs))}, {"X": 1}))          # terms reordered, ^ syntax
    # classic mistakes (score 0): one term replaced by a wrong transform
    order = sorted(range(len(specs)), key=lambda i: specs[i]["type"] != "left_exp")
    wrong = []
    for i in order:
        for label, bad in mistakes_for(specs[i]):
            expr = " + ".join(bad if j == i else zinv_form(s) for j, s in enumerate(specs))
            if differs(expr):
                wrong.append((label, expr))
    seen = set()
    for label, expr in wrong:
        if label in seen:
            continue
        seen.add(label)
        cases.append(({"X": expr}, {"X": 0}))
    if len(seen) < 2:
        cases.append(({"X": "2*(" + correct["X"] + ")"}, {"X": 0}))
    cases.append(({"X": "1/(1-n)"}, {"X": "invalid"}))
    # ROC: the correct option scores 1, each distractor 0
    for o in params["roc_choices"]:
        cases.append(({"roc": _key_for(params, o["text"])}, {"roc": 1 if o["correct"] else 0}))
    return cases
