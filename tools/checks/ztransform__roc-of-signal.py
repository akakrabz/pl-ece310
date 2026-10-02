"""Independent check of questions/ztransform/roc-of-signal.

Two independent decisions of the ROC, neither using ece310.zt:
  1. closed form from the samples: the radii come from the pieces that extend to +infinity (inner
     radius) / -infinity (outer radius), and z = 0 / z = infinity membership from the SUMMED signal's
     last / first nonzero index (so samples of different pieces that cancel at an end are handled);
  2. numerically: truncated sums sum_n x[n] z^-n over |n| <= 400 and |n| <= 800 (checklib.zsum) at
     radii just inside / just outside every boundary circle, between circles, near 0 and at large
     |z|: they must settle exactly where the closed form says the series converges.
Every option's displayed text is parsed back into a region and compared."""

import cmath
import html
import math
import re
from fractions import Fraction as F

from checklib import zsum


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
    raise ValueError(t)


def base_radius(spec):
    t = spec["type"]
    return abs(F(spec["a"])) if t in ("right_exp", "n_exp") else abs(F(spec["b"])) if t == "left_exp" else None


def term_region(spec):
    nz = [n for n in range(-70, 71) if sample(spec, n) != 0]
    first, last = min(nz), max(nz)
    to_right = any(sample(spec, n) != 0 for n in range(55, 71))
    to_left = any(sample(spec, n) != 0 for n in range(-70, -54))
    if to_right:
        return (base_radius(spec), None, False, first >= 0)
    if to_left:
        return (F(0), base_radius(spec), last <= 0, False)
    return (F(0), None, last <= 0, first >= 0)


def closed_form(specs):
    """ROC of the sum: (inner, outer|None, has0, hasinf) or 'EMPTY'."""
    x = lambda n: sum((sample(sp, n) for sp in specs), F(0))
    rights = [base_radius(sp) for sp in specs if term_region(sp)[1] is None and term_region(sp)[0] > 0]
    lefts = [base_radius(sp) for sp in specs if term_region(sp)[1] is not None]
    sum_right = any(x(n) != 0 for n in range(55, 71))
    sum_left = any(x(n) != 0 for n in range(-70, -54))
    if sum_right != bool(rights) or sum_left != bool(lefts):
        raise AssertionError("infinite tails of the pieces cancel: closed form does not apply")
    inner = max(rights, default=F(0))
    outer = min(lefts) if lefts else None
    if outer is not None and inner >= outer:
        return "EMPTY"
    has0 = not rights and all(x(n) == 0 for n in range(1, 71))
    hasinf = not lefts and all(x(n) == 0 for n in range(-70, 0))
    return (inner, outer, has0, hasinf)


def intersect(regions):
    inner = max(r[0] for r in regions)
    outs = [r[1] for r in regions if r[1] is not None]
    outer = min(outs) if outs else None
    if outer is not None and inner >= outer:
        return "EMPTY"
    return (inner, outer, all(r[2] for r in regions), all(r[3] for r in regions))


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
    if s == r"\text{all}z":
        return (F(0), None, True, True)
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


def converges(x, rho):
    zz = rho * cmath.exp(0.53j)
    try:
        s1 = zsum(x, zz, -400, 400)
        s2 = zsum(x, zz, -800, 800)
    except OverflowError:
        return False
    return abs(s2 - s1) <= 1e-7 * max(1.0, abs(s2))


def inside(reg, rho):
    if reg == "EMPTY":
        return False
    return reg[0] < rho and (reg[1] is None or rho < reg[1])


def check(params, correct):
    probs = []
    specs = params["specs"]
    if not 2 <= len(specs) <= 3:
        probs.append(f"{len(specs)} pieces")
    reg = closed_form(specs)
    if intersect([term_region(s) for s in specs]) != reg:
        probs.append(f"pieces cancel at an end of the signal: the intersection of the pieces' ROCs is not the ROC {reg}")
    if (reg == "EMPTY") != params["empty"]:
        probs.append(f"empty flag {params['empty']} but the closed-form ROC is {reg}")
    x = lambda n: sum((sample(s, n) for s in specs), F(0))
    radii = sorted({base_radius(s) for s in specs if base_radius(s) is not None})
    probe = [F(1, 20), F(20)] + [r * F(4, 5) for r in radii] + [r * F(5, 4) for r in radii]
    probe += [F(math.isqrt(int(1e6 * a * b)), 1000) for a, b in zip(radii, radii[1:])]
    for rho in probe:
        num = converges(x, float(rho))
        if num != inside(reg, rho):
            probs.append(f"|z| = {float(rho):.4g}: truncated sums {'converge' if num else 'diverge'} but the closed-form ROC {reg} says otherwise")
    opts = params["roc_choices"]
    if len(opts) < 4:
        probs.append(f"only {len(opts)} options")
    try:
        regions = [parse_roc(o["text"]) for o in opts]
    except ValueError as exc:
        return probs + [str(exc)]
    if len(set(regions)) != len(regions):
        probs.append("two options describe the same region")
    if regions.count("EMPTY") != 1:
        probs.append("the 'does not exist' option must appear exactly once")
    good = [r for r, o in zip(regions, opts) if o["correct"]]
    if len(good) != 1:
        probs.append(f"{len(good)} options marked correct")
    elif good[0] != reg:
        probs.append(f"marked {good[0]} but the ROC is {reg}")
    if reg != "EMPTY" and parse_roc(params["roc_tex"]) != reg:
        probs.append("answer panel ROC differs")
    for pc, s in zip(params["pieces"], specs):
        if parse_roc(pc["roc_tex"]) != term_region(s):
            probs.append(f"answer panel ROC of piece {pc['i']} is wrong")
    return probs


def _key_for(params, text):
    norm = lambda s: re.sub(r"\s+", "", html.unescape(s))
    for o in params["roc"]:
        if norm(o["html"]) == norm(text):
            return o["key"]
    raise KeyError(text)


def submissions(params, correct):
    cases = []
    for o in params["roc_choices"]:
        cases.append(({"roc": _key_for(params, o["text"])}, {"roc": 1 if o["correct"] else 0}))
    return cases
