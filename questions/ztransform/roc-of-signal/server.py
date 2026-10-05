"""ROC of a signal made of two or three standard pieces (Lectures 6-7 ROC rules; HW3 #1; exam family
"z-transform with ROC" and the True/False regulars about empty ROCs).

Pieces: right-sided exponentials A a^n u[n-k] (k may be negative: z = infinity excluded),
left-sided exponentials B b^n u[-n+m] (m may be positive: z = 0 excluded), n a^n u[n] and short
finite stretches c delta[n - n0]. In about 20 % of the variants the ROCs of the pieces do not
overlap and the correct answer is "the z-transform does not exist"; otherwise that option is a
distractor. Distractors toggle z = 0 / z = infinity, flip the side, or use a wrong radius."""

import random
from fractions import Fraction as F

from ece310 import fmt, mc, zt

MAGS = [F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(4, 3), F(3, 2), F(2), F(3)]
COEFS = [F(1), F(1), F(2), F(3), F(-1), F(-2), F(1, 2)]
LZ = r"\lvert z\rvert"
DNE = "The z-transform does not exist (the ROC is empty)."


def _base(lo=None, hi=None, allow_one=True, exclude=()):
    pool = [m for m in MAGS + ([F(1)] if allow_one else [])
            if (lo is None or m > lo) and (hi is None or m < hi) and m not in exclude]
    m = random.choice(pool)
    return -m if (m != 1 and random.random() < 0.3) else m


def _right(a=None, ks=(-3, -2, -1, 0, 0, 1, 2, 3)):
    a = _base() if a is None else a
    return zt.term_right_exp(random.choice(COEFS), a, random.choice(ks))


def _left(b=None, ms=(-2, -1, -1, 0, 1, 2)):
    b = _base() if b is None else b
    return zt.term_left_exp(random.choice(COEFS), b, random.choice(ms))


def _finite(starts=(-3, -2, -1, 1, 2, 3)):
    nz = [-3, -2, -1, 1, 2, 3]
    if random.random() < 0.5:
        return zt.term_finite([random.choice(nz)], random.choice(starts))
    vals = [random.choice(nz), random.randint(-3, 3), random.choice(nz)][: random.choice([2, 3])]
    vals[-1] = random.choice(nz)
    return zt.term_finite(vals, random.choice([-3, -2, -1, 0, 1, 2]))


def _note(t):
    s = t.spec
    if s["type"] == "right_exp":
        a, k = F(s["a"]), s["k"]
        tail = (rf"starts at $n = {k} \lt 0$, so $z = \infty$ is excluded" if k < 0
                else rf"starts at $n = {k} \ge 0$, so $z = \infty$ is included")
        return rf"right-sided, pole at $z = {fmt.tex_num(a)}$: outside the circle of radius ${fmt.tex_num(abs(a))}$; " + tail
    if s["type"] == "left_exp":
        b, m = F(s["b"]), s["m"]
        tail = (rf"ends at $n = {m} \gt 0$, so $z = 0$ is excluded" if m > 0
                else rf"ends at $n = {m} \le 0$, so $z = 0$ is included")
        return rf"left-sided, pole at $z = {fmt.tex_num(b)}$: inside the circle of radius ${fmt.tex_num(abs(b))}$; " + tail
    if s["type"] == "n_exp":
        a = F(s["a"])
        return rf"right-sided from $n = 0$, double pole at $z = {fmt.tex_num(a)}$: outside the circle of radius ${fmt.tex_num(abs(a))}$"
    first, last = s["start"], s["start"] + len(s["values"]) - 1
    parts = [r"samples at $n \gt 0$ exclude $z = 0$" if last > 0 else r"no samples at $n \gt 0$, so $z = 0$ is included",
             r"samples at $n \lt 0$ exclude $z = \infty$" if first < 0 else r"no samples at $n \lt 0$, so $z = \infty$ is included"]
    return "finite length, converges everywhere except possibly $z = 0$ or $z = \\infty$: " + "; ".join(parts)


# ----------------------------------------------------------------------------- ROC distractors
def roc_distractors(roc, radii):
    """Wrong ROCs ordered by priority: 0/infinity toggled, side flipped, then wrong radii."""
    R = zt.ROC
    first, rest = [], []
    k = roc.kind()
    if k == "right":
        r = roc.inner
        first.append(R(inner=r, outer=None, has0=False, hasinf=not roc.hasinf))
        first.append(R(inner=F(0), outer=r, has0=True, hasinf=False))
        for w in [1 / r] + list(radii):
            if w != r and w != 0:
                rest.append(R(inner=w, outer=None, has0=False, hasinf=roc.hasinf))
        rest.append(R(inner=F(0), outer=r, has0=False, hasinf=False))
        rest.append(R(inner=F(0), outer=None, has0=False, hasinf=roc.hasinf))
    elif k == "left":
        r = roc.outer
        first.append(R(inner=F(0), outer=r, has0=not roc.has0, hasinf=False))
        first.append(R(inner=r, outer=None, has0=False, hasinf=True))
        for w in [1 / r] + list(radii):
            if w != r and w != 0:
                rest.append(R(inner=F(0), outer=w, has0=roc.has0, hasinf=False))
        rest.append(R(inner=r, outer=None, has0=False, hasinf=False))
        rest.append(R(inner=F(0), outer=None, has0=roc.has0, hasinf=False))
    elif k == "two-sided":
        r1, r2 = roc.inner, roc.outer
        first.append(R(inner=r1, outer=None, has0=False, hasinf=True))
        first.append(R(inner=F(0), outer=r2, has0=True, hasinf=False))
        cands = sorted({1 / r1, 1 / r2, r1, r2} | set(radii))
        for lo in cands:
            for hi in cands:
                if 0 < lo < hi and (lo, hi) != (r1, r2):
                    rest.append(R(inner=lo, outer=hi, has0=False, hasinf=False))
        rest.append(R(inner=r1, outer=None, has0=False, hasinf=False))
        rest.append(R(inner=F(0), outer=r2, has0=False, hasinf=False))
    elif k == "finite":
        for h0, hi in [(not roc.has0, roc.hasinf), (roc.has0, not roc.hasinf), (not roc.has0, not roc.hasinf)]:
            first.append(R(inner=F(0), outer=None, has0=h0, hasinf=hi))
        for w in list(radii) + [F(1)]:
            rest.append(R(inner=w, outer=None, has0=False, hasinf=roc.hasinf))
    random.shuffle(rest)
    return first + rest


def empty_distractors(rights, lefts):
    """Plausible regions when the right-sided pieces need |z| > A and the left-sided ones |z| < B <= A."""
    R = zt.ROC
    A, B = rights, lefts
    out = []
    if B < A:
        out.append(R(inner=B, outer=A, has0=False, hasinf=False))     # the gap between the circles
    out.append(R(inner=A, outer=None, has0=False, hasinf=True))       # right-sided pieces only
    out.append(R(inner=F(0), outer=B, has0=True, hasinf=False))       # left-sided pieces only
    out.append(R(inner=A, outer=None, has0=False, hasinf=False))
    out.append(R(inner=F(0), outer=B, has0=False, hasinf=False))
    out.append(R(inner=F(0), outer=None, has0=False, hasinf=False))
    return out[:3] + random.sample(out[3:], len(out) - 3)


def _signal(want_empty):
    if want_empty:
        b = _base(hi=F(3))
        if random.random() < 0.15:
            a = -b if random.random() < 0.5 else b                       # equal radii: still no overlap
        else:
            a = _base(lo=abs(b))
        terms = [_right(a, ks=(-1, 0, 0, 1, 2)), _left(b, ms=(-1, -1, 0, 1))]
        if random.random() < 0.4:
            terms.append(_finite())
        return "empty", terms
    fam = random.choice(["RL", "RL", "RLD", "RD", "RD", "RRD", "LD", "LD", "LLD", "DD", "NL", "ND"])
    if fam in ("RL", "RLD", "NL"):
        b = _base(lo=F(1, 4))
        a = _base(hi=abs(b))
        first = _right(a) if fam != "NL" else zt.term_n_exp(random.choice([F(1), F(2), F(-1)]), a)
        terms = [first, _left(b)]
        if fam == "RLD":
            terms.append(_finite())
    elif fam == "RD":
        terms = [_right(), _finite()]
    elif fam == "RRD":
        a1 = _base()
        terms = [_right(a1), _right(_base(exclude=(abs(a1),))), _finite()]
    elif fam == "LD":
        terms = [_left(), _finite()]
    elif fam == "LLD":
        b1 = _base()
        terms = [_left(b1), _left(_base(exclude=(abs(b1),))), _finite()]
    elif fam == "DD":
        terms = [_finite(), _finite()]
    else:  # ND
        terms = [zt.term_n_exp(random.choice([F(1), F(2), F(-1)]), _base(allow_one=False)), _finite()]
    if fam in ("RRD", "LLD") and random.random() < 0.4:
        terms = terms[:2]
    return fam, terms


def ends_agree(terms, roc):
    """Linearity only promises 'ROC at least the intersection': if samples of different pieces cancel
    at an end of the summed signal, the true ROC can gain z = 0 or z = infinity. True when the summed
    signal's first/last nonzero index agrees with the intersection's 0/infinity membership."""
    if roc.empty:
        return True                     # infinite tails of opposite sides never cancel
    x = {n: sum((zt.x_value(t.spec, n) for t in terms), F(0)) for n in range(-20, 21)}
    if roc.inner == 0 and roc.has0 != all(x[n] == 0 for n in range(1, 21)):
        return False                    # no right-sided infinite piece: z = 0 in ROC iff no sample at n > 0
    if roc.outer is None and roc.hasinf != all(x[n] == 0 for n in range(-20, 0)):
        return False                    # no left-sided infinite piece: z = inf in ROC iff no sample at n < 0
    return True


def generate(data):
    want_empty = random.random() < 0.2
    while True:
        fam, terms = _signal(want_empty)
        random.shuffle(terms)
        if any(t.kind == "finite" for t in terms) and sum(t.kind == "finite" for t in terms) > 1:
            # two finite pieces must not overlap in index range (keeps the signal readable)
            fs = [t.spec for t in terms if t.kind == "finite"]
            ranges = [set(range(s["start"], s["start"] + len(s["values"]))) for s in fs]
            if ranges[0] & ranges[1]:
                continue
        x_tex, _, _, roc = zt.sum_terms(terms)
        if roc.empty != want_empty or not ends_agree(terms, roc):
            continue
        radii = sorted({zt.mag(p) for t in terms for p in t.poles})
        if want_empty:
            A = max(t.roc.inner for t in terms)
            B = min(t.roc.outer for t in terms if t.roc.outer is not None)
            regions = empty_distractors(A, B)
            ch = mc.choices(DNE, ["$" + r.tex() + "$" for r in regions], n=5)
        else:
            regions = roc_distractors(roc, radii)
            ch = mc.choices("$" + roc.tex() + "$", [DNE] + ["$" + r.tex() + "$" for r in regions], n=5)
        if len(ch) >= 4:
            break

    p = data["params"]
    p["fam"] = fam
    p["x_tex"] = x_tex
    p["specs"] = [t.spec for t in terms]
    p["pieces"] = [{"i": i + 1, "x_tex": t.x_tex, "roc_tex": t.roc.tex(), "note": _note(t)} for i, t in enumerate(terms)]
    p["empty"] = roc.empty
    p["roc_tex"] = roc.tex()
    if roc.empty:
        A = max(t.roc.inner for t in terms)
        B = min(t.roc.outer for t in terms if t.roc.outer is not None)
        p["why_tex"] = (rf"the right-sided pieces need ${LZ} \gt {fmt.tex_num(A)}$ and the left-sided pieces need "
                        rf"${LZ} \lt {fmt.tex_num(B)}$; since ${fmt.tex_num(B)} \le {fmt.tex_num(A)}$ no $z$ satisfies both")
    p["roc_choices"] = ch
