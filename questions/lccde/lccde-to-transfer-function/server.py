"""LCCDE <-> H(z) in both directions (Lecture 9 §1.1; HW4 #6(a); past exams FA2025 #6a, SP2025 #8a,
FA2023 #6b, FA2019 #10).

Direction "to_H": a causal LCCDE (Lecture 9 form, recursion form, or Lecture 9 form multiplied by an
integer so that a0 != 1) -> H(z) (pl-symbolic-input in z; any equivalent rational function scores 1).

Direction "to_ba": H(z) (ratio of polynomials in z^{-1}, possibly with a0 != 1; in positive powers
of z; or factored) -> coefficient vectors b = [b0 b1 ...] and a = [1 a1 ...] of the Lecture 9 form,
NORMALISED so that a[0] = 1 (stated in the prompt; unnormalised vectors score 0). b starts at b0
even when b0 = 0 (positive-power H with a lower-degree numerator). All direction-"to_ba"
coefficients are dyadic rationals, so they have short exact decimals (pl-matrix-input does not
accept fractions)."""

import math
import random
from fractions import Fraction as F

import numpy as np
import prairielearn as pl
from ece310 import fmt, lccde, poly

A_GEN = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(1, 6), F(-1, 6), F(3, 4), F(-3, 4),
         F(1), F(-1), F(2), F(-2), F(3, 2), F(-3, 2), F(2, 3), F(-2, 3), F(4, 3), F(-4, 3)]
B_GEN = [F(v) for v in (1, 2, 3, -1, -2, -3)] + [F(1, 2), F(-1, 2), F(3, 2)]
A_DY = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4), F(3, 4), F(-3, 4), F(1), F(-1), F(1, 8), F(-1, 8),
        F(3, 2), F(-3, 2), F(2), F(-2)]
B_DY = [F(v) for v in (1, 2, 3, -1, -2, -3)] + [F(1, 2), F(-1, 2), F(3, 2), F(-3, 2)]
POLES_DY = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4), F(3, 4), F(-3, 4), F(2), F(-2), F(3, 2), F(-3, 2)]
ZEROS_DY = [F(1), F(-1), F(2), F(-2), F(1, 2), F(-1, 2), F(3), F(-3), F(1, 4), F(-3, 2)]


def _poly_or_zero(choices, n):
    """n coefficients, first and last nonzero."""
    c = [random.choice(choices) for _ in range(n)]
    if n > 1 and random.random() < 0.3:
        c[1:-1] = [F(0)] * (n - 2)
    return c


def _sym_pos(b, a):
    """H in positive powers of z as a SymPy string (an equivalent form for the probes)."""
    top = max(len(b), len(a)) - 1

    def s(c):
        terms = [f"{fmt.sym_num(v)}*z**({top - k})" for k, v in enumerate(c) if v != 0]
        return "(" + " + ".join(terms) + ")"
    return f"{s(b)}/{s(a)}"


def _factors(roots, bare):
    """prod (1 - r z^{-1}); a single factor is shown without parentheses when bare=True."""
    if bare and len(roots) == 1:
        return fmt.tex_sum([(1, ""), (-roots[0], "z^{-1}")])
    return "".join(fmt.tex_factor(r) for r in roots)


def _vec(c):
    return pl.to_json(np.array([[float(v) for v in c]]))


def _plain_vec(c):
    return "[" + ", ".join(repr(float(v)) if float(v) != int(v) else str(int(v)) for v in c) + "]"


def _to_H(p, ca):
    N = random.choice([1, 2, 2])
    a = [F(1)] + [random.choice(A_GEN) for _ in range(N)]
    if N == 2 and random.random() < 0.25:
        a[1] = F(0)
    nb = random.choice([1, 2, 2, 3, 3])
    b = _poly_or_zero(B_GEN, nb)
    form = random.choice(["std", "rec", "scaled"])
    scale = 1
    if form == "scaled":
        scale = math.lcm(*[v.denominator for v in a + b])
        if scale == 1:
            scale = random.choice([2, 3])
    p["eq_tex"] = lccde.tex_lccde(b, a, "rec" if form == "rec" else "std", scale=scale)
    p["std_tex"] = lccde.tex_lccde(b, a, "std")
    p["scaled"] = form == "scaled"
    p["rec"] = form == "rec"
    p["scale"] = scale
    p["Y_tex"] = lccde.tex_poly_zinv(a)
    p["X_tex"] = lccde.tex_poly_zinv(b)
    p["H_tex"] = lccde.tex_H(b, a)
    p["Hpos_tex"] = lccde.tex_H(b, a, positive=True)
    p["b"], p["a"] = [str(v) for v in b], [str(v) for v in a]
    p["H_pos_sym"] = _sym_pos(b, a)
    ca["H"] = lccde.sym_H(b, a)


def _to_ba(p, ca):
    form = random.choice(["zinv", "zpos", "factored"])
    if form == "factored":
        while True:
            npole = random.choice([1, 2, 2])
            nzero = random.choice([0, 1, 1, 2]) if npole == 2 else random.choice([0, 1])
            poles = random.sample(POLES_DY, npole)
            zeros = random.sample(ZEROS_DY, nzero)
            if any(z == q for z in zeros for q in poles):
                continue
            if npole == 2 and abs(poles[0]) == abs(poles[1]) and random.random() < 0.7:
                continue
            break
        K = F(random.choice([1, 1, 2, 3, -1, -2]))
        b = [K * c for c in poly.from_roots(zeros)]
        a = list(poly.from_roots(poles))
        num_tex = (fmt.tex_num(K) if K != 1 or not zeros else "") + _factors(zeros, bare=(K == 1))
        if K == -1 and zeros:
            num_tex = "-" + _factors(zeros, bare=False)
        p["Hgiven_tex"] = fmt.tex_frac(num_tex or "1", _factors(poles, bare=True))
        c = 1
    else:
        N = random.choice([1, 2, 2])
        a = [F(1)] + [random.choice(A_DY) for _ in range(N)]
        if N == 2 and random.random() < 0.25:
            a[1] = F(0)
        if form == "zpos":
            m = random.choice(range(0, N)) if random.random() < 0.5 else N    # numerator degree in z
            num = _poly_or_zero(B_DY, m + 1)
            b = [F(0)] * (N - m) + num
            c = 1
            p["Hgiven_tex"] = lccde.tex_H(b, a, positive=True)
        else:
            nb = random.choice([1, 2, 2, 3])
            b = _poly_or_zero(B_DY, nb)
            c = 1
            if random.random() < 0.5:     # show H with a0 != 1: scale to integer coefficients
                c = math.lcm(*[v.denominator for v in a + b])
                c = c if c > 1 else random.choice([2, 3])
            p["Hgiven_tex"] = lccde.tex_H([c * v for v in b], [c * v for v in a])
        if not lccde.coprime(b, a):
            return _to_ba(p, ca)
    p["form_zpos"] = form == "zpos"
    p["form_factored"] = form == "factored"
    p["form_zinv"] = form == "zinv"
    p["c"] = c
    p["c_ne_1"] = c != 1
    p["lead_zeros"] = next(i for i, v in enumerate(b) if v != 0)
    p["has_lead_zeros"] = p["lead_zeros"] > 0
    p["zlead_tex"] = f"z^{{-{p['lead_zeros']}}}"
    p["N"] = len(a) - 1
    p["H_tex"] = lccde.tex_H(b, a)
    p["std_tex"] = lccde.tex_lccde(b, a, "std")
    p["b"], p["a"] = [str(v) for v in b], [str(v) for v in a]
    p["b_plain"], p["a_plain"] = _plain_vec(b), _plain_vec(a)
    p["b_tex"] = r"\left[" + r",\ ".join(fmt.tex_num(v) for v in b) + r"\right]"
    p["a_tex"] = r"\left[" + r",\ ".join(fmt.tex_num(v) for v in a) + r"\right]"
    ca["bvec"] = _vec(b)
    ca["avec"] = _vec(a)


def generate(data):
    p = data["params"]
    to_H = random.random() < 0.5
    p["to_H"] = to_H
    p["to_ba"] = not to_H
    # every template field exists in both directions (the unused ones stay empty)
    for k in ("eq_tex", "std_tex", "Y_tex", "X_tex", "H_tex", "Hpos_tex", "Hgiven_tex", "zlead_tex",
              "b_tex", "a_tex", "b_plain", "a_plain"):
        p[k] = ""
    for k in ("rec", "scaled", "form_zpos", "form_factored", "form_zinv", "c_ne_1", "has_lead_zeros"):
        p[k] = False
    p["scale"], p["c"], p["lead_zeros"] = 1, 1, 0
    if to_H:
        _to_H(p, data["correct_answers"])
    else:
        _to_ba(p, data["correct_answers"])
