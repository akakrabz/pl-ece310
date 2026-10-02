"""Self-tests for serverFilesCourse/ece310 (run: python3 tools/test_lib.py).

Every helper is checked against an independent computation: scipy.signal.residuez for partial
fractions, truncated z-transform sums at a point inside the ROC for the signal terms, numpy
convolution for seq.conv, and SymPy parsing of every generated SymPy string."""

from __future__ import annotations

import cmath
import math
import pathlib
import random
import sys
from fractions import Fraction as F

import numpy as np
from scipy import signal

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "serverFilesCourse"))
sys.path.insert(0, str(ROOT / "tools"))
import plsim  # noqa: E402,F401  (puts sympy on the path)
import sympy  # noqa: E402
from checklib import zsum  # noqa: E402

from ece310 import fmt, mc, poly, seq, zt  # noqa: E402

z, n = sympy.symbols("z n")
fails = 0


def ok(cond, msg):
    global fails
    if not cond:
        fails += 1
        print("FAIL", msg)


def S(expr: str):
    return sympy.sympify(expr, locals={"z": z, "n": n})


NICE = [F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(3, 2), F(2), F(3), F(-1, 4), F(-1, 3), F(-1, 2), F(-2, 3), F(-3, 2), F(-2), F(-3)]

# ------------------------------------------------------------------ fmt
ok(fmt.tex_num(F(-3, 4)) == r"-\tfrac{3}{4}", "tex_num")
ok(fmt.sym_num(F(-3, 4)) == "(-3/4)", "sym_num")
ok(fmt.tex_sum([(1, "x"), (-1, "y"), (F(1, 2), "z"), (0, "w"), (-2, "")]) == r"x - y + \tfrac{1}{2}z - 2", fmt.tex_sum([(1, "x"), (-1, "y"), (F(1, 2), "z"), (0, "w"), (-2, "")]))
ok(fmt.tex_seq([1, 2, 3], -1) == r"\{1,\ \underset{\uparrow}{2},\ 3\}", "tex_seq")
ok(fmt.tex_seq([5], 2) == r"\{\underset{\uparrow}{0},\ 0,\ 5\}", fmt.tex_seq([5], 2))
ok(fmt.tex_seq([5, 6], -4) == r"\{5,\ 6,\ 0,\ 0,\ \underset{\uparrow}{0}\}", fmt.tex_seq([5, 6], -4))
for c in ([1, F(-1, 2)], [2, 0, F(1, 4)], [F(3, 2), -1, 1]):
    e = S(fmt.sym_poly_zinv(c))
    ok(sympy.simplify(e - sum(sympy.Rational(x.numerator, x.denominator) if isinstance(x, F) else x * 1 for x in [0]) - sum((sympy.Rational(F(x).numerator, F(x).denominator)) * z ** (-k) for k, x in enumerate(c))) == 0, f"sym_poly_zinv {c}")
for b in NICE + [F(1), F(-1)]:
    e = S(fmt.sym_pow(b, "n"))
    ok(abs(complex(e.subs(n, 3)) - float(b) ** 3) < 1e-12, f"sym_pow {b}")

# ------------------------------------------------------------------ poly
for _ in range(200):
    a = [F(random.randint(-5, 5), random.randint(1, 4)) for _ in range(random.randint(1, 5))]
    b = [F(random.randint(1, 5))] + [F(random.randint(-5, 5), random.randint(1, 4)) for _ in range(random.randint(0, 3))]
    if poly.trim(b) == [0]:
        continue
    q, r = poly.pdivmod(a, b)
    back = poly.padd(poly.pmul(q, b), r)
    ok(poly.trim(back) == poly.trim(a), f"pdivmod {a} / {b}")
    ok(len(poly.trim(r)) < max(len(poly.trim(b)), 2) or poly.trim(r) == [0], f"deg r {r}")

# ------------------------------------------------------------------ pfe vs residuez
for _ in range(300):
    k = random.choice([1, 2, 3])
    poles = random.sample(NICE, k)
    num = [F(random.randint(-4, 4), random.choice([1, 2])) for _ in range(random.randint(1, k + 2))]
    if all(c == 0 for c in num):
        continue
    C, A = zt.pfe(num, poles)
    den = poly.from_roots(poles)
    r, p, kk = signal.residuez([float(c) for c in num], [float(c) for c in den])
    for pole, a in zip(poles, A):
        idx = int(np.argmin(np.abs(p - float(pole))))
        ok(abs(r[idx] - float(a)) < 1e-7, f"pfe residue pole {pole}: {a} vs {r[idx]}")
    kk = list(kk) if len(kk) else [0.0]
    Cf = [float(c) for c in C] + [0.0] * 5
    ok(all(abs(Cf[i] - (kk[i] if i < len(kk) else 0.0)) < 1e-7 for i in range(len(kk))), f"pfe direct terms {C} vs {kk}")
    # recombine at a random point
    zz = complex(1.7, 0.9)
    lhs = poly.peval(num, zz) / poly.peval(den, zz)
    rhs = poly.peval(C, zz) + sum(complex(a) / (1 - complex(pp) / zz) for a, pp in zip(A, poles))
    ok(abs(lhs - rhs) < 1e-9, "pfe recombination")

# complex conjugate poles
for th in (math.pi / 3, math.pi / 2, 2 * math.pi / 3):
    for r in (0.5, 1.0, 2.0):
        poles = [r * cmath.exp(1j * th), r * cmath.exp(-1j * th)]
        num = [1.0, -r * math.cos(th)]
        C, A = zt.pfe(num, poles)
        inv = zt.InverseZ(C, poles, A, ["R", "R"])
        for nn in range(8):
            ok(abs(inv.value(nn) - r**nn * math.cos(th * nn)) < 1e-9, f"cos inverse r={r} th={th} n={nn}")

# ------------------------------------------------------------------ inverse with every ROC
for _ in range(200):
    k = random.choice([2, 3])
    poles = random.sample(NICE, k)
    if len({abs(p) for p in poles}) < k:
        continue
    num = [F(random.randint(-4, 4)) for _ in range(random.randint(1, k))]
    if all(c == 0 for c in num):
        continue
    C, A = zt.pfe(num, poles)
    den = poly.from_roots(poles)
    for roc in zt.rocs_for_poles(poles):
        inv = zt.inverse_pfe(C, poles, A, roc)
        # sum_n x[n] z^-n at a point inside the ROC == X(z)
        rr = roc.test_point()
        zz = rr * cmath.exp(0.37j)
        tot = zsum(inv.value, zz)
        X = poly.peval(num, zz) / poly.peval(den, zz)
        ok(abs(tot - X) < 1e-6 * max(1, abs(X)), f"inverse ROC {roc.tex()} poles {poles}: {tot} vs {X}")
        # sympy pieces
        e_pos, e_neg = S(inv.sym_nonneg()), S(inv.sym_neg())
        for m in (0, 1, 4):
            ok(abs(complex(e_pos.subs(n, m)) - complex(inv.value(m))) < 1e-9, "sym_nonneg")
        for m in (-1, -3):
            ok(abs(complex(e_neg.subs(n, m)) - complex(inv.value(m))) < 1e-9, "sym_neg")
        # stability / causality flags
        ok(roc.contains_unit_circle() == (all((abs(p) < 1) == (s == "R") for p, s in zip(poles, inv.sides))), "stability flag")
        ok(roc.is_causal_system() == all(s == "R" for s in inv.sides), "causality flag")

# ------------------------------------------------------------------ signal terms: truncated sums
def trunc_sum(spec, zz, lo=-800, hi=800):
    return zsum(lambda m: zt.x_value(spec, m), zz, lo, hi)


def X_of(term, zz):
    return complex(S(term.X_sym).subs(z, zz))


for _ in range(300):
    kind = random.choice(["right", "left", "n", "finite", "cos"])
    if kind == "right":
        t = zt.term_right_exp(random.choice([1, 2, -1, F(1, 2), 3]), random.choice(NICE + [F(1)]), random.randint(-3, 3))
    elif kind == "left":
        t = zt.term_left_exp(random.choice([1, 2, -1, F(1, 2), 3]), random.choice(NICE + [F(1)]), random.randint(-3, 3))
    elif kind == "n":
        t = zt.term_n_exp(random.choice([1, 2, -1]), random.choice(NICE))
    elif kind == "finite":
        vals = [random.randint(-3, 3) for _ in range(random.randint(1, 5))]
        if all(v == 0 for v in vals):
            continue
        t = zt.term_finite(vals, random.randint(-3, 3))
    else:
        t = zt.term_cos(random.choice([1, 2, -1]), random.choice([F(1, 2), F(1), F(2, 3)]), random.choice(["pi/2", "pi/3", "2*pi/3", "pi"]))
    rr = t.roc.test_point()
    zz = rr * cmath.exp(0.61j)
    tot = trunc_sum(t.spec, zz)
    ok(abs(tot - X_of(t, zz)) < 1e-6 * max(1.0, abs(tot)), f"term {t.spec} roc {t.roc.tex()}: {tot} vs {X_of(t, zz)}")
    # ROC endpoints: 0 and infinity membership
    if t.roc.inner == 0:
        has_pos = any(zt.x_value(t.spec, m) != 0 for m in range(1, 60))
        ok(t.roc.has0 == (not has_pos), f"has0 {t.spec}")
    if t.roc.outer is None:
        has_neg = any(zt.x_value(t.spec, m) != 0 for m in range(-60, 0))
        ok(t.roc.hasinf == (not has_neg), f"hasinf {t.spec}")

# ROC intersection and tex
r1 = zt.ROC.right(F(1, 2), 0)
r2 = zt.ROC.left(F(3), -1)
ok(r1.intersect(r2).tex() == r"\tfrac{1}{2} \lt \lvert z\rvert \lt 3", r1.intersect(r2).tex())
ok(zt.ROC.right(3).intersect(zt.ROC.left(F(1, 2))).empty, "empty intersection")
ok(zt.ROC.finite(-2, 3).tex() == r"0 \lt \lvert z\rvert \lt \infty", zt.ROC.finite(-2, 3).tex())
ok(zt.ROC.finite(0, 3).tex() == r"\lvert z\rvert \gt 0", "finite causal")
ok(zt.ROC.finite(-3, 0).tex() == r"\lvert z\rvert \lt \infty", "finite anticausal")
ok(zt.ROC.right(2, -1).tex() == r"2 \lt \lvert z\rvert \lt \infty", "right with negative-time sample")
ok(zt.ROC.left(2, 1).tex() == r"0 \lt \lvert z\rvert \lt 2", "left with positive-time sample")
ok(zt.ROC.from_json(r1.to_json()) == r1, "ROC json roundtrip")

# ------------------------------------------------------------------ seq / mc
for _ in range(100):
    x = [random.randint(-3, 3) for _ in range(random.randint(1, 5))]
    h = [random.randint(-3, 3) for _ in range(random.randint(1, 4))]
    y, ny = seq.conv(x, -1, h, 2)
    ok([float(v) for v in y] == list(np.convolve(x, h).astype(float)) and ny == 1, "conv")
opts = mc.choices("A", ["B", "A", "C", "B", "D", "E"], n=4)
ok(len(opts) == 4 and sum(o["correct"] for o in opts) == 1 and len({o["text"] for o in opts}) == 4, "mc.choices")

print("ece310 library self-test:", "OK" if fails == 0 else f"{fails} failure(s)")
sys.exit(1 if fails else 0)
