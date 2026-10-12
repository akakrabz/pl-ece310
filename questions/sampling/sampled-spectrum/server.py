"""The "sketch X_d(omega)" exam problem answered with numbers (Lecture 17, ideal A/D conversion).

X_a(Omega) is a triangle, trapezoid, rectangle or V of height A and band edge Omega_B = 2 pi B, shown
as a figure; x[n] = x_a(nT). Asked: the band edge omega_B = Omega_B T of the copy centred at 0 (as a
multiple of pi), X_d(0), X_d(pi), X_d(omega_1) at one more frequency, and whether the copies are
separated, touch, or overlap (above / at / below the Nyquist rate).

Models: Lecture 17 Example 1 (box at T = 1/20, 1/10 s) and Example 2 (triangle at 150, 100, 75 Hz),
FA2024 MT2 #5(b) (triangle at T_N and 2T_N), SP2025 MT2 #5(a) (T = 1/300 s), SP2021 MT2 #8(b) (V at
5 kHz); about one variant in seven uses exactly those numbers, the others are random with b = 2BT a
simple fraction. Values are exact (ece310.w6_sampling); a query point is never a jump of X_d."""

import random
from fractions import Fraction as F

from ece310 import fmt, mc
from ece310 import w6_sampling as w

NAMED = [
    ("rect", 1, 5, 20, None, "Lecture 17, Example 1(a)"),
    ("rect", 1, 5, 10, None, "Lecture 17, Example 1(b)"),
    ("tri", 1, 50, 150, None, "Lecture 17, Example 2(a)"),
    ("tri", 1, 50, 100, None, "Lecture 17, Example 2(b)"),
    ("tri", 1, 50, 75, None, "Lecture 17, Example 2(c)"),
    ("tri", 1, 500, 1000, None, "FA2024 Midterm 2 #5(b) with $T = T_N$"),
    ("tri", 1, 500, 500, None, "FA2024 Midterm 2 #5(b) with $T = 2T_N$"),
    ("tri", 1, 100, 300, None, "SP2025 Midterm 2 #5(a)"),
    ("vee", 1, 5000, 5000, None, "SP2021 Midterm 2 #8(b)"),
]
B_CHOICES = [5, 10, 20, 25, 50, 100, 150, 200, 250, 300, 400, 500, 1000, 1500, 2000, 2500, 4000, 5000]
B_ABOVE = [F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(5, 6), F(2, 5), F(3, 5), F(1, 4)]
B_BELOW = [F(7, 6), F(5, 4), F(4, 3), F(3, 2), F(5, 3), F(2), F(5, 2), F(6, 5)]
U_CHOICES = sorted({F(k, 12) for k in range(1, 12)} | {F(k, 10) for k in range(1, 10)})
REGIMES = ["above", "at", "below"]
REGIME_TEXT = {
    "above": r"Above the Nyquist rate: a guard band separates the copies (no aliasing).",
    "at": r"Exactly at the Nyquist rate: the copies touch at $\omega = \pm\pi$ but do not overlap (no aliasing).",
    "below": r"Below the Nyquist rate: the copies overlap and add (aliasing).",
}


def _nice(v):
    return v.denominator in (1, 2, 3, 4, 5, 6, 8)


def _random_spec():
    while True:
        shape = random.choices(["tri", "trap", "rect", "vee"], weights=[4, 4, 0.4, 1.5])[0]
        reg = random.choices(REGIMES, weights=[4, 2, 4])[0]
        b = random.choice(B_ABOVE) if reg == "above" else (F(1) if reg == "at" else random.choice(B_BELOW))
        if shape in ("rect", "vee") and b > 2:
            continue
        B = F(random.choice(B_CHOICES))
        fs = 2 * B / b
        if fs.denominator != 1 or fs < 10:
            continue
        A = F(random.choice([1, 1, 2, 3, 4]))
        r = None
        if shape == "trap":
            r = random.choice([F(1, 4), F(1, 3), F(1, 2), F(2, 3)])
            if (2 * B * r).denominator != 1:
                continue
        return w.Spec(shape, A, B, fs, r), None


def _pick():
    while True:
        if random.random() < 0.14:
            shape, A, B, fs, r, src = random.choice(NAMED)
            sp = w.Spec(shape, A, B, fs, r)
        else:
            sp, src = _random_spec()
        X0, Xpi = sp.value(0), sp.value(1)
        if X0 is None or Xpi is None or not (_nice(X0) and _nice(Xpi)):
            continue
        cands = [u for u in U_CHOICES if sp.value(u) is not None and _nice(sp.value(u))]
        if src is None and all(sp.value(u) == X0 for u in cands) and Xpi == X0:
            continue                     # constant X_d: only the named exam/lecture instances keep it
        # X_d(omega_1) should differ from both X_d(0) and X_d(pi); prefer a point where the overlap matters.
        # A box sum takes only two values on [0, pi]: then at least differ from X_d(0).
        other = [u for u in cands if sp.value(u) not in (X0, Xpi)]
        pref = [u for u in other if sp.central_only(u) is not None and sp.central_only(u) != sp.value(u)]
        pool = (pref if pref and random.random() < 0.8 else other) or [u for u in cands if sp.value(u) != X0] or cands
        if not pool:
            continue
        return sp, src, random.choice(pool)


def generate(data):
    sp, src, u1 = _pick()
    X0, Xpi, X1 = sp.value(0), sp.value(1), sp.value(u1)
    reg = sp.regime()
    p = data["params"]
    p["spec"] = {"shape": sp.shape, "A": str(sp.A), "B": str(sp.B), "fs": str(sp.fs),
                 "r": str(sp.r) if sp.r is not None else None}
    p["u1"] = str(u1)
    p["source_html"] = src or ""
    p["named"] = src is not None
    p["fig_html"] = sp.figure()
    p["words_html"] = sp.words()
    p["xa_tex"] = sp.tex_xa()
    p["T_tex"] = sp.tex_T()
    p["u1_tex"] = w.tex_pi(u1)
    p["jumpy"] = sp.shape in ("rect", "vee")
    p["regime_choices"] = [{"text": REGIME_TEXT[k], "correct": k == reg} for k in REGIMES]

    # worked solution
    WB = sp.tex_omega()
    p["B_tex"] = w.tex_dec(sp.B)
    p["twoB_tex"] = w.tex_dec(2 * sp.B)
    p["fs_tex"] = w.tex_dec(sp.fs)
    p["cmp_tex"] = {"above": r"\gt", "at": "=", "below": r"\lt"}[reg]
    p["regime_html"] = {"above": "above the Nyquist rate: the copies are separated by a guard band",
                        "at": "exactly at the Nyquist rate: the copies touch at odd multiples of $\\pi$ without overlapping",
                        "below": "below the Nyquist rate: neighbouring copies overlap and must be added"}[reg]
    p["WB_tex"] = WB
    p["b_tex"] = w.tex_pi(sp.b)
    p["b_plain"] = fmt.plain_num(sp.b)
    p["b_work_tex"] = rf"{WB}\cdot\tfrac{{1}}{{{w.tex_dec(sp.fs)}}} = {w.tex_pi(sp.b)}"
    p["H_work_tex"] = rf"\tfrac{{1}}{{T}}\cdot{w.tex_dec(sp.A)} = {w.tex_dec(sp.fs)}\cdot{w.tex_dec(sp.A)} = {w.tex_dec(sp.H)}"
    p["central_tex"] = w.central_tex(sp)
    for key, u, nm in (("X0", F(0), r"X_d(0)"), ("Xpi", F(1), r"X_d(\pi)"),
                       ("X1", u1, rf"X_d\!\left({w.tex_pi(u1)}\right)")):
        p[f"{key}_line"] = w.value_tex(sp, u, nm)
        p[f"{key}_note_html"] = w.value_note_html(sp, u, nm)
    p["trap_note"] = sp.shape == "trap"
    if sp.shape == "trap":
        p["c_work_tex"] = rf"{sp.tex_omega(sp.r)}\cdot\tfrac{{1}}{{{w.tex_dec(sp.fs)}}} = {w.tex_pi(sp.c)}"

    c = data["correct_answers"]
    c["wB"] = float(sp.b)
    c["X0"] = float(X0)
    c["Xpi"] = float(Xpi)
    c["X1"] = float(X1)
