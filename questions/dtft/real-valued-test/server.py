"""Hermitian symmetry as a test for a real x[n] (Quiz 2, Oct 2026: "Hermitian symmetry as a test
for a real x[n]"; Lecture 14 symmetry properties; HW6).

(a) pl-checkbox: five spectra X_d(w) (given on -pi < w <= pi) drawn from a bank of look-alike
    families with random coefficients; select those whose x[n] is real, i.e. X_d(-w) = X_d^*(w).
    Real families and their non-real twins (same look, one symmetry broken):
        a cos kw + j b sin mw   vs  a cos kw + j b cos mw,  a cos kw + b sin mw
        c|w| e^{-jkw}           vs  c w e^{-jkw},  c e^{-jk|w|}
        j b sin kw              vs  b sin kw,  j b cos kw
        e^{-jkw}(a + b cos w)   vs  (1 + j) c cos kw
        j b sin w e^{-jkw}      vs  b sin w e^{-jkw}
        c w^2 e^{-jkw}, j c w   vs  c w
        1/(1 - a e^{-jw})       vs  1/(1 - j a e^{-jw})
        e^{jkw} + c             vs  e^{jkw} + j c
    1-4 of the five are real. The solution shows X_d(-w) and X_d^*(w) for every option.
(b) x real, X_d(w1) = a + jb (small nonzero integers, a < 0 in half the variants): Re/Im of
    X_d(-w1) = a - jb and of X_d(w2) for w2 in {2pi - w1, 2pi + w1, -2pi - w1}, and the principal phase
    of X_d(-w1) in radians (never +-pi since b != 0; quadrant II/III when a < 0)."""

import math
import random
from fractions import Fraction as Fr

from ece310 import fmt
from ece310.w1_signals import tex_pi

NZ3 = [-3, -2, -1, 1, 2, 3]
W1 = [Fr(1, 6), Fr(1, 4), Fr(1, 3), Fr(2, 5), Fr(1, 2), Fr(3, 5), Fr(2, 3), Fr(3, 4), Fr(5, 6)]


# ----------------------------------------------------------------------------- LaTeX pieces
def trig(fn, k):
    return rf"\{fn}\omega" if k == 1 else rf"\{fn} {k}\omega"


def ej(m, var=r"\omega"):
    """e^{j m var} (m a nonzero int)."""
    k = "" if abs(m) == 1 else str(abs(m))
    return "e^{" + ("-" if m < 0 else "") + "j" + k + var + "}"


def ts(*terms):
    return fmt.tex_sum(terms)


def cw(c, body):
    """c * body as a single term."""
    return fmt.tex_sum([(c, body)])


# ----------------------------------------------------------------------------- the bank
def _k():
    return random.choice([1, 2, 3])


def bank(kind):
    """(X, X(-w), X^*(w)) in LaTeX for a fresh random member of the family."""
    if kind == "cos_jsin":
        a, b, k, m = random.choice(NZ3), random.choice(NZ3), _k(), random.choice([1, 2])
        return (ts((a, trig("cos", k)), (b, "j" + trig("sin", m))),
                ts((a, trig("cos", k)), (-b, "j" + trig("sin", m))),
                ts((a, trig("cos", k)), (-b, "j" + trig("sin", m))))
    if kind == "cos_jcos":
        a, b, k, m = random.choice(NZ3), random.choice(NZ3), random.choice([1, 2]), _k()
        if k == m:
            m = k + 1
        return (ts((a, trig("cos", k)), (b, "j" + trig("cos", m))),
                ts((a, trig("cos", k)), (b, "j" + trig("cos", m))),
                ts((a, trig("cos", k)), (-b, "j" + trig("cos", m))))
    if kind == "cos_sin":
        a, b, k, m = random.choice(NZ3), random.choice(NZ3), _k(), random.choice([1, 2])
        return (ts((a, trig("cos", k)), (b, trig("sin", m))),
                ts((a, trig("cos", k)), (-b, trig("sin", m))),
                ts((a, trig("cos", k)), (b, trig("sin", m))))
    if kind == "abs_e":
        c, k = random.choice([1, 2, 3]), _k()
        return (cw(c, r"|\omega|\," + ej(-k)), cw(c, r"|\omega|\," + ej(k)), cw(c, r"|\omega|\," + ej(k)))
    if kind == "w_e":
        c, k = random.choice([1, 2, 3]), _k()
        return (cw(c, r"\omega\," + ej(-k)), cw(-c, r"\omega\," + ej(k)), cw(c, r"\omega\," + ej(k)))
    if kind == "e_jabs":
        c, k = random.choice([1, 2, 3]), _k()
        return (cw(c, ej(-k, r"|\omega|")), cw(c, ej(-k, r"|\omega|")), cw(c, ej(k, r"|\omega|")))
    if kind == "jsin":
        b, k = random.choice(NZ3), _k()
        return (cw(b, "j" + trig("sin", k)), cw(-b, "j" + trig("sin", k)), cw(-b, "j" + trig("sin", k)))
    if kind == "sin":
        b, k = random.choice(NZ3), _k()
        return (cw(b, trig("sin", k)), cw(-b, trig("sin", k)), cw(b, trig("sin", k)))
    if kind == "jcos":
        b, k = random.choice(NZ3), _k()
        return (cw(b, "j" + trig("cos", k)), cw(b, "j" + trig("cos", k)), cw(-b, "j" + trig("cos", k)))
    if kind == "e_cos":
        a, b, k = random.choice([1, 2, 3]), random.choice(NZ3), _k()
        inner = r"\left(" + ts((a, ""), (b, trig("cos", 1))) + r"\right)"
        return (ej(-k) + inner, ej(k) + inner, ej(k) + inner)
    if kind == "cplx_cos":
        c, k = random.choice([1, 2, 3]), _k()
        cc = "" if c == 1 else str(c)
        return (cc + r"(1 + j)" + trig("cos", k), cc + r"(1 + j)" + trig("cos", k), cc + r"(1 - j)" + trig("cos", k))
    if kind == "e_jsin":
        b, k = random.choice(NZ3), random.choice([1, 2])
        return (cw(b, "j" + trig("sin", 1) + r"\," + ej(-k)), cw(-b, "j" + trig("sin", 1) + r"\," + ej(k)),
                cw(-b, "j" + trig("sin", 1) + r"\," + ej(k)))
    if kind == "e_sin":
        b, k = random.choice(NZ3), random.choice([1, 2])
        return (cw(b, trig("sin", 1) + r"\," + ej(-k)), cw(-b, trig("sin", 1) + r"\," + ej(k)),
                cw(b, trig("sin", 1) + r"\," + ej(k)))
    if kind == "w2_e":
        c, k = random.choice([1, 2, 3]), _k()
        return (cw(c, r"\omega^2 " + ej(-k)), cw(c, r"\omega^2 " + ej(k)), cw(c, r"\omega^2 " + ej(k)))
    if kind == "jw":
        c = random.choice(NZ3)
        return (cw(c, r"j\omega"), cw(-c, r"j\omega"), cw(-c, r"j\omega"))
    if kind == "w":
        c = random.choice(NZ3)
        return (cw(c, r"\omega"), cw(-c, r"\omega"), cw(c, r"\omega"))
    if kind == "pole":
        a = random.choice([Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(-1, 3), Fr(2, 3)])
        f = lambda m: r"\dfrac{1}{" + ts((1, ""), (-a, ej(m))) + "}"  # noqa: E731
        return (f(-1), f(1), f(1))
    if kind == "pole_j":
        a = random.choice([Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(-1, 3), Fr(2, 3)])
        f = lambda s, m: (r"\dfrac{1}{1 " + ("+" if s * a > 0 else "-") + " j"  # noqa: E731
                          + fmt.tex_num(abs(a)) + r"\," + ej(m) + "}")
        return (f(-1, -1), f(-1, 1), f(1, 1))
    if kind == "e_plus":
        k, c = random.choice([-2, -1, 1, 2]), random.choice(NZ3)
        return (ts((1, ej(k)), (c, "")), ts((1, ej(-k)), (c, "")), ts((1, ej(-k)), (c, "")))
    if kind == "e_pj":
        k, c = random.choice([-2, -1, 1, 2]), random.choice(NZ3)
        return (ts((1, ej(k)), (c, "j")), ts((1, ej(-k)), (c, "j")), ts((1, ej(-k)), (-c, "j")))
    raise KeyError(kind)


REAL = ["cos_jsin", "abs_e", "jsin", "e_cos", "e_jsin", "w2_e", "jw", "pole", "e_plus"]
NOTREAL = ["cos_jcos", "cos_sin", "w_e", "e_jabs", "sin", "jcos", "cplx_cos", "e_sin", "w", "pole_j", "e_pj"]


# ----------------------------------------------------------------------------- (b)
def _cplx_tex(a, b):
    if b == 0:
        return str(a)
    jb = "j" if abs(b) == 1 else f"{abs(b)}j"
    return f"{a} {'-' if b < 0 else '+'} {jb}"


def _dec(v):
    return f"{v:.4f}"


def generate(data):
    p = data["params"]
    nreal = random.choice([1, 2, 2, 3, 3, 4])
    kinds = random.sample(REAL, nreal) + random.sample(NOTREAL, 5 - nreal)
    random.shuffle(kinds)
    choices, sol = [], []
    for kd in kinds:
        X, Xm, Xc = bank(kd)
        real = kd in REAL
        choices.append({"text": f"$X_d(\\omega) = {X}$", "correct": real})
        sol.append({"X_tex": X, "Xm_tex": Xm, "Xc_tex": Xc, "rel": "=" if real else r"\ne",
                    "verdict": "real" if real else "not real"})
    p["kinds"] = kinds
    p["real_choices"] = choices
    p["sol_rows"] = sol

    # (b)
    w1 = random.choice(W1)
    b = random.choice([-4, -3, -2, -1, 1, 2, 3, 4])
    a = random.choice([-4, -3, -2, -1]) if random.random() < 0.5 else random.choice([1, 2, 3, 4])
    shift = random.choices(["2pi-w", "2pi+w", "-2pi-w"], weights=[2, 1, 1])[0]
    w2 = {"2pi-w": 2 - w1, "2pi+w": 2 + w1, "-2pi-w": -2 - w1}[shift]
    conj2 = shift != "2pi+w"
    b2 = -b if conj2 else b
    ph = math.atan2(-b, a)
    p["a"], p["b"], p["w1"], p["w2"] = a, b, [w1.numerator, w1.denominator], [w2.numerator, w2.denominator]
    p["shift"] = shift
    p["w1_tex"] = tex_pi(w1)
    p["mw1_tex"] = tex_pi(-w1)
    p["w2_tex"] = tex_pi(w2)
    p["X1_tex"] = _cplx_tex(a, b)
    p["Xm_tex"] = _cplx_tex(a, -b)
    p["X2_tex"] = _cplx_tex(a, b2)
    red = {"2pi-w": r"- 2\pi", "2pi+w": r"- 2\pi", "-2pi-w": r"+ 2\pi"}[shift]
    p["w2_steps_tex"] = (rf"X_d\!\left({tex_pi(w2)}\right) = X_d\!\left({tex_pi(w2)} {red}\right) = "
                         + (rf"X_d\!\left({tex_pi(-w1)}\right) = X_d^*\!\left({tex_pi(w1)}\right)" if conj2
                            else rf"X_d\!\left({tex_pi(w1)}\right)"))
    ratio = fmt.tex_num(Fr(-b, a))
    q = Fr(ph / math.pi).limit_denominator(12)
    exact = rf"= {tex_pi(q)} " if abs(float(q) * math.pi - ph) < 1e-12 else ""
    base = rf"\arctan\!\left({ratio}\right)"
    if a > 0:
        p["ph_line_tex"] = rf"\mathrm{{Re}} \gt 0:\ \angle = {base} {exact}\approx {_dec(ph)}"
    else:
        pm = "+" if -b > 0 else "-"
        quad = "II" if -b > 0 else "III"
        p["ph_line_tex"] = (rf"\mathrm{{Re}} \lt 0 \text{{ (quadrant {quad})}}:\ \angle = {base} {pm} \pi {exact}\approx {_dec(ph)}")
    c = data["correct_answers"]
    c["rem"], c["imm"] = float(a), float(-b)
    c["re2"], c["im2"] = float(a), float(b2)
    c["ph"] = ph
