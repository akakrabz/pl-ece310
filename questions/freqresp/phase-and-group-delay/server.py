"""Magnitude, principal phase, group delay and the pi jumps of a frequency response given in closed
form (Lecture 16 sections 2, 4 and 7; SP2021 Midterm 2 #2, FA2023 Midterm 2 #4(a)).

Families (R = the real factor, so that H_d = (unit-modulus factor) * R):
  sin   g e^{-jMw} sin(Kw)          SP2021 #2 (an imaginary sequence); tau = M
  cos   g e^{-jMw} cos(Kw)          tau = M
  even  g (b + cos w) e^{js|w|}     FA2023 #4 (even phase, complex h); tau = -s for w > 0
  iir   g / (1 - a e^{-jw}), a = +-1/2: no real factor, no jumps;
        tau(w) = (a cos w - a^2) / (1 - 2a cos w + a^2)
Asked: |H_d| and the principal angle at w_a (R > 0) and at w_b (R < 0: the +-pi trap; for the IIR
w_b is -w_a or pi), the group delay at w_c, and which listed frequencies in (0, pi) carry a pi jump
(pl-checkbox, net-correct; the distractors are the 2pi wraps of the principal value first, then
other frequencies; "none of the listed frequencies" is an option). Angles are never exactly +-pi.
(The IIR uses a = +-1/2 only: for a = +-1/3 no frequency except 0 and pi has a pi-rational phase.)"""

import cmath
import math
import random
from fractions import Fraction as F

import numpy as np
import sympy as sp
from ece310 import fmt
from ece310.w5_freqresp import Polar, cos_pi, dec, plain, principal, rat, sin_pi, tex_ej, tex_pi, tex_r, tex_real

NICE = [F(k, 12) for k in range(1, 12) if F(k, 12).denominator in (2, 3, 4, 6)]   # pi/6 ... 5pi/6
NONE_TEXT = "none of the listed frequencies"


def _kw(k):
    return r"\omega" if k == 1 else f"{k}\\omega"


def _lin_tex(M):
    return r"e^{-j\omega}" if M == 1 else f"e^{{-j{M}\\omega}}"


def _gtex(g):
    return "" if g == 1 else fmt.tex_num(g)


# ----------------------------------------------------------------------------- families
def _draw_family():
    fam = random.choices(["sin", "cos", "even", "iir"], weights=[30, 22, 26, 22])[0]
    if fam in ("sin", "cos"):
        M = random.choice([1, 2, 3, 4])
        K = random.choice([1, 2, 3]) if fam == "sin" else random.choice([1, 2, 3])
        g = random.choice([1, 2, 3]) if fam == "sin" else random.choice([1, 2])
        fn = r"\sin" if fam == "sin" else r"\cos"
        H_tex = f"{_gtex(g)}{_lin_tex(M)}\\,{fn}({_kw(K)})"
        return dict(fam=fam, M=M, K=K, g=g, H_tex=H_tex, R_tex=f"{_gtex(g)}{fn}({_kw(K)})",
                    domain="")
    if fam == "even":
        b = random.choice([F(1, 2), F(-1, 2), F(0)])
        g = random.choice([1, 2]) if b != 0 else random.choice([1, 2, 3])
        s = random.choice([1, -1])
        R_tex = fmt.tex_sum([(g * b, ""), (g, r"\cos\omega")])
        ex = r"e^{j|\omega|}" if s == 1 else r"e^{-j|\omega|}"
        H_tex = (rf"\left({R_tex}\right){ex}" if b != 0 else f"{_gtex(g)}\\cos\\omega\\,{ex}")
        return dict(fam=fam, b=plain(b), g=g, s=s, H_tex=H_tex, R_tex=R_tex,
                    domain=r",\qquad -\pi \le \omega \le \pi")
    a = random.choice([F(1, 2), F(-1, 2)])
    g = random.choice([F(1), F(2), F(1, 2), F(3, 2)])
    den = fmt.tex_sum([(1, ""), (-a, r"e^{-j\omega}")])
    return dict(fam="iir", a=plain(a), g=plain(g), H_tex=fmt.tex_frac(fmt.tex_num(g), den), R_tex="",
                domain="")


def _R(sysd, f):
    """The real factor at w = f pi (exact), or None for the IIR."""
    fam = sysd["fam"]
    if fam == "sin":
        return sysd["g"] * sin_pi(sysd["K"] * f)
    if fam == "cos":
        return sysd["g"] * cos_pi(sysd["K"] * f)
    if fam == "even":
        return sysd["g"] * (rat(F(sysd["b"])) + cos_pi(f))
    return None


def _lin(sysd, f):
    """Angle (multiple of pi) of the unit-modulus factor at f."""
    fam = sysd["fam"]
    if fam in ("sin", "cos"):
        return -sysd["M"] * f
    if fam == "even":
        return sysd["s"] * abs(f)
    raise ValueError


def _H_num(sysd, w):
    fam = sysd["fam"]
    if fam == "sin":
        return sysd["g"] * cmath.exp(-1j * sysd["M"] * w) * math.sin(sysd["K"] * w)
    if fam == "cos":
        return sysd["g"] * cmath.exp(-1j * sysd["M"] * w) * math.cos(sysd["K"] * w)
    if fam == "even":
        return sysd["g"] * (float(F(sysd["b"])) + math.cos(w)) * cmath.exp(1j * sysd["s"] * abs(w))
    return float(F(sysd["g"])) / (1 - float(F(sysd["a"])) * cmath.exp(-1j * w))


def _polar(sysd, f):
    """(magnitude exact, principal angle Fraction, raw angle Fraction, R exact or None)."""
    if sysd["fam"] == "iir":
        a, g = rat(F(sysd["a"])), rat(F(sysd["g"]))
        d = 1 - a * cos_pi(f) + sp.I * a * sin_pi(f)
        mag = sp.nsimplify(sp.simplify(g / sp.sqrt(sp.expand((1 - a * cos_pi(f)) ** 2 + (a * sin_pi(f)) ** 2))))
        ang = -cmath.phase(complex(d)) / math.pi
        t = F(ang).limit_denominator(12)
        assert abs(float(t) - ang) < 1e-12
        return mag, principal(t), t, None
    R = sp.nsimplify(_R(sysd, f))
    raw = _lin(sysd, f) + (1 if float(R) < 0 else 0)
    return sp.Abs(R), principal(raw), raw, R


def _tau(sysd, f):
    fam = sysd["fam"]
    if fam in ("sin", "cos"):
        return F(sysd["M"])
    if fam == "even":
        return F(-sysd["s"]) if f > 0 else F(sysd["s"])
    a = F(sysd["a"])
    c = {F(0): F(1), F(1, 3): F(1, 2), F(1, 2): F(0), F(2, 3): F(-1, 2), F(1): F(-1)}[f]
    return (a * c - a * a) / (1 - 2 * a * c + a * a)


def _breaks(sysd):
    """(jumps, wraps): frequencies in (0, pi) (Fractions with denominator | 12) where the principal
    phase jumps by pi (sign change of R) or wraps by 2 pi."""
    fam = sysd["fam"]
    jumps, wraps = [], []
    if fam == "iir":
        return jumps, wraps
    for k in range(1, 12):
        f = F(k, 12)
        w = math.pi * float(f)
        eps = 1e-7
        lo, hi = _H_num(sysd, w - eps), _H_num(sysd, w + eps)
        if abs(lo) < 1e-12 or abs(hi) < 1e-12:
            continue
        d = cmath.phase(hi) - cmath.phase(lo)
        if abs(_R(sysd, f)) == 0 and float(_R(sysd, f - F(1, 1000))) * float(_R(sysd, f + F(1, 1000))) < 0:
            jumps.append(f)
        elif abs(abs(d) - 2 * math.pi) < 1e-3:
            wraps.append(f)
    return jumps, wraps


# ----------------------------------------------------------------------------- generate
def _choose(sysd):
    fam = sysd["fam"]
    if fam == "iir":
        fa = F(1, 3) if F(sysd["a"]) > 0 else F(2, 3)
        fb = random.choice([-fa, F(1)])
        return fa, fb
    pos, neg = [], []
    for f in NICE + [-f for f in NICE]:
        R = _R(sysd, f)
        if R == 0:
            continue
        mag, ph, raw, _ = _polar(sysd, f)
        if ph == 1:
            continue
        (pos if float(R) > 0 else neg).append(f)
    if not pos or not neg:
        return None
    pp = [f for f in pos if f > 0]
    nn = [f for f in neg if f > 0]
    fa = random.choice(pp if pp and random.random() < 0.8 else pos)
    fb = random.choice(nn if nn and random.random() < 0.8 else neg)
    if abs(fa) == abs(fb):
        return None
    if random.random() < 0.5:          # which part, (a) or (b), has the negative real factor
        fa, fb = fb, fa
    return fa, fb


def generate(data):
    while True:
        sysd = _draw_family()
        ch = _choose(sysd)
        if ch is None:
            continue
        fa, fb = ch
        jumps, wraps = _breaks(sysd)
        if sysd["fam"] == "iir":
            cands_c = [f for f in (F(0), F(1, 2), F(2, 3), F(1, 3), F(1)) if _tau(sysd, f) != 0]
        else:
            cands_c = [f for f in NICE if _R(sysd, f) != 0 and f not in (fa, fb) and f not in wraps]
        if not cands_c:
            continue
        fc = random.choice(cands_c)
        break

    p = data["params"]
    c = data["correct_answers"]
    p["sys"] = {k: v for k, v in sysd.items() if k not in ("H_tex", "R_tex", "domain")}
    p["fam"] = sysd["fam"]
    p["H_tex"] = sysd["H_tex"]
    p["domain_tex"] = sysd["domain"]
    p["is_iir"] = sysd["fam"] == "iir"
    p["is_lin"] = sysd["fam"] in ("sin", "cos")
    p["is_even"] = sysd["fam"] == "even"

    rows = []
    for tag, f in (("a", fa), ("b", fb)):
        mag, ph, raw, R = _polar(sysd, f)
        c[f"mag_{tag}"] = float(mag)
        p[f"mag_{tag}_tex"] = tex_r(mag)
        c[f"ph_{tag}"] = float(ph)
        p[f"f{tag}"] = plain(f)
        p[f"f{tag}_tex"] = tex_pi(f)
        p[f"ph_{tag}_plain"] = plain(ph)
        p[f"nojump_{tag}"] = plain(principal(raw - 1)) if (R is not None and float(R) < 0) else ""
        if sysd["fam"] == "iir":
            a = F(sysd["a"])
            d_tex = fmt.tex_sum([(1, ""), (-a, tex_ej(-f))])
            dmag = sp.sqrt(sp.nsimplify(sp.expand((1 - rat(a) * cos_pi(f)) ** 2 + (rat(a) * sin_pi(f)) ** 2)))
            dang = F(cmath.phase(1 - float(a) * cmath.exp(-1j * math.pi * float(f))) / math.pi).limit_denominator(12)
            dpol = Polar(dmag, dang)
            g_tex = fmt.tex_num(F(sysd["g"]))
            out = Polar(mag, ph).tex()
            line = (rf"H_d\!\left({tex_pi(f)}\right) = \dfrac{{{g_tex}}}{{{d_tex}}} = \dfrac{{{g_tex}}}{{{dpol.tex()}}}"
                    rf" = {out}")
            rows.append({"line": line, "mag_tex": tex_r(mag), "mag_dec": dec(mag), "irr": not mag.is_rational,
                         "ph_e_tex": tex_ej(ph),
                         "ph_tex": tex_pi(ph), "w_tex": tex_pi(f), "neg": False, "R_val": ""})
            continue
        neg = float(R) < 0
        R_val = tex_real(R)
        if neg:
            p["neg_letter"] = tag
            p["nojump_neg"] = plain(principal(raw - 1))
        lin = _lin(sysd, f)
        if neg:
            work = f"{tex_pi(lin)} + \\pi = {tex_pi(raw)}" + (f" \\to {tex_pi(ph)}" if raw != ph else "")
        else:
            work = tex_pi(lin) + (f" \\to {tex_pi(ph)}" if raw != ph else "")
        rows.append({"angle_work": work, "w_tex": tex_pi(f), "R_val": R_val, "neg": neg, "mag_tex": tex_r(mag),
                     "mag_dec": dec(mag), "irr": not mag.is_rational, "lin_tex": tex_pi(_lin(sysd, f)),
                     "raw_tex": tex_pi(raw), "ph_tex": tex_pi(ph), "wrapped": raw != ph, "line": ""})
    p["rows"] = rows

    tau = _tau(sysd, fc)
    c["tau"] = float(tau)
    p["fc"] = plain(fc)
    p["fc_tex"] = tex_pi(fc)
    p["tau_tex"] = fmt.tex_num(tau)
    p["tau_neg"] = tau < 0
    if sysd["fam"] == "iir":
        a = F(sysd["a"])
        cc = {F(0): F(1), F(1, 3): F(1, 2), F(1, 2): F(0), F(2, 3): F(-1, 2), F(1): F(-1)}[fc]
        p["tau_work_tex"] = (rf"\tau\!\left({tex_pi(fc)}\right) = \dfrac{{{fmt.tex_num(a, paren=True)}\cdot{fmt.tex_num(cc, paren=True)} - "
                             rf"{fmt.tex_num(a, paren=True)}^2}}{{1 - 2\cdot{fmt.tex_num(a, paren=True)}\cdot{fmt.tex_num(cc, paren=True)} + "
                             rf"{fmt.tex_num(a, paren=True)}^2}} = {fmt.tex_num(tau)}")
        p["a_tex"] = fmt.tex_num(a)
    elif sysd["fam"] == "even":
        p["s_word"] = "+" if sysd["s"] == 1 else "-"
        p["phase_lin_tex"] = r"|\omega|" if sysd["s"] == 1 else r"-|\omega|"
    else:
        p["M"] = sysd["M"]

    # checkbox: jumps (all of them) + distractors (wraps first) -> 4 frequencies + "none"
    distract = [f for f in wraps if f not in jumps]
    others = [f for f in NICE if f not in jumps and f not in distract]
    random.shuffle(others)
    need = max(0, 4 - len(jumps))
    picked = (distract[:need] + others)[:need]
    listed = sorted(jumps + picked)
    p["jump_choices"] = [{"text": f"$\\omega = {tex_pi(f)}$", "correct": f in jumps} for f in listed]
    p["jump_choices"].append({"text": NONE_TEXT, "correct": not jumps})
    p["listed"] = [plain(f) for f in listed]
    p["jump_freqs"] = [plain(f) for f in jumps]
    p["wrap_freqs"] = [plain(f) for f in wraps]
    p["jumps_tex"] = ",\\ ".join(tex_pi(f) for f in jumps)
    p["wraps_tex"] = ",\\ ".join(tex_pi(f) for f in wraps)
    p["has_jumps"] = bool(jumps)
    p["has_wraps"] = bool(wraps)
    p["none_text"] = NONE_TEXT
    if sysd["fam"] in ("sin", "cos", "even"):
        p["R_tex"] = sysd["R_tex"]
