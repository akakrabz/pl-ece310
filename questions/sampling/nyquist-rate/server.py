"""Bandwidth B (Hz), Nyquist rate 2B, Nyquist period T_N = 1/(2B) and "which sampling rates satisfy
the Nyquist criterion f_s > 2B" (Lecture 17, eqs. 17-21) for signals built from

  sum       sums of cosines / sines (B = the highest frequency)
  product   cos(W1 t) sin(W2 t) [+ a third tone]   (product-to-sum: B = f1 + f2)
  sinc      sin(Wc t)/(pi t)                       (a box: B = Wc/2pi)
  sinc2     its square                             (a triangle twice as wide: B doubles)
  sincmod   sin(Wc t)/(pi t) cos(W0 t)             (shifted copies: B = (Wc + W0)/2pi)
  sincprod  product of two sincs                   (bandwidths add)
  sincconv  convolution of two sincs               (the smaller bandwidth)
  scale     x_c(c t) for x_c bandlimited to B0     (B = c B0; SP2021 MT2 #1(b) x_c(2t), SP2023 MT2 #1(b) x_c(3t))
  square    x_c(t)^2                               (B = 2 B0)
  modulate  x_c(t) cos(W0 t)                       (B = B0 + W0/2pi)
  absprod   x_c(t) y_c(t)                          (B = B1 + B2)
  absconv   x_c(t) * y_c(t)                        (B = min(B1, B2))
Frequencies are written in rad/s (cos(300 pi t), Omega_max = 1000 pi rad/s) or in Hz (cos(2 pi 150 t),
"bandlimited to 500 Hz"); the prompt always says which. Exactly-at-Nyquist sampling does not satisfy
the strict criterion f_s > 2B (Lecture 17; a sine at B sampled at 2B gives sin(pi n) = 0)."""

import random
from fractions import Fraction as F

from ece310 import fmt

FREQS = [10, 20, 25, 30, 40, 50, 60, 75, 80, 100, 120, 125, 150, 200, 250, 300, 400, 500, 600, 750, 800,
         1000, 1200, 1500, 2000, 2500, 3000, 4000, 5000]
MULTS_BELOW = [F(1, 4), F(1, 3), F(2, 5), F(1, 2), F(3, 5), F(3, 4)]
MULTS_ABOVE = [F(5, 4), F(3, 2), F(2), F(5, 2), F(3)]


def _w_tex(f, form):
    """Angular frequency 2 pi f as LaTeX: '300\\pi' (rad/s form) or '2\\pi\\cdot 150' (Hz form)."""
    f = F(f)
    if form == "hz":
        return rf"2\pi\cdot {fmt.plain_num(f)}"
    c = 2 * f
    return (r"\pi" if c == 1 else f"{fmt.plain_num(c)}\\pi")


def _tone(A, kind, f, form, phase=None):
    arg = rf"{_w_tex(f, form)}\,t"
    if phase:
        arg += (" + " if phase > 0 else " - ") + _pi_tex(abs(phase))
    amp = "" if A == 1 else str(A)
    return rf"{amp}\{kind}\!\left({arg}\right)"


def _pi_tex(q):
    q = F(q)
    num = r"\pi" if q.numerator == 1 else f"{q.numerator}\\pi"
    return num if q.denominator == 1 else r"\tfrac{" + num + "}{" + str(q.denominator) + "}"


def _sinc_tex(f, form):
    return rf"\dfrac{{\sin\!\left({_w_tex(f, form)}\,t\right)}}{{\pi t}}"


def _band_text(f, unit):
    """'bandlimited to 500 Hz' or 'Omega_max = 1000 pi rad/s' (exactly: X is nonzero up to that edge)."""
    if unit == "hz":
        return rf"$B_0 = {fmt.plain_num(f)}$ Hz"
    return rf"$\Omega_{{\max}} = {_w_tex(f, 'rad')}$ rad/s"


def _dec(q):
    q = F(q)
    return fmt.plain_num(q) if q.denominator == 1 else f"{float(q):g}"


def _family():
    fam = random.choice(["sum", "product", "sinc", "sinc2", "sincmod", "sincprod", "sincconv",
                         "scale", "square", "modulate", "absprod", "absconv"])
    form = random.choice(["hz", "rad"])
    d = {"fam": fam, "form": form}
    if fam == "sum":
        k = random.choice([2, 3])
        fs = random.sample(FREQS[:-4], k)
        terms = []
        for f in fs:
            kind = random.choice(["cos", "sin"])
            ph = random.choice([None, None, F(1, 4), F(-1, 3), F(1, 6), F(1, 2)])
            terms.append(_tone(random.randint(1, 4), kind, f, form, ph))
        d["f"] = [str(f) for f in fs]
        B = F(max(fs))
        sig = " + ".join(terms)
        expl = rf"a sum of sinusoids at ${', '.join(fmt.plain_num(f) for f in fs)}$ Hz: the highest one sets the band edge"
        return d, rf"x_a(t) = {sig}", None, B, expl
    if fam == "product":
        f1, f2 = random.sample(FREQS[:-6], 2)
        sig = rf"\cos\!\left({_w_tex(f1, form)}\,t\right)\sin\!\left({_w_tex(f2, form)}\,t\right)"
        d["f"] = [str(f1), str(f2)]
        f3 = None
        if random.random() < 0.6:
            cands = [f for f in FREQS if max(f1, f2) < f < f1 + f2]
            if cands:
                f3 = random.choice(cands)
                sig += " + " + _tone(random.randint(1, 3), "cos", f3, form, random.choice([None, F(1, 3)]))
                d["f3"] = str(f3)
        B = F(f1 + f2)
        expl = (rf"$\cos a\,\sin b = \tfrac12\left[\sin(a+b) - \sin(a-b)\right]$: the product contains the tones "
                rf"${fmt.plain_num(f1 + f2)}$ Hz and ${fmt.plain_num(abs(f1 - f2))}$ Hz"
                + (rf", and the extra tone is at ${f3}$ Hz" if f3 else "")
                + r"; the highest frequency is the <b>sum</b> $f_1 + f_2$")
        return d, rf"x_a(t) = {sig}", None, B, expl
    if fam in ("sinc", "sinc2", "sincmod"):
        fc = random.choice(FREQS[:-3])
        d["f"] = [str(fc)]
        s = _sinc_tex(fc, form)
        if fam == "sinc":
            return d, rf"x_a(t) = {s}", None, F(fc), (
                rf"$\dfrac{{\sin(\Omega_c t)}}{{\pi t}} \leftrightarrow$ a box of height 1 on $\lvert\Omega\rvert \lt \Omega_c$, "
                rf"so the band edge is $\Omega_c = 2\pi\cdot{fc}$ rad/s")
        if fam == "sinc2":
            return d, rf"x_a(t) = \left({s}\right)^{{2}}", None, F(2 * fc), (
                rf"the square is a product, so its CTFT is $\tfrac{{1}}{{2\pi}}$ times the convolution of the box on "
                rf"$\lvert\Omega\rvert \lt 2\pi\cdot{fc}$ with itself: a triangle reaching twice as far, "
                rf"$2\cdot 2\pi\cdot{fc}$ rad/s")
        f0 = random.choice([f for f in FREQS if f > fc] or [2 * fc])
        d["f0"] = str(f0)
        return d, rf"x_a(t) = {s}\,\cos\!\left({_w_tex(f0, form)}\,t\right)", None, F(fc + f0), (
            rf"multiplying by $\cos(\Omega_0 t)$ moves the box to $\pm\Omega_0$ (modulation): it now reaches "
            rf"$\Omega_0 + \Omega_c = 2\pi({f0} + {fc})$ rad/s")
    if fam in ("sincprod", "sincconv"):
        f1, f2 = random.sample(FREQS[:-4], 2)
        d["f"] = [str(f1), str(f2)]
        s1, s2 = _sinc_tex(f1, form), _sinc_tex(f2, form)
        if fam == "sincprod":
            return d, rf"x_a(t) = {s1}\cdot{s2}", None, F(f1 + f2), (
                r"a product in time is a convolution of the spectra, and the convolution of two boxes reaches as far as the "
                rf"<b>sum</b> of their edges: $2\pi({f1} + {f2})$ rad/s")
        d["post"] = r"where $*$ denotes convolution"
        return d, rf"x_a(t) = {s1} * {s2}", None, F(min(f1, f2)), (
            r"a convolution in time is a product of the spectra, which is nonzero only where <b>both</b> boxes are: "
            rf"up to the <b>smaller</b> edge, $2\pi\cdot{min(f1, f2)}$ rad/s")
    # abstract families: x_c(t) with a given band edge
    B0 = random.choice(FREQS[:-3])
    unit = random.choice(["hz", "rad"])
    d["f"] = [str(B0)]
    d["unit"] = unit
    given = rf"$x_c(t)$ has highest frequency {_band_text(B0, unit)} ($X_c(\Omega) \ne 0$ up to that edge, $0$ beyond)"
    if fam == "scale":
        c = random.choice([F(2), F(3), F(4), F(1, 2), F(3, 2)])
        d["c"] = str(c)
        arg = fmt.plain_num(c) + "t" if c.denominator == 1 else rf"\tfrac{{{c.numerator}}}{{{c.denominator}}}t"
        return d, rf"y(t) = x_c({arg})", given, c * B0, (
            rf"$x_c(at) \leftrightarrow \tfrac{{1}}{{\lvert a\rvert}}X_c\!\left(\tfrac{{\Omega}}{{a}}\right)$: scaling time by "
            rf"$a = {fmt.tex_num(c)}$ scales every frequency, and the band edge, by $a$: multiply by ${fmt.tex_num(c)}$")
    if fam == "square":
        return d, r"y(t) = x_c^2(t)", given, F(2 * B0), (
            r"$x_c^2$ is a product, so its spectrum is a convolution of $X_c$ with itself, which reaches twice as far")
    if fam == "modulate":
        f0 = random.choice(FREQS[:-3])
        d["f0"] = str(f0)
        return d, rf"y(t) = x_c(t)\cos\!\left({_w_tex(f0, form)}\,t\right)", given, F(B0 + f0), (
            r"modulation by $\cos(\Omega_0 t)$ shifts $X_c$ to $\pm\Omega_0$: the band edge moves to $\Omega_0 + \Omega_{\max}$")
    B1 = random.choice([f for f in FREQS[:-3] if f != B0])
    unit1 = "rad" if unit == "hz" else "hz"
    d["f1"] = str(B1)
    given += rf"; $v_c(t)$ has highest frequency {_band_text(B1, unit1)}"
    if fam == "absprod":
        return d, r"y(t) = x_c(t)\,v_c(t)", given, F(B0 + B1), (
            r"a product in time is a convolution of the spectra: the band edges <b>add</b>")
    d["post"] = r"where $*$ denotes convolution"
    return d, r"y(t) = x_c(t) * v_c(t)", given, F(min(B0, B1)), (
        r"a convolution in time is a product of the spectra, nonzero only where both are: the <b>smaller</b> band edge")


def generate(data):
    d, sig_tex, given_html, B, expl = _family()
    p = data["params"]
    p["data"] = d
    p["sig_tex"] = sig_tex
    p["has_given"] = given_html is not None
    p["post_html"] = d.pop("post", "")
    p["given_html"] = given_html or ""
    p["name_tex"] = "x_a(t)" if sig_tex.startswith("x_a") else "y(t)"
    p["expl_html"] = expl[0].upper() + expl[1:] if expl[0].isalpha() else expl
    p["B_str"] = str(B)
    p["B_tex"] = _dec(B)
    p["WB_tex"] = rf"2\pi\cdot{_dec(B)} = {_dec(2 * B)}\pi"
    p["fN_tex"] = fmt.plain_num(2 * B)
    p["TN_tex"] = rf"\tfrac{{1}}{{{fmt.plain_num(2 * B)}}}"
    TN = 1 / (2 * B)
    p["TN_dec"] = f"{float(TN):.6g}"
    # checkbox: exactly 2B (the strict-inequality trap) and four other multiples of 2B, at least one
    # below and one above; only rates with at most one decimal place
    def ok(m):
        return (2 * B * m * 10).denominator == 1
    below = [m for m in MULTS_BELOW if ok(m)]
    above = [m for m in MULTS_ABOVE if ok(m)]
    while True:
        nb = random.randint(1, min(3, len(below)))
        ms = random.sample(below, nb) + random.sample(above, 4 - nb)
        if len(ms) == 4:
            break
    rates = sorted([F(1)] + ms)
    opts = []
    for m in rates:
        f = 2 * B * m
        txt = fmt.plain_num(f) if f.denominator == 1 else f"{float(f):g}"
        opts.append({"text": rf"$f_s = {txt}$ Hz", "correct": bool(m > 1), "m": str(m)})
    p["rate_choices"] = opts
    p["rate_lines"] = [{"line_html": o["text"] + (" &mdash; satisfies $f_s \\gt 2B$" if o["correct"] else
                                              (" &mdash; exactly the Nyquist rate: does <b>not</b> satisfy the strict criterion "
                                               "$f_s \\gt 2B$ (the copies only touch; a tone exactly at $B$ can be lost)"
                                               if F(o["m"]) == 1 else " &mdash; below the Nyquist rate"))}
                       for o in opts]
    c = data["correct_answers"]
    c["B"] = float(B)
    c["fN"] = float(2 * B)
    c["TN"] = float(TN)
