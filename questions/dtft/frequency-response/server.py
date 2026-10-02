"""Frequency response of a causal, stable LTI filter given by its LCCDE (Lecture 9 form):
|H_d(0)|, |H_d(pi)|, the filter type (lowpass / highpass / bandpass / bandstop / allpass) and the
output amplitude for x[n] = A cos(pi n / 2 + theta) applied for all n (eigenfunction property).
Lecture 14 (H_d(w) = H(z) on the unit circle), HW5; course summary "sinusoidal response of an LSI system".

Families (all with rational |H_d| at w = 0, pi/2, pi and a magnitude response that is monotone or
unimodal on [0, pi], so the type is unambiguous; the checker verifies the shape with freqz):
  fir2  y = p x[n] + q x[n-1], (|p|, |q|) a Pythagorean pair     -> LP (pq > 0) or HP (pq < 0)
  fir3  y = c x[n] + d x[n-1] + c x[n-2], |d| >= 2|c|             -> LP (cd > 0) or HP (cd < 0)
        y = c (x[n] + x[n-2])  -> bandstop;  y = c (x[n] - x[n-2])  -> bandpass
  iir1  y[n] - a y[n-1] = g x[n], a = +-3/4 (|1 + ja| = 5/4)       -> LP (a > 0) or HP (a < 0)
  iir2  y[n] + r y[n-2] = g (x[n] -+ x[n-2]), 0 < r < 1           -> bandpass / bandstop
  ap    y[n] - a y[n-1] = g (-a x[n] + x[n-1])  or  y[n] = g x[n-k]  -> allpass
"""

import random
from fractions import Fraction as F

from ece310 import fmt

TYPES = ["lowpass", "highpass", "bandpass", "bandstop", "allpass"]
THETAS = [("", "0"), (r" + \tfrac{\pi}{4}", "pi/4"), (r" - \tfrac{\pi}{3}", "-pi/3"), (r" + \tfrac{\pi}{6}", "pi/6")]


# ----------------------------------------------------------------------------- exact helpers
def _sqrt_frac(q):
    q = F(q)
    if q < 0:
        return None
    rn, rd = round(q.numerator ** 0.5), round(q.denominator ** 0.5)
    return F(rn, rd) if rn * rn == q.numerator and rd * rd == q.denominator else None


def _eval(c, w0):
    """sum_k c[k] e^{-j w0 k} as an exact Gaussian rational (re, im), w0 in {'0', 'pi/2', 'pi'}."""
    unit = {"0": lambda k: (1, 0), "pi": lambda k: ((-1) ** k, 0),
            "pi/2": lambda k: [(1, 0), (0, -1), (-1, 0), (0, 1)][k % 4]}[w0]
    re = sum((F(v) * unit(k)[0] for k, v in enumerate(c)), F(0))
    im = sum((F(v) * unit(k)[1] for k, v in enumerate(c)), F(0))
    return re, im


def _mag(b, a, w0):
    """|H_d(w0)| = sqrt(|B|^2 / |A|^2) if that is rational, else None."""
    br, bi = _eval(b, w0)
    ar, ai = _eval(a, w0)
    den = ar * ar + ai * ai
    if den == 0:
        return None
    return _sqrt_frac((br * br + bi * bi) / den)


def _cplx_tex(re, im):
    if im == 0:
        return fmt.tex_num(re)
    jt = ("" if abs(im) == 1 else fmt.tex_num(abs(im))) + "j"
    if re == 0:
        return ("-" if im < 0 else "") + jt
    return f"{fmt.tex_num(re)} {'-' if im < 0 else '+'} {jt}"


def _ejw_tex(k):
    return "" if k == 0 else ("e^{-j\\omega}" if k == 1 else f"e^{{-j{k}\\omega}}")


def _poly_ejw_tex(c):
    return fmt.tex_sum((v, _ejw_tex(k)) for k, v in enumerate(c))


# ----------------------------------------------------------------------------- families
def _draw():
    fam = random.choices(["fir2", "fir3", "fir3s", "iir1", "iir2", "ap"], weights=[16, 16, 14, 20, 16, 18])[0]
    if fam == "fir2":
        p, q = random.choice([(3, 4), (4, 3), (6, 8), (8, 6), (5, 12), (12, 5),
                              (F(3, 5), F(4, 5)), (F(4, 5), F(3, 5)), (F(3, 2), 2), (2, F(3, 2))])
        s, t = random.choice([1, -1]), random.choice([1, -1])
        b, a = [F(s * p), F(s * t * q)], [F(1)]
        kind = "lowpass" if t > 0 else "highpass"
        why = "the two taps have the same sign, so they add at $\\omega = 0$ and partly cancel at $\\omega = \\pi$" if t > 0 \
            else "the two taps have opposite signs, so they partly cancel at $\\omega = 0$ and add at $\\omega = \\pi$"
    elif fam == "fir3":
        c = random.choice([-2, -1, 1, 2])
        d = random.choice([1, -1]) * random.randint(2 * abs(c), 2 * abs(c) + 3)
        b, a = [F(c), F(d), F(c)], [F(1)]
        kind = "lowpass" if c * d > 0 else "highpass"
        amp = fmt.tex_sum([(d, ""), (2 * c, r"\cos\omega")])
        why = (f"$H_d(\\omega) = e^{{-j\\omega}}\\left({amp}\\right)$ and the bracket "
               "never changes sign, so the magnitude changes monotonically with $\\cos\\omega$")
    elif fam == "fir3s":
        c = random.choice([-3, -2, -1, 1, 2, 3, F(-1, 2), F(1, 2)])
        sgn = random.choice([1, -1])
        b, a = [F(c), F(0), F(sgn * c)], [F(1)]
        kind = "bandstop" if sgn > 0 else "bandpass"
        why = (f"$|H_d(\\omega)| = {2 * abs(c)}\\,|\\cos\\omega|$: zero at $\\omega = \\tfrac{{\\pi}}{{2}}$, largest at $0$ and $\\pi$" if sgn > 0
               else f"$|H_d(\\omega)| = {2 * abs(c)}\\,|\\sin\\omega|$: zero at $\\omega = 0$ and $\\pi$, largest at $\\tfrac{{\\pi}}{{2}}$")
    elif fam == "iir1":
        av = random.choice([F(3, 4), F(-3, 4)])
        g = random.choice([F(1, 4), F(1, 2), F(3, 4), F(1), F(5, 4), F(3, 2), F(7, 4), F(2)]) * random.choice([1, -1])
        b, a = [g], [F(1), -av]
        kind = "lowpass" if av > 0 else "highpass"
        why = (f"the pole at $z = {fmt.tex_num(av)}$ lies on the positive real axis, closer to $z = 1$ ($\\omega = 0$) than to "
               f"$z = -1$ ($\\omega = \\pi$): the denominator $|1 - {fmt.tex_num(av)}\\,e^{{-j\\omega}}|$ grows steadily from "
               "$\\omega = 0$ to $\\omega = \\pi$, so the gain falls" if av > 0 else
               f"the pole at $z = {fmt.tex_num(av)}$ lies on the negative real axis, closer to $z = -1$ ($\\omega = \\pi$) than to "
               f"$z = 1$ ($\\omega = 0$): the denominator $|1 + {fmt.tex_num(-av)}\\,e^{{-j\\omega}}|$ shrinks steadily from "
               "$\\omega = 0$ to $\\omega = \\pi$, so the gain rises")
    elif fam == "iir2":
        r = random.choice([F(1, 4), F(1, 2), F(1, 9), F(4, 9), F(1, 3)])
        g = random.choice([F(1, 2), F(1), F(2), F(1, 4)]) * random.choice([1, -1])
        sgn = random.choice([1, -1])
        b, a = [g, F(0), sgn * g], [F(1), F(0), r]
        kind = "bandstop" if sgn > 0 else "bandpass"
        sr = _sqrt_frac(r)
        poles = (f"$z = \\pm {fmt.tex_num(sr)}\\,j$" if sr is not None else f"$z = \\pm j\\sqrt{{{fmt.tex_num(r)}}}$")
        why = (f"the zeros at $z = \\pm j$ null $\\omega = \\tfrac{{\\pi}}{{2}}$; the poles at {poles} sit at the same angles "
               "and only sharpen the notch" if sgn > 0 else
               f"the zeros at $z = \\pm 1$ null $\\omega = 0$ and $\\omega = \\pi$, and the poles at {poles} boost $\\omega = \\tfrac{{\\pi}}{{2}}$")
    else:
        if random.random() < 0.75:
            av = random.choice([F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(2, 3), F(-2, 3), F(3, 4), F(-3, 4),
                                F(1, 4), F(-1, 4)])
            g = random.choice([F(1), F(2), F(1, 2), F(3)]) * random.choice([1, -1])
            b, a = [-av * g, g], [F(1), -av]
            why = ("the numerator $-a + e^{-j\\omega}$ equals $e^{-j\\omega}$ times the conjugate of the denominator "
                   "$1 - a e^{-j\\omega}$, so their magnitudes are equal at every $\\omega$")
        else:
            k = random.randint(1, 3)
            g = random.choice([F(1), F(2), F(3), F(1, 2)]) * random.choice([1, -1])
            b, a = [F(0)] * k + [g], [F(1)]
            why = "a pure delay only adds the linear phase $e^{-jk\\omega}$; $|e^{-jk\\omega}| = 1$"
        kind = "allpass"
    return fam, b, a, kind, why


def _lccde_tex(b, a):
    lhs = fmt.tex_sum((v, "y[n]" if k == 0 else f"y[n-{k}]") for k, v in enumerate(a))
    rhs = fmt.tex_sum((v, "x[n]" if k == 0 else f"x[n-{k}]") for k, v in enumerate(b))
    return f"{lhs} = {rhs}"


def generate(data):
    while True:
        fam, b, a, kind, why = _draw()
        mags = {w0: _mag(b, a, w0) for w0 in ("0", "pi/2", "pi")}
        if any(m is None for m in mags.values()):
            continue
        # input amplitude A in 2..10 with an integer output amplitude B = A |H(pi/2)| (else redraw)
        good = [A for A in range(2, 11) if (A * mags["pi/2"]).denominator == 1]
        if good:
            break
    m0, m2, mpi = mags["0"], mags["pi/2"], mags["pi"]
    A = random.choice(good[:4])
    B = A * m2
    th_tex, th_name = random.choice(THETAS)

    p = data["params"]
    p["b"] = [fmt.plain_num(v) for v in b]
    p["a"] = [fmt.plain_num(v) for v in a]
    p["family"], p["kind"], p["A"], p["theta"] = fam, kind, A, th_name
    p["lccde_tex"] = _lccde_tex(b, a)
    p["is_fir"] = len(a) == 1
    if len(a) > 1:                       # factor the gain out of the numerator: g * B'(z) / A(z)
        g = next(v for v in reversed(b) if v != 0) if fam == "ap" else b[0]
        bn = [v / g for v in b]
        gt = "" if g == 1 else ("-" if g == -1 else fmt.tex_num(g) + r"\,")
        p["H_tex"] = gt + fmt.tex_frac(_poly_ejw_tex(bn), _poly_ejw_tex(a))
        p["Hz_tex"] = gt + fmt.tex_frac(fmt.tex_poly_zinv(bn), fmt.tex_poly_zinv(a))
    else:
        p["H_tex"], p["Hz_tex"] = _poly_ejw_tex(b), fmt.tex_poly_zinv(b)
    p["x_tex"] = rf"{A}\cos\!\left(\tfrac{{\pi}}{{2}}n{th_tex}\right)"
    p["theta_tex"] = th_tex if th_tex else ""
    p["why_html"] = why

    def at(w0, zt):
        br, bi = _eval(b, w0)
        ar, ai = _eval(a, w0)
        val = mags[w0]
        nb, na = _sqrt_frac(br * br + bi * bi), _sqrt_frac(ar * ar + ai * ai)
        num = _cplx_tex(br, bi)
        if len(a) == 1:
            return rf"H(z)\big|_{{z={zt}}} = {num}" + (rf", \quad |H_d| = {fmt.tex_num(val)}" if (bi != 0 or br < 0) else "")
        den = _cplx_tex(ar, ai)
        head = rf"H(z)\big|_{{z={zt}}} = \left({num}\right)\big/\left({den}\right), \quad "
        if nb is not None and na is not None:
            return head + rf"|H_d| = {fmt.tex_num(nb)}\big/{fmt.tex_num(na, paren=True)} = {fmt.tex_num(val)}"
        return head + (rf"|H_d| = \sqrt{{{fmt.tex_num(br * br + bi * bi)}\big/{fmt.tex_num(ar * ar + ai * ai, paren=True)}}}"
                       rf" = {fmt.tex_num(val)}")

    p["at0_tex"] = at("0", "1")
    p["at2_tex"] = at("pi/2", "j")
    p["atpi_tex"] = at("pi", "-1")
    p["m0_tex"], p["m2_tex"], p["mpi_tex"] = fmt.tex_num(m0), fmt.tex_num(m2), fmt.tex_num(mpi)
    p["B_tex"] = fmt.tex_num(B)
    p["type_choices"] = [{"text": t, "correct": t == kind} for t in TYPES]   # fixed natural order

    c = data["correct_answers"]
    c["H0"] = float(m0)
    c["Hpi"] = float(mpi)
    c["B"] = float(B)
