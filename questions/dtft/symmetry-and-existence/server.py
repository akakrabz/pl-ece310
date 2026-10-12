"""Symmetry, periodicity and existence of the DTFT (Lecture 13 §6.1-6.2, Lecture 14 §1 and §4; the
site's T/F bank §1; SP2023 MT2 #2-#3, FA2024 MT2 #4(a), SP2023 MT2 #6, Lecture 14 slides 6 and 12-13,
HW6 #4). Every variant: one numeric part and three multiple-choice items.

(a) X_d at one frequency, Re and Im (numbers):
      herm  x real, X_d = F(w) given on [0, pi]; w0 in [-pi, 0): X_d(w0) = F(-w0)^*      (SP2023 #2)
      per   x arbitrary, X_d = F(w) given on [-pi, pi]; w0 in [2pi, 3pi] or [-3pi, -2pi]:
            X_d(w0) = F(w0 -+ 2pi)                                                      (SP2023 #3)
      both  x real, F on [0, pi]; w0 in (pi, 2pi): X_d(w0) = X_d(w0 - 2pi) = F(2pi - w0)^*  (slide 13 (e))
(b) is h[n] real-valued for a closed-form H_d(w)? (H_d(-w) = H_d(w)^* for all w <=> h real)
(c), (d) which statement describes the DTFT of a given signal: absolutely summable and X_d = X(e^{jw});
    not summable but a DTFT with Dirac impulses (u[n], cosines, e^{jw0 n}); not summable but an ordinary
    finite-energy function (sin(wc n)/(pi n): Lecture 14 Table 1; study-site Lecture 13 page §6.2, §8); no DTFT at all (growing)."""

import cmath
import math
import random
from fractions import Fraction as F

from ece310 import fmt, mc
from ece310.w4_dtft import tex_pi, tex_pi_n

T_VALUES = [F(1, 6), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(5, 6)]


def ej(angle):
    """e^{j angle pi} in LaTeX (angle a Fraction; '1' for 0)."""
    if angle == 0:
        return "1"
    return "e^{" + ("-" if angle < 0 else "") + "j" + tex_pi(abs(angle)) + "}"


def wrap(t):
    """tex of t*pi, parenthesised when negative."""
    s = tex_pi(t)
    return rf"\left({s}\right)" if t < 0 else s


# id: (LaTeX of F(w), LaTeX of F at t*pi, numerical F)
FBANK = {
    "jw": (r"j\omega", lambda t: rf"j{wrap(t)}", lambda w: 1j * w),
    "1pw_e1": (r"(1 + \omega)\,e^{-j\omega}", lambda t: rf"\left(1 + {wrap(t)}\right){ej(-t)}",
               lambda w: (1 + w) * cmath.exp(-1j * w)),
    "w2_e2": (r"\omega^2 e^{-j2\omega}", lambda t: rf"\left({tex_pi(t)}\right)^2 {ej(-2 * t)}",
              lambda w: w * w * cmath.exp(-2j * w)),
    "pimw_e1": (r"(\pi - \omega)\,e^{j\omega}", lambda t: rf"\left(\pi - {wrap(t)}\right){ej(t)}",
                lambda w: (math.pi - w) * cmath.exp(1j * w)),
    "w_e3": (r"\omega\,e^{-j3\omega}", lambda t: rf"{wrap(t)}\,{ej(-3 * t)}", lambda w: w * cmath.exp(-3j * w)),
    "w_jw2": (r"\omega + j\omega^2", lambda t: rf"{tex_pi(t)} + j\left({tex_pi(t)}\right)^2", lambda w: w + 1j * w * w),
    "j1pw": (r"j(1 + \omega)", lambda t: rf"j\left(1 + {wrap(t)}\right)", lambda w: 1j * (1 + w)),
    "jw_pimw": (r"j\omega(\pi - \omega)", lambda t: rf"j{wrap(t)}\left(\pi - {wrap(t)}\right)",
                lambda w: 1j * w * (math.pi - w)),
    "1pw2_e1": (r"(1 + \omega^2)\,e^{-j\omega}", lambda t: rf"\left(1 + \left({tex_pi(t)}\right)^2\right){ej(-t)}",
                lambda w: (1 + w * w) * cmath.exp(-1j * w)),
    "w2_ep1": (r"\omega^2 e^{j\omega}", lambda t: rf"\left({tex_pi(t)}\right)^2 {ej(t)}",
               lambda w: w * w * cmath.exp(1j * w)),
}
# A real x needs X_d(0) and X_d(pi) real (X_d(0) = sum x[n], X_d(pi) = sum (-1)^n x[n]): only these
# formulas are offered with a real x; jw, w_jw2 and j1pw appear only with an arbitrary x (SP2023 #3).
REAL_OK = ["1pw_e1", "w2_e2", "pimw_e1", "w_e3", "jw_pimw", "1pw2_e1", "w2_ep1"]

# id: (LaTeX of H_d(w), h real?)
HBANK = {
    "absw_sin": (r"|\omega|\,e^{-j\pi\sin\omega}", True),
    "w_cos": (r"\omega\,e^{j\pi\cos\omega}", False),
    "w2_cos": (r"\omega^2 e^{j\cos\omega}", False),
    "jw_e2": (r"j\omega\,e^{-j2\omega}", True),
    "w_e1": (r"\omega\,e^{-j\omega}", False),
    "cos_jsin2": (r"\cos\omega + j\sin 2\omega", True),
    "cos_jcos2": (r"\cos\omega + j\cos 2\omega", False),
    "e_jw2": (r"e^{-j\omega^2}", False),
    "absw_e3": (r"|\omega|\,e^{-j3\omega}", True),
    "jsin_e1": (r"j\sin(\omega)\,e^{-j\omega}", True),
    "sin_e1": (r"\sin(\omega)\,e^{-j\omega}", False),
    "w2_jsin": (r"\omega^2 e^{-j\sin\omega}", True),
    "pole_real": (r"\dfrac{1}{1 - \frac{1}{2}e^{-j\omega}}", True),
    "pole_j": (r"\dfrac{1}{1 - \frac{j}{2}e^{-j\omega}}", False),
    "e_jabs": (r"e^{j|\omega|}", False),
    "pimabs_e1": (r"(\pi - |\omega|)\,e^{j\omega}", True),
}
H_WHY = {
    True: r"$H_d(-\omega) = H_d^*(\omega)$ for every $\omega$ (even magnitude, odd phase), so $h[n]$ is real.",
    False: r"$H_d(-\omega) \ne H_d^*(\omega)$, so $h[n]$ is complex.",
}

EXIST_OPTS = [
    r"It is absolutely summable, and its DTFT is $X(z)$ evaluated at $z = e^{j\omega}$",
    r"It is not absolutely summable, but it has a DTFT that contains Dirac impulses $\delta(\omega - \omega_k)$",
    r"It is not absolutely summable, but it has a DTFT that is an ordinary function (no impulses)",
    r"It has no DTFT at all",
]


def _signal():
    """(id, params, LaTeX, category 0..3, reason)."""
    kind = random.choices(["geo", "nageo", "absgeo", "leftgrow_ok", "finite", "geo_big", "ngrow", "left_bad",
                           "step", "alt_step", "cos", "cexp", "cexp_step", "sinc"],
                          weights=[8, 6, 5, 6, 4, 8, 4, 6, 7, 5, 7, 5, 5, 14])[0]
    if kind == "geo":
        a = random.choice([F(1, 2), F(-1, 3), F(3, 4), F(-4, 5), F(2, 3)])
        return kind, [str(a)], rf"{fmt.tex_pow(a)}u[n]", 0, rf"$\sum_n |x[n]| = \frac{{1}}{{1 - {fmt.tex_num(abs(a))}}}$; ROC $|z| \gt {fmt.tex_num(abs(a))}$ contains the unit circle."
    if kind == "nageo":
        a = random.choice([F(1, 2), F(-1, 2), F(2, 3), F(-3, 4)])
        return kind, [str(a)], rf"n\,{fmt.tex_pow(a)}u[n]", 0, rf"$|a| = {fmt.tex_num(abs(a))} \lt 1$ beats the factor $n$; ROC $|z| \gt {fmt.tex_num(abs(a))}$ contains $|z| = 1$."
    if kind == "absgeo":
        a = random.choice([F(1, 2), F(1, 3), F(-1, 2), F(3, 4)])
        return kind, [str(a)], fmt.tex_pow(a, r"|n|"), 0, rf"Two decaying geometric tails; the ROC ${fmt.tex_num(abs(a))} \lt |z| \lt {fmt.tex_num(1 / abs(a))}$ contains the unit circle."
    if kind == "leftgrow_ok":
        b = random.choice([2, 3, -2])
        return kind, [str(b)], rf"{fmt.tex_pow(b)}u[-n]", 0, rf"For $n \le 0$, $|{b}^n| = {abs(b)}^{{-|n|}}$ decays; ROC $|z| \lt {abs(b)}$ contains the unit circle."
    if kind == "finite":
        L = random.choice([3, 4, 5, 8])
        return kind, [str(L)], rf"u[n] - u[n - {L}]", 0, "A finite-length sequence is always absolutely summable."
    if kind == "geo_big":
        a = random.choice([F(2), F(-3, 2), F(5, 4), F(3)])
        return kind, [str(a)], rf"{fmt.tex_pow(a)}u[n]", 3, rf"$|x[n]| = {fmt.tex_pow(abs(a))}$ grows; the ROC $|z| \gt {fmt.tex_num(abs(a))}$ misses the unit circle (FA2024 MT2 1(a))."
    if kind == "ngrow":
        a = random.choice([F(2), F(-3, 2), F(3)])
        return kind, [str(a)], rf"n\,{fmt.tex_pow(a)}u[n]", 3, rf"$|x[n]| = n\,{fmt.tex_pow(abs(a))}$ grows exponentially; the ROC $|z| \gt {fmt.tex_num(abs(a))}$ (double pole) misses the unit circle."
    if kind == "left_bad":
        a = random.choice([F(1, 2), F(1, 3), F(-1, 2)])
        return kind, [str(a)], rf"{fmt.tex_pow(a)}u[-n]", 3, rf"For $n \to -\infty$, $|x[n]| = {fmt.tex_num(1 / abs(a))}^{{|n|}}$ grows; the ROC $|z| \lt {fmt.tex_num(abs(a))}$ misses the unit circle."
    if kind == "step":
        return kind, [], r"u[n]", 1, r"Bounded but not summable (pole on $|z| = 1$): $\frac{1}{1 - e^{-j\omega}} + \pi\delta(\omega)$, which is not $X(e^{j\omega})$."
    if kind == "alt_step":
        return kind, [], r"(-1)^{n}u[n]", 1, r"Pole at $z = -1$ on the unit circle: $\frac{1}{1 + e^{-j\omega}} + \pi\delta(\omega - \pi)$ (FA2021 MT2 #2)."
    if kind == "cos":
        w = random.choice([F(1, 4), F(1, 3), F(2, 3), F(1, 5)])
        return kind, [str(w)], rf"\cos\left({tex_pi_n(w)}\right)\ \text{{for all }} n", 1, rf"Periodic, not summable, no z-transform: $\pi[\delta(\omega - {tex_pi(w)}) + \delta(\omega + {tex_pi(w)})]$ on $[-\pi, \pi]$."
    if kind == "cexp":
        w = random.choice([F(1, 4), F(1, 2), F(-1, 3)])
        sgn, dsg = ("-", "+") if w < 0 else ("", "-")
        return kind, [str(w)], rf"e^{{{sgn}j{tex_pi_n(abs(w))}}}\ \text{{for all }} n", 1, rf"$|x[n]| = 1$ forever: $2\pi\delta(\omega {dsg} {tex_pi(abs(w))})$ on $[-\pi, \pi]$; its z-transform does not exist (empty ROC)."
    if kind == "cexp_step":
        w = random.choice([F(1, 4), F(1, 2), F(1, 3)])
        return kind, [str(w)], rf"e^{{j{tex_pi_n(w)}}}u[n]", 1, rf"Pole $e^{{j{tex_pi(w)}}}$ on the unit circle: the DTFT needs the impulse $\pi\delta(\omega - {tex_pi(w)})$ (HW6 #4(b))."
    wc = random.choice([F(1, 4), F(1, 3), F(1, 2), F(2, 3)])
    return kind, [str(wc)], rf"\dfrac{{\sin\left({tex_pi_n(wc)}\right)}}{{\pi n}}", 2, rf"Decays only like $1/n$: $\sum|x[n]| = \infty$, but the energy is finite, and the DTFT is the rectangle $1$ on $|\omega| \le {tex_pi(wc)}$ (Lecture 14 Table 1 and §1: not absolutely summable, not stable; study-site Lecture 13 page §8)."


def _n_mistakes(fam, f, t):
    """How many classic wrong values (formula outside its range, conjugate without the sign flip,
    sign flip without the conjugate, ...) differ clearly from the true X_d(w0)."""
    w = float(t) * math.pi
    if fam == "per":
        true, wrong = f(w), [f(w + 2 * math.pi), f(-w).conjugate(), f(-w)]
    else:
        true = f(w).conjugate()
        wrong = [f(w), f(-w), f(-w).conjugate(), -f(w)]
    tol = lambda v: 2 * (1e-4 + 1e-3 * abs(v))  # noqa: E731
    return sum(1 for v in wrong if abs(v.real - true.real) > tol(true.real) or abs(v.imag - true.imag) > tol(true.imag))


def generate(data):
    p = data["params"]
    while True:
        fam = random.choice(["herm", "herm", "per", "both"])
        fid = random.choice(REAL_OK if fam in ("herm", "both") else list(FBANK))
        t = random.choice(T_VALUES) if fam != "per" else random.choice(T_VALUES + [-v for v in T_VALUES])
        if _n_mistakes(fam, FBANK[fid][2], t) >= 2:
            break
    F_tex, at_tex, fnum = FBANK[fid]
    if fam == "herm":
        w0 = -t
        val = fnum(float(t) * math.pi).conjugate()
        steps = rf"X_d\!\left({tex_pi(w0)}\right) = X_d^*\!\left({tex_pi(t)}\right) = \left[{at_tex(t)}\right]^*"
        why = r"$x[n]$ is real, so $X_d(-\omega) = X_d^*(\omega)$: conjugate <b>and</b> flip the sign of $\omega$."
        rng = r"0 \le \omega \le \pi"
    elif fam == "per":
        k = random.choice([1, 1, -1]) if t > 0 else random.choice([-1, -1, 1])
        w0 = t + 2 * k
        val = fnum(float(t) * math.pi)
        steps = rf"X_d\!\left({tex_pi(w0)}\right) = X_d\!\left({tex_pi(w0)} {'-' if k > 0 else '+'} 2\pi\right) = X_d\!\left({tex_pi(t)}\right) = {at_tex(t)}"
        why = r"Nothing is known about symmetry, but every DTFT is $2\pi$-periodic: shift the argument into $[-\pi, \pi]$."
        rng = r"-\pi \le \omega \le \pi"
    else:
        w0 = 2 - t
        val = fnum(float(t) * math.pi).conjugate()
        steps = (rf"X_d\!\left({tex_pi(w0)}\right) = X_d\!\left({tex_pi(w0)} - 2\pi\right) = X_d\!\left({tex_pi(-t)}\right)"
                 rf" = X_d^*\!\left({tex_pi(t)}\right) = \left[{at_tex(t)}\right]^*")
        why = r"Use both: periodicity brings $\omega_0$ into $[-\pi, 0]$, then $x$ real gives $X_d(-\omega) = X_d^*(\omega)$."
        rng = r"0 \le \omega \le \pi"
    p["family"] = fam
    p["fid"] = fid
    p["F_tex"] = F_tex
    p["range_tex"] = rng
    p["is_real"] = fam in ("herm", "both")
    p["w0"] = [w0.numerator, w0.denominator]
    p["t"] = [t.numerator, t.denominator]
    p["w0_tex"] = tex_pi(w0)
    p["steps_tex"] = steps
    p["why_html"] = why
    re_t, im_t = f"{val.real:.4f}", f"{abs(val.imag):.4f}"
    if abs(val.imag) < 5e-5:
        p["val_tex"] = re_t
    elif abs(val.real) < 5e-5:
        p["val_tex"] = ("-" if val.imag < 0 else "") + rf"j\,{im_t}"
    else:
        p["val_tex"] = f"{re_t} {'-' if val.imag < 0 else '+'} j\\,{im_t}"
    data["correct_answers"]["re"] = val.real if abs(val.real) > 1e-12 else 0.0
    data["correct_answers"]["im"] = val.imag if abs(val.imag) > 1e-12 else 0.0

    hid = random.choice(list(HBANK))
    p["hid"] = hid
    p["H_tex"] = HBANK[hid][0]
    p["h_choices"] = mc.yes_no(HBANK[hid][1])
    p["h_why_html"] = H_WHY[HBANK[hid][1]]

    while True:
        s1, s2 = _signal(), _signal()
        if s1[0] != s2[0]:
            break
    for name, s in (("s1", s1), ("s2", s2)):
        p[name] = {"kind": s[0], "args": s[1], "tex": s[2], "cat": s[3], "why_html": s[4],
                   "verdict": ["Absolutely summable, $X_d = X(e^{j\\omega})$", "A DTFT with Dirac impulses",
                               "An ordinary DTFT, but not absolutely summable", "No DTFT"][s[3]]}
        p[name + "_choices"] = [{"text": txt, "correct": i == s[3]} for i, txt in enumerate(EXIST_OPTS)]
