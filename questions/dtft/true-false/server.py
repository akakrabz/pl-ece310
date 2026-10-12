"""Six True/False statements on the DTFT and the frequency response (exam family "True/False",
problem 1 on 7 of 7 past Midterm 2 exams; Lectures 13-16).

Fixed statements are the Unit 3 statements of the study site's Midterm 2 True/False bank (official
answers, Q1-Q16) plus a few facts stated in Lectures 13-16. Parameterized statements are randomized
and their answer is recomputed, e.g. FA2024 #1(a) "x[n] = 2^n u[n] has X(z) = 1/(1 - 2z^-1), so its
DTFT is 1/(1 - 2e^-jw)" becomes "x[n] = a^n u[n] ..." (True iff |a| < 1). One statement per topic
group, six groups out of nine; about half of the statements are True on average."""

import math
import random
from fractions import Fraction as Fr

from ece310 import fmt, mc


def T(x):
    return fmt.tex_num(x)


def pj(x):
    x = Fr(x)
    return [x.numerator, x.denominator]


def pi_tex(c):
    """c*pi for a rational c: 0, \\pi, -\\pi, 2\\pi, \\tfrac{\\pi}{3}, -\\tfrac{2\\pi}{3}."""
    c = Fr(c)
    if c == 0:
        return "0"
    s = "-" if c < 0 else ""
    p, q = abs(c).numerator, abs(c).denominator
    num = r"\pi" if p == 1 else rf"{p}\pi"
    return s + (num if q == 1 else rf"\tfrac{{{num}}}{{{q}}}")


def ejw(k):
    """e^{-j k w} for a rational k: '' (k = 0), e^{-j\\omega}, e^{j2\\omega}, e^{-j\\omega/2}, e^{-j3\\omega/2}."""
    k = Fr(k)
    if k == 0:
        return ""
    p, q = abs(k).numerator, abs(k).denominator
    body = (r"\omega" if p == 1 else rf"{p}\omega") + ("" if q == 1 else f"/{q}")
    return rf"e^{{{'-' if k > 0 else ''}j{body}}}"


def shift_tex(name, k):
    """x[n - k] in LaTeX."""
    return f"{name}[n]" if k == 0 else f"{name}[n {'-' if k > 0 else '+'} {abs(k)}]"


def xd_at(w):
    return rf"X_d\left({w}\right)"


# ----------------------------------------------------------------------------- fixed statements
# sid: (source, statement, answer, reason)
FIXED = {
    "ztrans_always": (
        "FA2023 #1(a)",
        r"A discrete-time signal $x[n]$ with z-transform $X(z)$ always has a finite DTFT $X_d(\omega)$ given by $X(z)$ "
        r"evaluated at $z = e^{j\omega}$, i.e. evaluate $X(z)$ along the unit circle.",
        False,
        r"Only when the ROC contains the unit circle: $x[n] = 2^n u[n]$ has $X(z) = \dfrac{1}{1 - 2z^{-1}}$, "
        r"$\lvert z\rvert \gt 2$, but $\sum_n\lvert x[n]\rvert = \infty$ and no DTFT."),
    "finite_periodic": (
        "FA2024 #1(b)",
        r"The DTFT $X_d(\omega)$ of a finite-length signal $x[n]$ is $2\pi$-periodic.",
        True,
        r"Every DTFT is $2\pi$-periodic: $e^{-j(\omega + 2\pi)n} = e^{-j\omega n}$ for every integer $n$. "
        r"Finite length only guarantees that the DTFT exists."),
    "periodic_only_infinite": (
        "FA2019 #1(e)",
        r"The DTFT $X_d(\omega)$ of a discrete-time signal $x[n]$ is periodic only when $x[n]$ has infinite length.",
        False,
        r"Periodicity comes from $n$ being an integer, not from the length of $x$: "
        r"$\delta[n]$ has $X_d(\omega) = 1$ and $\delta[n-1]$ has $e^{-j\omega}$, both $2\pi$-periodic."),
    "uniqueness": (
        "Lecture 13, DTFT existence and uniqueness",
        r"Two different absolutely summable sequences can have the same DTFT.",
        False,
        r"The DTFT is unique: the inverse DTFT $x[n] = \frac{1}{2\pi}\int_{-\pi}^{\pi}X_d(\omega)e^{j\omega n}\,d\omega$ "
        r"recovers every sample, so equal DTFTs mean equal sequences."),
    "real_hermitian": (
        "FA2023 #1(c)",
        r"The DTFT of a real-valued signal is always Hermitian symmetric.",
        True,
        r"For real $x[n]$, $X_d(-\omega) = \sum_n x[n]e^{j\omega n} = X_d^*(\omega)$: even magnitude, odd phase (Lecture 14)."),
    "real_X0_Xpi": (
        "Lecture 14, Hermitian symmetry",
        r"If $x[n]$ is real-valued and absolutely summable, then $X_d(0)$ and $X_d(\pi)$ are real numbers.",
        True,
        r"$X_d(0) = \sum_n x[n]$ and $X_d(\pi) = \sum_n (-1)^n x[n]$ are sums of real numbers "
        r"(equivalently $X_d(\pi) = X_d(-\pi) = X_d^*(\pi)$)."),
    "real_Xd_real": (
        "Lecture 14, Hermitian symmetry",
        r"If $x[n]$ is real-valued, then its DTFT $X_d(\omega)$ is real for every $\omega$.",
        False,
        r"A real signal has a Hermitian DTFT, not a real one: $x[n] = \delta[n-1]$ has $X_d(\omega) = e^{-j\omega}$. "
        r"$X_d$ is real when $x$ is real and even."),
    "Xd_real_x_real": (
        "Lecture 14, Hermitian symmetry",
        r"If $X_d(\omega)$ is real for every $\omega$, then $x[n]$ is real-valued.",
        False,
        r"$x[n] = j\delta[n+1] - j\delta[n-1]$ is purely imaginary, yet $X_d(\omega) = je^{j\omega} - je^{-j\omega} = -2\sin\omega$ "
        r"is real. A real $X_d$ means $x[-n] = x^*[n]$, not a real $x$."),
    "conv_product": (
        "Lecture 14, convolution property",
        r"If $y[n] = x[n] * h[n]$, where $x[n]$ and $h[n]$ are absolutely summable, then $Y_d(\omega) = X_d(\omega)H_d(\omega)$.",
        True,
        r"Convolution in time is multiplication in frequency: it is $Y(z) = X(z)H(z)$ on the unit circle, "
        r"which lies in both ROCs."),
    "product_product": (
        "Lecture 14, windowing property",
        r"If $y[n] = x[n]\,h[n]$, where $x[n]$ and $h[n]$ are absolutely summable, then $Y_d(\omega) = X_d(\omega)H_d(\omega)$.",
        False,
        r"A product in time is a periodic convolution in frequency, $Y_d(\omega) = \frac{1}{2\pi}\int_{-\pi}^{\pi}"
        r"X_d(\theta)H_d(\omega - \theta)\,d\theta$. Counterexample: $x = h = \delta[n] + \delta[n-1]$ gives $y = x$, but "
        r"$X_dH_d = (1 + e^{-j\omega})^2 \ne 1 + e^{-j\omega}$."),
    "zp_can_change": (
        "FA2019 #1(a)",
        r"Zero padding a finite-length signal can change its DTFT.",
        False,
        r"The appended zeros add only terms $0\cdot e^{-j\omega n}$ to $\sum_n x[n]e^{-j\omega n}$. Zero padding changes the "
        r"DFT (a finer grid of samples of the same $X_d$), never the DTFT."),
    "zp_equal": (
        "FA2023 #1(f), SP2023 #1(c)",
        r"Let $\{x[n]\}_{n=0}^{N-1}$ be a length-$N$ signal with DTFT $X_d(\omega)$ and $\{y[n]\}_{n=0}^{M+N-1}$ be $x[n]$ "
        r"zero-padded with $M$ zeros. The DTFT of $y[n]$, $Y_d(\omega)$, is equal to $X_d(\omega)$.",
        True,
        r"$Y_d(\omega) = \sum_{n=0}^{N-1}x[n]e^{-j\omega n} + \sum_{n=N}^{N+M-1}0\cdot e^{-j\omega n} = X_d(\omega)$."),
    "cos_eigen_every": (
        "FA2021 #1(a), reworded with “every”",
        r"For a fixed frequency $0 \lt \omega_0 \lt \pi$, $x[n] = \cos(\omega_0 n)$ is an eigenfunction of every stable LTI system.",
        False,
        r"Only $e^{j\omega_0 n}$ is an eigenfunction of every LTI system. The stable one-sample delay $y[n] = x[n-1]$ maps "
        r"$\cos(\omega_0 n)$ to $\cos(\omega_0 n - \omega_0)$, which is not a multiple of $\cos(\omega_0 n)$ when $0 \lt \omega_0 \lt \pi$."),
    "exp_eigen": (
        "Lecture 15, eigenfunctions",
        r"For every stable LTI system with frequency response $H_d(\omega)$, the input $x[n] = e^{j\omega_0 n}$ (for all $n$) "
        r"produces the output $y[n] = H_d(\omega_0)\,e^{j\omega_0 n}$.",
        True,
        r"$y[n] = \sum_k h[k]e^{j\omega_0(n-k)} = e^{j\omega_0 n}\sum_k h[k]e^{-j\omega_0 k} = H_d(\omega_0)e^{j\omega_0 n}$."),
    "real_cos_formula": (
        "FA2019 #1(d)",
        r"An LSI system has transfer function $H(z)$. Suppose that the output of the system to the input "
        r"$x[n] = \cos(\omega_0 n + \phi)$ is $y[n] = \lvert H(e^{j\omega_0})\rvert\cos\left(\omega_0 n + \phi + \angle H(e^{j\omega_0})\right)$ "
        r"for any $\omega_0$, $\phi$, where $H(e^{j\omega}) = H(z)\big|_{z=e^{j\omega}}$. Then the system $h[n]$ can be either real or complex.",
        False,
        r"The output is $\tfrac12 H(e^{j\omega_0})e^{j(\omega_0 n + \phi)} + \tfrac12 H(e^{-j\omega_0})e^{-j(\omega_0 n + \phi)}$; it equals the "
        r"stated formula for every $\omega_0$, $\phi$ only if $H(e^{-j\omega}) = H^*(e^{j\omega})$, i.e. only if $h[n]$ is real."),
}


# ----------------------------------------------------------------------------- G1: existence, ROC
def roc_substitute():
    side = random.choice(["right", "right", "left"])
    if side == "right":
        a = random.choice([Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(3, 4), Fr(2), Fr(3), Fr(-2), Fr(3, 2), Fr(1), Fr(-1)])
    else:
        a = random.choice([Fr(2), Fr(3), Fr(-2), Fr(4, 3), Fr(1, 2), Fr(1, 3), Fr(-1, 2)])
    Xz = rf"\dfrac{{1}}{{{fmt.tex_sum([(1, ''), (-a, 'z^{-1}')])}}}"
    Xw = rf"\dfrac{{1}}{{{fmt.tex_sum([(1, ''), (-a, r'e^{-j\omega}')])}}}"
    if side == "right":
        xt = "u[n]" if a == 1 else rf"{fmt.tex_pow(a, 'n')}u[n]"
        roc = rf"\lvert z\rvert \gt {T(abs(a))}"
        ans = abs(a) < 1
    else:
        base = fmt.tex_pow(a, 'n')
        xt = (rf"-\left({base}\right)u[-n-1]" if (a > 0 and a.denominator == 1) else rf"-{base}u[-n-1]")
        roc = rf"\lvert z\rvert \lt {T(abs(a))}"
        ans = abs(a) > 1
    text = (rf"The z-transform of $x[n] = {xt}$ is $X(z) = {Xz}$, ${roc}$. "
            rf"Therefore, the DTFT of $x[n]$ is $X_d(\omega) = {Xw}$.")
    sub = r"$X_d(\omega) = X(z)\big|_{z=e^{j\omega}}$"
    if ans:
        why = (rf"The ROC ${roc}$ contains the unit circle, so $x[n]$ is absolutely summable and {sub} (Lecture 14).")
    elif abs(a) == 1:
        pair = (r"$u[n] \leftrightarrow \dfrac{1}{1 - e^{-j\omega}} + \pi\delta(\omega)$" if a == 1 else
                r"$(-1)^n u[n] \leftrightarrow \dfrac{1}{1 + e^{-j\omega}} + \pi\delta(\omega - \pi)$ on $-\pi \lt \omega \le \pi$")
        why = (rf"The pole $z = {T(a)}$ lies on the unit circle, so the ROC ${roc}$ does not contain $\lvert z\rvert = 1$ and "
               rf"$\sum_n\lvert x[n]\rvert = \infty$. The DTFT exists only with an impulse, {pair} (Lecture 14), "
               "which the formula misses.")
    elif side == "right":
        why = (rf"The ROC ${roc}$ misses the unit circle: $\lvert x[n]\rvert = {fmt.tex_pow(abs(a), 'n')}$ grows, "
               r"$\sum_n\lvert x[n]\rvert = \infty$, and $x[n]$ has no DTFT. Substituting $z = e^{j\omega}$ outside the ROC "
               "gives numbers that mean nothing.")
    else:
        why = (rf"The ROC ${roc}$ lies inside the unit circle: for $n \le -1$, $\lvert x[n]\rvert = "
               rf"{fmt.tex_pow(1 / abs(a), r'\lvert n\rvert')}$ grows as $n \to -\infty$, so $\sum_n\lvert x[n]\rvert = \infty$ "
               "and there is no DTFT.")
    return "FA2024 #1(a), FA2023 #1(a)", text, ans, why, {"a": pj(a), "side": side}


def abs_summable():
    form = random.choice(["right", "anti", "two"])
    a = random.choice([Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(3, 4), Fr(2), Fr(-2), Fr(3), Fr(4, 3)])
    r = abs(a)
    if form == "right":
        xt, ans = rf"{fmt.tex_pow(a, 'n')}u[n]", r < 1
    elif form == "anti":
        xt, ans = rf"{fmt.tex_pow(a, 'n')}u[-n]", r > 1
    else:
        xt, ans = fmt.tex_pow(a, r"\lvert n\rvert"), r < 1
    text = (rf"The DTFT sum $\sum_{{n=-\infty}}^{{\infty}} x[n]\,e^{{-j\omega n}}$ of $x[n] = {xt}$ converges absolutely "
            r"for every $\omega$.")
    if ans:
        q = r if form != "anti" else 1 / r
        total = {"right": 1 / (1 - q), "anti": 1 / (1 - q), "two": (1 + q) / (1 - q)}[form]
        why = (rf"$\sum_n\lvert x[n]\rvert$ is a geometric series with ratio ${T(q)} \lt 1$: it equals ${T(total)} \lt \infty$, "
               r"and $\lvert X_d(\omega)\rvert \le \sum_n\lvert x[n]\rvert$ for every $\omega$ (Lecture 13).")
    else:
        grow = {"right": rf"$\lvert x[n]\rvert = {fmt.tex_pow(r, 'n')}$ grows as $n \to \infty$",
                "anti": rf"for $n \le 0$, $\lvert x[n]\rvert = {fmt.tex_pow(1 / r, r'\lvert n\rvert')}$ grows as $n \to -\infty$",
                "two": rf"$\lvert x[n]\rvert = {fmt.tex_pow(r, r'\lvert n\rvert')}$ grows in both directions"}[form]
        why = (rf"The samples do not decay: {grow}. So $\sum_n\lvert x[n]\rvert = \infty$; the terms of the DTFT sum do not "
               "even tend to zero.")
    return "Lecture 13, existence of the DTFT", text, ans, why, {"a": pj(a), "form": form}


# ----------------------------------------------------------------------------- G2: periodicity
def bandlimited_zero():
    wc = random.choice([Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(1, 6)])
    kind = random.choice(["beyond", "gap", "gap", "upper"])
    W = pi_tex(wc)
    if kind == "beyond":
        lo, hi, rng, ans = wc, None, rf"\lvert\omega\rvert \gt {W}", False
    elif kind == "gap":
        k = random.choice([0, 0, 1])
        lo, hi = 2 * k + wc, 2 * (k + 1) - wc
        rng, ans = rf"{pi_tex(lo)} \lt \omega \lt {pi_tex(hi)}", True
    else:
        lo, hi, rng, ans = Fr(1), Fr(2), r"\pi \lt \omega \lt 2\pi", False
    text = (rf"Suppose that $x[n]$ is bandlimited to ${W}$, that is, $X_d(\omega) = 0$ for ${W} \lt \lvert\omega\rvert \le \pi$. "
            rf"Then $X_d(\omega) = 0$ for ${rng}$.")
    lp = rf"$x[n] = \dfrac{{\sin\left({W}n\right)}}{{\pi n}}$"
    if kind == "beyond":
        why = (rf"Every DTFT is $2\pi$-periodic, so the band $\lvert\omega\rvert \lt {W}$ reappears around $\pm2\pi, \pm4\pi, \dots$: "
               rf"$X_d(2\pi) = X_d(0)$, which is $1$ for the ideal low-pass {lp}.")
    elif kind == "gap":
        why = (rf"The interval lies between the copy of the band centred at ${pi_tex(2 * k)}$ and the one centred at "
               rf"${pi_tex(2 * k + 2)}$: shifting by a multiple of $2\pi$ moves each of its points into "
               rf"${W} \lt \lvert\omega\rvert \le \pi$, where $X_d = 0$.")
    else:
        why = (rf"$X_d$ on $\pi \lt \omega \lt 2\pi$ repeats $X_d$ on $-\pi \lt \omega \lt 0$, which contains the band "
               rf"$-{W} \lt \omega \lt 0$: for the ideal low-pass {lp}, "
               rf"$X_d\left(2\pi - {pi_tex(wc / 2)}\right) = X_d\left(-{pi_tex(wc / 2)}\right) = 1$.")
    return "FA2019 #1(b)", text, ans, why, {"wc": pj(wc), "lo": pj(lo), "hi": None if hi is None else pj(hi)}


# ----------------------------------------------------------------------------- G3: impulses
def finite_cos_impulses():
    w0 = random.choice([Fr(1, 32), Fr(1, 16), Fr(1, 8), Fr(1, 4), Fr(1, 3)])
    N = random.choice([8, 16, 32])
    kind = random.choice(["finite", "all"])
    W = pi_tex(w0)
    if kind == "finite":
        sig = (rf"the finite-length signal $x[n] = \cos\left({W}n\right)$, $n = 0, 1, \dots, {N - 1}$ "
               r"(and $x[n] = 0$ otherwise),")
    else:
        sig = rf"the signal $x[n] = \cos\left({W}n\right)$, defined for all $n$,"
    text = (rf"The DTFT of {sig} is $X_d(\omega) = \pi\left[\delta\left(\omega - {W}\right) + \delta\left(\omega + {W}\right)\right]$ "
            r"for $-\pi \le \omega \lt \pi$.")
    if kind == "finite":
        why = (rf"A finite-length signal is absolutely summable, so its DTFT is an ordinary bounded function, "
               rf"$X_d(\omega) = \sum_{{n=0}}^{{{N - 1}}}\cos\left({W}n\right)e^{{-j\omega n}}$ (two shifted Dirichlet kernels) with "
               rf"$\lvert X_d(\omega)\rvert \le {N}$: no impulses. Impulses need the cosine to last forever.")
    else:
        why = (rf"$\cos\left({W}n\right) = \tfrac12 e^{{j{W}n}} + \tfrac12 e^{{-j{W}n}}$ and $e^{{j\omega_0 n}} \leftrightarrow "
               r"2\pi\delta(\omega - \omega_0)$ on one period (Lectures 13-14): two impulses of area $\pi$, repeated every $2\pi$.")
    return "FA2019 #1(c)", text, kind == "all", why, {"w0": pj(w0), "N": N, "kind": kind}


def coef_tex(c):
    """c*pi for a Gaussian-integer coefficient c = (re, im) with one part zero."""
    re, im = c
    if im == 0:
        return {1: r"\pi", -1: r"-\pi"}.get(re, rf"{re}\pi")
    return {1: r"j\pi", -1: r"-j\pi"}.get(im, rf"{im}j\pi")


def impulses_tex(imps):
    """Sum of impulses c_i*pi*delta(w - loc_i*pi); two impulses of equal size are factored as c[d(...) +- d(...)]."""
    def d(loc):
        if loc == 0:
            return r"\delta(\omega)"
        return rf"\delta\left(\omega {'-' if loc > 0 else '+'} {pi_tex(abs(loc))}\right)"
    if len(imps) == 1:
        (re, im), loc = imps[0][:2], imps[0][2]
        return coef_tex((re, im)) + d(loc)
    (r1, i1, l1), (r2, i2, l2) = imps
    sign = "+" if (r1, i1) == (r2, i2) else "-"
    return rf"{coef_tex((r1, i1))}\left[{d(l1)} {sign} {d(l2)}\right]"


def impulse_pair():
    w0 = random.choice([Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4)])
    sig = random.choice(["exp", "cos", "sin", "const"])
    good = random.random() < 0.5
    A = random.choice([2, 3, 5])
    W = pi_tex(w0)
    sig_tex = {"exp": rf"e^{{j{W}n}}", "cos": rf"\cos\left({W}n\right)", "sin": rf"\sin\left({W}n\right)", "const": str(A)}[sig]
    right = {"exp": [(2, 0, w0)], "cos": [(1, 0, w0), (1, 0, -w0)], "sin": [(0, -1, w0), (0, 1, -w0)],
             "const": [(2 * A, 0, 0)]}[sig]
    # (impulses, what they are really the DTFT of)
    wrong = {"exp": [([(1, 0, w0)], rf"\tfrac12 e^{{j{W}n}}"), ([(2, 0, -w0)], rf"e^{{-j{W}n}}")],
             "cos": [([(2, 0, w0), (2, 0, -w0)], rf"2\cos\left({W}n\right)"),
                     ([(1, 0, w0), (-1, 0, -w0)], rf"j\sin\left({W}n\right)")],
             "sin": [([(0, 1, w0), (0, -1, -w0)], rf"-\sin\left({W}n\right)"),
                     ([(1, 0, w0), (-1, 0, -w0)], rf"j\sin\left({W}n\right)")],
             "const": [([(A, 0, 0)], rf"\tfrac{{{A}}}{{2}}")]}[sig]
    if good:
        imps, really = right, None
    else:
        imps, really = random.choice(wrong)
    text = (rf"For $-\pi \le \omega \lt \pi$, the DTFT of $x[n] = {sig_tex}$ (for all $n$) is $X_d(\omega) = {impulses_tex(imps)}$.")
    pair = rf"$X_d(\omega) = {impulses_tex(right)}$"
    base = {"exp": r"$e^{j\omega_0 n} \leftrightarrow 2\pi\delta(\omega - \omega_0)$ (sifting: $\frac{1}{2\pi}\int 2\pi\delta(\omega - \omega_0)e^{j\omega n}d\omega = e^{j\omega_0 n}$)",
            "cos": r"Euler, $\cos\theta = \tfrac12(e^{j\theta} + e^{-j\theta})$, and $e^{j\omega_0 n} \leftrightarrow 2\pi\delta(\omega - \omega_0)$",
            "sin": r"Euler, $\sin\theta = \tfrac{1}{2j}(e^{j\theta} - e^{-j\theta})$, and $e^{j\omega_0 n} \leftrightarrow 2\pi\delta(\omega - \omega_0)$; $\tfrac{2\pi}{2j} = -j\pi$",
            "const": rf"${A} = {A}e^{{j0n}}$ and $e^{{j\omega_0 n}} \leftrightarrow 2\pi\delta(\omega - \omega_0)$ with $\omega_0 = 0$"}[sig]
    if good:
        why = rf"{base} (Lectures 13-14)."
    else:
        why = (rf"The pair is {pair}, from {base}. The stated impulses are the DTFT of ${really}$ (invert them by sifting).")
    prm = {"sig": sig, "w0": pj(w0), "A": A, "imps": [[re, im, pj(loc)] for re, im, loc in imps]}
    return "Lectures 13-14, impulse pairs", text, good, why, prm


# ----------------------------------------------------------------------------- G4: uniqueness, inverse
def inverse_finite():
    n1 = random.choice([0, 0, 0, -1, -2])
    L = random.choice([3, 4])
    n2 = n1 + L - 1
    coefs = [random.choice([1, 2, 3, 4, 5, -1, -2, -3]) for _ in range(L)]
    M = random.choice([n2, n2, n2 - 1]) if n1 == 0 else n2
    ans = n1 >= 0 and n2 <= M
    X = fmt.tex_sum([(c, ejw(k)) for k, c in zip(range(n1, n2 + 1), coefs)])
    text = (rf"If $x[n]$ is the inverse DTFT of $X_d(\omega) = {X}$, then $x[n]$ must be zero for $n \lt 0$ and $n \gt {M}$.")
    seq = fmt.tex_seq(coefs, n1)
    why = rf"The DTFT is unique and $\delta[n - n_0] \leftrightarrow e^{{-j\omega n_0}}$, so $x[n] = {seq}$ and nothing else"
    if ans:
        why += rf": zero outside $0 \le n \le {M}$."
    elif n1 < 0:
        why += rf". The term with $e^{{j{'' if n1 == -1 else -n1}\omega}}$ sits at $n = {n1}$: $x[{n1}] = {coefs[0]} \ne 0$."
    else:
        why += rf", and $x[{n2}] = {coefs[-1]} \ne 0$."
    return "SP2025 #1(a)", text, ans, why, {"coefs": coefs, "n1": n1, "M": M}


def fractional_delay():
    d = random.choice([Fr(1, 2), Fr(1, 2), Fr(3, 2), Fr(1, 3), Fr(5, 2), Fr(-1, 2), Fr(1), Fr(2), Fr(3), Fr(-1), Fr(4), Fr(-2)])
    ans = d.denominator == 1
    sg = "-" if d > 0 else "+"
    dt = rf"\delta\left[n {sg} {T(abs(d))}\right]"
    text = (rf"Let $x[n]$ be a discrete-time signal with DTFT $X_d(\omega) = {ejw(d)}$, $-\pi \le \omega \le \pi$. "
            rf"The corresponding signal is $x[n] = {dt}$.")
    if ans:
        why = (rf"For an integer delay, ${ejw(d)}$ is the DTFT of ${dt}$ (time shift: $x[n-k] \leftrightarrow e^{{-j\omega k}}X_d(\omega)$, "
               r"with $\delta[n] \leftrightarrow 1$), and the DTFT is unique.")
    else:
        x0 = math.sin(math.pi * float(d)) / (math.pi * float(d))
        arg = rf"\left(n {sg} {T(abs(d))}\right)"
        why = (rf"$n$ is an integer, so ${dt}$ is not a signal. The inverse DTFT is $x[n] = \dfrac{{1}}{{2\pi}}\int_{{-\pi}}^{{\pi}}"
               rf"e^{{j\omega{arg}}}\,d\omega = \dfrac{{\sin\left(\pi{arg}\right)}}{{\pi{arg}}}$, a sinc sampled off its zeros: "
               rf"nonzero for every $n$ (e.g. $x[0] \approx {x0:.3f}$).")
    return "FA2023 #1(b)", text, ans, why, {"d": pj(d)}


# ----------------------------------------------------------------------------- G5: symmetry
def ej(c):
    """e^{j c pi} with the sign of c folded into the exponent: 1, e^{j\\tfrac{\\pi}{3}}, e^{-j\\tfrac{\\pi}{2}}."""
    c = Fr(c)
    if c == 0:
        return "1"
    return rf"e^{{{'-' if c < 0 else ''}j{pi_tex(abs(c))}}}"


def hermitian_values():
    rel = random.choice(["eq", "conj"])
    if random.random() < 0.2:
        a = random.choice([Fr(0), Fr(1)])
        b = a + random.choice([Fr(2), Fr(-2), Fr(4), Fr(6), Fr(1)])
    else:
        a = random.choice([Fr(1, 3), Fr(1, 4), Fr(2, 3), Fr(1, 2), Fr(3, 4), Fr(1, 6), Fr(5, 6)])
        b = random.choice([a + 2, -a, 2 - a, a + 1, a - 2, 4 - a])

    def even_int(q):
        return q.denominator == 1 and q.numerator % 2 == 0

    a_int = a.denominator == 1
    if rel == "eq":
        ans = even_int(b - a) or (even_int(b + a) and a_int)
    else:
        ans = even_int(b + a) or (even_int(b - a) and a_int)
    star = "^{*}" if rel == "conj" else ""
    A, B = pi_tex(a), pi_tex(b)
    text = rf"If $x[n]$ is real-valued, then ${xd_at(A)} = X_d{star}\left({B}\right)$."
    if ans and rel == "eq":
        why = rf"${B}$ and ${A}$ differ by a multiple of $2\pi$, and every DTFT is $2\pi$-periodic."
    elif ans and a_int:                                      # a = 0 or pi: X_d(a) is real, b = a modulo 2pi
        what = r"X_d(0) = \sum_n x[n]" if a == 0 else r"X_d(\pi) = \sum_n(-1)^n x[n]"
        why = (rf"${B}$ and ${A}$ differ by a multiple of $2\pi$, so $X_d\left({B}\right) = {xd_at(A)}$ (periodicity), and "
               rf"${what}$ is real for real $x$, so it equals its conjugate.")
    elif ans and b == -a:
        why = (rf"This is Hermitian symmetry itself: $X_d(-\omega) = X_d^*(\omega)$ at $\omega = {A}$ gives "
               rf"$X_d\left({B}\right) = X_d^*\left({A}\right)$; conjugate both sides.")
    elif ans:
        why = (rf"${B} = -{A} + {pi_tex(b + a)}$, so $X_d\left({B}\right) = X_d\left(-{A}\right) = X_d^*\left({A}\right)$ "
               "(periodicity, then Hermitian symmetry); conjugate both sides.")
    else:
        rhs = ej(b) if rel == "conj" else ej(-b)
        rule = (r"Periodicity relates $\omega$ and $\omega + 2\pi$; Hermitian symmetry relates $X_d(-\omega)$ and $X_d^*(\omega)$; "
                "neither gives this identity. ")
        why = (rule + rf"Counterexample: the real $x[n] = \delta[n-1]$ has $X_d(\omega) = e^{{-j\omega}}$, so the left side is "
               rf"${ej(-a)}$ and the right side is ${rhs}$, which differ.")
    return "Lecture 14, slide 13", text, ans, why, {"rel": rel, "a": pj(a), "b": pj(b)}


# id: (X_d tex, real?, X_d(-w) tex, X_d^*(w) tex)
XD_REAL = {
    "w2_ecos": (r"\omega^2 e^{j\cos\omega}", False, r"\omega^2 e^{j\cos\omega}", r"\omega^2 e^{-j\cos\omega}"),
    "w2_e2w": (r"\omega^2 e^{-j2\omega}", True, r"\omega^2 e^{j2\omega}", r"\omega^2 e^{j2\omega}"),
    "jw": (r"j\omega", True, r"-j\omega", r"-j\omega"),
    "w": (r"\omega", False, r"-\omega", r"\omega"),
    "absw_ew2": (r"\lvert\omega\rvert e^{j\omega^2}", False, r"\lvert\omega\rvert e^{j\omega^2}", r"\lvert\omega\rvert e^{-j\omega^2}"),
    "twocos_ew": (r"(2 + \cos\omega)\,e^{-j\omega}", True, r"(2 + \cos\omega)\,e^{j\omega}", r"(2 + \cos\omega)\,e^{j\omega}"),
    "twocos_eabs": (r"(2 + \cos\omega)\,e^{j\lvert\omega\rvert}", False, r"(2 + \cos\omega)\,e^{j\lvert\omega\rvert}",
                    r"(2 + \cos\omega)\,e^{-j\lvert\omega\rvert}"),
    "cos2_esin": (r"\cos^2(\omega)\,e^{j\sin\omega}", True, r"\cos^2(\omega)\,e^{-j\sin\omega}", r"\cos^2(\omega)\,e^{-j\sin\omega}"),
}
XD_REAL_SRC = {"w2_ecos": "Lecture 14, slide 12", "cos2_esin": "Lecture 15, slide 10", "jw": "Lecture 14, Hermitian symmetry",
               "w": "HW6 #7 (purely imaginary signal)"}


def real_from_Xd():
    sid = random.choice(sorted(XD_REAL))
    tex, ans, neg, conj = XD_REAL[sid]
    text = rf"The signal $x[n]$ whose DTFT is $X_d(\omega) = {tex}$ for $-\pi \le \omega \le \pi$ is real-valued."
    if ans:
        why = (rf"$X_d(-\omega) = {neg} = X_d^*(\omega)$ for every $\omega$: the DTFT is Hermitian (even magnitude, odd phase), "
               r"so $x[n]$ is real (Lecture 14: real $\iff$ Hermitian).")
    else:
        why = (rf"$X_d(-\omega) = {neg}$ but $X_d^*(\omega) = {conj}$: not Hermitian, so $x[n]$ is complex "
               r"(real $\iff$ $X_d(-\omega) = X_d^*(\omega)$, Lecture 14).")
        if sid == "w":
            why += r" Here $X_d(-\omega) = -X_d^*(\omega)$: $x[n]$ is purely imaginary."
    return XD_REAL_SRC.get(sid, "Lecture 14, Hermitian symmetry"), text, ans, why, {"id": sid}


# ----------------------------------------------------------------------------- G6: properties
def prop_shift():
    k = random.choice([1, 2, 3, -1, -2])
    form = random.choice(["ok", "ok", "sign", "freq"])
    ok = rf"{ejw(k)}X_d(\omega)"
    claim = {"ok": ok, "sign": rf"{ejw(-k)}X_d(\omega)",
             "freq": rf"X_d(\omega {'-' if k > 0 else '+'} {abs(k)})"}[form]
    text = rf"If $x[n] \leftrightarrow X_d(\omega)$, then the DTFT of $y[n] = {shift_tex('x', k)}$ is $Y_d(\omega) = {claim}$."
    if form == "ok":
        why = (rf"Substitute $m = n {'-' if k > 0 else '+'} {abs(k)}$: $\sum_n x[n {'-' if k > 0 else '+'} {abs(k)}]e^{{-j\omega n}} = "
               rf"{ejw(k)}\sum_m x[m]e^{{-j\omega m}}$ (time-shift property, Lecture 14).")
    elif form == "sign":
        why = (rf"Time shift: $Y_d(\omega) = {ok}$. The stated ${claim}$ belongs to ${shift_tex('x', -k)}$, the shift in the "
               "other direction.")
    else:
        why = (rf"A shift in time multiplies by a linear phase: $Y_d(\omega) = {ok}$. Shifting $X_d$ along $\omega$ is what "
               r"modulation by $e^{j\omega_0 n}$ does.")
    return "Lecture 14, time-shift property", text, form == "ok", why, {"k": k, "form": form}


def prop_modulation():
    w0 = random.choice([Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3)])
    kind = random.choice(["cos", "exp", "alt"])
    good = random.random() < 0.5
    W = pi_tex(w0)
    mm, pp = rf"X_d\left(\omega - {W}\right)", rf"X_d\left(\omega + {W}\right)"
    if kind == "cos":
        sig = rf"x[n]\cos\left({W}n\right)"
        ok, bad = rf"\tfrac12 {mm} + \tfrac12 {pp}", rf"{mm} + {pp}"
        reason = (rf"$\cos\left({W}n\right) = \tfrac12 e^{{j{W}n}} + \tfrac12 e^{{-j{W}n}}$, and each exponential shifts the spectrum: "
                  "two half-height copies.")
    elif kind == "exp":
        sig = rf"e^{{j{W}n}}x[n]"
        ok, bad = mm, pp
        reason = rf"$\sum_n x[n]e^{{j{W}n}}e^{{-j\omega n}} = \sum_n x[n]e^{{-j(\omega - {W})n}}$: the spectrum moves to the right."
    else:
        sig = r"(-1)^n x[n]"
        ok, bad = r"X_d(\omega - \pi)", r"-X_d(\omega)"
        reason = r"$(-1)^n = e^{j\pi n}$ is a modulation by $\pi$, not a constant factor: $Y_d(\omega) = X_d(\omega - \pi)$."
    claim = ok if good else bad
    text = rf"If $x[n] \leftrightarrow X_d(\omega)$, then the DTFT of $y[n] = {sig}$ is $Y_d(\omega) = {claim}$."
    why = reason + ("" if good else rf" So $Y_d(\omega) = {ok}$, not ${bad}$.")
    return "Lecture 14, frequency shift and modulation", text, good, why, {"kind": kind, "w0": pj(w0), "good": good}


def parseval():
    L = random.choice([3, 4])
    n0 = random.choice([-1, 0, -2]) if L == 3 else random.choice([-1, -2, 0])
    while True:
        x = [random.choice([1, 2, 3, -1, -2, -3]) for _ in range(L)]
        E, S, x0 = sum(v * v for v in x), sum(x), x[-n0]
        if S != x0:
            break
    kind = random.choice(["energy", "area"])
    good = random.random() < 0.5
    seq = fmt.tex_seq(x, n0)
    if kind == "energy":
        lhs = r"\int_{-\pi}^{\pi}\lvert X_d(\omega)\rvert^2\,d\omega"
        if good:
            K, kp = 2 * E, True
        else:
            K, kp = random.choice([(E, True), (E, False)])
        why = (rf"Parseval: ${lhs} = 2\pi\sum_n\lvert x[n]\rvert^2 = 2\pi\cdot{E} = {pi_tex(2 * E)}$.")
    else:
        lhs = r"\int_{-\pi}^{\pi}X_d(\omega)\,d\omega"
        K, kp = (2 * x0, True) if good else (2 * S, True)
        why = (rf"Setting $n = 0$ in the inverse DTFT: ${lhs} = 2\pi x[0] = 2\pi\cdot({x0}) = {pi_tex(2 * x0)}$"
               + ("." if good else rf"; $2\pi\sum_n x[n] = {pi_tex(2 * S)}$ would be $2\pi X_d(0)$."))
    claim = pi_tex(K) if kp else str(K)
    text = rf"For $x[n] = {seq}$, $\displaystyle {lhs} = {claim}$."
    if not good and kind == "energy":
        why += rf" The stated ${claim}$ drops a factor of {'2' if kp else r'$2\pi$'}."
    return "Lecture 14, Parseval; HW5 #2", text, good, why, {"x": x, "n0": n0, "kind": kind, "K": K, "kpi": kp}


# ----------------------------------------------------------------------------- G7: zero padding
def zp_scaled():
    c = random.choice([2, 3, 4])
    form = random.choice(["same", "same", "exam", "upsample"])
    claim = {"same": r"Y_d(\omega) = X_d(\omega)", "exam": rf"X_d(\omega) = Y_d({c}\omega)",
             "upsample": rf"Y_d(\omega) = X_d({c}\omega)"}[form]
    text = (rf"Let $x[n]$ be a length-$N$ signal with DTFT $X_d(\omega)$. We zero-pad $x[n]$ to length ${c}N$ to obtain $y[n]$ "
            rf"with DTFT $Y_d(\omega)$. We may say that ${claim}$.")
    base = r"Appending zeros adds only zero terms to the DTFT sum, so $Y_d(\omega) = X_d(\omega)$"
    if form == "same":
        why = base + "."
    else:
        why = (base + rf"; the stated relation would need $X_d(\omega) = X_d({c}\omega)$, false in general (e.g. $x = \delta[n-1]$). "
               rf"Inserting ${c - 1}$ zero{'s' if c > 2 else ''} between the samples (upsampling) is what gives $X_d({c}\omega)$.")
    return "SP2025 #1(d)", text, form == "same", why, {"c": c, "form": form}


# ----------------------------------------------------------------------------- G8: eigenfunctions, sinusoids
TWOCOS = {Fr(1, 4): r"\sqrt{2}", Fr(1, 3): "1", Fr(1, 2): "0", Fr(2, 3): "-1", Fr(3, 4): r"-\sqrt{2}"}


def cos_multiple():
    while True:
        w0 = random.choice([Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4)])
        kind = random.choice(["sym", "sym", "delay", "delay", "fir2"])
        c = random.choice([1, 2, 3])
        k = random.choice([1, 2, 3])
        if kind == "sym" and c == 1 and w0 == Fr(2, 3):
            continue                      # H_d(w0) = 0
        if kind == "fir2" and w0 == Fr(1, 2):
            continue                      # H_d(w0) = 0
        break
    W = pi_tex(w0)
    if kind == "sym":
        htex = rf"\delta[n+1] + {'' if c == 1 else c}\delta[n] + \delta[n-1]"
        ans = True
        tc = TWOCOS[w0]
        if tc in ("1", "0", "-1"):
            v = c + int(tc)
            coef = {1: "", -1: "-"}.get(v, str(v))
        else:
            coef = rf"\left({c} {'+' if tc[0] != '-' else '-'} \sqrt{{2}}\right)"
        why = (rf"$h$ is real and even, so $H_d(\omega) = {c} + 2\cos\omega$ is real: "
               rf"$y[n] = H_d\left({W}\right)\cos\left({W}n\right) = {coef}\cos\left({W}n\right)$.")
        prm = {"kind": kind, "c": c}
    elif kind == "delay":
        htex = rf"\delta[n-{k}]"
        kw = k * w0
        ans = kw.denominator == 1
        why = (rf"$y[n] = x[n-{k}] = \cos\left({W}n - {pi_tex(kw)}\right)$" +
               (rf", and ${pi_tex(kw)}$ is a multiple of $\pi$, so $y[n] = {'-' if kw.numerator % 2 else ''}x[n]$." if ans else
                rf": the phase shift ${pi_tex(kw)}$ is not a multiple of $\pi$, so $y$ is not a multiple of $x$ "
                r"(compare the ratios $y[n]/x[n]$ at two values of $n$ where $x[n] \ne 0$)."))
        prm = {"kind": kind, "k": k}
    else:
        htex = r"\delta[n] + \delta[n-2]"
        ans = False
        why = (rf"$H_d(\omega) = 1 + e^{{-j2\omega}} = 2e^{{-j\omega}}\cos\omega$, so $y[n] = 2\cos\left({W}\right)"
               rf"\cos\left({W}n - {W}\right)$: a phase shift of ${W}$, not a multiple of $\pi$, so $y$ is not a multiple of $x$.")
        prm = {"kind": kind}
    prm["w0"] = pj(w0)
    text = (rf"The LTI system with impulse response $h[n] = {htex}$ maps $x[n] = \cos\left({W}n\right)$ (for all $n$) to an output "
            r"that is a constant multiple of $x[n]$.")
    return "FA2021 #1(a); Lectures 15-16", text, ans, why, prm


# id: (H_d tex, real h?, H_d(-w) tex, H_d^*(w) tex, source, frequencies to avoid)
SYS_H = {
    "absw_esin": (r"\lvert\omega\rvert e^{-j\pi\sin\omega}", True, r"\lvert\omega\rvert e^{j\pi\sin\omega}",
                  r"\lvert\omega\rvert e^{j\pi\sin\omega}", "FA2024 #4", []),
    "cos2_esin": (r"\cos^2(\omega)\,e^{j\sin\omega}", True, r"\cos^2(\omega)\,e^{-j\sin\omega}",
                  r"\cos^2(\omega)\,e^{-j\sin\omega}", "Lecture 15, slide 10", [Fr(1, 2)]),
    "jw_esin": (r"j\omega\,e^{j\pi\sin\omega}", True, r"-j\omega\,e^{-j\pi\sin\omega}", r"-j\omega\,e^{-j\pi\sin\omega}",
                "SP2025 #4", []),
    "w_ecos": (r"\omega\,e^{j\pi\cos\omega}", False, r"-\omega\,e^{j\pi\cos\omega}", r"\omega\,e^{-j\pi\cos\omega}",
               "SP2023 #6", [Fr(1, 3), Fr(2, 3)]),
    "half_eabs": (r"\left(\tfrac12 + \cos\omega\right)e^{j\lvert\omega\rvert}", False,
                  r"\left(\tfrac12 + \cos\omega\right)e^{j\lvert\omega\rvert}", r"\left(\tfrac12 + \cos\omega\right)e^{-j\lvert\omega\rvert}",
                  "FA2023 #4", [Fr(2, 3)]),
    "w_esin": (r"\omega\,e^{j\sin\omega}", False, r"-\omega\,e^{-j\sin\omega}", r"\omega\,e^{-j\sin\omega}", "HW6 #7", []),
}


def sinusoid_formula():
    sid = random.choice(sorted(SYS_H))
    tex, ans, neg, conj, src, avoid = SYS_H[sid]
    w0 = random.choice([w for w in (Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4)) if w not in avoid])
    W = pi_tex(w0)
    text = (rf"The LTI system with frequency response $H_d(\omega) = {tex}$, $-\pi \le \omega \le \pi$, maps "
            rf"$x[n] = \cos\left({W}n\right)$ to $y[n] = \left\lvert H_d\left({W}\right)\right\rvert\cos\left({W}n + "
            rf"\angle H_d\left({W}\right)\right)$.")
    if ans:
        why = (rf"$H_d(-\omega) = {neg} = H_d^*(\omega)$: Hermitian, so $h[n]$ is real, and for a real system the two "
               "exponentials of the cosine recombine into exactly this formula (Lecture 15).")
    else:
        why = (rf"$H_d(-\omega) = {neg}$ but $H_d^*(\omega) = {conj}$, so $h[n]$ is complex. Split the cosine: "
               rf"$y[n] = \tfrac12 H_d\left({W}\right)e^{{j{W}n}} + \tfrac12 H_d\left(-{W}\right)e^{{-j{W}n}}$, and here "
               rf"$H_d\left(-{W}\right) \ne H_d^*\left({W}\right)$, so $y[n]$ is not of the stated form.")
    return src, text, ans, why, {"id": sid, "w0": pj(w0)}


# ----------------------------------------------------------------------------- G9: magnitude, phase, group delay
GD_TAPS = {
    "sym": [[1, 2, 1], [1, 1, 1], [2, 1, 2], [1, 3, 3, 1], [1, 0, 1], [1, 2, 2, 1], [1, 0, 0, 1], [1, 2, 3, 2, 1]],
    "anti": [[1, 0, -1], [1, -1], [1, 2, -2, -1], [2, 1, 0, -1, -2]],
    "asym": [[3, 2, 1], [1, 2], [2, 1, 1], [1, 1, 2], [3, 1]],
}


def group_delay():
    kind = random.choice(["sym", "sym", "anti", "asym"])
    h = random.choice(GD_TAPS[kind])
    L = len(h)
    M = Fr(L - 1, 2)
    if kind == "asym":
        g0 = Fr(sum(n * v for n, v in enumerate(h)), sum(h))
        gpi = Fr(sum(n * v * (-1) ** n for n, v in enumerate(h)), sum(v * (-1) ** n for n, v in enumerate(h)))
        tau = g0
        ans = False
        why = (r"The taps are not symmetric about any point, so the phase is not linear: "
               rf"$\tau_{{gd}}(0) = \dfrac{{\sum_n n\,h[n]}}{{\sum_n h[n]}} = {T(g0)}$ but "
               rf"$\tau_{{gd}}(\pi) = \dfrac{{\sum_n n(-1)^n h[n]}}{{\sum_n(-1)^n h[n]}} = {T(gpi)}$ (compare Lecture 16, Example 3).")
    else:
        tau = M if random.random() < 0.6 else Fr(L - 1)
        ans = tau == M
        form = (rf"e^{{-j{T(M) if M != 1 else ''}\omega}}A(\omega)" if kind == "sym" else
                rf"je^{{-j{T(M) if M != 1 else ''}\omega}}A(\omega)")
        why = (rf"The taps are {'symmetric' if kind == 'sym' else 'antisymmetric'} about $M = {T(M)}$, so $H_d(\omega) = {form}$ "
               rf"with $A$ real: the phase is a line of slope $-{T(M)}$ (plus $\pm\pi$ jumps and $2\pi$ wraps), and "
               rf"$\tau_{{gd}} = -\dfrac{{d\angle H_d}}{{d\omega}} = {T(M)}$" + ("." if ans else rf", not ${T(tau)}$."))
    text = (rf"The FIR filter $h[n] = {fmt.tex_seq(h, 0)}$ has a constant group delay $\tau_{{gd}}(\omega) = {T(tau)}$ "
            rf"sample{'' if tau == 1 else 's'} "
            r"at every frequency where $H_d(\omega) \ne 0$.")
    return "Lecture 16, group delay", text, ans, why, {"h": h, "tau": pj(tau), "kind": kind}


# taps, M, amplitude terms [(coef, k)] meaning sum coef*cos(k w), nonneg?, (point where negative, value) or reason
MAG_TAPS = [
    ([1, 2, 1], Fr(1), [(2, 0), (2, 1)], True, r"it equals $4\cos^2\tfrac{\omega}{2}$"),
    ([1, 0, 1], Fr(1), [(2, 1)], False, (r"\pi", "-2")),
    ([1, 1, 1], Fr(1), [(1, 0), (2, 1)], False, (r"\pi", "-1")),
    ([1, 4, 1], Fr(1), [(4, 0), (2, 1)], True, r"it is at least $2$"),
    ([1, 3, 1], Fr(1), [(3, 0), (2, 1)], True, r"it is at least $1$"),
    ([2, 1, 2], Fr(1), [(1, 0), (4, 1)], False, (r"\pi", "-3")),
    ([1, 2, 3, 2, 1], Fr(2), [(3, 0), (4, 1), (2, 2)], True, r"it equals $(1 + 2\cos\omega)^2$"),
    ([1, 0, 2, 0, 1], Fr(2), [(2, 0), (2, 2)], True, r"it equals $4\cos^2\omega$"),
    ([1, 1, 1, 1, 1], Fr(2), [(1, 0), (2, 1), (2, 2)], False, (r"\tfrac{2\pi}{3}", "-1")),
    ([1, 0, 0, 1], Fr(3, 2), [(2, Fr(3, 2))], False, (r"\tfrac{2\pi}{3}", "-2")),
    ([1, 3, 3, 1], Fr(3, 2), [(6, Fr(1, 2)), (2, Fr(3, 2))], True, r"it equals $8\cos^3\tfrac{\omega}{2}$ and $\cos\tfrac{\omega}{2} \ge 0$ for $\lvert\omega\rvert \le \pi$"),
]


def cos_tex(k):
    k = Fr(k)
    if k == 0:
        return ""
    if k.denominator == 1:
        return r"\cos\omega" if k == 1 else rf"\cos {k.numerator}\omega"
    p = k.numerator
    return rf"\cos\tfrac{{{'' if p == 1 else p}\omega}}{{{k.denominator}}}"


def mag_nonneg():
    h, M, terms, ans, info = random.choice(MAG_TAPS)
    A = fmt.tex_sum([(c, cos_tex(k)) for c, k in terms])
    lin = rf"e^{{-j{T(M) if M != 1 else ''}\omega}}"
    text = (rf"The FIR filter $h[n] = {fmt.tex_seq(h, 0)}$ has $H_d(\omega) = {lin}\left({A}\right)$. Its magnitude response is "
            rf"$\lvert H_d(\omega)\rvert = {A}$ for $-\pi \le \omega \le \pi$.")
    if ans:
        why = (rf"The real amplitude ${A}$ is never negative on $[-\pi, \pi]$ ({info}), and $\lvert {lin}\rvert = 1$, so it is "
               "the magnitude.")
    else:
        wn, val = info
        why = (rf"${A}$ is negative for some $\omega$ (at $\omega = {wn}$ it equals ${val}$), and a magnitude is never negative: "
               rf"$\lvert H_d(\omega)\rvert = \left\lvert {A}\right\rvert$, and the sign goes into the phase as $\pm\pi$ (Lecture 16).")
    return "Lecture 16; FA2024 #3", text, ans, why, {"h": h, "M": pj(M), "terms": [[c, pj(k)] for c, k in terms]}


# id: (H_d tex, has a pi jump?, reason, source)
PHASE_JUMP = {
    "cos4": (r"\cos(4\omega)", True, r"$\cos 4\omega$ is real and changes sign at $\omega = \pm\tfrac{\pi}{8}, \pm\tfrac{3\pi}{8}, \dots$: the phase switches between $0$ and $\pm\pi$.", "Lecture 16, slide 6"),
    "H2": (r"e^{-j2\omega}\left(2 + \cos 2\omega\right)", False, r"$2 + \cos 2\omega \ge 1 \gt 0$, so the phase is $-2\omega$ with $2\pi$ wraps at $\pm\tfrac{\pi}{2}$ only.", "Lecture 16, slide 6"),
    "H3": (r"je^{-j3\omega}\sin(3\omega)", True, r"The real factor $\sin 3\omega$ changes sign at $\omega = 0, \pm\tfrac{\pi}{3}, \pm\tfrac{2\pi}{3}$.", "Lecture 16, slide 6"),
    "H4": (r"e^{-j2\omega}\left(1 - 2\cos\omega - \cos 2\omega\right)", True, r"The real factor $-2(\cos^2\omega + \cos\omega - 1)$ is $-2$ at $\omega = 0$ and $2$ at $\pm\pi$: it changes sign where $\cos\omega = \tfrac{\sqrt5 - 1}{2}$.", "Lecture 16, slide 6"),
    "H5": (r"e^{-j2\omega}\left(3 - 2\cos\omega - \cos 2\omega\right)", False, r"The real factor equals $2(1 - \cos\omega)(2 + \cos\omega) \ge 0$: it touches $0$ at $\omega = 0$ without changing sign, so the phase is just $-2\omega$ (wrapped).", "Lecture 16, slide 6"),
    "sp2023": (r"2e^{-j3\omega}\cos(3\omega)", True, r"$\cos 3\omega$ changes sign at $\omega = \pm\tfrac{\pi}{6}, \pm\tfrac{\pi}{2}, \pm\tfrac{5\pi}{6}$.", "SP2023 #5"),
    "sp2025": (r"e^{-j3\omega}\left(2 + 2\cos 3\omega\right)", False, r"$2 + 2\cos 3\omega \ge 0$ (double zeros at $\pm\tfrac{\pi}{3}$ and $\pm\pi$), so the phase is $-3\omega$ with $2\pi$ wraps only.", "SP2025 #3"),
    "delay": (r"e^{-j3\omega}", False, r"$\lvert H_d\rvert = 1$ and the phase is the line $-3\omega$, wrapped into $[-\pi, \pi]$: $2\pi$ wraps only.", "Lecture 16, slides 10-11"),
}


def phase_jump():
    sid = random.choice(sorted(PHASE_JUMP))
    tex, ans, why, src = PHASE_JUMP[sid]
    text = (rf"The phase response $\angle H_d(\omega)$ of the system with $H_d(\omega) = {tex}$ has a jump of $\pm\pi$ somewhere in "
            r"$-\pi \lt \omega \lt \pi$ ($2\pi$ wraps do not count).")
    return src, text, ans, why, {"id": sid}


def principal_phase():
    while True:
        k = random.choice([1, 2, 3])
        w0 = random.choice([Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4), Fr(1, 6), Fr(5, 6)])
        kw = k * w0
        if not (kw.denominator == 1 and kw.numerator % 2 == 1):      # avoid an angle of exactly +-pi
            break
    claim = -kw
    ans = claim > -1
    W = pi_tex(w0)
    text = (rf"The phase response of the delay $h[n] = \delta[n-{k}]$ at $\omega = {W}$ is $\angle H_d\left({W}\right) = {pi_tex(claim)}$ "
            r"(phase responses are principal angles in $[-\pi, \pi]$).")
    if ans:
        why = rf"$H_d\left({W}\right) = e^{{-j{pi_tex(kw)}}}$, and ${pi_tex(claim)}$ already lies in $[-\pi, \pi]$."
    else:
        m = 0
        while claim + 2 * m <= -1:
            m += 1
        why = (rf"$H_d\left({W}\right) = e^{{-j{pi_tex(kw)}}}$, but ${pi_tex(claim)}$ lies outside $[-\pi, \pi]$: adding "
               rf"${pi_tex(2 * m)}$ gives the principal angle ${pi_tex(claim + 2 * m)}$ (the same complex number; the unwrapped "
               rf"line $-{'' if k == 1 else k}\omega$ still says the delay is ${k}$ sample{'' if k == 1 else 's'}).")
    return "Lecture 16, slides 10-11", text, ans, why, {"k": k, "w0": pj(w0)}


# ----------------------------------------------------------------------------- assembly
GROUPS = {
    "exist": ["ztrans_always", "roc_substitute", "abs_summable"],
    "periodic": ["finite_periodic", "periodic_only_infinite", "bandlimited_zero"],
    "impulses": ["finite_cos_impulses", "impulse_pair"],
    "unique": ["inverse_finite", "fractional_delay", "uniqueness"],
    "symmetry": ["real_hermitian", "real_X0_Xpi", "real_Xd_real", "Xd_real_x_real", "hermitian_values", "real_from_Xd"],
    "props": ["prop_shift", "prop_modulation", "conv_product", "product_product", "parseval"],
    "zeropad": ["zp_can_change", "zp_equal", "zp_scaled"],
    "eigen": ["cos_eigen_every", "exp_eigen", "real_cos_formula", "cos_multiple", "sinusoid_formula"],
    "phase": ["group_delay", "mag_nonneg", "phase_jump", "principal_phase"],
}

PARAM = {f.__name__: f for f in (roc_substitute, abs_summable, bandlimited_zero, finite_cos_impulses, impulse_pair,
                                 inverse_finite, fractional_delay, hermitian_values, real_from_Xd, prop_shift,
                                 prop_modulation, parseval, zp_scaled, cos_multiple, sinusoid_formula, group_delay,
                                 mag_nonneg, phase_jump, principal_phase)}


def make(sid):
    if sid in FIXED:
        src, text, ans, why = FIXED[sid]
        return src, text, ans, why, {}
    return PARAM[sid]()


def generate(data):
    groups = random.sample(list(GROUPS), 6)
    items = []
    for i, g in enumerate(groups):
        sid = random.choice(GROUPS[g])
        src, text, ans, why, prm = make(sid)
        items.append({"name": f"tf{i + 1}", "label": "abcdef"[i], "sid": sid, "stmt_html": text, "src": src,
                      "param": prm, "ans": "True" if ans else "False", "why": why,
                      "choices": mc.true_false(ans)})
    data["params"]["statements"] = items
