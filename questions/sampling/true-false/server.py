"""Six True/False statements on sampling, reconstruction and the DFT (exam family "True/False",
problem 1 on 7 of 7 past Midterm 2 exams; Lecture 17 and the material after it).

Fixed statements are the later-topic statements of the study site's Midterm 2 True/False bank
(official answers, Q17-Q43; Q42 is left out because its answer depends on a convention) plus facts
stated in the 2019 course notes (Singer & Munson: ideal and zero-order-hold D/A, A/D -> H_d -> D/A,
DFT properties, zero padding, windows, FFT). Parameterized statements are randomized and their
answer recomputed, e.g. SP2023 #1(b) "the Nyquist rate of x(3t) is F/3" becomes "... of x(ct) is ..."
(True iff the claim is cF). One statement per topic group, six groups out of nine."""

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


def pin(c):
    """c*pi*n in LaTeX, with a space after a bare \\pi: \\tfrac{\\pi}{3}n, \\pi n, 2\\pi n."""
    t = pi_tex(c)
    return t + (" n" if t.endswith(r"\pi") else "n")


def ej(c):
    """e^{j c pi} with the sign of c folded into the exponent; 1 for c = 0."""
    c = Fr(c)
    if c == 0:
        return "1"
    return rf"e^{{{'-' if c < 0 else ''}j{pi_tex(abs(c))}}}"


def omega_t(f):
    """2*pi*f*t for f in Hz, as LaTeX: 100\\pi t."""
    w = 2 * Fr(f)
    if w.denominator == 1:
        return rf"{w.numerator}\pi t"
    return rf"{float(w):g}\pi t"


def seq_tex(vals):
    return r"\{" + r",\ ".join(T(v) for v in vals) + r"\}"


def big(n):
    """Integer with thin spaces: 1\\,048\\,576."""
    return f"{n:,}".replace(",", r"\,") if n >= 10000 else str(n)


def even_int(q):
    q = Fr(q)
    return q.denominator == 1 and q.numerator % 2 == 0


def fold(q):
    """Digital frequency in units of pi, folded into [0, 1], of a tone at q*fs (q = f0/fs)."""
    r = Fr(q) % 1
    return 2 * min(r, 1 - r)


# ----------------------------------------------------------------------------- fixed statements
# sid: (source, statement, answer, reason)
FIXED = {
    "bandlimited_def": (
        "SP2023 #1(e)",
        r"If $x_a(t)$ is a bandlimited continuous-time signal, then there must exist a finite $\Omega_{\max}$ such that "
        r"$X_a(\Omega) = 0$ for $\lvert\Omega\rvert \gt \Omega_{\max}$, where $X_a(\Omega)$ is the CTFT of $x_a(t)$.",
        True,
        r"That is the definition of bandlimited (Lecture 17: $X_a(\Omega) = 0$ for $\lvert\Omega\rvert \gt 2\pi B$)."),
    "every_signal_nyquist": (
        "FA2023 #1(d)",
        r"Every continuous-time signal has a minimum sampling frequency to avoid aliasing during analog-to-digital conversion.",
        False,
        r"Only bandlimited signals have a Nyquist rate. A signal with energy at arbitrarily high frequencies (a rectangular "
        r"pulse, whose CTFT is a sinc) aliases at every sampling rate."),
    "larger_T": (
        "FA2021 #1(e), reworded",
        r"A bandlimited signal $x_a(t)$ is sampled without aliasing. Increasing the sampling period $T$ (while still avoiding "
        r"aliasing) moves the band edge of the DTFT $X_d(\omega)$ of the samples closer to $\omega = 0$.",
        False,
        r"An analog frequency $\Omega$ lands at $\omega = \Omega T$, so the band edge $2\pi BT$ moves outward, toward $\pi$: "
        r"a larger $T$ stretches $X_d$ in $\omega$ (only its height $\tfrac1T$ shrinks)."),
    "recon_causal": (
        "FA2021 #1(c)",
        r"The ideal reconstruction filter is causal.",
        False,
        r"Its impulse response $\dfrac{\sin(\pi t/T)}{\pi t/T}$ is nonzero for $t \lt 0$."),
    "da_through_samples": (
        "2019 notes, ideal reconstruction formula",
        r"The output of an ideal D/A converter with sampling period $T$ passes through the samples: $x_a(nT) = x[n]$ for every "
        r"integer $n$.",
        True,
        r"$x_a(t) = \sum_m x[m]\dfrac{\sin(\pi(t - mT)/T)}{\pi(t - mT)/T}$; at $t = nT$ the sinc centred at $nT$ equals $1$ and "
        r"every other one is $\sin(\pi(n - m))/(\pi(n - m)) = 0$."),
    "zoh_flat_lpf": (
        "2019 notes, zero-order hold",
        r"For every input sequence, following a zero-order-hold D/A converter by an ideal analog low-pass filter with constant "
        r"gain on $\lvert\Omega\rvert \lt \pi/T$ reproduces the output of the ideal D/A converter exactly.",
        False,
        r"The zero-order hold multiplies the spectrum by $Te^{-j\Omega T/2}\,\dfrac{\sin(\Omega T/2)}{\Omega T/2}$: a half-sample "
        r"delay and a droop toward $\Omega = \pm\pi/T$. Matching the ideal D/A needs a filter that undoes both, "
        r"$e^{j\Omega T/2}\big/\tfrac{\sin(\Omega T/2)}{\Omega T/2}$ in the band, which boosts the high frequencies; a flat gain does not."),
    "adda_lti_bandlimited": (
        "2019 notes, DT processing of CT signals",
        r"If $x_a(t)$ is bandlimited to less than $f_s/2$ Hz and ideal A/D and D/A converters with rate $f_s$ are used, then "
        r"A/D $\to H_d(\omega) \to$ D/A behaves as an LTI system with frequency response $H_d(\Omega T)$ for "
        r"$\lvert\Omega\rvert \lt \pi/T$ and $0$ beyond.",
        True,
        r"Without aliasing only the $n = 0$ copy of $X_a$ lies in $\lvert\Omega\rvert \lt \pi/T$: $Y_d(\Omega T) = H_d(\Omega T)"
        r"\tfrac1T X_a(\Omega)$, and the ideal D/A multiplies by $T$ in that band and removes the rest, so "
        r"$Y_a(\Omega) = H_d(\Omega T)X_a(\Omega)$."),
    "adda_always_ti": (
        "2019 notes, DT processing of CT signals",
        r"The system A/D $\to H_d(\omega) \to$ D/A (ideal converters, sampling period $T$) is time-invariant for every input "
        r"$x_a(t)$, bandlimited or not.",
        False,
        r"With aliasing it is not. Take $H_d = 1$, $T = 1$ s, $x_a(t) = \cos(1.5\pi t)$: the output is $\cos(0.5\pi t)$. Delaying the "
        r"input by $0.5$ s gives samples $\cos(1.5\pi n - 0.75\pi) = \cos(0.5\pi n + 0.75\pi)$ and the output "
        r"$\cos(0.5\pi t + 0.75\pi)$, not the delayed output $\cos(0.5\pi t - 0.25\pi)$. (The notes: linear, but in general not "
        r"time-invariant.)"),
    "dft_same": (
        "FA2021 #1(b)",
        r"Two different sequences of the same length can have the same DFT.",
        False,
        r"The DFT is invertible: the inverse DFT recovers $x[n]$, $0 \le n \le N - 1$, from $X[k]$, so equal DFTs mean equal sequences."),
    "zp_mainlobe": (
        "2019 notes, zero-padding in DFT spectral analysis",
        r"Zero-padding a length-$N$ record of a sinusoid before computing the DFT narrows the main lobe of the sinusoid's peak, "
        r"so closely spaced sinusoids become easier to resolve.",
        False,
        r"Zero padding gives a denser set of samples of exactly the same DTFT; the main lobe keeps its width $4\pi/N$. "
        r"Resolution needs a longer observation $NT$, not more zeros."),
    "more_samples_narrower": (
        "2019 notes, DFT spectral analysis of sinusoids",
        r"Recording a sinusoid for longer (more samples at the same sampling rate) narrows the main lobe of its peak in the DTFT "
        r"of the recorded samples.",
        True,
        r"The $N$ samples of a sinusoid give two periodic sincs whose main lobes have width $4\pi/N$: a larger $N$ makes them "
        r"narrower."),
    "fft_is_dft": (
        "FA2019 #1(j)",
        r"FFT is just an efficient algorithm to evaluate DFT.",
        True,
        r"The FFT returns exactly the DFT values, with about $N\log_2 N$ multiply-accumulates instead of $N^2$ (radix-2)."),
    "fft_power2": (
        "2019 notes, zero-padding in DFT spectral analysis",
        r"A common reason to zero-pad a sequence before computing its DFT is to reach a length, such as a power of 2, for "
        r"which a radix-2 FFT can be used.",
        True,
        r"Zero padding does not change the DTFT, so padding to the next power of 2 costs nothing and allows a fast radix-2 FFT "
        r"(the notes give this as one of the two reasons to zero-pad)."),
}


# ----------------------------------------------------------------------------- G1: Nyquist rate
def mul_F(r):
    r = Fr(r)
    if r == 1:
        return "F"
    if r.denominator == 1:
        return f"{r.numerator}F"
    return rf"\tfrac{{{'' if r.numerator == 1 else r.numerator}F}}{{{r.denominator}}}"


def scaled_arg(c):
    c = Fr(c)
    return f"{c.numerator}t" if c.denominator == 1 else f"t/{c.denominator}"


def nyq_scale():
    c = random.choice([Fr(2), Fr(3), Fr(4), Fr(1, 2), Fr(1, 3)])
    form = random.choice(["ok", "ok", "inv", "same"])
    claim = {"ok": c, "inv": 1 / c, "same": Fr(1)}[form]
    text = (rf"If the Nyquist rate of $x_a(t)$ is $F$ Hz, then the Nyquist rate of $y_a(t) = x_a({scaled_arg(c)})$ is "
            rf"${mul_F(claim)}$ Hz.")
    why = (rf"Time scaling, $x_a(ct) \leftrightarrow \frac{{1}}{{\lvert c\rvert}}X_a(\Omega/c)$ with $c = {T(c)}$: every frequency is "
           rf"multiplied by ${T(c)}$, so the highest frequency and the Nyquist rate become ${mul_F(c)}$ Hz"
           + ("." if form == "ok" else rf", not ${mul_F(claim)}$ Hz."))
    return "SP2021 #1(b), SP2023 #1(b)", text, form == "ok", why, {"c": pj(c), "claim": pj(claim)}


def nyq_product():
    f1, f2 = random.sample([50, 100, 150, 200, 250, 300, 400], 2)
    form = random.choice(["prod", "prod", "square", "sum"])
    good = random.random() < 0.5
    if form == "prod":
        expr = rf"\cos({omega_t(f1)})\sin({omega_t(f2)})"
        R_true, R_slip = 2 * (f1 + f2), 2 * max(f1, f2)
        why = (rf"$\cos a\sin b = \tfrac12[\sin(a + b) - \sin(a - b)]$: the product contains tones at ${f1 + f2}$ Hz and "
               rf"${abs(f1 - f2)}$ Hz, so the highest frequency is ${f1 + f2}$ Hz and the Nyquist rate is ${R_true}$ Hz.")
    elif form == "square":
        expr = rf"\cos^2({omega_t(f1)})"
        R_true, R_slip = 4 * f1, 2 * f1
        why = (rf"$\cos^2\theta = \tfrac12 + \tfrac12\cos 2\theta$: squaring doubles the frequency to ${2 * f1}$ Hz, so the "
               rf"Nyquist rate is ${R_true}$ Hz.")
    else:
        expr = rf"\cos({omega_t(f1)}) + 3\sin({omega_t(f2)})"
        R_true, R_slip = 2 * max(f1, f2), 2 * (f1 + f2)
        why = (rf"A sum contains only the frequencies of its terms, ${f1}$ Hz and ${f2}$ Hz: the Nyquist rate is "
               rf"$2\cdot{max(f1, f2)} = {R_true}$ Hz.")
    R = R_true if good else R_slip
    if not good:
        why += rf" Not ${R}$ Hz."
    text = rf"The Nyquist rate of $x_a(t) = {expr}$ is ${R}$ Hz."
    return "Lecture 17; sampling and the Nyquist rate, Practice 1", text, good, why, {"form": form, "f1": f1, "f2": f2, "R": R}


# ----------------------------------------------------------------------------- G2: sampled sinusoids, aliasing
def cos_to_sequence():
    f0 = random.choice([25, 50, 60, 100])
    kind = random.choice(["cos", "cos", "sin"])
    p = random.choice([Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(1, 3), Fr(2, 3)])
    form = random.choice(["principal", "plus", "minus", "half", "half"])
    r = {"principal": p, "plus": p + 2, "minus": 2 - p, "half": p / 2}[form]     # Omega0*T in units of pi
    Tp = r / (2 * f0)
    ans = even_int(r - p) or (kind == "cos" and even_int(r + p))
    P, R = pi_tex(p), pi_tex(r)
    text = (rf"Sampling $x_a(t) = \{kind}({omega_t(f0)})$ with sampling period $T = {T(Tp)}$ s gives "
            rf"$x[n] = \{kind}\left({P}n\right)$.")
    why = rf"$x[n] = x_a(nT) = \{kind}\left({2 * f0}\pi\cdot{T(Tp)}\,n\right) = \{kind}\left({R}n\right)$"
    if form == "principal":
        why += "."
    elif form == "plus":
        why += rf", and ${R} = {P} + 2\pi$: frequencies $2\pi$ apart give the same sequence."
    elif form == "minus" and kind == "cos":
        why += (rf", and ${R} = 2\pi - {P}$: $\cos\left({R}n\right) = \cos\left(-{P}n\right) = \cos\left({P}n\right)$ "
                "(the cosine is even).")
    elif form == "minus":
        why += (rf", and ${R} = 2\pi - {P}$: $\sin\left({R}n\right) = \sin\left(-{P}n\right) = -\sin\left({P}n\right)$: the sign "
                "flips (the sine is odd).")
    else:
        why += rf": a different digital frequency, since ${R}$ is not $\pm{P}$ modulo $2\pi$."
    return "SP2021 #1(a), SP2023 #1(a), FA2023 #3", text, ans, why, {"f0": f0, "kind": kind, "p": pj(p), "T": pj(Tp)}


QS = [Fr(1, 10), Fr(1, 5), Fr(3, 10), Fr(2, 5), Fr(3, 5), Fr(7, 10), Fr(4, 5), Fr(9, 10), Fr(11, 10), Fr(6, 5), Fr(13, 10)]


def alias_reconstruct():
    fs = random.choice([200, 400, 500, 1000])
    q = random.choice(QS)
    f0 = int(fs * q)
    w = fold(q)                                   # digital frequency / pi in [0, 1]
    f1 = int(w * fs / 2)
    aliased = q > Fr(1, 2)
    good = random.random() < 0.5
    fc = f1 if good else (f0 if aliased else fs - f0)
    text = (rf"The signal $x_a(t) = \cos({omega_t(f0)})$ is sampled at $f_s = {fs}$ Hz, and the samples are passed to an ideal "
            rf"D/A converter with the same sampling rate. The output is $y_a(t) = \cos({omega_t(fc)})$.")
    why = rf"$x[n] = \cos\left({pin(2 * q)}\right)"
    why += rf" = \cos\left({pin(w)}\right)$ (aliased: ${f0}$ Hz is above $f_s/2 = {fs // 2}$ Hz)" if 2 * q != w else "$"
    why += (rf". The ideal D/A returns the one signal bandlimited to $f_s/2$ with these samples, $\cos({omega_t(f1)})$"
            + ("." if good else (rf"; it does not give back the input above $f_s/2$." if aliased else
                                 rf"; the copy at $f_s - {f0} = {fs - f0}$ Hz is removed by the ideal reconstruction filter.")))
    return "Lecture 17, sampling a sinusoid; ideal D/A", text, good, why, {"fs": fs, "f0": f0, "fc": fc}


# ----------------------------------------------------------------------------- G3: the sampled spectrum
def guard_band():
    B = random.choice([100, 200, 250, 500])
    ratio = random.choice([Fr(5, 2), Fr(3), Fr(4), Fr(2), Fr(3, 2), Fr(2)])
    fs = int(B * ratio)
    e = Fr(2 * B, fs)                             # band edge / pi
    ans = fs > 2 * B
    text = (rf"$x_a(t)$ has the triangular CTFT $X_a(\Omega) = 1 - \dfrac{{\lvert\Omega\rvert}}{{{2 * B}\pi}}$ for "
            rf"$\lvert\Omega\rvert \le {2 * B}\pi$ rad/s and $X_a(\Omega) = 0$ otherwise. It is sampled at $f_s = {fs}$ Hz. Then the "
            r"DTFT $X_d(\omega)$ of the samples is zero on some interval inside $0 \lt \omega \lt \pi$ (a guard band).")
    if ans:
        why = (rf"$B = {B}$ Hz and $f_s \gt 2B = {2 * B}$ Hz: the band edge lands at $\omega = 2\pi B/f_s = {pi_tex(e)} \lt \pi$ and "
               rf"the copy centred at $2\pi$ starts at ${pi_tex(2 - e)} \gt \pi$, so $X_d = 0$ for ${pi_tex(e)} \lt \lvert\omega\rvert \le \pi$.")
    elif fs == 2 * B:
        why = (rf"$f_s = 2B = {2 * B}$ Hz exactly: the band edge lands at $\omega = \pi$, the copies touch there, and "
               r"$X_d(\omega) \gt 0$ for every $\lvert\omega\rvert \lt \pi$: no guard band.")
    else:
        why = (rf"$f_s \lt 2B = {2 * B}$ Hz: the band edge $2\pi B/f_s = {pi_tex(e)}$ lies beyond $\pi$, the copies overlap "
               r"(aliasing) and $X_d(\omega) \gt 0$ everywhere.")
    return "Lecture 17, Nyquist criterion", text, ans, why, {"B": B, "fs": fs}


def height_scaling():
    A = random.choice([1, 2, 3])
    B = random.choice([100, 200, 250])
    fs = random.choice([3 * B, 4 * B, 1000 if 1000 > 2 * B else 5 * B])
    form = random.choice(["ok", "ok", "noT", "T"])
    V = {"ok": Fr(A * fs), "noT": Fr(A), "T": Fr(A, fs)}[form]
    amp = "" if A == 1 else str(A)
    Xa = (rf"1 - \dfrac{{\lvert\Omega\rvert}}{{{2 * B}\pi}}" if A == 1 else
          rf"{amp}\left(1 - \dfrac{{\lvert\Omega\rvert}}{{{2 * B}\pi}}\right)")
    text = (rf"$x_a(t)$ has CTFT $X_a(\Omega) = {Xa}$ for $\lvert\Omega\rvert \le {2 * B}\pi$ rad/s and $0$ otherwise, and "
            rf"$x[n] = x_a(nT)$ with $T = \tfrac{{1}}{{{fs}}}$ s. Then the DTFT of $x[n]$ satisfies $X_d(0) = {T(V)}$.")
    why = (rf"$X_d(\omega) = \frac{{1}}{{T}}\sum_k X_a\left(\frac{{\omega - 2\pi k}}{{T}}\right)$. Since $f_s = {fs} \gt 2B = {2 * B}$ Hz "
           rf"the copies do not overlap, and at $\omega = 0$ only $k = 0$ contributes: $X_d(0) = f_s X_a(0) = {fs}\cdot{A} = {A * fs}$ "
           r"(every copy is scaled by $1/T$)" + ("." if form == "ok" else rf", not ${T(V)}$."))
    return "Lecture 17, the sampled spectrum", text, form == "ok", why, {"A": A, "B": B, "fs": fs, "V": pj(V)}


def sifting():
    while True:
        c = random.choice([1, 2, 3, -1, -2])
        t0 = random.choice([-1, 0, 2, 3, 5])
        a, b = random.choice([(3, None), (0, 4), (-2, 1), (1, 6), (-3, 2), (4, None)])
        if t0 not in (a, b) and t0 + c != 0:
            break
    inside = a < t0 and (b is None or t0 < b)
    true_val = t0 + c if inside else 0
    good = random.random() < 0.5
    V = true_val if good else (0 if inside else t0 + c)
    b_tex = r"\infty" if b is None else str(b)
    tc = f"t {'+' if c > 0 else '-'} {abs(c)}"
    dl = r"\delta(t)" if t0 == 0 else rf"\delta(t {'-' if t0 > 0 else '+'} {abs(t0)})"
    text = rf"The value of $\displaystyle\int_{{{a}}}^{{{b_tex}}}({tc})\,{dl}\,dt$ is ${V}$."
    if inside:
        why = rf"The impulse sits at $t = {t0}$, inside the limits: sifting gives $({tc})\big|_{{t={t0}}} = {t0 + c}$."
    else:
        why = (rf"The impulse sits at $t = {t0}$, outside $[{a}, {b_tex}{')' if b is None else ']'}$: the integrand is zero on the whole interval, so the "
               r"integral is $0$ (sifting only fires inside the limits).")
    return "SP2021 #1(d)", text, good, why, {"c": c, "t0": t0, "a": a, "b": b, "V": V}


# ----------------------------------------------------------------------------- G4: D/A conversion
def sinc_tex(A, k):
    arg = rf"\pi(t - {tT(k)})/T"
    return rf"{A}\,\dfrac{{\sin\left({arg}\right)}}{{{arg}}}"


def da_impulse():
    A = random.choice([2, 3, 4, 5])
    k = random.choice([1, 2, 3, 5])
    form = random.choice(["sinc", "sinc", "sinc", "impulse", "impulse", "rect"])
    claim = {"sinc": sinc_tex(A, k), "impulse": rf"{A}\delta(t - {tT(k)})",
             "rect": rf"{A}\left[u(t - {tT(k)}) - u(t - {tT(k + 1)})\right]"}[form]
    text = (rf"A discrete-time signal $x[n] = {A}\delta[n - {k}]$ is passed to an ideal D/A converter with sampling period $T$ "
            rf"to recover $x_a(t)$. The recovered signal is $x_a(t) = {claim}$.")
    why = (r"An ideal D/A converter turns each sample into a sinc, $x_a(t) = \sum_n x[n]\,\dfrac{\sin(\pi(t - nT)/T)}{\pi(t - nT)/T}$, "
           rf"so $x_a(t) = {sinc_tex(A, k)}$: it equals ${A}$ at $t = {tT(k)}$, $0$ at the other sampling instants, and is nonzero "
           "in between.")
    if form == "impulse":
        why += " A sample never becomes an impulse."
    elif form == "rect":
        why += " The rectangle is what a zero-order hold produces."
    return "FA2023 #1(e), FA2024 #1(c)", text, form == "sinc", why, {"A": A, "k": k, "form": form}


def da_value():
    while True:
        a = random.choice([1, 2, 3])
        b = random.choice([1, 2, 3, -1, -2])
        n1 = random.choice([1, 2])
        n2 = n1 + random.choice([2, 3])
        form = random.choice(["zero_int", "zero_int", "sample", "zero_at_sample", "zero_half", "zero_half"])
        if form == "zero_int":
            t0 = Fr(random.choice([m for m in range(0, n2 + 3) if m not in (n1, n2)]))
            V = 0
        elif form == "sample":
            t0 = Fr(random.choice([n1, n2]))
            V = a if t0 == n1 else b
        elif form == "zero_at_sample":
            t0 = Fr(random.choice([n1, n2]))
            V = 0
        else:
            t0 = Fr(random.choice([n1, n2 - 1, n2])) + Fr(1, 2)
            V = 0
        xa = sum(v * (1.0 if float(t0) == m else math.sin(math.pi * (float(t0) - m)) / (math.pi * (float(t0) - m)))
                 for v, m in ((a, n1), (b, n2)))
        if form != "zero_half" or abs(xa) > 0.05:
            break
    ans = form in ("zero_int", "sample")
    bt = rf"{'+' if b > 0 else '-'} {'' if abs(b) == 1 else abs(b)}\delta[n - {n2}]"
    at = rf"{'' if a == 1 else a}\delta[n - {n1}]"
    t0t = str(t0.numerator) if t0.denominator == 1 else f"{float(t0):g}"
    text = (rf"Let $x[n] = {at} {bt}$. We recover $x_a(t)$ by ideal digital-to-analog conversion with sampling period "
            rf"$T = 1$ s. Then $x_a(t) = {V}$ at $t = {t0t}$ s.")

    def term(v, m, first):
        s = ("" if v > 0 else "-") if first else (" + " if v > 0 else " - ")
        return rf"{s}{'' if abs(v) == 1 else abs(v)}\dfrac{{\sin(\pi(t - {m}))}}{{\pi(t - {m})}}"
    why = rf"$x_a(t) = {term(a, n1, True)}{term(b, n2, False)}$. "
    if t0.denominator == 1:
        xv = a if t0 == n1 else (b if t0 == n2 else 0)
        why += (rf"At the integer $t = {t0t}$ every sinc vanishes except one centred there (if any), so $x_a({t0t}) = x[{t0t}] = {xv}$: "
                "ideal interpolation passes through the samples.")
    else:
        why += (rf"At $t = {t0t}$ s, halfway between two samples, neither sinc vanishes: $x_a({t0t}) \approx {xa:.3f} \ne 0$; "
                "only at the sampling instants does $x_a$ reduce to a sample value.")
    return "SP2025 #1(c)", text, ans, why, {"a": a, "b": b, "n1": n1, "n2": n2, "t0": pj(t0), "V": V}


def tT(m):
    return {0: "0", 1: "T"}.get(m, f"{m}T")


def zoh_value():
    while True:
        y = [random.choice([1, 2, 3, 4, 5, -1, -2, -3]) for _ in range(4)]
        m = random.choice([0, 1, 2])
        if y[m] != y[m + 1]:
            break
    Tms = random.choice([1, 2, 4])
    f = random.choice([Fr(1, 4), Fr(1, 2), Fr(3, 4)])
    t0 = (m + f) * Tms
    form = random.choice(["hold", "hold", "next", "lin"])
    V = {"hold": Fr(y[m]), "next": Fr(y[m + 1]), "lin": y[m] + f * (y[m + 1] - y[m])}[form]
    t0t = str(t0.numerator) if t0.denominator == 1 else f"{float(t0):g}"
    text = (rf"A zero-order-hold D/A converter with sampling period $T = {Tms}$ ms converts $y[n] = {fmt.tex_seq(y, 0)}$ "
            rf"(and $y[n] = 0$ for other $n$) into $y_a(t)$. Then $y_a(t) = {T(V)}$ at $t = {t0t}$ ms.")
    why = (rf"A zero-order hold keeps each sample for one sampling period, $y_a(t) = y[n]$ for $nT \le t \lt (n + 1)T$: a staircase. "
           rf"$t = {t0t}$ ms lies in $[{tT(m)}, {tT(m + 1)})$, so $y_a(t) = y[{m}] = {y[m]}$")
    why += {"hold": ".", "next": rf"; ${y[m + 1]}$ is the next step, which starts only at $t = {tT(m + 1)}$.",
            "lin": rf"; ${T(V)}$ would be a straight-line interpolation between $y[{m}]$ and $y[{m + 1}]$."}[form]
    return "2019 notes, zero-order hold", text, form == "hold", why, {"y": y, "T": Tms, "t0": pj(t0), "V": pj(V)}


# ----------------------------------------------------------------------------- G5: A/D -> H_d -> D/A
def adda_lowpass():
    while True:
        fs = random.choice([400, 500, 1000])
        q = random.choice(QS)
        wc = random.choice([Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(1, 3), Fr(2, 3)])
        w = fold(q)
        if w != wc:
            break
    f0 = int(fs * q)
    aliased = q > Fr(1, 2)
    passes = w < wc
    fa = int(w * fs / 2)
    claim = random.choice(["x", "zero", "alias"] if aliased else ["x", "x", "zero", "image"])
    ctex = {"x": rf"\cos({omega_t(f0)})", "zero": "0", "alias": rf"\cos({omega_t(fa)})",
            "image": rf"\cos({omega_t(fs - f0)})"}[claim]
    ans = {"x": passes and not aliased, "zero": not passes, "alias": passes, "image": False}[claim]
    text = (rf"$x_a(t) = \cos({omega_t(f0)})$ is sampled at $f_s = {fs}$ Hz, filtered by an ideal digital low-pass filter with "
            rf"cutoff $\omega_c = {pi_tex(wc)}$ ($H_d(\omega) = 1$ for $\lvert\omega\rvert \lt \omega_c$, $0$ for "
            r"$\omega_c \lt \lvert\omega\rvert \le \pi$), and converted back by an ideal D/A converter with the same $f_s$. "
            rf"The output is $y_a(t) = {ctex}$.")
    why = rf"$x[n] = \cos\left({pin(2 * q)}\right)"
    if aliased:
        why += rf" = \cos\left({pin(w)}\right)$ (aliased: ${f0}$ Hz is above $f_s/2 = {fs // 2}$ Hz)"
    else:
        why += "$"
    if passes:
        why += (rf". The digital frequency ${pi_tex(w)}$ is below the cutoff ${pi_tex(wc)}$, so the filter passes it, and the D/A "
                rf"maps it to $\Omega = \omega f_s = {2 * fa}\pi$ rad/s: $y_a(t) = \cos({omega_t(fa)})$.")
    else:
        why += rf". The digital frequency ${pi_tex(w)}$ is above the cutoff ${pi_tex(wc)}$, so the filter removes it: $y_a(t) = 0$."
    fc = {"x": f0, "zero": None, "alias": fa, "image": fs - f0}[claim]
    return "2019 notes, DT processing of CT signals; FA2019 #9", text, ans, why, {"fs": fs, "f0": f0, "wc": pj(wc), "fc": fc}


def adda_cutoff():
    fs = random.choice([1200, 2400, 4800])
    wc = random.choice([Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4)])
    true_fc = int(wc / 2 * fs)
    form = random.choice(["ok", "ok", "double", "half"])
    fc = {"ok": true_fc, "double": int(wc * fs), "half": true_fc // 2}[form]
    text = (rf"Inputs bandlimited to $f_s/2$ are sampled at $f_s = {fs}$ Hz, filtered by an ideal digital low-pass filter with "
            rf"cutoff $\omega_c = {pi_tex(wc)}$, and converted back by an ideal D/A converter. The overall system acts as an ideal "
            rf"analog low-pass filter with cutoff ${fc}$ Hz.")
    why = (r"With ideal converters and no aliasing the overall response is $H_d(\Omega T)$ for $\lvert\Omega\rvert \lt \pi/T$: "
           r"the cutoff $\omega_c$ maps to $\Omega_c = \omega_c/T = \omega_c f_s$, i.e. "
           rf"$f_c = \dfrac{{\omega_c}}{{2\pi}}f_s = {true_fc}$ Hz" + ("." if form == "ok" else rf", not ${fc}$ Hz."))
    return "2019 notes, DT processing of CT signals", text, form == "ok", why, {"fs": fs, "wc": pj(wc), "fc": fc}


# ----------------------------------------------------------------------------- G6: DFT basics
def dft_X0():
    N = random.choice([4, 5, 6, 8])
    kind = random.choice(["X0", "X0", "XN2"]) if N % 2 == 0 else "X0"
    while True:
        x = [random.choice([1, 2, 3, 4, 5, 6, 7, 9, -1, -2, -3, -4, -6, -7]) for _ in range(N)]
        S, A = sum(x), sum(v * (-1) ** n for n, v in enumerate(x))
        if S != 0 and A != 0 and S != A and x[0] != S:
            break
    form = random.choice(["ok", "ok", "zero", "other"])
    if kind == "X0":
        k, true_val = 0, S
        V = {"ok": S, "zero": 0, "other": x[0]}[form]
        why = rf"$X[0] = \sum_{{n=0}}^{{{N - 1}}} x[n]e^{{0}} = {S}$" + ("." if form == "ok" else rf", not ${V}$.")
    else:
        k, true_val = N // 2, A
        V = {"ok": A, "zero": 0, "other": S}[form]
        why = (rf"$X[{N // 2}] = \sum_n x[n]e^{{-j2\pi\cdot{N // 2}n/{N}}} = \sum_n(-1)^n x[n] = {A}$"
               + {"ok": ".", "zero": ", not $0$.", "other": rf", not ${S}$, which is $X[0] = \sum_n x[n]$."}[form])
    text = (rf"Let $\{{x[n]\}}_{{n=0}}^{{{N - 1}}} = {seq_tex(x)}$. Consider the corresponding ${N}$-point DFT "
            rf"$\{{X[k]\}}_{{k=0}}^{{{N - 1}}}$. Then $X[{k}] = {V}$.")
    return "FA2019 #1(i), SP2021 #1(e)", text, V == true_val, why, {"x": x, "k": k, "V": V}


def dft_real_k():
    N = random.choice([5, 6, 7, 8])
    kind = random.choice(["real", "conj"])
    good = random.random() < 0.5
    if kind == "real":
        trues = [0] + ([N // 2] if N % 2 == 0 else [])
        k1 = random.choice(trues) if good else random.choice([k for k in range(1, N) if k not in trues])
        k2 = k1
        claim = rf"$X[{k1}]$ is real-valued"
        ans = (2 * k1) % N == 0
    else:
        k1 = random.choice(range(1, N))
        if good:
            k2 = N - k1
        else:
            k2 = random.choice([k for k in range(N) if (k1 + k) % N != 0])
        claim = rf"$X[{k1}] = X^*[{k2}]$"
        ans = (k1 + k2) % N == 0
    text = (rf"If $\{{x[n]\}}_{{n=0}}^{{{N - 1}}}$ is real-valued and $\{{X[k]\}}_{{k=0}}^{{{N - 1}}}$ is its DFT, then {claim}.")
    if ans:
        why = (rf"For real $x$, $X[k] = X^*[\langle N - k\rangle_N]$ (2019 notes, Table 2.5)"
               + (rf": $X[{k1}]$ equals its own conjugate, so it is real (e.g. $X[0] = \sum_n x[n]$)." if kind == "real" else
                  rf", and ${k2} = {N} - {k1}$."))
    else:
        rhs = rf"X[{k1}]^*" if kind == "real" else rf"X^*[{k2}]"
        why = (rf"Real $x$ only gives $X[k] = X^*[\langle N - k\rangle_N]$. Counterexample: $x[n] = \delta[n-1]$ has "
               rf"$X[k] = e^{{-j2\pi k/{N}}}$, so $X[{k1}] = {ej(Fr(-2 * k1, N))}$ differs from ${rhs} = {ej(Fr(2 * k2, N))}$.")
    return "SP2021 #1(c); 2019 notes, DFT symmetry", text, ans, why, {"N": N, "kind": kind, "k1": k1, "k2": k2}


def idft_periodic():
    N = 4
    X = random.choice([[1, 2, 3, 4], [4, 3, 2, 1], [2, 0, 1, 1], [3, 1, -1, 1], [1, 1, 2, 4], [6, 2, 0, 2]])
    x0 = Fr(sum(X), N)
    form = random.choice(["zero", "periodic", "value", "value_zero"])
    claim = {"zero": rf"$x[n]$ must be zero for $n \lt 0$ and $n \gt {N - 1}$",
             "periodic": rf"$x[n + {N}] = x[n]$ for every integer $n$",
             "value": rf"$x[{N}] = {T(x0)}$", "value_zero": rf"$x[{N}] = 0$"}[form]
    text = (rf"If $x[n]$ is the inverse DFT of $X[k] = {seq_tex(X)}$, $0 \le k \le {N - 1}$, that is "
            rf"$x[n] = \frac{{1}}{{{N}}}\sum_{{k=0}}^{{{N - 1}}}X[k]e^{{j2\pi kn/{N}}}$, then {claim}.")
    why = (rf"The inverse-DFT formula is ${N}$-periodic in $n$, since $e^{{j2\pi k(n + {N})/{N}}} = e^{{j2\pi kn/{N}}}$: "
           rf"$x[{N}] = x[0] = \frac{{1}}{{{N}}}\sum_k X[k] = {T(x0)}$" + (r", so $x[n]$ is not zero outside $0 \le n \le 3$."
                                                                       if form == "zero" else "."))
    ans = form in ("periodic", "value")
    V = {"value": x0, "value_zero": Fr(0)}.get(form)
    return "SP2023 #1(d), SP2025 #1(b)", text, ans, why, {"X": X, "form": form, "V": None if V is None else pj(V)}


# ----------------------------------------------------------------------------- G7: DFT and DTFT, zero padding
def dft_dtft_sample():
    N = random.choice([8, 10, 12, 14, 16])
    k0 = random.choice(range(N // 2 + 1, N))
    form = random.choice(["neg", "neg", "pos", "mirror", "half"])
    th = {"neg": Fr(2 * k0, N) - 2, "pos": Fr(2 * k0, N), "mirror": Fr(2 * (N - k0), N), "half": Fr(k0, N)}[form]
    text = (rf"The DTFT and DFT of a length-${N}$ sequence $x[n]$ are given by $X_d(\omega)$ and $X[k]$, respectively. We know "
            rf"then that $X[{k0}] = X_d\left({pi_tex(th)}\right)$.")
    why = (rf"$X[k] = X_d\left(\frac{{2\pi k}}{{{N}}}\right)$, so $X[{k0}] = X_d\left({pi_tex(Fr(2 * k0, N))}\right)$")
    if form == "neg":
        why += rf" $= X_d\left({pi_tex(th)}\right)$ by $2\pi$-periodicity."
    elif form == "pos":
        why += "."
    elif form == "mirror":
        why += (rf"; ${pi_tex(th)}$ is the frequency of $X[{N - k0}]$, which differs from $X[{k0}]$ in general (the two are "
                "conjugates only when $x$ is real).")
    else:
        why += rf"; ${pi_tex(th)}$ is a different frequency (the DFT samples are $\frac{{2\pi}}{{{N}}}$ apart)."
    ans = form in ("neg", "pos")
    return "FA2023 #1(h)", text, ans, why, {"N": N, "k0": k0, "th": pj(th)}


ZP_PAIRS = [(20, 32), (12, 16), (10, 25), (6, 9), (8, 12), (16, 24), (10, 16)]


def zp_index():
    L, M = random.choice(ZP_PAIRS)
    g = math.gcd(L, M)
    a, b = L // g, M // g
    j = random.randint(1, g - 1)
    form = random.choice(["ok", "ok", "same", "off"])
    k0, m0 = {"ok": (a * j, b * j), "same": (a * j, a * j), "off": (a * j, b * j + 1)}[form]
    ans = Fr(k0, L) == Fr(m0, M)
    text = (rf"Assume $x[n]$ is a finite-duration sequence of length ${L}$, and $y[n]$ is obtained by zero-padding $x[n]$ to "
            rf"length ${M}$: $y[n] = x[n]$ for $n = 0, \dots, {L - 1}$ and $y[n] = 0$ for $n = {L}, \dots, {M - 1}$. Let "
            rf"$\{{X[k]\}}_{{k=0}}^{{{L - 1}}}$ and $\{{Y[k]\}}_{{k=0}}^{{{M - 1}}}$ be the DFTs of $x$ and $y$. Then $X[{k0}] = Y[{m0}]$.")
    why = (rf"Both DFTs sample the same DTFT: $X[k] = X_d\left(\frac{{2\pi k}}{{{L}}}\right)$ and $Y[m] = X_d\left(\frac{{2\pi m}}{{{M}}}\right)$. "
           rf"Here $\frac{{2\pi\cdot{k0}}}{{{L}}} = {pi_tex(Fr(2 * k0, L))}$ and $\frac{{2\pi\cdot{m0}}}{{{M}}} = {pi_tex(Fr(2 * m0, M))}$"
           + (": the same frequency." if ans else ": different frequencies, so the values differ in general."))
    return "FA2021 #1(d), FA2024 #1(d)", text, ans, why, {"L": L, "M": M, "k0": k0, "m0": m0}


# ----------------------------------------------------------------------------- G8: spectral analysis, windows
WORD = {1: "one", 2: "two"}


def cos_bins():
    kind = random.choice(["digital", "digital", "analog"])
    if kind == "digital":
        N = random.choice([16, 20, 24, 32, 64])
        c = random.choice([Fr(1), Fr(2), Fr(3), Fr(3, 2), Fr(5, 2), Fr(N // 2)])
        f0 = fs = None
    else:
        fs = random.choice([40, 50, 64, 100])
        N = fs                                              # one second of samples
        f0 = random.choice([Fr(5), Fr(10), Fr(fs // 2), Fr(15, 2), Fr(25, 2)])
        c = f0                                              # periods in one second
    w0 = 2 * c / N                                          # digital frequency / pi
    if c.denominator == 1:
        count = 1 if c == N // 2 and N % 2 == 0 else 2
    else:
        count = N
    if count <= 2 and random.random() < 0.75:
        m = count
    else:
        m = random.choice([1, 2]) if count > 2 else 3 - count
    ans = m == count
    s = "s" if m > 1 else ""
    if kind == "digital":
        text = (rf"The ${N}$-point DFT of the finite-length signal $x[n] = \cos\left({pin(w0)}\right)$, $0 \le n \le {N - 1}$, "
                rf"has exactly {WORD[m]} nonzero value{s}.")
    else:
        text = (rf"An analog signal $x_a(t) = \cos({omega_t(f0)})$ is sampled at $f_s = {fs}$ Hz for 1 second to obtain "
                rf"$\{{x[n]\}}_{{n=0}}^{{{N - 1}}}$. Its DFT $X[k]$ is nonzero for exactly {WORD[m]} value{s} of $k$, $0 \le k \le {N - 1}$.")
    if count == 2:
        why = (rf"$x[n] = \cos\left({pin(w0)}\right)$ with ${pi_tex(w0)} = \frac{{2\pi\cdot{c}}}{{{N}}}$: exactly ${c}$ periods in "
               rf"${N}$ samples, so $x[n] = \tfrac12 e^{{j2\pi\cdot{c}n/{N}}} + \tfrac12 e^{{-j2\pi\cdot{c}n/{N}}}$ falls on the bins "
               rf"$k = {c}$ and $k = {N - c}$: $X[{c}] = X[{N - c}] = {N // 2}$ and every other bin is $0$ (two nonzero values).")
    elif count == 1:
        why = (rf"$x[n] = \cos(\pi n) = (-1)^n = e^{{j2\pi\cdot{N // 2}n/{N}}}$: one exponential exactly on bin ${N // 2}$, so "
               rf"$X[{N // 2}] = {N}$ and all other bins are $0$ (one nonzero value).")
    else:
        why = (rf"$x[n] = \cos\left({pin(w0)}\right)$ and $\frac{{{N}\omega_0}}{{2\pi}} = {T(c)}$ is not an integer: the record holds "
               rf"${T(c)}$ periods, the cosine falls between the bins and leaks into all ${N}$ of them.")
    prm = {"kind": kind, "N": N, "c": pj(c), "m": m}
    return "FA2019 #1(f), FA2023 #1(g), SP2025 #1(e)", text, ans, why, prm


def hamming():
    good = random.random() < 0.5
    if good:
        side, main = "lower", "wider"
    else:
        side, main = random.choice([("lower", "narrower"), ("lower", "narrower"), ("higher", "wider"), ("higher", "narrower")])
    text = (rf"Applying the Hamming window (instead of the rectangular window) for spectral analysis leads to {side} sidelobes "
            rf"and {main} main lobes in the spectral components.")
    why = ("It is a trade: the Hamming window gives much lower sidelobes but a main lobe about twice as wide (2019 notes, "
           "zero-padding and windows)")
    if good:
        why += "."
    else:
        wrong = [w for w, ok in (("the sidelobes", side == "lower"), ("the main lobes", main == "wider")) if not ok]
        why += f"; the statement gets {' and '.join(wrong)} wrong."
    return "SP2025 #1(f)", text, good, why, {"side": side, "main": main}


def rect_mainlobe():
    N = random.choice([8, 10, 16, 20, 25, 32])
    form = random.choice(["ok", "ok", "half", "double"])
    c = {"ok": Fr(2, N), "half": Fr(1, N), "double": Fr(4, N)}[form]
    text = (rf"The DTFT of the length-${N}$ rectangular window $w[n] = 1$ for $0 \le n \le {N - 1}$ (and $0$ otherwise) first "
            rf"reaches zero, for $\omega \gt 0$, at $\omega = {pi_tex(c)}$.")
    why = (rf"$W_d(\omega) = e^{{-j{T(Fr(N - 1, 2))}\omega}}\,\dfrac{{\sin({N}\omega/2)}}{{\sin(\omega/2)}}$ first vanishes where "
           rf"${N}\omega/2 = \pi$, i.e. at $\omega = \frac{{2\pi}}{{{N}}} = {pi_tex(Fr(2, N))}$ (a main lobe of width ${pi_tex(Fr(4, N))}$)"
           + ("." if form == "ok" else rf", not ${pi_tex(c)}$."))
    return "2019 notes, periodic sinc (Fig. 2.20)", text, form == "ok", why, {"N": N, "c": pj(c)}


# ----------------------------------------------------------------------------- G9: circular operations, FFT
def circ_shift():
    N = random.choice([4, 5])
    x = random.sample([1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3], N)
    m = random.choice([1, 2]) if N == 5 else 1          # N = 4, m = 2 would make both rotations equal
    form = random.choice(["ok", "ok", "left", "linear"])
    ok = [x[(n - m) % N] for n in range(N)]
    claim = {"ok": ok, "left": [x[(n + m) % N] for n in range(N)], "linear": [0] * m + x[:N - m]}[form]
    ex = rf"{2 * m}\pi k/{N}"                               # 2*pi*m*k/N with the factors multiplied out
    text = (rf"Let $\{{x[n]\}}_{{n=0}}^{{{N - 1}}} = {seq_tex(x)}$ with ${N}$-point DFT $X[k]$. The sequence whose ${N}$-point DFT "
            rf"is $e^{{-j{ex}}}X[k]$ is ${seq_tex(claim)}$ ($n = 0, \dots, {N - 1}$).")
    why = (rf"Circular time shift: $x[\langle n - {m}\rangle_{{{N}}}] \leftrightarrow e^{{-j{ex}}}X[k]$. The samples rotate "
           rf"right by ${m}$, and the last {'sample wraps' if m == 1 else 'two samples wrap'} around to the front: ${seq_tex(ok)}$")
    why += {"ok": ".", "left": rf"; the stated sequence is the rotation to the left, which belongs to $e^{{j{ex}}}X[k]$.",
            "linear": "; the stated sequence is a linear shift, which drops samples instead of wrapping them around."}[form]
    return "2019 notes, circular time shift", text, claim == ok, why, {"x": x, "m": m, "claim": claim}


def circ_conv():
    N = 4
    while True:
        x = [random.choice([0, 1, 2, 3, -1]) for _ in range(N)]
        h = [random.choice([1, 2, -1])] + [random.choice([0, 1, 2, -1]) for _ in range(N - 1)]
        lin = [sum(x[i] * h[n - i] for i in range(N) if 0 <= n - i < N) for n in range(2 * N - 1)]
        circ = [lin[n] + (lin[n + N] if n + N < len(lin) else 0) for n in range(N)]
        if circ != lin[:N] and x[-1] != 0:
            break
    good = random.random() < 0.5
    claim = circ if good else lin[:N]
    text = (rf"The ${N}$-point circular convolution of ${seq_tex(x)}$ and ${seq_tex(h)}$ (both indexed $n = 0, \dots, {N - 1}$) "
            rf"is ${seq_tex(claim)}$.")
    why = (rf"The linear convolution is ${seq_tex(lin)}$ (length ${2 * N - 1}$); the ${N}$-point circular convolution adds the samples "
           rf"at $n \ge {N}$ onto $n - {N}$: ${seq_tex(circ)}$"
           + ("." if good else rf". The stated sequence is the first ${N}$ samples of the linear convolution, without the wrap-around."))
    return "2019 notes, circular convolution", text, good, why, {"x": x, "h": h, "claim": claim}


def circ_equals_linear():
    Lx, Lh = random.choice([3, 4, 5, 6, 8]), random.choice([2, 3, 4, 5])
    L = Lx + Lh - 1
    good = random.random() < 0.55
    if good:
        N = random.choice([L, L + 1, 1 << (L - 1).bit_length()])
    else:
        N = random.choice(range(max(Lx, Lh), L))
    text = (rf"To compute the linear convolution of a length-${Lx}$ sequence and a length-${Lh}$ sequence with the DFT, it is enough "
            rf"to zero-pad both to length ${N}$, multiply their ${N}$-point DFTs and take the inverse ${N}$-point DFT.")
    why = (rf"The linear convolution has length ${Lx} + {Lh} - 1 = {L}$. Multiplying ${N}$-point DFTs gives the ${N}$-point circular "
           rf"convolution, which equals the linear one only if $N \ge {L}$ (otherwise the samples beyond $n = N - 1$ wrap around "
           rf"and add onto the first ones): $N = {N}$ is " + ("enough." if good else "too short."))
    return "FA2019 #1(g)", text, good, why, {"Lx": Lx, "Lh": Lh, "N": N}


def fft_count():
    m = random.choice([6, 7, 8, 9, 10, 11, 12])
    N = 2 ** m
    form = random.choice(["ok", "ok", "swap", "fullN", "double"])
    s, b = {"ok": (m, N // 2), "swap": (N // 2, m), "fullN": (m, N), "double": (2 * m, N // 2)}[form]
    text = rf"A radix-2 decimation-in-time FFT of length $N = {big(N)}$ consists of ${big(s)}$ stages of ${big(b)}$ butterflies each."
    why = (rf"Each stage halves the DFT length: $N = 2^{{{m}}}$ gives $\log_2 N = {m}$ stages of $N/2 = {big(N // 2)}$ butterflies, "
           rf"about $N\log_2 N = {big(N * m)}$ multiply-accumulates instead of $N^2 = {big(N * N)}$ for the direct DFT"
           + ("." if form == "ok" else rf", not ${big(s)}$ stages of ${big(b)}$."))
    return "2019 notes, radix-2 FFT", text, form == "ok", why, {"m": m, "s": s, "b": b}


# ----------------------------------------------------------------------------- assembly
GROUPS = {
    "nyquist": ["nyq_scale", "nyq_product", "bandlimited_def", "every_signal_nyquist"],
    "sinusoid": ["cos_to_sequence", "alias_reconstruct"],
    "spectrum": ["guard_band", "height_scaling", "larger_T", "sifting"],
    "da": ["da_impulse", "da_value", "recon_causal", "da_through_samples", "zoh_value", "zoh_flat_lpf"],
    "processing": ["adda_lowpass", "adda_cutoff", "adda_lti_bandlimited", "adda_always_ti"],
    "dft_basic": ["dft_X0", "dft_real_k", "dft_same", "idft_periodic"],
    "dft_dtft": ["dft_dtft_sample", "zp_index", "zp_mainlobe", "more_samples_narrower"],
    "spectral": ["cos_bins", "hamming", "rect_mainlobe"],
    "circ_fft": ["circ_shift", "circ_conv", "circ_equals_linear", "fft_is_dft", "fft_count", "fft_power2"],
}

PARAM = {f.__name__: f for f in (nyq_scale, nyq_product, cos_to_sequence, alias_reconstruct, guard_band, height_scaling,
                                 sifting, da_impulse, da_value, zoh_value, adda_lowpass, adda_cutoff, dft_X0, dft_real_k,
                                 idft_periodic, dft_dtft_sample, zp_index, cos_bins, hamming, rect_mainlobe, circ_shift,
                                 circ_conv, circ_equals_linear, fft_count)}


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
