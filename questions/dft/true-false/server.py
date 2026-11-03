"""Five True/False statements on the DFT and spectral analysis (exam family "True/False", problem 1 on
every past Midterm 2; Lectures 20-23 and HW9).

Most of the pool re-randomizes the past exam statements of the study site's Midterm 2 True/False bank
(sections 6-7 and the DFT halves of section 3): FA2019 #1(i) / SP2021 #1(e) (X[0] = 0?), SP2021 #1(c)
(X[0] real), FA2025 #1(a) (DFT from X_d = 3 + cos 2w), FA2021 #1(b) (same DFT), SP2023 #1(d) /
SP2025 #1(b) / FA2025 #1(b) (the inverse DFT is periodic), FA2023 #1(h) and FA2025 #6(a) (bin <-> DTFT
frequency), FA2024 #1(d) and FA2021 #1(d) (zero padding and bins), FA2019 #1(f) / FA2023 #1(g) /
SP2025 #1(e) (nonzero bins of a cosine), SP2025 #1(f) (Hamming), SP2021 #6 (CTFT, DTFT or DFT?),
FA2021 #2(c) (is H_d(w) = H(e^{jw})?), FA2025 #7 (a repeated DFT), FA2025 #8 (four DFT plots),
FA2019 #4 / SP2025 #8 (phase ramp), SP2021 #4 / SP2023 #9(b) / HW9 #2 (alternating signs), HW9 #5 and
FA2025 #8(b) (peak bin -> frequency). The rest are Lecture 22-23 facts: an on-bin cosine's peak AN/2,
full- and half-lobe separation (N > 4pi/dw), zero padding does not improve resolution, recording time
2/df does not depend on f_s, windowing is linear but not time-invariant and convolves the spectra,
main-lobe widths 4pi/N and 8pi/N. Circular-convolution computations and the FFT are left out (not on
the Fall 2026 midterm).

One statement from each of five of the eight groups; each statement is True about half the time; a
variant never has all five answers equal."""

import random
from fractions import Fraction as Fr
from math import gcd

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


def pin(c, var="n"):
    """c*pi*var: \\tfrac{\\pi}{8}n, \\pi n, 2\\pi n."""
    t = pi_tex(c)
    return t + (f" {var}" if t.endswith(r"\pi") else var)


def seq(vals):
    return r"\{" + r",\ ".join(T(v) for v in vals) + r"\}"


def dec(v):
    """A short decimal: 0.55, 0.5, 0.125."""
    return f"{float(v):g}"


def eqdec(v):
    """'= 0.0625' or, for longer decimals, '\approx 0.0156'."""
    t = f"{float(v):g}"
    if "." in t and len(t.split(".")[1]) > 4:
        return rf"\approx {float(v):.4f}"
    return "= " + t


def signed_sum(vals, alt=False):
    """1 - 2 + 3 (alt: x[0] - x[1] + x[2] with parentheses around negative entries)."""
    out = ""
    for n, v in enumerate(vals):
        if alt:
            t = T(v) if v >= 0 else rf"({T(v)})"
            out += t if n == 0 else ((" - " if n % 2 else " + ") + t)
        else:
            out += T(v) if n == 0 else ((" - " if v < 0 else " + ") + T(abs(v)))
    return out


# ============================================================================ G1: DFT values
def x0_zero():
    N = random.choice([6, 8])
    kind = random.choice(["X0", "X0", "XN2"])
    zero = random.random() < 0.45
    vals = [1, 2, 3, 4, 5, 6, 7, 9, -1, -2, -3, -4, -6, -7]
    sg = [1] * N if kind == "X0" else [(-1) ** n for n in range(N)]
    while True:
        x = [random.choice(vals) for _ in range(N)]
        if zero:
            last = -sum(s * v for s, v in zip(sg[:-1], x[:-1])) * sg[-1]
            if last == 0 or abs(last) > 9:
                continue
            x[-1] = last
        val = sum(s * v for s, v in zip(sg, x))
        if (val == 0) == zero:
            break
    k = 0 if kind == "X0" else N // 2
    text = (rf"Let $\{{x[n]\}}_{{n=0}}^{{{N - 1}}} = {seq(x)}$. Consider the corresponding ${N}$-point DFT "
            rf"$\{{X[k]\}}_{{k=0}}^{{{N - 1}}}$. Then, $X[{k}] = 0$.")
    if kind == "X0":
        why = rf"$X[0] = \sum_n x[n]e^{{0}} = \sum_n x[n] = {signed_sum(x)} = {val}$."
    else:
        why = (rf"$X[{k}] = \sum_n x[n]e^{{-j\pi n}} = \sum_n(-1)^n x[n] = {signed_sum(x, alt=True)} = {val}$.")
    return "FA2019 #1(i), SP2021 #1(e)", text, val == 0, why, {"x": x, "k": k}


def real_bins():
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
        k1 = random.choice([k for k in range(1, N) if 2 * k != N])
        k2 = N - k1 if good else random.choice([v for v in (N - 1 - k1, k1) if (k1 + v) % N and v != 0])
        claim = rf"$X[{k1}] = X^*[{k2}]$"
        ans = (k1 + k2) % N == 0
    text = rf"If $\{{x[n]\}}_{{n=0}}^{{{N - 1}}}$ is real-valued and $\{{X[k]\}}_{{k=0}}^{{{N - 1}}}$ is its DFT, then {claim}."
    if ans and kind == "real":
        what = r"X[0] = \sum_n x[n]" if k1 == 0 else rf"X[{k1}] = \sum_n (-1)^n x[n]"
        why = rf"${what}$ is a sum of real numbers" + (
            "." if k1 == 0 else rf" (equivalently $X[k] = X^*[\langle -k\rangle_N]$ and $\langle -{k1}\rangle_{{{N}}} = {k1}$).")
    elif ans:
        why = rf"For real $x$, $X[k] = X^*[\langle -k\rangle_N] = X^*[N - k]$, and ${k2} = {N} - {k1}$ (conjugate symmetry)."
    else:
        rhs = rf"X^*[{k1}]" if kind == "real" else rf"X^*[{k2}]"
        why = (rf"Real $x$ only gives $X[k] = X^*[N - k]$; the partner of bin ${k1}$ is bin ${N - k1}$. Counterexample: "
               rf"$x[n] = \delta[n-1]$ has $X[k] = e^{{-j2\pi k/{N}}}$, so $X[{k1}] = e^{{-j{pi_tex(Fr(2 * k1, N))}}}$ but "
               rf"${rhs} = e^{{j{pi_tex(Fr(2 * k2, N))}}}$.")
    return "SP2021 #1(c); Lecture 21, conjugate symmetry", text, ans, why, {"N": N, "kind": kind, "k1": k1, "k2": k2}


def fa25_dtft():
    b = random.choice([1, 2, 2, 3])
    L = 2 * b + 1
    a = random.choice([2, 3, 4, 5])
    c = random.choice([1, 1, 2, -1])
    forms = ["ok", "ok", "ok", "nosub", "halfgrid"] + (["drop"] if b > 1 else [])
    form = random.choice(forms)
    cw = r"\cos\omega" if b == 1 else rf"\cos({b}\omega)"
    coef = {1: "+ ", 2: "+ 2", -1: "- "}[c]
    Xd = rf"{a} {coef}{cw}"

    def ck(t):
        return rf"\cos\left({pin(t, 'k')}\right)"

    arg = {"ok": Fr(2 * b, L), "drop": Fr(2, L), "halfgrid": Fr(b, L)}.get(form)
    ckt = ck(arg) if arg is not None else (r"\cos k" if b == 1 else rf"\cos({b}k)")
    claim = rf"{a} {coef}{ckt}"
    xs = fmt.tex_sum([(Fr(c, 2), rf"\delta[n + {b}]"), (a, r"\delta[n]"), (Fr(c, 2), rf"\delta[n - {b}]")])
    text = (rf"The DTFT of a length-${L}$ signal $x[n]$ is given by $X_d(\omega) = {Xd}$. The DFT of $x[n]$ is then "
            rf"$X[k] = {claim}$.")
    why = (rf"The DFT evaluates the DTFT at $\omega_k = \frac{{2\pi k}}{{{L}}}$: $X[k] = X_d\left(\frac{{2\pi k}}{{{L}}}\right) = "
           rf"{a} {coef}{ck(Fr(2 * b, L))}$. (Here $x[n] = {xs}$; "
           rf"on $0 \le n \le {L - 1}$ it is its circular shift, which has the same DFT.)")
    if form != "ok":
        why += " The stated formula is a different function of $k$ (compare $k = 1$)."
    return "FA2025 #1(a)", text, form == "ok", why, {"a": a, "b": b, "c": c, "form": form}


def parseval_dft():
    N = random.choice([4, 5, 6])
    while True:
        x = [random.choice([1, 2, 3, -1, -2, -3, 0]) for _ in range(N)]
        if x[0] != 0 and x[-1] != 0 and sum(v != 0 for v in x) >= N - 1:
            break
    E = sum(v * v for v in x)
    form = random.choice(["ok", "ok", "ok", "noN", "N2", "overN"])
    V = {"ok": N * E, "noN": E, "N2": N * N * E, "overN": Fr(E, N)}[form]
    text = (rf"Let $\{{x[n]\}}_{{n=0}}^{{{N - 1}}} = {seq(x)}$ with ${N}$-point DFT $X[k]$. Then "
            rf"$\displaystyle\sum_{{k=0}}^{{{N - 1}}}\lvert X[k]\rvert^2 = {T(V)}$.")
    sq = " + ".join(rf"{T(v)}^2" if v >= 0 else rf"({T(v)})^2" for v in x)
    why = (rf"Parseval for the DFT: $\sum_k\lvert X[k]\rvert^2 = N\sum_n\lvert x[n]\rvert^2 = {N}\left({sq}\right) = {N}\cdot{E} = {N * E}$"
           + ("." if form == "ok" else rf", not ${T(V)}$."))
    return "Lectures 20-21, Parseval for the DFT", text, form == "ok", why, {"x": x, "V": pj(V)}


# ============================================================================ G2: inverse DFT, uniqueness, repeats
def idft_periodic():
    kind = random.choice(["list", "list", "fa25", "flat"])
    if kind == "list":
        N = 4
        X = random.sample([1, 2, 3, 4, 5, 6], 4)
        form = random.choice(["zero", "zero", "periodic", "value", "value_zero"])
        x0 = Fr(sum(X), N)
        claim = {"zero": rf"$x[n]$ must be zero for $n \lt 0$ and $n \gt {N - 1}$",
                 "periodic": rf"$x[n + {N}] = x[n]$ for every integer $n$",
                 "value": rf"$x[{N}] = {T(x0)}$", "value_zero": rf"$x[{N}] = 0$"}[form]
        text = (rf"If $x[n]$ is the inverse DFT of $X[k] = {seq(X)}$, $0 \le k \le {N - 1}$, then {claim}.")
        why = (rf"The inverse DFT $x[n] = \frac{{1}}{{{N}}}\sum_{{k=0}}^{{{N - 1}}}X[k]e^{{j2\pi kn/{N}}}$ is ${N}$-periodic in $n$, since "
               rf"$e^{{j2\pi k(n + {N})/{N}}} = e^{{j2\pi kn/{N}}}$: $x[{N}] = x[0] = \frac{{1}}{{{N}}}\sum_k X[k] = {T(x0)}$"
               + (rf", so $x[n]$ is not zero outside $0 \le n \le {N - 1}$." if form == "zero" else "."))
        ans = form in ("periodic", "value")
        prm = {"kind": kind, "X": X, "form": form}
    elif kind == "fa25":
        N = random.choice([3, 4, 5])
        cc = random.choice([1, 2, 3])
        X = [cc * N] + [0] * (N - 1)
        form = random.choice(["window", "window", "const"])
        if form == "window":
            items = ["0", "0", rf"\underset{{\uparrow}}{{{cc}}}"] + [str(cc)] * (N - 1) + ["0", "0"]
            body = r",\ ".join([r"\cdots"] + items + [r"\cdots"])
            text = (rf"Suppose we have a length-${N}$ DFT $\{{X[k]\}}_{{k=0}}^{{{N - 1}}} = {seq(X)}$. The signal "
                    rf"$\{{x[n]\}} = \mathrm{{DFT}}^{{-1}}\{{X[k]\}} = \{{{body}\}}$.")
        else:
            text = (rf"Suppose we have a length-${N}$ DFT $\{{X[k]\}}_{{k=0}}^{{{N - 1}}} = {seq(X)}$. Its inverse DFT, "
                    rf"$x[n] = \frac{{1}}{{{N}}}\sum_{{k=0}}^{{{N - 1}}}X[k]e^{{j2\pi kn/{N}}}$, equals ${cc}$ for every integer $n$.")
        why = (rf"$x[n] = \frac{{1}}{{{N}}}\cdot {cc * N}\cdot e^{{0}} = {cc}$ for <b>every</b> $n$: the inverse-DFT formula is "
               rf"${N}$-periodic in $n$, not zero outside $0 \le n \le {N - 1}$.")
        ans = form == "const"
        prm = {"kind": kind, "X": X, "form": form}
    else:
        N = random.choice([3, 4, 5])
        B = random.choice([2, 3, 4, 5])
        X = [B] * N
        form = random.choice(["delta_one", "delta_one", "spread", "delta_all"])
        claim = {"delta_one": rf"$x[0] = {B}$ and $x[n] = 0$ for $1 \le n \le {N - 1}$",
                 "spread": rf"$x[n] = {T(Fr(B, N))}$ for $0 \le n \le {N - 1}$",
                 "delta_all": rf"$x[n] = {B}\delta[n]$ for every integer $n$"}[form]
        text = (rf"If $x[n]$ is the inverse DFT of $X[k] = {seq(X)}$, $0 \le k \le {N - 1}$, that is "
                rf"$x[n] = \frac{{1}}{{{N}}}\sum_{{k=0}}^{{{N - 1}}}X[k]e^{{j2\pi kn/{N}}}$, then {claim}.")
        why = (rf"$x[n] = \frac{{{B}}}{{{N}}}\sum_{{k=0}}^{{{N - 1}}}e^{{j2\pi kn/{N}}}$ is ${B}$ when $n$ is a multiple of ${N}$ and $0$ "
               rf"otherwise (a geometric sum over the ${N}$ roots of unity $e^{{j2\pi k/{N}}}$): ${B}\delta[n]$ on $0 \le n \le {N - 1}$, repeated every ${N}$ "
               rf"samples, so $x[{N}] = {B}$ too.")
        ans = form == "delta_one"
        prm = {"kind": kind, "X": X, "form": form}
    return "SP2023 #1(d), SP2025 #1(b), FA2025 #1(b)", text, ans, why, prm


def dft_same():
    form = random.choice(["same", "unique", "unique", "zp_same"])
    if form == "same":
        text = r"Two different sequences of the same length can have the same DFT."
        why = (r"The DFT is invertible: the inverse DFT recovers $x[n]$, $0 \le n \le N-1$, from $X[k]$, so equal DFTs mean "
               r"equal sequences.")
        ans = False
    elif form == "unique":
        text = (r"If two length-$N$ sequences have the same $N$-point DFT, then they are equal for $0 \le n \le N-1$.")
        why = r"The inverse DFT recovers every sample from the $N$ DFT values: the DFT is one-to-one."
        ans = True
    else:
        text = (r"A length-$N$ sequence and the same sequence zero-padded to length $2N$ have the same DFT.")
        why = (r"They have the same <b>DTFT</b>, but the DFTs differ: $N$ values at $\frac{2\pi k}{N}$ versus $2N$ values at "
               r"$\frac{2\pi k}{2N}$ (the padded DFT evaluates the same DTFT twice as densely).")
        ans = False
    return "FA2021 #1(b)", text, ans, why, {"form": form}


def dup_rule():
    form = random.choice(["dup_ok", "dup_copies", "inter_ok", "inter_noh"])
    if form.startswith("dup"):
        rule = (r"$\{Y[k]\}_{k=0}^{2N-1} = \{X[0], \dots, X[N-1], X[0], \dots, X[N-1]\}$ (the $N$-point DFT of $x$ written twice)")
        claim = (r"$y[2n] = x[n]$ and $y[2n+1] = 0$ for $0 \le n \le N-1$" if form == "dup_ok" else
                 r"$y[n] = x[\langle n\rangle_N]$ for $0 \le n \le 2N-1$ (two copies of $x$)")
        why = (r"Split the length-$2N$ inverse DFT at $k = N$: $y[n] = \frac{1 + (-1)^n}{2N}\sum_{k=0}^{N-1}X[k]e^{j\frac{\pi k}{N}n}$, "
               r"which is $x[n/2]$ for even $n$ and $0$ for odd $n$ (FA2025 #7): repeating the DFT upsamples the signal.")
    else:
        rule = r"$Y[2k] = X[k]$ and $Y[2k+1] = 0$, $0 \le k \le N-1$ (a zero after each DFT value)"
        claim = (r"$y[n] = \tfrac12\,x[\langle n\rangle_N]$ for $0 \le n \le 2N-1$" if form == "inter_ok" else
                 r"$y[n] = x[\langle n\rangle_N]$ for $0 \le n \le 2N-1$ (two copies of $x$)")
        why = (r"Only the even bins survive: $y[n] = \frac{1}{2N}\sum_{k=0}^{N-1}X[k]e^{j\frac{2\pi (2k)}{2N}n} = "
               r"\frac12\cdot\frac1N\sum_k X[k]e^{j\frac{2\pi k}{N}n} = \tfrac12 x[\langle n\rangle_N]$: two copies at half height.")
    text = (rf"Let $x[n]$ be a length-$N$ signal with $N$-point DFT $X[k]$. If the length-$2N$ DFT {rule} is inverted with the "
            rf"length-$2N$ inverse DFT, the result is {claim}.")
    return "FA2025 #7", text, form in ("dup_ok", "inter_ok"), why, {"form": form}


# ============================================================================ G3: the DFT evaluates the DTFT
def bin_freq():
    N = random.choice([8, 10, 12, 14, 16, 32])
    k0 = random.choice(range(N // 2 + 1, N))
    form = random.choice(["neg", "pos", "conjneg", "conjpos", "mirror", "half"])
    th = {"neg": Fr(2 * k0, N) - 2, "pos": Fr(2 * k0, N), "conjneg": Fr(2 * (N - k0), N),
          "conjpos": Fr(2 * k0, N) - 2, "mirror": Fr(2 * (N - k0), N), "half": Fr(k0, N)}[form]
    conj = form.startswith("conj")
    real = " real-valued" if conj else ""
    text = (rf"The DTFT and DFT of a{real} length-${N}$ sequence $x[n]$ are given by $X_d(\omega)$ and $X[k]$, respectively. "
            rf"We know then that $X[{k0}] = X_d{'^*' if conj else ''}\left({pi_tex(th)}\right)$.")
    w = Fr(2 * k0, N)
    why = rf"$X[k] = X_d\left(\frac{{2\pi k}}{{{N}}}\right)$, so $X[{k0}] = X_d\left({pi_tex(w)}\right) = X_d\left({pi_tex(w - 2)}\right)$ ($2\pi$-periodicity)"
    if form in ("neg", "pos"):
        why += "."
    elif form == "conjneg":
        why += (rf"; for real $x$, $X_d^*\left({pi_tex(th)}\right) = X_d\left({pi_tex(-th)}\right)$, the same frequency.")
    elif form == "conjpos":
        why += (rf"; but $X_d^*\left({pi_tex(th)}\right) = X_d\left({pi_tex(-th)}\right) = X^*[{k0}]$, the conjugate, which differs "
                rf"unless $X[{k0}]$ is real.")
    elif form == "mirror":
        why += (rf"; ${pi_tex(th)}$ is the frequency of $X[{N - k0}]$, which differs from $X[{k0}]$ in general (they are "
                "conjugates only when $x$ is real).")
    else:
        why += rf"; ${pi_tex(th)}$ is a different frequency (the bins are $\frac{{2\pi}}{{{N}}}$ apart)."
    return "FA2023 #1(h), FA2025 #6(a)", text, form in ("neg", "pos", "conjneg"), why, {"N": N, "k0": k0, "th": pj(th), "conj": conj}


def zp_scale():
    c = random.choice([2, 3, 4])
    form = random.choice(["ok", "ok", "ok", "inv", "scaled", "same"])
    zeros = "$N$ zeros" if c == 2 else rf"${c - 1}N$ zeros"
    claim = {"ok": rf"$X[k] = Y[{c}k]$ and $X_d(\omega) = Y_d(\omega)$",
             "inv": rf"$X[{c}k] = Y[k]$ and $X_d(\omega) = Y_d(\omega)$",
             "scaled": rf"$X[k] = Y[{c}k]$ and $X_d(\omega) = Y_d({c}\omega)$",
             "same": r"$X[k] = Y[k]$ and $X_d(\omega) = Y_d(\omega)$"}[form]
    text = (rf"A length-$N$ signal $x[n]$ with DTFT $X_d(\omega)$ and DFT $X[k]$ is zero-padded with {zeros} to length ${c}N$ to "
            rf"form $y[n]$ (DTFT $Y_d$, ${c}N$-point DFT $Y[k]$). Then, {claim}.")
    why = (rf"$Y_d = X_d$ (the zeros add nothing to the DTFT sum), and $Y[{c}k] = X_d\left(\frac{{2\pi\cdot {c}k}}{{{c}N}}\right) = "
           rf"X_d\left(\frac{{2\pi k}}{{N}}\right) = X[k]$")
    why += {"ok": ".", "inv": r"; $Y[k]$ sits at $\frac{2\pi k}{" + str(c) + r"N}$, a different frequency from $X[" + str(c) + "k]$.",
            "scaled": rf"; $Y_d(\omega) = X_d({c}\omega)$ would be upsampling (zeros <i>between</i> samples), not padding.",
            "same": rf"; $Y[k]$ for $k \ne 0$ sits at $\frac{{2\pi k}}{{{c}N}}$, not at $\frac{{2\pi k}}{{N}}$."}[form]
    return "FA2024 #1(d), SP2025 #1(d)", text, form == "ok", why, {"c": c, "form": form}


ZP_PAIRS = [(20, 32), (12, 16), (10, 25), (6, 9), (8, 12), (16, 24), (10, 16), (6, 30), (8, 60), (18, 30)]


def zp_index():
    L, M = random.choice(ZP_PAIRS)
    g = gcd(L, M)
    a, b = L // g, M // g
    j = random.randint(1, g - 1)
    form = random.choice(["ok", "ok", "same", "off", "conj_ok", "conj_bad"])
    k0, m0 = a * j, b * j
    if form == "same":
        m0 = k0
    elif form == "off":
        m0 = m0 + 1
    elif form == "conj_ok":
        m0 = (M - m0) % M
    conj = form.startswith("conj")
    real = "real-valued " if conj else ""
    text = (rf"Assume $x[n]$ is a {real}finite-duration sequence of length ${L}$, and $y[n]$ is obtained by zero-padding $x[n]$ to "
            rf"length ${M}$. Let $\{{X[k]\}}_{{k=0}}^{{{L - 1}}}$ and $\{{Y[k]\}}_{{k=0}}^{{{M - 1}}}$ be the DFTs of $x$ and $y$. Then "
            rf"$X[{k0}] = Y{'^*' if conj else ''}[{m0}]$.")
    wx, wy = Fr(2 * k0, L), Fr(2 * m0, M)
    why = (rf"Both DFTs evaluate the same DTFT: $X[{k0}] = X_d\left({pi_tex(wx)}\right)$ and $Y[{m0}] = X_d\left({pi_tex(wy)}\right)$")
    if conj:
        ok = (wx + wy) % 2 == 0
        why += (rf"; for real $x$, $Y^*[{m0}] = X_d\left({pi_tex(-wy)}\right) = X_d\left({pi_tex((-wy) % 2)}\right)$"
                + (": the same frequency." if ok else ": a different frequency."))
    else:
        ok = (wx - wy) % 2 == 0
        why += ": the same frequency." if ok else ": different frequencies."
    return "FA2021 #1(d), SP2025 #7(c)", text, ok, why, {"L": L, "M": M, "k0": k0, "m0": m0, "conj": conj}


# ============================================================================ G4: which transform (SP2021 #6), H_d vs H(z)
WHICH = {
    "dft_cont": (r"The DFT $X[k]$ of a finite-length sequence is a function of a continuous frequency variable.", False,
                 r"The DFT is $N$ numbers, one per integer $k$ ($\omega_k = \frac{2\pi k}{N}$); the <b>DTFT</b> is the function of the "
                 r"continuous $\omega$."),
    "dtft_periodic": (r"The DTFT of a finite-length sequence is periodic in $\omega$ with period $2\pi$.", True,
                      r"$e^{-j(\omega + 2\pi)n} = e^{-j\omega n}$ for every integer $n$, so every DTFT is $2\pi$-periodic."),
    "dft_computer": (r"The DFT can be applied to recorded data of finite length using a computer.", True,
                     r"$N$ recorded samples in, $N$ numbers out, each a finite sum (SP2021 #6(b): DFT)."),
    "ctft_periodic": (r"The CTFT $X_a(\Omega)$ of a continuous-time signal is periodic in $\Omega$.", False,
                      r"The CTFT is defined for every real $\Omega$ and is typically not periodic (e.g. $e^{-t}u(t) \leftrightarrow "
                      r"\frac{1}{1 + j\Omega}$); periodicity in frequency comes from sampling."),
    "dtft_from_z": (r"If $x[n]$ is absolutely summable, its DTFT can be computed from its z-transform by setting $z = e^{j\omega}$.", True,
                    r"Absolute summability means the ROC contains the unit circle, so $X_d(\omega) = X(z)\big|_{z = e^{j\omega}}$ "
                    r"(SP2021 #6(a): DTFT)."),
    "dtft_only_N": (r"The DTFT of a length-$N$ sequence is defined only at the $N$ frequencies $\omega_k = \frac{2\pi k}{N}$.", False,
                    r"The DTFT is defined for every real $\omega$; the DFT evaluates it at the $N$ frequencies $\frac{2\pi k}{N}$."),
    "dft_periodic_k": (r"The DFT formula $X[k] = \sum_{n=0}^{N-1}x[n]e^{-j2\pi kn/N}$, evaluated at every integer $k$, is periodic "
                       r"in $k$ with period $N$.", True,
                       r"$e^{-j2\pi(k + N)n/N} = e^{-j2\pi kn/N}e^{-j2\pi n} = e^{-j2\pi kn/N}$, so $X[k + N] = X[k]$."),
}


def which_transform():
    sid = random.choice(sorted(WHICH))
    text, ans, why = WHICH[sid]
    return "SP2021 #6 (CTFT, DTFT or DFT?)", text, ans, why, {"id": sid}


def hd_substitute():
    form = random.choice(["sys", "sys", "sys", "fa21", "every", "stable"])
    if form == "sys":
        a = random.choice([Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(3, 4), Fr(-3, 4), Fr(2), Fr(-2), Fr(1), Fr(-1), Fr(3, 2)])
        den = fmt.tex_sum([(1, ""), (-a, r"e^{-j\omega}")])
        sysx = fmt.tex_sum([(a, "y[n-1]"), (1, "x[n]")])
        text = (rf"The causal LTI system $y[n] = {sysx}$ has frequency response "
                rf"$H_d(\omega) = H(z)\big|_{{z = e^{{j\omega}}}} = \dfrac{{1}}{{{den}}}$.")
        ans = abs(a) < 1
        hn = rf"{fmt.tex_pow(a, 'n')}u[n]" if a != 1 else "u[n]"
        if ans:
            why = (rf"$h[n] = {hn}$ and the ROC $\lvert z\rvert \gt {T(abs(a))}$ contains the unit circle, so $H_d$ is $H(z)$ on "
                   r"$\lvert z\rvert = 1$.")
        elif abs(a) == 1:
            why = (rf"$h[n] = {hn}$: the pole $z = {T(a)}$ is on the unit circle, the ROC $\lvert z\rvert \gt 1$ misses it, and the "
                   r"DTFT exists only with an impulse (FA2021 #2: $(-1)^n u[n] \leftrightarrow \frac{1}{1 + e^{-j\omega}} + \pi\delta(\omega - \pi)$ "
                   r"on one period), which the formula misses.")
        else:
            why = (rf"$h[n] = {hn}$ grows: the ROC $\lvert z\rvert \gt {T(abs(a))}$ misses the unit circle, the system is unstable and "
                   r"has no frequency response; substituting $z = e^{j\omega}$ outside the ROC means nothing.")
        prm = {"form": form, "a": pj(a)}
    elif form == "fa21":
        text = (r"The causal LTI system $y[n] = y[n-2] + x[n] - x[n-1]$ has frequency response "
                r"$H_d(\omega) = H(z)\big|_{z = e^{j\omega}} = \dfrac{1}{1 + e^{-j\omega}}$.")
        ans = False
        why = (r"$H(z) = \frac{1 - z^{-1}}{1 - z^{-2}} = \frac{1}{1 + z^{-1}}$, $\lvert z\rvert \gt 1$, so $h[n] = (-1)^n u[n]$: the ROC "
               r"misses the unit circle and $H_d(\omega) = \frac{1}{1 + e^{-j\omega}} + \pi\delta(\omega - \pi)$ on $(-\pi, \pi]$ (FA2021 #2).")
        prm = {"form": form}
    elif form == "every":
        text = r"For every LTI system with transfer function $H(z)$, the frequency response is $H_d(\omega) = H(z)\big|_{z = e^{j\omega}}$."
        ans = False
        why = (r"Only when the ROC contains the unit circle. Counterexample: $h[n] = 2^n u[n]$, $H(z) = \frac{1}{1 - 2z^{-1}}$, "
               r"$\lvert z\rvert \gt 2$: $h$ is not summable and there is no $H_d$.")
        prm = {"form": form}
    else:
        text = (r"For every BIBO-stable LTI system with transfer function $H(z)$, the frequency response is "
                r"$H_d(\omega) = H(z)\big|_{z = e^{j\omega}}$.")
        ans = True
        why = (r"Stable means $\sum_n\lvert h[n]\rvert \lt \infty$, i.e. the ROC contains $\lvert z\rvert = 1$, where $H(z)$ and the "
               r"DTFT sum coincide.")
        prm = {"form": form}
    return "FA2021 #2(c)", text, ans, why, prm


# ============================================================================ G5: properties
def half_shift():
    N = random.choice([4, 6, 8])
    h = N // 2
    form = random.choice(["alt_shift", "alt_mod", "kshift_mod", "kshift_shift"])
    Ytex = r"(-1)^k X[k]" if form.startswith("alt") else rf"X[\langle k - {h}\rangle_{{{N}}}]"
    ytex = (rf"x[\langle n - {h}\rangle_{{{N}}}]" if form in ("alt_shift", "kshift_shift") else r"(-1)^n x[n]")
    text = (rf"Let $x[n]$ be a length-${N}$ signal with ${N}$-point DFT $X[k]$. The length-${N}$ signal whose DFT is ${Ytex}$ is "
            rf"$y[n] = {ytex}$.")
    ans = form in ("alt_shift", "kshift_mod")
    if form.startswith("alt"):
        why = (rf"$(-1)^k = e^{{-j\pi k}} = e^{{-j\frac{{2\pi}}{{{N}}}k\cdot {h}}}$: a linear phase, i.e. a circular shift by $\frac{{N}}{{2}} = {h}$: "
               rf"$y[n] = x[\langle n - {h}\rangle_{{{N}}}]$" + ("." if ans else r", not a sign alternation of $x$ (SP2023 #9(b))."))
    else:
        why = (rf"Frequency shift by $\ell = {h}$: $e^{{j\frac{{2\pi}}{{{N}}}{h}n}}x[n] = (-1)^n x[n] \leftrightarrow X[\langle k - {h}\rangle_{{{N}}}]$"
               + ("." if ans else r"; a shift in $k$ modulates $x$, it does not shift it."))
    return "SP2021 #4, SP2023 #9(b), HW9 #2", text, ans, why, {"N": N, "form": form}


def phase_ramp():
    if random.random() < 0.4:
        x = list(range(1, 9))
    else:
        x = random.sample([1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3], 8)
    m = random.choice([1, 2, 3, 5, 6, 7])
    adv = random.random() < 0.5                       # e^{+j..}: an advance by m
    t = Fr(2 * m, 8)
    form = random.choice(["ok", "ok", "dir", "linear"])
    if form == "linear" and m > 3:
        form = "dir"                                  # a long linear shift would leave an obviously wrong list of zeros
    shift = -m if adv else m                          # circular delay
    ok = [x[(n - shift) % 8] for n in range(8)]
    claim = {"ok": ok, "dir": [x[(n + shift) % 8] for n in range(8)],
             "linear": [x[n - shift] if 0 <= n - shift < 8 else 0 for n in range(8)]}[form]
    if claim == ok:
        form = "ok"
    ph = rf"e^{{{'' if adv else '-'}j{pin(t, 'k')}}}"
    text = (rf"Let $x[n] = {seq(x)}$, $0 \le n \le 7$, with 8-point DFT $X[k]$. The sequence whose DFT is $Y[k] = {ph}X[k]$ is "
            rf"$y[n] = {seq(claim)}$.")
    word = rf"an <b>advance</b> by ${m}$" if adv else rf"a circular delay by ${m}$"
    sh_tex = "" if shift == 1 else (rf"\cdot {shift}" if shift > 0 else rf"\cdot({shift})")
    why = (rf"${ph} = e^{{-j\frac{{2\pi}}{{8}}k{sh_tex}}}$: {word}, $y[n] = x[\langle n {'+' if adv else '-'} {m}\rangle_8] = {seq(ok)}$")
    why += {"ok": ".", "dir": ": the stated list is the shift in the other direction.",
            "linear": ": the stated list is a linear shift, which drops samples instead of wrapping them around."}[form]
    return "FA2019 #4, SP2025 #8", text, form == "ok", why, {"x": x, "m": m, "adv": adv, "claim": claim}


# ============================================================================ G6: DFT of a sinusoid
WORD = {1: "one", 2: "two"}


def cos_bins():
    kind = random.choice(["digital", "digital", "analog"])
    if kind == "digital":
        N = random.choice([16, 20, 24, 32, 64])
        c = random.choice([Fr(1), Fr(2), Fr(3), Fr(3, 2), Fr(5, 2), Fr(N // 2)])
        fs = f0 = None
    else:
        fs = random.choice([40, 50, 64, 100])
        N = fs
        f0 = random.choice([Fr(5), Fr(10), Fr(fs // 2), Fr(15, 2), Fr(25, 2)])
        c = f0
    w0 = 2 * c / N
    count = (1 if c == Fr(N, 2) else 2) if c.denominator == 1 else N
    m = count if (count <= 2 and random.random() < 0.6) else (random.choice([1, 2]) if count > 2 else 3 - count)
    ans = m == count
    s = "s" if m > 1 else ""
    if kind == "digital":
        text = (rf"The ${N}$-point DFT of the finite-length signal $x[n] = \cos\left({pin(w0)}\right)$, $0 \le n \le {N - 1}$, "
                rf"has only {WORD[m]} nonzero element{s}.")
    else:
        text = (rf"An analog signal $x_a(t) = \cos({T(2 * f0)}\pi t)$ is sampled at $f_s = {fs}$ Hz for 1 second to obtain "
                rf"$\{{x[n]\}}_{{n=0}}^{{{N - 1}}}$. Its DFT $X[k]$ is nonzero for only {WORD[m]} value{s} of $k$, $0 \le k \le {N - 1}$.")
    if count == 2:
        why = (rf"$\omega_0 = {pi_tex(w0)} = \frac{{2\pi\cdot{c}}}{{{N}}}$: exactly ${c}$ periods in ${N}$ samples, so the two exponentials of "
               rf"the cosine sit on bins $k = {c}$ and $k = {N - c}$: $X[{c}] = X[{N - c}] = {N // 2}$, every other bin is $0$.")
    elif count == 1:
        why = (rf"$x[n] = \cos(\pi n) = (-1)^n = e^{{j\frac{{2\pi\cdot{N // 2}}}{{{N}}}n}}$: one exponential exactly on bin ${N // 2}$, "
               rf"so $X[{N // 2}] = {N}$ and all other bins are $0$.")
    else:
        why = (rf"$\omega_0 = {pi_tex(w0)}$ and $\frac{{{N}\omega_0}}{{2\pi}} = {T(c)}$ is not an integer: the record holds ${T(c)}$ periods, "
               rf"the cosine falls between bins and leaks into all ${N}$ of them.")
    return "FA2019 #1(f), FA2023 #1(g), SP2025 #1(e)", text, ans, why, {"kind": kind, "N": N, "c": pj(c), "m": m}


def peak_height():
    N = random.choice([16, 32, 64])
    m = random.randint(1, N // 2 - 1)
    A = random.choice([Fr(1), Fr(2), Fr(3), Fr(1, 2), Fr(3, 4), Fr(4)])
    trig = random.choice(["cos", "cos", "sin"])
    form = random.choice(["ok", "ok", "ok", "full", "half", "quarter"])
    V = {"ok": A * N / 2, "full": A * N, "half": A / 2, "quarter": A * N / 4}[form]
    At = "" if A == 1 else T(A)
    text = (rf"The ${N}$-point DFT of $x[n] = {At}\{trig}\left({pin(Fr(2 * m, N))}\right)$, $0 \le n \le {N - 1}$, satisfies "
            rf"$\lvert X[{m}]\rvert = {T(V)}$.")
    why = (rf"$\omega_0 = {pi_tex(Fr(2 * m, N))} = \frac{{2\pi\cdot {m}}}{{{N}}}$ is a bin frequency. By Euler the {'cosine' if trig == 'cos' else 'sine'} "
           r"is two exponentials of magnitude $\frac{A}{2}$ at $\pm\omega_0$; the one at $+\omega_0$ adds up to $\frac{A}{2}\cdot N$ in bin "
           rf"${m}$ and the other one contributes $0$ there, so $\lvert X[{m}]\rvert = \frac{{AN}}{{2}} = {T(A * N / 2)}$"
           + ("." if form == "ok" else rf", not ${T(V)}$."))
    return "FA2025 #8(b), HW9 #4; Lecture 22", text, form == "ok", why, {"N": N, "m": m, "A": pj(A), "trig": trig, "V": pj(V)}


def peak_freq():
    while True:
        N = random.choice([32, 64, 100, 128, 200, 256])
        df = random.choice([2, 4, 5, 10, 25])
        fs = N * df
        if fs <= 10000:
            break
    k0 = random.randint(3, N // 2 - 3)
    form = random.choice(["ok", "ok", "rad", "half", "mirror"])
    f0 = k0 * df
    claim = {"ok": rf"$f_0 = {f0}$ Hz", "rad": rf"$\Omega_0 = {2 * f0}\pi$ rad/s", "half": rf"$f_0 = {T(Fr(f0, 2))}$ Hz",
             "mirror": rf"$f_0 = {(N - k0) * df}$ Hz"}[form]
    text = (rf"A tone $x_a(t) = \cos(2\pi f_0 t)$ with $0 \lt f_0 \lt f_s/2$ is sampled at $f_s = {fs}$ Hz and the ${N}$-point DFT of "
            rf"${N}$ samples is computed. $\lvert X[k]\rvert$ is nonzero only at $k = {k0}$ and $k = {N - k0}$. Then {claim}.")
    why = (rf"Bin $k$ sits at $\omega_k = \frac{{2\pi k}}{{N}}$, i.e. at $f = \frac{{k f_s}}{{N}} = k\cdot{df}$ Hz; bin ${k0}$ gives "
           rf"$f_0 = {f0}$ Hz ($\Omega_0 = 2\pi f_0 = {2 * f0}\pi$ rad/s), and bin ${N - k0}$ is its negative-frequency mirror")
    why += {"ok": ".", "rad": ".", "half": r" (the factor $\frac12$ comes from mixing up $\frac{\omega}{2\pi}$ and $\frac{\omega}{\pi}$).",
            "mirror": rf", at $-{f0}$ Hz, not ${(N - k0) * df}$ Hz."}[form]
    return "HW9 #5, FA2025 #8(b)", text, form in ("ok", "rad"), why, {"N": N, "fs": fs, "k0": k0, "form": form}


# ============================================================================ G7: resolution
DELTA = {32: [Fr(2, 25), Fr(1, 10)], 64: [Fr(1, 25), Fr(1, 20)], 128: [Fr(1, 50), Fr(1, 40)]}


def zp_resolution():
    N = random.choice([32, 64, 128])
    d = random.choice(DELTA[N])                       # 2pi/N < dw < 4pi/N: main lobes overlap
    w1 = random.choice([Fr(3, 10), Fr(2, 5), Fr(1, 2), Fr(3, 5)])
    c = random.choice([2, 4, 8])
    form = random.choice(["pad", "pad", "more", "grid"])
    sig = (rf"Two equal-amplitude sinusoids at $\omega_1 = {dec(w1)}\pi$ and $\omega_2 = {dec(w1 + d)}\pi$ are recorded for $N = {N}$ "
           r"samples (rectangular window), so their main lobes overlap.")
    claim = {"pad": rf"Zero padding the record to length ${c * N}$ before the DFT makes the two main lobes full-lobe separated.",
             "more": rf"Recording ${c * N}$ samples instead (same sampling rate) makes the two main lobes full-lobe separated.",
             "grid": (rf"Zero padding the record to length ${c * N}$ evaluates the same DTFT at ${c}$ times as many frequencies; each "
                      rf"main lobe keeps the width $\frac{{4\pi}}{{{N}}}$.")}[form]
    text = sig + " " + claim
    why = (rf"Full-lobe separation needs $\Delta\omega \gt \frac{{4\pi}}{{N}}$ for the number $N$ of <b>recorded</b> samples; here "
           rf"$\Delta\omega = {dec(d)}\pi$ and $\frac{{4\pi}}{{{N}}} {eqdec(Fr(4, N))}\pi$.")
    if form == "pad":
        why += " Zero padding does not change the DTFT, so the lobes keep their width and still overlap: it refines the grid, not the resolution."
    elif form == "more":
        why += rf" With ${c * N}$ samples the threshold is $\frac{{4\pi}}{{{c * N}}} {eqdec(Fr(4, c * N))}\pi \lt {dec(d)}\pi$."
    else:
        why += r" Padding adds zero terms only: $Y_d(\omega) = X_d(\omega)$, sampled more densely."
    return "Lectures 22-23; FA2025 #8(a)", text, form in ("more", "grid"), why, {"N": N, "d": pj(d), "c": c, "form": form}


def full_lobe():
    N = random.choice([16, 20, 32, 50, 64, 100, 128])
    kind = random.choice(["full", "full", "half"])
    thr = Fr(4 if kind == "full" else 2, N)            # in units of pi
    cands = [Fr(v, 1000) for v in range(5, 400, 5)]
    above = random.random() < 0.5
    cands = [d for d in cands if (thr < d < 2 * thr) if above] or [d for d in cands if Fr(1, 2) * thr < d < thr]
    d = random.choice(cands)
    w1 = random.choice([Fr(1, 5), Fr(3, 10), Fr(2, 5), Fr(1, 2)])
    crit = "full-lobe separated (their main lobes do not overlap)" if kind == "full" else \
        "half-lobe separated (at most half of each main lobe overlaps)"
    text = (rf"Two equal-amplitude sinusoids at $\omega_1 = {dec(w1)}\pi$ and $\omega_2 = {dec(w1 + d)}\pi$ observed for $N = {N}$ "
            rf"samples with a rectangular window are {crit}.")
    ans = d > thr
    why = (rf"The rectangular window's main lobe has width $\frac{{4\pi}}{{N}}$. {'Full' if kind == 'full' else 'Half'}-lobe separation "
           rf"needs $\Delta\omega \gt \frac{{{4 if kind == 'full' else 2}\pi}}{{{N}}} {eqdec(thr)}\pi$; here $\Delta\omega = {dec(d)}\pi$, which is "
           + ("larger." if ans else "smaller."))
    return "Lecture 23, separating main lobes", text, ans, why, {"N": N, "d": pj(d), "kind": kind}


def record_time():
    form = random.choice(["count", "count", "time", "time", "double"])
    fs = random.choice([4000, 8000, 16000])
    f1 = random.choice([82, 98, 110, 131, 220])
    df = random.choice([4, 5, 8, 10])
    f2 = f1 + df
    if form == "count":
        sub = random.choice(["ok", "ok", "half", "double"])
        V = {"ok": Fr(2 * fs, df), "half": Fr(fs, df), "double": Fr(4 * fs, df)}[sub]
        text = (rf"Two tones of similar amplitude at $f_1 = {f1}$ Hz and $f_2 = {f2}$ Hz are recorded at $f_s = {fs}$ Hz. With a "
                rf"rectangular window, full-lobe separation of their peaks needs more than $N = {V}$ samples.")
        why = (rf"$\Delta\omega = \frac{{2\pi\cdot{df}}}{{{fs}}}$ and full-lobe separation needs $\Delta\omega \gt \frac{{4\pi}}{{N}}$, i.e. "
               rf"$N \gt \frac{{2f_s}}{{\Delta f}} = \frac{{2\cdot{fs}}}{{{df}}} = {2 * fs // df}$ samples ($\frac{{2}}{{\Delta f}} = {dec(Fr(2, df))}$ s; "
               r"Lecture 23's guitar example)" + ("." if sub == "ok" else rf", not ${V}$."))
        ans = sub == "ok"
        prm = {"form": form, "fs": fs, "df": df, "V": int(V)}
    elif form == "double":
        text = (rf"To full-lobe separate two tones $\Delta f = {df}$ Hz apart (rectangular window), doubling the sampling rate from "
                rf"${fs}$ Hz to ${2 * fs}$ Hz halves the recording time needed.")
        why = (rf"$N \gt \frac{{2f_s}}{{\Delta f}}$ samples take $\frac{{N}}{{f_s}} \gt \frac{{2}}{{\Delta f}} = {dec(Fr(2, df))}$ s at any $f_s$: "
               r"doubling $f_s$ doubles the samples needed (the tones are half as far apart in $\omega$) but not the time.")
        ans = False
        prm = {"form": form, "fs": fs, "df": df}
    else:
        sub = random.choice(["ok", "half"])
        tsec = Fr(2, df) if sub == "ok" else Fr(1, df)
        text = (rf"To full-lobe separate two tones $\Delta f = {df}$ Hz apart (rectangular window), the recording must last longer than "
                rf"${dec(tsec)}$ s, whatever the sampling rate.")
        why = (rf"$\Delta\omega = \frac{{2\pi\Delta f}}{{f_s}} \gt \frac{{4\pi}}{{N}}$ gives $N \gt \frac{{2f_s}}{{\Delta f}}$, a duration "
               rf"$\frac{{N}}{{f_s}} \gt \frac{{2}}{{\Delta f}} = {dec(Fr(2, df))}$ s, independent of $f_s$"
               + ("." if sub == "ok" else rf"; ${dec(tsec)}$ s is the half-lobe figure."))
        ans = sub == "ok"
        prm = {"form": form, "df": df, "t": pj(tsec)}
    return "Lecture 23, Exercise 1", text, ans, why, prm


def fa25_plots():
    N = random.choice([16, 32])
    form = random.choice(["zp_same", "more_taller", "zp_taller", "ham_narrow"])
    lead = (rf"A sum of two sinusoids is sampled; consider the DFT of $N = {N}$ samples. ")
    claim = {
        "zp_same": (rf"Zero padding the ${N}$ samples to length ${2 * N}$ gives a DFT whose main lobes have the same heights and widths "
                    rf"as in the ${N}$-point DFT, evaluated at twice as many frequencies."),
        "more_taller": (rf"Taking ${2 * N}$ samples instead (same sampling rate) gives main lobes about twice as tall and half as wide."),
        "zp_taller": (rf"Zero padding the ${N}$ samples to length ${2 * N}$ gives main lobes about twice as tall and half as wide, just like "
                      rf"taking ${2 * N}$ samples."),
        "ham_narrow": (rf"Applying a Hamming window to ${2 * N}$ samples gives narrower main lobes than the rectangular window on the same "
                       rf"${2 * N}$ samples."),
    }[form]
    text = lead + claim
    why = {
        "zp_same": r"Padding leaves the DTFT unchanged (height $\approx\frac{AN}{2}$, width $\frac{4\pi}{N}$); only the grid gets finer (FA2025 #8: DFT II vs DFT I).",
        "more_taller": rf"With $2N$ samples each peak is $\approx\frac{{A\cdot 2N}}{{2}}$ high and the main lobe $\frac{{4\pi}}{{2N}}$ wide (FA2025 #8: DFT III).",
        "zp_taller": r"Zero padding does not change the DTFT: heights and widths stay those of $N$ samples; only more data narrows the lobes.",
        "ham_narrow": r"The Hamming window lowers the side lobes but doubles the main-lobe width ($\approx\frac{8\pi}{N}$ instead of $\frac{4\pi}{N}$) (FA2025 #8: DFT IV).",
    }[form]
    return "FA2025 #8(a)", text, form in ("zp_same", "more_taller"), why, {"N": N, "form": form}


# ============================================================================ G8: windows
def hamming():
    good = random.random() < 0.5
    if good:
        side, main = "lower", "wider"
    else:
        side, main = random.choice([("lower", "narrower"), ("lower", "narrower"), ("higher", "wider"), ("higher", "narrower")])
    text = (rf"Applying the Hamming window for spectral analysis (instead of the rectangular window) leads to {side} sidelobes and "
            rf"{main} main lobes in the spectral components.")
    why = ("It is a trade-off: the Hamming window's sidelobes are about $-42$ dB instead of $-13$ dB, but its main lobe is about twice "
           r"as wide ($\frac{8\pi}{N}$ instead of $\frac{4\pi}{N}$)")
    if good:
        why += "."
    else:
        wrong = [w for w, ok in (("the sidelobes", side == "lower"), ("the main lobes", main == "wider")) if not ok]
        why += f"; the statement gets {' and '.join(wrong)} wrong."
    return "SP2025 #1(f)", text, good, why, {"side": side, "main": main}


def window_lti():
    form = random.choice(["linear", "lti", "ti", "notti"])
    claim = {"linear": "is a linear system", "lti": "is a linear time-invariant (LTI) system",
             "ti": "is a time-invariant system", "notti": "is linear but not time-invariant"}[form]
    text = (rf"Windowing a signal, $y[n] = w[n]\,x[n]$ with a fixed Hamming window $w[n]$ of length $N$, {claim}.")
    why = (r"Linear: $w[n](a x_1[n] + b x_2[n]) = a y_1[n] + b y_2[n]$. Not time-invariant: delaying the input by one sample gives "
           r"$w[n]x[n-1]$, not $y[n-1] = w[n-1]x[n-1]$ (the window stays put). Windowing is multiplication, not convolution "
           "(Lecture 23).")
    return "Lecture 23, window functions", text, form in ("linear", "notti"), why, {"form": form}


def window_product():
    form = random.choice(["conv", "product"])
    claim = (r"$Y_d(\omega) = \dfrac{1}{2\pi}\displaystyle\int_{-\pi}^{\pi}X_d(\theta)\,W_d(\omega - \theta)\,d\theta$" if form == "conv"
             else r"$Y_d(\omega) = X_d(\omega)\,W_d(\omega)$")
    text = rf"If $y[n] = x[n]\,w[n]$ (windowing), then {claim}."
    why = (r"Windowing property: a product in time is a periodic convolution of the DTFTs, divided by $2\pi$. That is why every sinusoid "
           r"shows up as a copy of $W_d$ (a main lobe with side lobes) centred at its frequency. A product of DTFTs would be a "
           r"convolution in time.")
    return "Lecture 23, eq. (18)", text, form == "conv", why, {"form": form}


def mainlobe_width():
    N = random.choice([8, 16, 20, 32, 64])
    form = random.choice(["rect_ok", "rect_ok", "rect_half", "rect_height", "ham_ok", "ham_same"])
    if form.startswith("rect"):
        wdt = Fr(4, N) if form != "rect_half" else Fr(2, N)
        hgt = N if form != "rect_height" else Fr(N, 2)
        text = (rf"The DTFT of the length-${N}$ rectangular window $w[n] = 1$, $0 \le n \le {N - 1}$, has a main lobe of width "
                rf"${pi_tex(wdt)}$ (between its first zeros) and height ${T(hgt)}$ at $\omega = 0$.")
        why = (rf"$\lvert W_d(\omega)\rvert = \left\lvert\dfrac{{\sin({N}\omega/2)}}{{\sin(\omega/2)}}\right\rvert$: first zeros at "
               rf"$\pm\frac{{2\pi}}{{{N}}}$, so the width is ${pi_tex(Fr(4, N))}$, and $W_d(0) = \sum_n 1 = {N}$.")
        ans = form == "rect_ok"
    else:
        N = random.choice([16, 20, 32, 64])
        wdt = Fr(8, N) if form == "ham_ok" else Fr(4, N)
        text = (rf"The DTFT of the length-${N}$ Hamming window has a main lobe of width about ${pi_tex(wdt)}$ (between its first zeros).")
        why = rf"The Hamming main lobe is about twice as wide as the rectangular one: about $\frac{{8\pi}}{{N}} = {pi_tex(Fr(8, N))}$ (Lecture 23)."
        ans = form == "ham_ok"
    return "Lectures 22-23", text, ans, why, {"N": N, "form": form}


# ============================================================================ assembly
GROUPS = {
    "values": ["x0_zero", "real_bins", "fa25_dtft", "parseval_dft"],
    "inverse": ["idft_periodic", "idft_periodic", "dft_same", "dup_rule"],
    "dtft": ["bin_freq", "zp_scale", "zp_index"],
    "transforms": ["which_transform", "hd_substitute"],
    "props": ["half_shift", "phase_ramp"],
    "sinusoid": ["cos_bins", "peak_height", "peak_freq"],
    "resolution": ["zp_resolution", "full_lobe", "record_time", "fa25_plots"],
    "window": ["hamming", "window_lti", "window_product", "mainlobe_width"],
}
FUNCS = {f.__name__: f for f in (x0_zero, real_bins, fa25_dtft, parseval_dft, idft_periodic, dft_same, dup_rule, bin_freq,
                                 zp_scale, zp_index, which_transform, hd_substitute, half_shift, phase_ramp, cos_bins,
                                 peak_height, peak_freq, zp_resolution, full_lobe, record_time, fa25_plots, hamming,
                                 window_lti, window_product, mainlobe_width)}


def generate(data):
    while True:
        groups = random.sample(list(GROUPS), 5)
        items = []
        for i, g in enumerate(groups):
            sid = random.choice(GROUPS[g])
            src, text, ans, why, prm = FUNCS[sid]()
            items.append({"name": f"tf{i + 1}", "label": "abcde"[i], "sid": sid, "stmt_html": text, "src": src,
                          "param": prm, "ans": "True" if ans else "False", "why_html": why, "choices": mc.true_false(ans)})
        n_true = sum(it["ans"] == "True" for it in items)
        if 1 <= n_true <= 4:
            break
    data["params"]["statements"] = items
