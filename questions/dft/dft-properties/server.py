"""DFT properties: find y[n] from a DFT Y[k] defined by one or two operations on X[k] (Midterm 2 family
"DFT computation and properties", 7 of 7 past exams), plus a conceptual question on the conjugate
symmetry of a real signal's DFT.

Models: FA2019 MT2 #4 (Y = e^{-j pi k/4} X), SP2021 #4 (Y = e^{-j(6 pi k/8 + pi)} X), SP2025 #8
(Y = e^{j 3 pi k/4} X: an advance), FA2024 #8 (circular reversal + shift), FA2021 #3(b)
(x[<3 - n>_8] j^n), SP2025 #7(b) and SP2023 #9(b) (a permuted, sign-alternated list of the X_k),
FA2021 #6 (products of DFTs are circular convolutions), SP2021 #5 (is x real?). Properties from the
2019 course notes, Sec. 2.5.1 (Table 2.5).

x: length N in {6, 8}, small integers, no shorter period, not circularly even. Families:
  shift    Y = e^{-j 2 pi k m/N} X (delay m, or written as an advance), optionally times e^{-j pi} = -1
  reverse  Y = X[<-k>_N]  or  Y = X^*[k]                         -> y = x[<-n>_N]
  revshift Y = e^{-j 2 pi k m/N} X[<-k>_N]                        -> y = x[<m - n>_N]
  half     Y = (-1)^k X  or  Y = X[<k - N/2>_N]                   -> x[<n - N/2>_N] or (-1)^n x[n]
  list     Y = s sigma^k X[<k - m>_N] shown as a list of +-X_i      -> s (-1)^m e^{j 2 pi m n/N} x[<n - N/2>]
           (sigma = -1) or e^{j 2 pi m n/N} x[n] (sigma = 1); m with e^{j 2 pi m/N} in {+-1, +-j}, so
           y has integer real and imaginary parts
  square   Y = X^2                                                -> circular self-convolution
  evenodd  Y = Re{X} or j Im{X}                                   -> (x[n] +- x[<-n>_N]) / 2 (integers)
(b) four candidate DFT lists, exactly one Hermitian (X[0], X[N/2] real, X[k] = X^*[<-k>_N]); the
others break it in classic ways (X[0] complex, equal instead of conjugate pairs, pairs matched by a
plain flip, X[N/2] complex / a negated pair)."""

import random
from collections import Counter
from fractions import Fraction as Fr

import numpy as np
import prairielearn as pl
from ece310 import fmt, w8_dft as W

FAMILIES = ["shift", "shift", "shift", "reverse", "revshift", "revshift", "half", "list", "list", "listf", "listf",
            "square", "square", "evenodd"]
PHASE = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}      # j^q as (re, im)


# ----------------------------------------------------------------------------- helpers
def _frac_pi(t, var="k"):
    """t*pi*var inside an exponent: \\frac{3\\pi}{4}k, \\frac{\\pi}{2}k, \\pi k."""
    t = Fr(t)
    p, q = abs(t.numerator), t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return (num if q == 1 else r"\frac{" + num + "}{" + str(q) + "}") + (" " if q == 1 else "") + var


def _idx_list_tex(N, idx, coef_tex=None):
    """{x[i0], x[i1], ...} with optional coefficient LaTeX per entry."""
    items = []
    for n, i in enumerate(idx):
        c = coef_tex[n] if coef_tex else ""
        items.append(f"{c}x[{i}]")
    return W.tex_list(items)


def _cvals_tex(vals):
    return W.tex_list([W.tex_cplx(re, im) for re, im in vals])


def _rand_x(N, lo=-3, hi=5):
    while True:
        x = [random.randint(lo, hi) for _ in range(N)]
        if sum(1 for v in x if v == 0) > 2 or len(set(x)) < N - 2 or x[0] == 0 or x[-1] == 0:
            continue
        shifts = {tuple(x[(n - m) % N] for n in range(N)) for m in range(N)}
        if len(shifts) < N:
            continue
        if [x[(-n) % N] for n in range(N)] == x or x[::-1] == x:
            continue
        if all(x[(n - N // 2) % N] == -x[n] for n in range(N)):
            continue                                      # anti-periodic: shift by N/2 = negation
        return x


def _circ(x, f):
    """[x[<f(n)>_N] for n]"""
    N = len(x)
    return [x[f(n) % N] for n in range(N)]


# ----------------------------------------------------------------------------- families
def _fam_shift(p, N):
    x = _rand_x(N)
    m = random.choice([v for v in range(1, N) if 2 * v != N])
    neg = random.random() < 0.3
    advance = m > N // 2 and random.random() < 0.6
    t = Fr(2 * (N - m), N) if advance else Fr(2 * m, N)          # displayed exponent coefficient (times pi)
    sgn = "" if advance else "-"
    if neg:
        expo = sgn + r"j\left(" + _frac_pi(t) + r" + \pi\right)"
    else:
        expo = sgn + "j" + _frac_pi(t)
    p["Y_def_tex"] = rf"Y[k] = e^{{{expo}}}\,X[k]"
    p["disp"] = {"t": [t.numerator, t.denominator], "advance": advance, "neg": neg}
    y = [(-1 if neg else 1) * v for v in _circ(x, lambda n: n - m)]
    idx = [(n - m) % N for n in range(N)]
    mm = N - m
    steps = []
    if advance:
        steps.append(rf"$e^{{j{_frac_pi(t)}}} = e^{{-j\frac{{2\pi}}{{{N}}}k\cdot(-{mm})}}$: a circular <b>advance</b> by ${mm}$, "
                     rf"i.e. a circular delay by $m = {N} - {mm} = {m}$.")
    else:
        steps.append(rf"$e^{{-j{_frac_pi(t)}}} = e^{{-j\frac{{2\pi}}{{{N}}}k\cdot {m}}}$: a circular delay by $m = {m}$, "
                     rf"$x[{W.tex_mod('n - m', N)}] \leftrightarrow e^{{-j\frac{{2\pi}}{{N}}km}}X[k]$.")
    if neg:
        steps.append(r"The constant $e^{\mp j\pi} = -1$ is not part of the shift: it negates every sample.")
    p["steps"] = [{"step_html": s} for s in steps]
    ytex = ("-" if neg else "") + f"x[{W.tex_mod(f'n - {m}', N)}]"
    p["y_formula_tex"] = rf"y[n] = {ytex}"
    p["y_index_tex"] = _idx_list_tex(N, idx, ["-" if neg else ""] * N)
    wrong = [(-1 if neg else 1) * v for v in _circ(x, lambda n: n + m)]
    lin = [(-1 if neg else 1) * (x[n - m] if n - m >= 0 else 0) for n in range(N)]
    return x, [(v, 0) for v in y], {"wrong_dir": wrong, "linear": lin}


def _fam_reverse(p, N):
    x = _rand_x(N)
    conj = random.random() < 0.5
    p["Y_def_tex"] = r"Y[k] = X^*[k]" if conj else rf"Y[k] = X[{W.tex_mod('-k', N)}]"
    y = _circ(x, lambda n: -n)
    if conj:
        steps = [rf"Conjugation property: $x^*[{W.tex_mod('-n', N)}] \leftrightarrow X^*[k]$ "
                 rf"(conjugating the definition turns $e^{{-j\frac{{2\pi}}{{N}}kn}}$ into $e^{{+j\frac{{2\pi}}{{N}}kn}}$, a reversal of $n$).",
                 r"$x$ is real, so $x^*[n] = x[n]$; the conjugate of a real signal's DFT is the DFT of the circularly reversed signal."]
    else:
        steps = [rf"Time reversal: $x[{W.tex_mod('-n', N)}] \leftrightarrow X[{W.tex_mod('-k', N)}]$."]
    steps.append(rf"Circular reversal keeps $x[0]$ in place: $x[{W.tex_mod('-n', N)}] = \{{x[0], x[{N - 1}], \dots, x[1]\}}$.")
    p["steps"] = [{"step_html": s} for s in steps]
    p["y_formula_tex"] = rf"y[n] = x[{W.tex_mod('-n', N)}]"
    p["y_index_tex"] = _idx_list_tex(N, [(-n) % N for n in range(N)])
    p["disp"] = {"conj": conj}
    return x, [(v, 0) for v in y], {"plain_flip": x[::-1], "identity": list(x)}


def _fam_revshift(p, N):
    x = _rand_x(N)
    m = random.choice([v for v in range(1, N) if 2 * v != N])
    t = Fr(2 * m, N)
    p["Y_def_tex"] = rf"Y[k] = e^{{-j{_frac_pi(t)}}}\,X[{W.tex_mod('-k', N)}]"
    p["disp"] = {"t": [t.numerator, t.denominator]}
    y = _circ(x, lambda n: m - n)
    steps = [rf"Reversal: $z[n] = x[{W.tex_mod('-n', N)}] \leftrightarrow Z[k] = X[{W.tex_mod('-k', N)}]$.",
             rf"$e^{{-j{_frac_pi(t)}}} = e^{{-j\frac{{2\pi}}{{{N}}}k\cdot {m}}}$ delays $z$ circularly by ${m}$: "
             rf"$y[n] = z[{W.tex_mod(f'n - {m}', N)}] = x[{W.tex_mod(f'-(n - {m})', N)}] = x[{W.tex_mod(f'{m} - n', N)}]$."]
    p["steps"] = [{"step_html": s} for s in steps]
    p["y_formula_tex"] = rf"y[n] = x[{W.tex_mod(f'{m} - n', N)}]"
    p["y_index_tex"] = _idx_list_tex(N, [(m - n) % N for n in range(N)])
    wrong = _circ(x, lambda n: -n - m)
    noflip = _circ(x, lambda n: n - m)
    return x, [(v, 0) for v in y], {"wrong_dir": wrong, "no_reversal": noflip}


def _fam_half(p, N):
    x = _rand_x(N)
    h = N // 2
    alt = random.random() < 0.5
    sh = _circ(x, lambda n: n - h)
    mo = [v * (-1) ** n for n, v in enumerate(x)]
    if alt:
        p["Y_def_tex"] = r"Y[k] = (-1)^k\,X[k]"
        y = sh
        steps = [rf"$(-1)^k = e^{{-j\pi k}} = e^{{-j\frac{{2\pi}}{{{N}}}k\cdot {h}}}$: a circular shift by $\frac{{N}}{{2}} = {h}$, not a sign change of $y$."]
        p["y_formula_tex"] = rf"y[n] = x[{W.tex_mod(f'n - {h}', N)}]"
        p["y_index_tex"] = _idx_list_tex(N, [(n - h) % N for n in range(N)])
        other = mo
    else:
        p["Y_def_tex"] = rf"Y[k] = X[{W.tex_mod(f'k - {h}', N)}]"
        y = mo
        steps = [rf"Frequency shift: $e^{{j\frac{{2\pi}}{{N}}\ell n}}x[n] \leftrightarrow X[{W.tex_mod(r'k - \ell', N)}]$ with $\ell = {h} = \frac{{N}}{{2}}$, "
                 rf"and $e^{{j\frac{{2\pi}}{{{N}}}{h}n}} = e^{{j\pi n}} = (-1)^n$."]
        p["y_formula_tex"] = r"y[n] = (-1)^n\,x[n]"
        p["y_index_tex"] = _idx_list_tex(N, list(range(N)), ["" if n % 2 == 0 else "-" for n in range(N)])
        other = sh
    p["steps"] = [{"step_html": s} for s in steps]
    p["disp"] = {"alt": alt}
    return x, [(v, 0) for v in y], {"swapped_property": other, "negated": [-v for v in x]}


def _fam_list(p, N):
    x = _rand_x(N)
    h = N // 2
    if N == 6:
        sigma, m = -1, 3
    else:
        sigma = random.choice([-1, -1, 1])
        m = random.choice([2, 4, 6] if sigma == -1 else [2, 6])
    s = random.choice([1, -1]) if sigma == -1 else 1
    q = (4 * m // N) % 4                       # e^{j 2 pi m n / N} = j^{q n}
    idx = [(k - m) % N for k in range(N)]
    sgn = [s * (sigma ** k) for k in range(N)]
    p["list_idx"], p["list_sgn"] = idx, sgn
    items = [("-" if g < 0 else "") + f"X_{{{i}}}" for i, g in zip(idx, sgn)]
    p["Y_def_tex"] = r"\{Y[k]\}_{k=0}^{" + str(N - 1) + "} = " + W.tex_list(items)
    pre = (s * (-1) ** m) if sigma == -1 else 1
    base = _circ(x, lambda n: n - h) if sigma == -1 else list(x)
    y = []
    for n, v in enumerate(base):
        re, im = PHASE[(q * n) % 4]
        y.append((pre * re * v, pre * im * v))
    phase_tex = {0: "", 1: "j^n", 2: "(-1)^n", 3: "(-j)^n"}[q]
    ell_tex = W.tex_mod(f"k - {m}", N)
    steps = []
    if sigma == -1:
        steps.append(rf"Read the list as $Y[k] = {'-' if s < 0 else ''}(-1)^k\,X[{ell_tex}]$ (check a few entries: $Y[0] = {items[0]}$, $Y[1] = {items[1]}$).")
        steps.append(rf"Step 1: $V[k] = (-1)^k X[k]$ is a circular shift by $\frac{{N}}{{2}} = {h}$: $v[n] = x[{W.tex_mod(f'n - {h}', N)}]$.")
        steps.append(rf"Step 2: $(-1)^k = (-1)^{{{m}}}(-1)^{{k - {m}}}$, so $Y[k] = {'-' if s * (-1) ** m < 0 else ''}V[{ell_tex}]$, "
                     rf"a frequency shift by $\ell = {m}$: multiply $v[n]$ by $e^{{j\frac{{2\pi}}{{{N}}}{m}n}} = {phase_tex}$.")
        coef = "" if pre == 1 else "-"
        p["y_formula_tex"] = rf"y[n] = {coef}{phase_tex}\,x[{W.tex_mod(f'n - {h}', N)}]"
    else:
        steps.append(rf"Read the list as $Y[k] = X[{ell_tex}]$ (check: $Y[0] = {items[0]}$, $Y[1] = {items[1]}$).")
        steps.append(rf"Frequency shift: $e^{{j\frac{{2\pi}}{{N}}\ell n}}x[n] \leftrightarrow X[{W.tex_mod(r'k - \ell', N)}]$ with $\ell = {m}$, "
                     rf"and $e^{{j\frac{{2\pi}}{{{N}}}{m}n}} = {phase_tex}$.")
        p["y_formula_tex"] = rf"y[n] = {phase_tex}\,x[n]"
    p["steps"] = [{"step_html": t} for t in steps]
    bidx = [(n - h) % N for n in range(N)] if sigma == -1 else list(range(N))
    coefs = []
    for n in range(N):
        re, im = PHASE[(q * n) % 4]
        c = W.tex_cplx(pre * re, pre * im)
        coefs.append("" if c == "1" else ("-" if c == "-1" else c + r"\,"))
    p["y_index_tex"] = _idx_list_tex(N, bidx, coefs)
    p["disp"] = {"sigma": sigma, "m": m, "s": s}
    mistakes = {}
    if q in (1, 3):                              # modulation in the wrong direction (conjugate phase)
        mistakes["wrong_dir"] = [(re, -im) for re, im in y]
    if sigma == -1:
        mistakes["linear"] = [((pre * PHASE[(q * n) % 4][0]) * (x[n - h] if n >= h else 0),
                               (pre * PHASE[(q * n) % 4][1]) * (x[n - h] if n >= h else 0)) for n in range(N)]
        mistakes["no_shift"] = [(pre * PHASE[(q * n) % 4][0] * x[n], pre * PHASE[(q * n) % 4][1] * x[n]) for n in range(N)]
        if (-1) ** m == -1:
            mistakes["sign"] = [(-re, -im) for re, im in y]
    else:
        mistakes["shift_instead"] = [(v, 0) for v in _circ(x, lambda n: n - m)]
    return x, y, mistakes


def _fam_listf(p, N):
    """SP2025 #7(b) / SP2023 #9(b) shape with any l (e^{j 2 pi l n/N} may be irrational): y[n] = c e^{j 2 pi l n/N} x[<n - n0>_N]."""
    h = N // 2
    nn = np.arange(N)

    def form(x, c, l, n0):
        return c * np.exp(2j * np.pi * l * nn / N) * np.asarray(x, float)[(nn - n0) % N]

    triples = [(c, l, n0) for c in (1, -1) for l in range(N) for n0 in range(N)]
    while True:                                           # every (c, l, n0) must give a different y for this x
        x = _rand_x(N)
        if len({tuple(np.round(form(x, *t), 6).tolist()) for t in triples}) == len(triples):
            break
    sigma = random.choice([-1, -1, 1])
    m = random.choice([1, 3, 5, 7, 3, 2, 6] if N == 8 else [1, 2, 4, 5, 2])
    s = random.choice([1, -1])
    idx = [(k - m) % N for k in range(N)]
    sgn = [s * (sigma ** k) for k in range(N)]
    p["list_idx"], p["list_sgn"] = idx, sgn
    items = [("-" if g < 0 else "") + f"X_{{{i}}}" for i, g in zip(idx, sgn)]
    p["Y_def_tex"] = r"\{Y[k]\}_{k=0}^{" + str(N - 1) + "} = " + W.tex_list(items)
    cval = s * (-1) ** m if sigma == -1 else s
    n0 = h if sigma == -1 else 0
    ell_tex = W.tex_mod(f"k - {m}", N)
    ph = rf"e^{{j{_frac_pi(Fr(2 * m, N), 'n')}}}"
    steps = []
    if sigma == -1:
        steps.append(rf"Read the list as $Y[k] = {'-' if s < 0 else ''}(-1)^k\,X[{ell_tex}]$ (check a few entries: $Y[0] = {items[0]}$, $Y[1] = {items[1]}$).")
        steps.append(rf"Step 1: $V[k] = (-1)^k X[k]$ is a circular shift by $\frac{{N}}{{2}} = {h}$: $v[n] = x[{W.tex_mod(f'n - {h}', N)}]$.")
        steps.append(rf"Step 2: $(-1)^k = (-1)^{{{m}}}(-1)^{{k - {m}}}$, so $Y[k] = {'-' if cval < 0 else ''}V[{ell_tex}]$, "
                     rf"a frequency shift by $\ell = {m}$: multiply $v[n]$ by $e^{{j\frac{{2\pi}}{{{N}}}\ell n}} = {ph}$.")
        p["y_formula_tex"] = rf"y[n] = {'-' if cval < 0 else ''}{ph}\,x[{W.tex_mod(f'n - {h}', N)}]"
    else:
        steps.append(rf"Read the list as $Y[k] = {'-' if s < 0 else ''}X[{ell_tex}]$ (check: $Y[0] = {items[0]}$, $Y[1] = {items[1]}$).")
        steps.append(rf"Frequency shift: $e^{{j\frac{{2\pi}}{{N}}\ell n}}x[n] \leftrightarrow X[{W.tex_mod(r'k - \ell', N)}]$ with $\ell = {m}$: "
                     rf"the factor is ${ph}$, and there is no time shift.")
        p["y_formula_tex"] = rf"y[n] = {'-' if cval < 0 else ''}{ph}\,x[n]"
    p["steps"] = [{"step_html": t} for t in steps]
    p["disp"] = {"sigma": sigma, "m": m, "s": s}
    p["formula"] = {"c": cval, "l": m, "n0": n0}
    return x, cval, m, n0


def _fam_square(p, N):
    x = _rand_x(N, -2, 3)
    lin = np.convolve(x, x).astype(int).tolist()               # 2N - 1 samples
    y = [lin[n] + (lin[n + N] if n + N < len(lin) else 0) for n in range(N)]
    p["Y_def_tex"] = r"Y[k] = \big(X[k]\big)^2"
    steps = [rf"$Y[k] = X[k]\,X[k]$, and a product of $N$-point DFTs is an $N$-point <b>circular</b> convolution: "
             rf"$y[n] = \sum_{{m=0}}^{{{N - 1}}} x[m]\,x[{W.tex_mod('n - m', N)}]$.",
             rf"Shortcut: compute the linear convolution $c = x * x$ (length $2N - 1 = {2 * N - 1}$) and wrap it: $y[n] = c[n] + c[n + {N}]$ for $0 \le n \le {N - 1}$."]
    p["steps"] = [{"step_html": s} for s in steps]
    p["lin_tex"] = W.tex_list([str(v) for v in lin])
    p["wrap_tex"] = W.tex_list([f"{lin[n]} {'+' if n + N >= len(lin) or lin[n + N] >= 0 else '-'} {abs(lin[n + N]) if n + N < len(lin) else 0}"
                                for n in range(N)])
    p["y_formula_tex"] = rf"y[n] = c[n] + c[n + {N}]"
    p["y_index_tex"] = p["wrap_tex"]
    p["is_square"] = True
    return x, [(v, 0) for v in y], {"truncated": lin[:N], "elementwise": [v * v for v in x]}


def _fam_evenodd(p, N):
    even = random.random() < 0.6
    while True:
        x = _rand_x(N)
        for n in range(1, (N + 1) // 2):
            if (x[n] - x[N - n]) % 2:
                x[N - n] += 1 if x[N - n] < 5 else -1
        if len({tuple(x[(n - m) % N] for n in range(N)) for m in range(N)}) < N or [x[(-n) % N] for n in range(N)] == x:
            continue
        r = _circ(x, lambda n: -n)
        part = [(a + b) // 2 if even else (a - b) // 2 for a, b in zip(x, r)]
        if N == 8 and sum(1 for v in part if v) >= 3 or N == 6 and part[1] != 0 and part[2] != 0:
            break                                         # no mostly-zero answers
    if even:
        p["Y_def_tex"] = r"Y[k] = \mathrm{Re}\{X[k]\}"
        y = [(a + b) // 2 for a, b in zip(x, r)]
        steps = [r"$\mathrm{Re}\{X[k]\} = \frac12\big(X[k] + X^*[k]\big)$, and for real $x$ the conjugate $X^*[k]$ is the DFT of "
                 rf"$x[{W.tex_mod('-n', N)}]$ (conjugation property).",
                 rf"So $y[n] = \frac12\big(x[n] + x[{W.tex_mod('-n', N)}]\big)$, the circular even part of $x$."]
        p["y_formula_tex"] = rf"y[n] = \tfrac12\left(x[n] + x[{W.tex_mod('-n', N)}]\right)"
        p["y_index_tex"] = W.tex_list([rf"\tfrac12({a} {'+' if b >= 0 else '-'} {abs(b)})" for a, b in zip(x, r)])
        nohalf = [a + b for a, b in zip(x, r)]
        flip = [(a + b) / 2 for a, b in zip(x, x[::-1])]
    else:
        p["Y_def_tex"] = r"Y[k] = j\,\mathrm{Im}\{X[k]\}"
        y = [(a - b) // 2 for a, b in zip(x, r)]
        steps = [r"$j\,\mathrm{Im}\{X[k]\} = \frac12\big(X[k] - X^*[k]\big)$, and for real $x$ the conjugate $X^*[k]$ is the DFT of "
                 rf"$x[{W.tex_mod('-n', N)}]$ (conjugation property).",
                 rf"So $y[n] = \frac12\big(x[n] - x[{W.tex_mod('-n', N)}]\big)$, the circular odd part of $x$ (always $y[0] = 0$)."]
        p["y_formula_tex"] = rf"y[n] = \tfrac12\left(x[n] - x[{W.tex_mod('-n', N)}]\right)"
        p["y_index_tex"] = W.tex_list([rf"\tfrac12({a} {'-' if b >= 0 else '+'} {abs(b)})" for a, b in zip(x, r)])
        nohalf = [a - b for a, b in zip(x, r)]
        flip = [(a - b) / 2 for a, b in zip(x, x[::-1])]
    p["steps"] = [{"step_html": s} for s in steps]
    p["disp"] = {"even": even}
    return x, [(v, 0) for v in y], {"no_half": nohalf, "plain_flip": flip}


# ----------------------------------------------------------------------------- (b) which list is real
def _gauss(nonreal=False):
    while True:
        a, b = random.randint(-3, 4), random.randint(-3, 3)
        if (b != 0 or not nonreal) and (a, b) != (0, 0):
            return a, b


def _hermitian(M):
    """A random Hermitian list of Gaussian integers: (X0, pairs, mid, full list)."""
    while True:
        X0 = random.randint(-3, 6)
        pairs = [_gauss(nonreal=(i == 0)) for i in range((M - 1) // 2)]
        mid = [random.randint(-3, 4)] if M % 2 == 0 else []
        full = [(X0, 0)] + pairs + [(v, 0) for v in mid] + [(a, -b) for a, b in reversed(pairs)]
        if X0 != 0 and len(set(full)) >= M - 1:
            return X0, pairs, mid, full


def _real_lists():
    """Four candidate lists, each built from its OWN random Hermitian base and then (for the three
    distractors) broken in one classic way, so comparing the options entry by entry gives nothing away."""
    M = random.choice([5, 6, 7, 8])
    lists = {"good": _hermitian(M)[3]}
    X0, _, _, full = _hermitian(M)
    lists["x0_complex"] = [(X0, random.choice([1, -1, 2, -2]))] + full[1:]
    X0, pairs, mid, _ = _hermitian(M)
    lists["equal_pairs"] = [(X0, 0)] + pairs + [(v, 0) for v in mid] + list(reversed(pairs))
    L = M // 2                                            # plain-flip pairing: entry i with M-1-i
    first = [(random.randint(-3, 6), random.choice([1, -1, 2]))] + [_gauss(nonreal=True) for _ in range(L - 1)]
    middle = [(random.randint(-3, 4), 0)] if M % 2 else []
    lists["plain_flip"] = first + middle + [(a, -b) for a, b in reversed(first)]
    _, _, _, full = _hermitian(M)
    bad = list(full)
    if M % 2 == 0:
        bad[M // 2] = (bad[M // 2][0], random.choice([1, -1, 2, -2]))
        lists["mid_complex"] = bad
    else:
        a1, b1 = bad[1]
        bad[M - 1] = (-a1, b1)                            # -X^*[1] instead of X^*[1]
        lists["neg_pair"] = bad
    _, _, _, full = _hermitian(M)
    bad = list(full)
    a1, b1 = bad[M - 1]
    bad[M - 1] = (a1 + random.choice([1, -1, 2, -2]), b1)   # Im stays odd, Re no longer even
    lists["re_broken"] = bad
    return M, lists


def odd_imag(lst):
    M = len(lst)
    return all(lst[k][1] == -lst[(-k) % M][1] for k in range(M))


def even_real(lst):
    M = len(lst)
    return all(lst[k][0] == lst[(-k) % M][0] for k in range(M))


def majority(lists):
    """Per-entry majority vote over the options (None where the top count is tied)."""
    out = []
    for i in range(len(lists[0])):
        cnt = Counter(tuple(l[i]) for l in lists).most_common()
        out.append(None if len(cnt) > 1 and cnt[0][1] == cnt[1][1] else cnt[0][0])
    return out


def agreement(lists, i):
    """Number of entries option i shares with the other options (position by position)."""
    return sum(tuple(a) == tuple(b) for j, l in enumerate(lists) if j != i for a, b in zip(lists[i], l))


REASON = {
    "good": r"$Z[0]${mid} real and $Z[k] = Z^*[{M} - k]$ for every pair: Hermitian, so this signal is real.",
    "x0_complex": r"$Z[0] = \sum_n z[n]$ must be real for a real signal, but here it is complex.",
    "equal_pairs": r"The pairs are equal, $Z[k] = Z[{M} - k]$, not conjugate: e.g. $Z[1] \ne Z^*[{Mm1}]$.",
    "plain_flip": r"Conjugate pairs matched by a plain flip ($k \leftrightarrow {Mm1} - k$) instead of the circular one ($k \leftrightarrow {M} - k$); $Z[0]$ is not even real.",
    "mid_complex": r"$Z[{h}]$ pairs with itself ($\langle -{h}\rangle_{M} = {h}$), so it must be real; here it is complex.",
    "neg_pair": r"$Z[{Mm1}]$ is not the conjugate of $Z[1]$.",
    "re_broken": r"The imaginary parts are odd, but the real parts must be even too: $\mathrm{Re}\,Z[1] \ne \mathrm{Re}\,Z[{Mm1}]$.",
}


def generate(data):
    p, c = data["params"], data["correct_answers"]
    N = random.choice([6, 8])
    fam = random.choice(FAMILIES)
    p["is_square"] = False
    p["disp"] = {}
    p["family"], p["N"], p["Nm1"] = fam, N, N - 1
    p["ask_formula"], p["ask_vec"] = fam == "listf", fam != "listf"
    if fam == "listf":
        x, cval, ell, n0 = _fam_listf(p, N)
        c["coef"], c["ell"], c["n0"] = cval, ell, n0
        p["y_complex"] = True
    else:
        x, y, mistakes = {"shift": _fam_shift, "reverse": _fam_reverse, "revshift": _fam_revshift, "half": _fam_half,
                          "list": _fam_list, "square": _fam_square, "evenodd": _fam_evenodd}[fam](p, N)
        del mistakes                               # probes are recomputed by the checker
        p["y"] = [list(t) for t in y]
        p["y_complex"] = any(im != 0 for _, im in y)
        p["y_tex"] = _cvals_tex(y)
        if p["y_complex"]:
            c["y"] = pl.to_json(np.array([[complex(a, b) for a, b in y]]))
        else:
            c["y"] = pl.to_json(np.array([[float(a) for a, _ in y]]))
    p["x"] = x
    p["x_tex"] = W.tex_list([fmt.tex_num(v) for v in x])
    p["X_names_tex"] = W.tex_list([f"X_{{{k}}}" for k in range(N)])

    while True:
        M, lists = _real_lists()
        keys = ["good", "re_broken"] + random.sample([k for k in lists if k not in ("good", "re_broken")], 2)
        cand = [lists[k] for k in keys]
        if (len({tuple(v) for v in cand}) == 4 and majority(cand) != cand[0]
                and agreement(cand, 0) <= min(agreement(cand, i) for i in (1, 2, 3))
                and sum(map(odd_imag, cand)) >= 2 and sum(map(even_real, cand)) >= 2):
            break                                         # neither half of the symmetry alone finds the answer
    opts = []
    for k in keys:
        opts.append({"text": "$" + _cvals_tex(lists[k]) + "$", "correct": k == "good", "kind": k})
    random.shuffle(opts)
    p["M"], p["Mm1"], p["M_odd"], p["Mh"] = M, M - 1, M % 2 == 1, M // 2
    p["real_list_choices"] = [{"text": o["text"], "correct": o["correct"]} for o in opts]
    p["real_lists"] = [{"kind": o["kind"], "vals": [list(t) for t in lists[o["kind"]]]} for o in opts]
    mid = rf", $Z[{M // 2}]$" if M % 2 == 0 else ""
    p["real_reasons"] = [{"list_html": o["text"],
                          "verdict": "real" if o["correct"] else "not real",
                          "why_html": REASON[o["kind"]].replace("{mid}", mid).replace("{M}", str(M))
                          .replace("{Mm1}", str(M - 1)).replace("{h}", str(M // 2))}
                         for o in opts]
