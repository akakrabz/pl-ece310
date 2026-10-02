"""Systems y[n] = T{x}[n] with known Linear / Time-invariant / Causal / BIBO-stable flags.

Two sources (used by questions/systems/system-properties):

* ``BANK`` — every row of the study site's system-property bank whose answers are official (the
  seven past Midterm 1 tables, HW1 #5-6, HW2 #1(c), the Lecture 3 practice list) plus the three
  practice-problem systems of "Classifying system properties". Each entry: id, LaTeX of the right
  side, flags "LTCS" as Y/N, one reason per property, source.
* ``generate_system()`` — a compositional generator. A system is a sum of one or two terms
  (plus possibly a constant). A basic term is c * g[n] * f(x[phi(n)]) with a gain g, a value map f
  and an index map phi; a special term is a sum over input samples (moving sum, accumulator,
  geometric sum). Flags combine by rules that are true for the restricted combinations allowed here:
    - L  iff every term is linear and there is no constant (at most ONE nonlinear term, and a constant
         only when no term has f(0) != 0, so the zero-input / scaling counterexamples are exact);
    - TI iff every term is TI (gain constant and phi a pure shift);
    - C  iff every term is causal (never a gain that vanishes where the term looks into the future, never
         a plain sample x[n-k] inside a moving-sum window, never x[|n|] together with x[-n] under the same
         f: those could cancel the future samples, see _cancels);
    - S  iff every term is stable (at most ONE unstable term, so nothing can cancel).
  Gains are applied OUTSIDE f, never |.| or squares of +-1 gains (those would silently become TI).
  The independent checker (tools/checks/systems__system-properties.py) re-evaluates every system
  numerically and verifies the four flags.

All randomness uses the module ``random`` (PrairieLearn seeds it). Pure-JSON output.
"""

from __future__ import annotations

import random
from fractions import Fraction as Fr

from . import fmt

PROPS = ["Linear", "Time-invariant", "Causal", "BIBO stable"]

# ----------------------------------------------------------------------------- the bank
# (id, LaTeX of the right side of y[n] = ..., flags LTCS, [reason L, TI, C, S], source)
BANK = [
    ("abs_n_gain", r"\lvert n\rvert\,x[n]", "YNYN",
     [r"A fixed gain $\lvert n\rvert$ times $x[n]$ is linear.",
      r"$x = \delta[n]$ gives $y = 0$ (the gain is $0$ at $n = 0$), but $x = \delta[n-1]$ gives $\delta[n-1] \ne 0$.",
      r"$y[n]$ uses only the present sample $x[n]$.",
      r"The bounded input $x[n] = 1$ gives $y[n] = \lvert n\rvert$, which grows without bound."], "FA2025 #2"),
    ("abs_diff", r"\lvert x[n] - x[n-1]\rvert", "NYYY",
     [r"$x$ and $-x$ give the same output, so $T\{-x\} \ne -T\{x\}$.",
      r"No $n$ outside the brackets: delaying $x$ delays $y$.",
      r"Uses $x[n]$ and $x[n-1]$ only.",
      r"$\lvert x\rvert \le B \Rightarrow \lvert y\rvert \le 2B$."], "FA2025 #2"),
    ("conv_2n_uneg", r"x[n] * 2^{n}u[-n]", "YYNY",
     [r"Convolution with a fixed $h[n]$ is linear.", r"Convolution with a fixed $h[n]$ is time-invariant.",
      r"$h[-1] = \tfrac12 \ne 0$: $y[n]$ uses the future samples $x[n+1], x[n+2], \dots$",
      r"$\sum_n \lvert h[n]\rvert = \sum_{n \le 0} 2^{n} = 2 \lt \infty$."], "FA2025 #2"),
    ("conv_jn", r"x[n] * j^{n}u[n]", "YYYN",
     [r"Convolution with a fixed $h[n]$ is linear.", r"Convolution with a fixed $h[n]$ is time-invariant.",
      r"$h[n] = j^n u[n]$ is zero for $n \lt 0$.",
      r"$\lvert h[n]\rvert = 1$ for all $n \ge 0$, so $\sum\lvert h\rvert = \infty$: $x = j^n u[n]$ gives $(n+1)j^n u[n]$."], "SP2025 #2"),
    ("abs_idx_aff", r"2x[\lvert n\rvert] + 10", "NNNY",
     [r"$x = 0$ gives $y = 10 \ne 0$.",
      r"$x[\lvert n\rvert]$ folds the time axis at $n = 0$; a delayed input is not folded the same way.",
      r"$y[-3] = 2x[3] + 10$ needs a future sample.", r"$\lvert y\rvert \le 2B + 10$."], "SP2025 #2"),
    ("exp_plus1", r"e^{x[n]+1}", "NYYY",
     [r"$x = 0$ gives $y = e \ne 0$.", r"The same memoryless map at every $n$.", r"Uses only $x[n]$.",
      r"$\lvert x\rvert \le B \Rightarrow \lvert y\rvert \le e^{B+1}$."], "SP2025 #2"),
    ("prod_next", r"x[n]\,x[n+1]", "NYNY",
     [r"A product of input samples: doubling $x$ quadruples $y$.", r"No explicit $n$: shifting $x$ shifts $y$.",
      r"Needs $x[n+1]$, a future sample.", r"$\lvert y\rvert \le B^2$."], "FA2024 #2"),
    ("inv_abs_gain", r"\dfrac{x[n]}{\lvert n\rvert + 1}", "YNYY",
     [r"A fixed gain $\tfrac{1}{\lvert n\rvert+1}$ times $x[n]$ is linear.", r"The gain changes with $n$: $\delta[n] \mapsto \delta[n]$ but $\delta[n-1] \mapsto \tfrac12\delta[n-1]$.",
      r"Uses only $x[n]$.", r"$\lvert n\rvert + 1 \ge 1$, so $\lvert y\rvert \le \lvert x\rvert \le B$."], "FA2024 #2"),
    ("sin_plus_x0", r"\sin(x[n]) + x[0]", "NNNY",
     [r"$\sin$ is not linear: $\sin(2v) \ne 2\sin v$.", r"$x[0]$ is a fixed-time sample: delay $x$ and $y[n]$ still reads time $0$.",
      r"For $n \lt 0$, $y[n]$ needs $x[0]$, a future sample.", r"$\lvert y\rvert \le 1 + B$."], "FA2024 #2"),
    ("log_gain", r"\log(\lvert n\rvert + 1)\,x[n]", "YNYN",
     [r"A fixed gain times $x[n]$ is linear.", r"The gain $\log(\lvert n\rvert+1)$ changes with $n$.", r"Uses only $x[n]$.",
      r"$x = 1$ gives $y = \log(\lvert n\rvert + 1) \to \infty$ (slowly, but without bound)."], "FA2023 #2"),
    ("conv_u_np1", r"x[n] * u[n+1]", "YYNN",
     [r"Convolution with a fixed $h[n]$ is linear.", r"Convolution with a fixed $h[n]$ is time-invariant.",
      r"$h[-1] = 1 \ne 0$: $y[n]$ uses $x[n+1]$.", r"$\sum\lvert h\rvert = \sum_{n \ge -1} 1 = \infty$; $x = u[n]$ gives an output growing like $n$."], "FA2023 #2"),
    ("plus3", r"x[n] + 3", "NYYY",
     [r"Affine, not linear: $x = 0$ gives $y = 3 \ne 0$.", r"No explicit $n$.", r"Uses only $x[n]$.", r"$\lvert y\rvert \le B + 3$."], "FA2023 #2"),
    ("conv_m1n", r"x[n] * (-1)^{n}u[n]", "YYYN",
     [r"Convolution with a fixed $h[n]$ is linear.", r"Convolution with a fixed $h[n]$ is time-invariant.",
      r"$h[n] = (-1)^n u[n]$ is zero for $n \lt 0$.",
      r"$\lvert h[n]\rvert = 1$ for $n \ge 0$: $x = (-1)^n u[n]$ gives $(n+1)(-1)^n u[n]$."], "SP2023 #2"),
    ("div_x2", r"\dfrac{x[n]}{x[2]}", "NNNN",
     [r"Scaling $x$ by $2$ leaves $y$ unchanged, so $T\{2x\} \ne 2T\{x\}$.", r"$x[2]$ is a fixed-time sample.",
      r"For $n \lt 2$, $y[n]$ needs the future sample $x[2]$.",
      r"The bounded input $x = \delta[n]$ has $x[2] = 0$: division by zero."], "SP2023 #2"),
    ("cos2_gain", r"\cos^{2}\!\left(\tfrac{\pi}{2}n\right)x[n]", "YNYY",
     [r"A fixed gain times $x[n]$ is linear.", r"The gain is $1, 0, 1, 0, \dots$: $\delta[n]$ passes, $\delta[n-1]$ is blocked.",
      r"Uses only $x[n]$.", r"$\lvert y\rvert \le \lvert x\rvert \le B$."], "SP2023 #2"),
    ("cos_shift_gain", r"x[n]\cos\!\left(\tfrac{\pi(n-2)}{3}\right)", "YNYY",
     [r"A fixed gain times $x[n]$ is linear.", r"The $-2$ is inside the cosine: an $n$-dependent gain, not a delay.",
      r"Uses only $x[n]$.", r"$\lvert\cos\rvert \le 1$, so $\lvert y\rvert \le B$."], "SP2021 #3"),
    ("x3_prod", r"x[3]\,x[n]", "NNNY",
     [r"Doubling $x$ quadruples $y$.", r"$x[3]$ is a fixed-time sample.", r"For $n \lt 3$, $y[n]$ needs $x[3]$, a future sample.",
      r"$\lvert y\rvert \le B^2$."], "SP2021 #3"),
    ("cgain", r"(0.8 + 0.8j)^{n}\,x[n]", "YNYN",
     [r"A fixed (complex) gain times $x[n]$ is linear.", r"The gain depends on $n$.", r"Uses only $x[n]$.",
      r"$\lvert 0.8+0.8j\rvert = 0.8\sqrt2 \approx 1.13 \gt 1$: $x = 1$ gives $\lvert y[n]\rvert = 1.13^n \to \infty$."], "SP2021 #3"),
    ("idx_absn_plus_n", r"x[\lvert n\rvert + n]", "YNNY",
     [r"It only re-indexes $x$, so it is linear.", r"The index map ($0$ for $n \le 0$, $2n$ for $n \gt 0$) is not a plain shift.",
      r"$y[1] = x[2]$ needs a future sample.", r"Re-indexing cannot make a bounded signal unbounded."], "FA2019 #2(a)"),
    ("clip", r"\min\{\max\{x[n], -2\},\ 2\}\quad\text{(clipping to } [-2, 2]\text{)}", "NYYY",
     [r"$x = 3$ gives $2$, but $x = 6$ gives $2$, not $4$.", r"The clipping rule is the same at every $n$.",
      r"Uses only $x[n]$.", r"$-2 \le y \le 2$ always."], "HW1 #5"),
    ("window", r"x[n]\big(u[n] - u[n-4]\big)", "YNYY",
     [r"Multiplication by a fixed window is linear.", r"The window stays put: $\delta[n-3]$ passes but $\delta[n-4]$ is cut off.",
      r"Uses only $x[n]$.", r"$\lvert y\rvert \le \lvert x\rvert$."], "HW1 #6 (N = 4)"),
    ("half_abs_gain", r"\left(\tfrac{1}{2}\right)^{\lvert n\rvert}x[n]", "YNYY",
     [r"A fixed gain times $x[n]$ is linear.", r"$\delta[n] \mapsto \delta[n]$ but $\delta[n-1] \mapsto \tfrac12\delta[n-1]$.",
      r"Uses only $x[n]$.", r"$\lvert y\rvert \le \lvert x\rvert \le B$."], "HW2 #1(c)"),
    ("avg2", r"\tfrac{1}{2}x[n] + \tfrac{1}{2}x[n-1]", "YYYY",
     [r"A weighted sum of input samples is linear.", r"Constant weights, no explicit $n$.", r"Present and past samples only.",
      r"$\lvert y\rvert \le \tfrac12 B + \tfrac12 B = B$."], "Lecture 3"),
    ("absx", r"\lvert x[n]\rvert", "NYYY",
     [r"$\lvert -x\rvert = \lvert x\rvert$, not $-\lvert x\rvert$.", r"No explicit $n$.", r"Uses only $x[n]$.", r"$\lvert y\rvert = \lvert x\rvert \le B$."], "Lecture 3"),
    ("expx", r"e^{x[n]}", "NYYY",
     [r"$x = 0$ gives $y = 1 \ne 0$.", r"No explicit $n$.", r"Uses only $x[n]$.", r"$\lvert y\rvert \le e^{B}$."], "Lecture 3"),
    ("n_diff", r"n\big(x[n] - x[n-1]\big)", "YNYN",
     [r"A fixed gain $n$ times a linear combination: linear.", r"The gain $n$ depends on time.", r"Present and past samples only.",
      r"$x = (-1)^n$ gives $y = 2n(-1)^n$, unbounded."], "Lecture 3"),
    ("max0", r"\max\{0,\ x[n]\}", "NYYY",
     [r"$x = -1$ gives $0$, but $-x = 1$ gives $1 \ne -0$.", r"No explicit $n$.", r"Uses only $x[n]$.", r"$0 \le y \le B$."], "Lecture 3"),
    ("x0_prod", r"x[n]\,x[0]", "NNNY",
     [r"Doubling $x$ quadruples $y$.", r"$x[0]$ is a fixed-time sample.", r"For $n \lt 0$, $y[n]$ needs $x[0]$, a future sample.",
      r"$\lvert y\rvert \le B^2$."], "Lecture 3 (annotated slides)"),
    ("diff", r"x[n] - x[n-1]", "YYYY",
     [r"A difference of input samples is linear.", r"No explicit $n$.", r"Present and past samples only.", r"$\lvert y\rvert \le 2B$."], "Lecture 3"),
    ("down2", r"x[2n]", "YNNY",
     [r"Re-indexing is linear.", r"$T\{x[n-n_0]\} = x[2n-n_0]$ but $y[n-n_0] = x[2n-2n_0]$.", r"$y[1] = x[2]$ needs a future sample.",
      r"Re-indexing keeps $\lvert y\rvert \le B$."], "Lecture 3"),
    ("fold", r"x[\lvert n\rvert]", "YNNY",
     [r"Re-indexing is linear.", r"Folding at $n = 0$ does not commute with shifts.", r"$y[-2] = x[2]$ needs a future sample.",
      r"Re-indexing keeps $\lvert y\rvert \le B$."], "Lecture 3"),
    ("avg3", r"\tfrac{1}{3}\big(x[n] + x[n-1] + x[n-2]\big)", "YYYY",
     [r"A weighted sum of input samples is linear.", r"Constant weights, no explicit $n$.", r"Present and past samples only.",
      r"An average of bounded samples: $\lvert y\rvert \le B$."], "Lecture 3"),
    ("idx_absn_minus_n", r"x[\lvert n\rvert - n]", "YNNY",
     [r"Re-indexing is linear.", r"The index map $\lvert n\rvert - n$ is not a plain shift.", r"$y[-1] = x[2]$ is a future sample.",
      r"Re-indexing keeps $\lvert y\rvert \le B$."], "Lecture 3"),
    ("np1_diff", r"(n+1)\big(x[n] - x[n-1]\big)", "YNYN",
     [r"A fixed gain $(n+1)$ times a linear combination: linear.", r"The gain depends on $n$.", r"Present and past samples only.",
      r"$x = (-1)^n$ gives $y = 2(n+1)(-1)^n$, unbounded."], "Lecture 3"),
    ("sin2", r"\sin^{2}(x[n])", "NYYY",
     [r"$\sin^2$ is not linear: $x = \tfrac{\pi}{2}$ gives $1$, $x = \pi$ gives $0 \ne 2$.", r"No explicit $n$.", r"Uses only $x[n]$.",
      r"$0 \le y \le 1$ whatever the input."], "Lecture 3"),
    ("ln_abs", r"\ln\lvert x[n]\rvert", "NYYN",
     [r"$x = 1$ gives $0$, $x = 2$ gives $\ln 2 \ne 0$.", r"No explicit $n$.", r"Uses only $x[n]$.",
      r"$x[n] = 0$ is a bounded input, but $\ln 0 = -\infty$."], "Lecture 3"),
    ("median", r"\operatorname{median}\{x[n],\ x[n-1],\ x[n-2]\}", "NYYY",
     [r"Windows $(1,0,0)$ and $(0,1,0)$ have median $0$, their sum $(1,1,0)$ has median $1$.", r"The same rule at every $n$.",
      r"Present and past samples only.", r"The median is one of the samples, so $\lvert y\rvert \le B$."], "Lecture 3"),
    ("square", r"x^{2}[n]", "NYYY",
     [r"Doubling $x$ quadruples $y$.", r"No explicit $n$.", r"Uses only $x[n]$.", r"$\lvert y\rvert \le B^2$."], "Lecture 3 notes, Ex. 2"),
    ("n_x3n", r"n\,x[3n]", "YNNN",
     [r"Gain times re-indexing: linear.", r"Both the gain $n$ and the index $3n$ break time-invariance.", r"$y[1] = x[3]$ needs a future sample.",
      r"$x = 1$ gives $y = n$, unbounded."], "Lecture 3 notes, Ex. 5"),
    ("fir3", r"x[n] + x[n-2] - x[n-4]", "YYYY",
     [r"A combination of input samples is linear.", r"Constant coefficients, no explicit $n$.", r"Present and past samples only.",
      r"$\lvert y\rvert \le 3B$."], "Lecture 3 notes, Ex. 6"),
    ("x10_exp", r"x^{10}[n] + e^{x[n]}", "NYYY",
     [r"Powers and exponentials are not linear ($x = 0$ gives $y = 1$).", r"No explicit $n$.", r"Uses only $x[n]$.",
      r"$\lvert y\rvert \le B^{10} + e^{B}$."], "Lecture 3 notes, Ex. 7"),
    ("conv_half_np2", r"x[n] * \left(\tfrac{1}{2}\right)^{n}u[n+2]", "YYNY",
     [r"Convolution with a fixed $h[n]$ is linear.", r"Convolution with a fixed $h[n]$ is time-invariant.",
      r"$h[-2] = 4$, $h[-1] = 2 \ne 0$: $y[n]$ uses $x[n+1]$ and $x[n+2]$.",
      r"$\sum\lvert h\rvert = 4 + 2 + \sum_{n\ge0}(\tfrac12)^n = 8$."], "Classifying system properties, Practice 1(a)"),
    ("x_plus_xrev", r"x[n] + x[-n]", "YNNY",
     [r"A sum of re-indexed samples is linear.", r"The index $-n$ is not a shift.", r"At $n = -1$ it reads $x[1]$, a future sample.",
      r"$\lvert y\rvert \le 2B$."], "Classifying system properties, Practice 1(b)"),
    ("accum", r"\sum_{k=-\infty}^{n} x[k]", "YYYN",
     [r"A sum of input samples is linear.", r"It is the convolution with $h[n] = u[n]$, hence time-invariant.",
      r"Only $x[k]$ with $k \le n$.", r"The bounded input $u[n]$ gives $(n+1)u[n]$."], "Classifying system properties, Practice 1(c)"),
    ("cos_times_prev", r"\cos(x[n])\,x[n-1]", "NYYY",
     [r"$x = u[n]$ gives $y[1] = \cos 1$, but $x = 2u[n]$ gives $y[1] = 2\cos 2 \ne 2\cos 1$.", r"No explicit $n$.",
      r"Present and past samples only.", r"$\lvert y\rvert \le \lvert x[n-1]\rvert \le B$."], "Classifying system properties, Practice 1(d)"),
    ("half_n_prev", r"x[n] + \left(\tfrac{1}{2}\right)^{n}x[n-1]", "YNYN",
     [r"A sum of fixed gains times samples: linear.", r"$\delta[n-1] \mapsto \delta[n-1] + \tfrac14\delta[n-2]$ but $\delta[n-2] \mapsto \delta[n-2] + \tfrac18\delta[n-3]$.",
      r"Uses $x[n]$ and $x[n-1]$ only.", r"The gain $(\tfrac12)^n = 2^{-n}$ explodes as $n \to -\infty$: $x = u[-n]$ gives $y[-10] = 1025$."],
     "Classifying system properties, Practice 2"),
    ("future_sum", r"\sum_{k=n-2}^{\infty}\left(\tfrac{1}{3}\right)^{k-n}x[k]", "YYNY",
     [r"A sum of fixed weights times samples: linear.", r"It is the convolution with $h[m] = 3^{m}u[2-m]$.",
      r"$h[-1] = \tfrac13 \ne 0$: $y[n]$ reads $x[k]$ for $k$ up to $+\infty$.", r"$\sum\lvert h\rvert = 9 + 3 + 1 + \tfrac12 = \tfrac{27}{2}$."],
     "Classifying system properties, Practice 3"),
]

# ----------------------------------------------------------------------------- generator atoms
# gain: (LaTeX, time-varying, bounded)
GAINS = {
    "one": ("", False, True),
    "n": ("n", True, False),
    "half_n": (r"\left(\tfrac{1}{2}\right)^{n}", True, False),
    "m1n": (r"(-1)^{n}", True, True),
    "cos3": (r"\cos\!\left(\tfrac{\pi n}{3}\right)", True, True),
    "un": ("u[n]", True, True),
    "half_abs": (r"\left(\tfrac{1}{2}\right)^{\lvert n\rvert}", True, True),
}
# value map f: (linear, f(0) != 0, scale factor for 2x (nonlinear f(0)=0 only), bound LaTeX)
FUNCS = {
    "id": (True, False, 2, "B"),
    "sq": (False, False, 4, "B^{2}"),
    "cube": (False, False, 8, "B^{3}"),
    "exp": (False, True, None, "e^{B}"),
    "cos": (False, True, None, ""),
}
# index maps: ("shift", k) is x[n-k]; others below
MAPS_NONSHIFT = {"rev": "-n", "down2": "2n", "fold": r"\lvert n\rvert"}


def _idx_tex(phi):
    kind, k = phi
    if kind == "shift":
        return "n" if k == 0 else (f"n - {k}" if k > 0 else f"n + {-k}")
    return MAPS_NONSHIFT[kind]


def _f_tex(f, phi):
    xs = f"x[{_idx_tex(phi)}]"
    return {"id": xs, "sq": f"x^{{2}}[{_idx_tex(phi)}]", "cube": f"x^{{3}}[{_idx_tex(phi)}]",
            "exp": f"e^{{{xs}}}", "cos": rf"\cos\!\left({xs}\right)"}[f]


def _basic_body(g, f, phi):
    gt = GAINS[g][0]
    ft = _f_tex(f, phi)
    if not gt:
        return ft
    sep = r"\," if (g in ("n",) and f in ("id", "sq", "cube")) or g == "un" else r"\,"
    return gt + sep + ft


# special linear terms: (LaTeX builder, TI, causal, stable)
def _special(kind, a=0, b=0):
    if kind == "mavg":
        lo = "n" if a == 0 else f"n-{a}"
        hi = "n" if b == 0 else (f"n+{b}" if b > 0 else f"n-{-b}")
        return rf"\sum_{{k={lo}}}^{{{hi}}} x[k]", True, b <= 0, True
    if kind == "accum":
        hi = "n" if a == 0 else f"n-{a}"
        return rf"\sum_{{k=-\infty}}^{{{hi}}} x[k]", True, True, False
    if kind == "geo":
        return r"\sum_{k=0}^{\infty}\left(\tfrac{1}{2}\right)^{k}x[n-k]", True, True, True
    if kind == "geo_future":
        return r"\sum_{k=0}^{\infty}\left(\tfrac{1}{2}\right)^{k}x[n+k]", True, False, True
    if kind == "accum_future":
        return r"\sum_{k=n}^{\infty} x[k]", True, False, False
    raise ValueError(kind)


COEFS = [Fr(1), Fr(1), Fr(2), Fr(3), Fr(1, 2), Fr(-1), Fr(-2)]


GAIN_W = {"one": 5, "n": 1, "half_n": 1, "m1n": 1, "cos3": 1, "un": 1, "half_abs": 1}
FUNC_W = {"id": 5, "sq": 1, "cube": 1, "exp": 1, "cos": 1}


def _wchoice(w):
    keys = list(w)
    return random.choices(keys, weights=[w[k] for k in keys])[0]


def _rand_basic():
    g = _wchoice(GAIN_W)
    f = _wchoice(FUNC_W)
    if g == "un" or random.random() < 0.75:
        phi = ("shift", random.choice([-2, -1, 0, 0, 1, 1, 2, 3]))
    else:
        phi = (random.choice(list(MAPS_NONSHIFT)), 0)
    return {"type": "basic", "g": g, "f": f, "phi": list(phi)}


def _rand_special():
    kind = random.choice(["mavg", "mavg", "accum", "geo", "geo_future", "accum_future"])
    a = b = 0
    if kind == "mavg":
        a, b = random.choice([(2, 1), (1, 1), (2, 0), (3, 0), (0, 2), (1, 2)])
    elif kind == "accum":
        a = random.choice([0, 1])
    return {"type": "special", "kind": kind, "a": a, "b": b}


def term_flags(t):
    """(linear, TI, causal, stable, f0_nonzero) of one term."""
    if t["type"] == "special":
        _, ti, c, s = _special(t["kind"], t["a"], t["b"])
        return True, ti, c, s, False
    g, f, phi = t["g"], t["f"], tuple(t["phi"])
    lin = FUNCS[f][0]
    ti = (not GAINS[g][1]) and phi[0] == "shift"
    causal = phi[0] == "shift" and phi[1] >= 0
    stable = GAINS[g][2]
    return lin, ti, causal, stable, FUNCS[f][1]


def term_tex(t):
    if t["type"] == "special":
        return _special(t["kind"], t["a"], t["b"])[0]
    return _basic_body(t["g"], t["f"], tuple(t["phi"]))


def _key(t):
    return (t["type"], t.get("g"), t.get("f"), tuple(t.get("phi", ())), t.get("kind"), t.get("a"), t.get("b"))


def _indices_read(t):
    """LaTeX list of the input samples a causal term reads."""
    if t["type"] == "special":
        k = t["kind"]
        if k == "mavg":
            return [f"x[{_idx_tex(('shift', j))}]" for j in range(-t["b"], t["a"] + 1)]
        if k == "accum":
            return [r"x[k],\ k \le " + ("n" if t["a"] == 0 else f"n-{t['a']}")]
        return [r"x[n-k],\ k \ge 0"]
    return [f"x[{_idx_tex(tuple(t['phi']))}]"]


def _sys_tex(terms, d):
    out = ""
    for t in terms:
        c = Fr(t["c"][0], t["c"][1])
        body = term_tex(t)
        mag = "" if abs(c) == 1 else fmt.tex_num(abs(c))
        if mag and (body[0].isdigit() or body.startswith((r"\left", "(", r"\tfrac"))):
            mag += r"\,"
        piece = mag + body
        if not out:
            out = ("-" if c < 0 else "") + piece
        else:
            out += (" - " if c < 0 else " + ") + piece
    if d != 0:
        out += (" - " if d < 0 else " + ") + fmt.tex_num(abs(d))
    return out


def _reasons(terms, d, flags):
    L, TI, C, S = flags
    tf = [term_flags(t) for t in terms]
    r = []
    # ---- linearity
    if L:
        r.append(r"Every term is a fixed factor times samples of $x$ (no constant, no function of $x$), so superposition holds.")
    elif d != 0:
        r.append(rf"$x[n] = 0$ gives $y[n] = {fmt.tex_num(d)} \ne 0$, but a linear system maps $0$ to $0$.")
    else:
        i = next(k for k, f in enumerate(tf) if not f[0])
        t = terms[i]
        c = Fr(t["c"][0], t["c"][1])
        if FUNCS[t["f"]][1]:                       # exp or cos: f(0) = 1
            g = GAINS[t["g"]][0]
            z = fmt.tex_num(c) if not g else ((("-" if c < 0 else "") + ("" if abs(c) == 1 else fmt.tex_num(abs(c)) + r"\,")) + g)
            r.append(rf"$x[n] = 0$ gives $y[n] = {z}$, not identically $0$ (the term ${term_tex(t)}$ maps $0$ to "
                     rf"{'$e^0 = 1$' if t['f'] == 'exp' else r'$\cos 0 = 1$'}).")
        else:
            fac = FUNCS[t["f"]][2]
            r.append(rf"Doubling the input multiplies the term ${term_tex(t)}$ by ${fac}$ but every other term only by $2$, "
                     rf"so $T\{{2x\}} \ne 2\,T\{{x\}}$.")
    # ---- time invariance
    if TI:
        r.append(r"Every sample is read as $x[n-k]$ with a constant $k$ and $n$ appears nowhere else (summation limits move with $n$), "
                 r"so delaying $x$ delays $y$.")
    else:
        t = next(t for t, f in zip(terms, tf) if not f[1])
        g, phi = t["g"], tuple(t["phi"])
        if GAINS[g][1]:
            shifted = {"n": "(n-n_0)", "half_n": r"\left(\tfrac{1}{2}\right)^{n-n_0}", "m1n": r"(-1)^{n-n_0}",
                       "cos3": r"\cos\!\left(\tfrac{\pi (n-n_0)}{3}\right)", "un": "u[n-n_0]",
                       "half_abs": r"\left(\tfrac{1}{2}\right)^{\lvert n-n_0\rvert}"}[g]
            r.append(rf"The factor ${GAINS[g][0]}$ (in ${term_tex(t)}$) depends on $n$: for the input $x[n-n_0]$ that term "
                     rf"still has ${GAINS[g][0]}$, but $y[n-n_0]$ has ${shifted}$ instead, and these differ.")
        else:
            m = {"rev": (r"x[-n-n_0]", r"x[-n+n_0]"), "down2": (r"x[2n-n_0]", r"x[2n-2n_0]"),
                 "fold": (r"x[\lvert n\rvert-n_0]", r"x[\lvert n-n_0\rvert]")}[phi[0]]
            r.append(rf"In ${term_tex(t)}$ the index is not a plain shift: the input $x[n-n_0]$ produces ${m[0]}$ there, "
                     rf"but $y[n-n_0]$ contains ${m[1]}$.")
    # ---- causality
    if C:
        reads = []
        for t in terms:
            for s in _indices_read(t):
                if s not in reads:
                    reads.append(s)
        sep = r",\ "
        r.append(rf"$y[n]$ reads only ${sep.join(reads)}$: present and past samples.")
    else:
        t = next(t for t, f in zip(terms, tf) if not f[2])
        if t["type"] == "special":
            if t["kind"] == "mavg":
                r.append(rf"$y[n]$ reads $x[n+{t['b']}]$, a future sample.")
            else:
                r.append(r"$y[n]$ reads $x[n+1], x[n+2], \dots$: future samples.")
        else:
            kind, k = t["phi"]
            msg = {"rev": r"$y[-1]$ reads $x[1]$", "down2": r"$y[1]$ reads $x[2]$", "fold": r"$y[-1]$ reads $x[1]$"}
            r.append((msg[kind] if kind != "shift" else rf"$y[n]$ reads $x[n+{-k}]$") + ", a future sample.")
    # ---- stability
    if S:
        bound = []
        for t in terms:
            c = abs(Fr(t["c"][0], t["c"][1]))
            if t["type"] == "special":
                if t["kind"] == "mavg":
                    bound.append((c * (t["a"] + t["b"] + 1), "B"))
                else:                               # geo, geo_future: sum (1/2)^k = 2
                    bound.append((2 * c, "B"))
            else:
                bound.append((c, FUNCS[t["f"]][3]))
        merged = {}
        for c, body in bound:
            merged[body] = merged.get(body, Fr(0)) + c
        if d != 0:
            merged[""] = merged.get("", Fr(0)) + abs(Fr(d))
        r.append(rf"If $\lvert x[n]\rvert \le B$ then $\lvert y[n]\rvert \le {fmt.tex_sum([(c, body) for body, c in merged.items()])}$ "
                 r"(every gain here is bounded by $1$).")
    else:
        t = next(t for t, f in zip(terms, tf) if not f[3])
        if t["type"] == "special":
            if t["kind"] == "accum":
                r.append(r"The bounded input $x[n] = u[n]$ makes the running sum grow like $n$: unbounded.")
            else:
                r.append(r"The bounded input $x[n] = 1$ makes $\sum_{k \ge n} x[k]$ diverge.")
        else:
            inp = "0" if FUNCS[t["f"]][1] else "1"
            how = r"grow like $\lvert n\rvert$" if t["g"] == "n" else r"grow like $2^{-n}$ as $n \to -\infty$"
            r.append(rf"The bounded input $x[n] = {inp}$ makes the term ${term_tex(t)}$ {how}: the output is unbounded.")
    return r


def _cancels(terms):
    """Pairs that could cancel the very samples that decide a flag are not used:
    * a plain sample x[n-k] (gain 1, f = id) inside a moving-sum window
      (sum_{k=n-1}^{n+1} x[k] - x[n+1] = x[n-1] + x[n] is causal);
    * x[|n|] and x[-n] under the same f: the two index maps agree on the whole half-line n <= 0, where
      both look into the future (x[|n|] - x[-n] = 0 for n <= 0 is causal);
    * the same f and index map with gains that agree on a half-line ((1/2)^n and (1/2)^|n|, 1 and u[n]).
    Every other pair of distinct terms has weights that coincide only at isolated n, or different f."""
    basics = [t for t in terms if t["type"] == "basic"]
    for a in basics:
        for b in basics:
            if a is b or a["f"] != b["f"]:
                continue
            if {a["phi"][0], b["phi"][0]} == {"fold", "rev"}:
                return True
            # gains that agree on a half-line ((1/2)^n = (1/2)^|n| and 1 = u[n] for n >= 0) could cancel there:
            # (1/2)^n x[2n] - (1/2)^|n| x[2n] = 0 for n >= 0, exactly where x[2n] looks ahead
            if a["phi"] == b["phi"] and {a["g"], b["g"]} in ({"half_n", "half_abs"}, {"one", "un"}):
                return True
    win = [t for t in terms if t["type"] == "special" and t["kind"] == "mavg"]
    for m in win:
        for t in terms:
            if (t["type"] == "basic" and t["g"] == "one" and t["f"] == "id" and t["phi"][0] == "shift"
                    and -m["b"] <= t["phi"][1] <= m["a"]):
                return True
    return False


def generate_system():
    """Random compositional system -> dict(tex, flags [L, TI, C, S], reasons, struct)."""
    while True:
        n_terms = random.choice([1, 2, 2, 2])
        terms = []
        for _ in range(n_terms):
            t = _rand_special() if random.random() < 0.25 else _rand_basic()
            terms.append(t)
        keys = [_key(t) for t in terms]
        if len(set(keys)) != len(keys):
            continue
        tf = [term_flags(t) for t in terms]
        if sum(not f[0] for f in tf) > 1 or sum(not f[3] for f in tf) > 1:
            continue
        if n_terms == 2 and all(t["type"] == "special" for t in terms):
            continue
        if _cancels(terms):
            continue
        d = Fr(0)
        if random.random() < 0.25 and not any(f[4] for f in tf):
            d = Fr(random.choice([1, 2, 3, 5, -1, -2]))
        for i, t in enumerate(terms):
            c = random.choice(COEFS)
            t["c"] = [c.numerator, c.denominator]
        flags = [all(f[0] for f in tf) and d == 0, all(f[1] for f in tf), all(f[2] for f in tf), all(f[3] for f in tf)]
        if all(flags):                 # keep the fully LTI-causal-stable case for the bank
            continue
        if not any(flags) and random.random() < 0.7:
            continue
        return {"tex": _sys_tex(terms, d), "flags": flags, "reasons": _reasons(terms, d, flags),
                "struct": {"terms": terms, "d": [d.numerator, d.denominator]}}


def bank_system(i=None):
    e = BANK[random.randrange(len(BANK)) if i is None else i]
    sid, tex, fl, why, src = e
    return {"tex": tex, "flags": [c == "Y" for c in fl], "reasons": list(why), "id": sid, "source": src}
