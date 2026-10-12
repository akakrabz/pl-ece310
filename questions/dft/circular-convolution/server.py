"""Linear vs circular convolution (Midterm 2 family "linear vs circular convolution, the FFT":
FA2021 MT2 #6, FA2019 MT2 #11 and the true/false items FA2019 #1(g), (h); taught after Lecture 17,
tag `ahead`). DFT property from the 2019 course notes, section 2.5.1.8: the inverse N-point DFT of
X[k]H[k] is the circular convolution sum_{m=0}^{N-1} x[m] h[<n - m>_N].

x (length L1 = 3..5) and h (length L2 = 2..4) are small integers with nonzero end samples, both
starting at n = 0; L = L1 + L2 - 1 is the length of the linear convolution y = x * h.
  (a) the N-point circular convolution y_N for max(L1, L2) <= N < L: the linear result with the
      overflow y[N], ..., y[L-1] wrapped onto n = 0, ..., L-N-1 (time aliasing, y_N[n] = sum_r y[n + rN]);
  (b) the smallest N for which the circular convolution equals the linear one: L;
  (c) checkbox: the indices n where y_N[n] differs from y[n]: the indices 0..L-N-1 that receive a nonzero
      wrapped sample. In 30% of the variants two or more samples wrap and exactly one wrapped
      sample equals 0 (never the last one, y[L-1] = x[L1-1] h[L2-1] != 0), so the answer has a gap,
      e.g. {0, 2}, and "the first L - N indices" is not enough;
  (d) FA2019 #11: two sequences of lengths P and Q convolved with a radix-2 FFT: the minimum number of
      zeros appended to each = (next power of two >= P + Q - 1) - P and - Q; P + Q - 1 is never a power
      of two, so the plain DFT answer (P + Q - 1 - P, ...) is wrong.
"""

import random

import numpy as np
import prairielearn as pl
from ece310 import fmt, seq

NZ = [-3, -2, -1, 1, 2, 3]


def _rand(L):
    v = [random.randint(-3, 3) for _ in range(L)]
    v[0], v[-1] = random.choice(NZ), random.choice(NZ)
    return v


def _seq_tex(name, v):
    return rf"\{{{name}[n]\}}_{{n=0}}^{{{len(v) - 1}}} = \{{" + r",\ ".join(str(a) for a in v) + r"\}"


def _vec_tex(v):
    return r"\{" + r",\ ".join(str(a) for a in v) + r"\}"


def _table_html(x, h, y):
    """Shift-and-add table: one row x[k] h[n-k] per nonzero x[k], then the column sums."""
    L = len(y)
    head = "".join(f"<th>$n={n}$</th>" for n in range(L))
    rows = []
    for k, xv in enumerate(x):
        if xv == 0:
            continue
        cells = "".join(f"<td>${fmt.tex_num(xv * h[n - k])}$</td>" if 0 <= n - k < len(h) else "<td></td>"
                        for n in range(L))
        lab = r"x[0]\,h[n]" if k == 0 else rf"x[{k}]\,h[n-{k}]"
        rows.append(f"<tr><th>${lab}$</th>{cells}</tr>")
    tot = "".join(f"<td><b>${fmt.tex_num(v)}$</b></td>" for v in y)
    return ('<table class="table table-sm table-bordered text-center" style="width:auto">'
            f"<tr><th></th>{head}</tr>{''.join(rows)}<tr><th>$y[n]$</th>{tot}</tr></table>")


def generate(data):
    want_gap = random.random() < 0.3
    while True:
        L1, L2 = random.choice([(3, 2), (3, 3), (4, 2), (4, 3), (4, 4), (5, 2), (5, 3), (5, 4), (3, 4)])
        x, h = _rand(L1), _rand(L2)
        y, _ = seq.conv(x, 0, h, 0)
        y = [int(v) for v in y]
        L = L1 + L2 - 1
        N = L - random.choice([1, 2, 2, 3])            # wrap 1-3 samples (FA2021 #6 wraps 2)
        if N < max(L1, L2):
            continue
        tail = y[N:]
        zeros = [i for i, v in enumerate(tail) if v == 0]
        if (want_gap and (len(tail) < 2 or len(zeros) != 1)) or (not want_gap and zeros):
            continue                                   # gap: exactly one wrapped zero; else none
        yN = [y[n] + (y[n + N] if n + N < L else 0) for n in range(N)]
        break
    while True:                                        # FA2019 #11: radix-2 FFT padding
        P, Q = random.randint(6, 40), random.randint(3, 20)
        M = P + Q - 1
        if M & (M - 1) == 0:
            continue                                   # M itself a power of two: DFT and FFT answers coincide
        N2 = 1 << (M - 1).bit_length()
        break
    W = L - N
    wrap_terms = []
    for n in range(N):
        if n + N < L:
            wrap_terms.append(rf"{y[n]} {'-' if y[n + N] < 0 else '+'} {abs(y[n + N])}")
        else:
            wrap_terms.append(str(y[n]))
    p = data["params"]
    p.update({
        "x": x, "h": h, "Nc": N, "L1": L1, "L2": L2, "L": L, "Lm1": L - 1, "W": W, "Wm1": W - 1,
        "Nm1": N - 1, "x_tex": _seq_tex("x", x), "h_tex": _seq_tex("h", h),
        "y_tex": _vec_tex(y), "yN_tex": _vec_tex(yN), "wrap_tex": r"\{" + r",\ ".join(wrap_terms) + r"\}",
        "table_html": _table_html(x, h, y),
        "tail_tex": ", ".join(rf"y[{N + i}] = {v}" for i, v in enumerate(tail)),
        "idx_list_tex": ", ".join(str(i) for i in range(W) if tail[i] != 0),
        "c_expl_html": (rf"Exactly the indices that receive a wrapped sample: $n = {', '.join(str(i) for i in range(W))}$. "
                        r"The remaining samples agree with $y[n]$." if not zeros else
                        rf"The wrapped samples land on $n = {', '.join(str(i) for i in range(W))}$, but "
                        rf"$y[{N + zeros[0]}] = 0$, so $y_{{{N}}}[{zeros[0]}] = y[{zeros[0]}] + 0 = y[{zeros[0]}]$: only "
                        rf"$n = {', '.join(str(i) for i in range(W) if tail[i] != 0)}$ differ."),
        "wrap_phrase_html": (rf"the sample beyond $n = {N - 1}$ (${', '.join(rf'y[{N + i}] = {v}' for i, v in enumerate(tail))}$) "
                             rf"wrapped onto $n = 0$") if W == 1 else
                            (rf"the ${W}$ samples beyond $n = {N - 1}$ (${', '.join(rf'y[{N + i}] = {v}' for i, v in enumerate(tail))}$) "
                             rf"wrapped onto $n = {', '.join(str(i) for i in range(W))}$"),
        "sumcheck_tex": rf"\textstyle\sum x \cdot \sum h = {sum(x)}\cdot{fmt.tex_num(sum(h), paren=True)} = {sum(y)}",
        "P": P, "Q": Q, "Pm1": P - 1, "Qm1": Q - 1, "M": M, "N2": N2,
        "diff_choices": [{"text": f"$n = {n}$", "correct": n < W and tail[n] != 0} for n in range(N)],
    })
    c = data["correct_answers"]
    c["yN"] = pl.to_json(np.array([[float(v) for v in yN]]))
    c["Nmin"] = L
    c["zu"] = N2 - P
    c["zv"] = N2 - Q
