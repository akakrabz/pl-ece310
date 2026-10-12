"""Matched filter (Lecture 12, "Convolution as template matching", notes §1.1 eqs. 1-2, slides 3-5).

A template t (3 small integers, not symmetric, start index a in {-1, 0, 1}) is hidden, scaled by
c in {1, 2}, in an integer input x of length 8-10 (start 0), plus small integer noise (usually
also on the copy itself). With the matched filter h[n] = t[-n],
    y[n] = x[n] * h[n] = sum_k x[k] t[k - n] = sum_{m} x[n + m] t[m],
so y peaks at n0 = (start of the copy) - a. Asked: h (first index + row vector), the peak index
and value, y[n1] at one other index, and the normalized score at the peak
    rho[n0] = y[n0] / ( sqrt(sum t^2) sqrt(sum of x^2 over the window under the template) ) <= 1
(Cauchy-Schwarz; the course site's Lecture 12 page, section 3, gives this normalization - the
notes motivate it with the "loud is not similar" effect but do not define it, so the question
states the formula). Variants are redrawn until the maximum of y is unique (margin >= 2) and
sits on the copy."""

import math
import random

import numpy as np
import prairielearn as pl
from ece310 import fmt, svg

NZ = [-2, -1, 1, 2, 3]


def _draw():
    while True:
        t = [random.choice(NZ), random.choice([-2, -1, 0, 1, 2, 3]), random.choice(NZ)]
        if t == t[::-1] or sum(v * v for v in t) < 5:
            continue
        a = random.choice([0, 0, -1, 1])
        L = random.choice([8, 9, 10])
        c = random.choice([1, 1, 2])
        p = random.randint(1, L - 4)                      # copy occupies p..p+2 (never the very first sample)
        x = [random.choice([0, 0, 0, 1, -1, 1, -1, 2, -2]) for _ in range(L)]
        for m in range(3):
            x[p + m] = c * t[m]
        if random.random() < 0.8:                         # noise on the copy as well
            m = random.randrange(3)
            x[p + m] += random.choice([-1, 1])
        if x[0] == 0 or x[-1] == 0 or any(abs(v) > 6 for v in x):
            continue
        # y[n] = sum_m x[n + a + m] t[m]; x occupies 0..L-1, so y lives on -(a+2)..L-1-a
        ns = list(range(-(a + 2), L - a))
        y = [sum((x[n + a + m] if 0 <= n + a + m < L else 0) * t[m] for m in range(3)) for n in ns]
        n0 = p - a
        y0 = y[ns.index(n0)]
        others = [v for n, v in zip(ns, y) if n != n0]
        if not y0 - max(others) >= 2:
            continue
        ew = sum(x[p + m] ** 2 for m in range(3))
        et = sum(v * v for v in t)
        if y0 <= 0:
            continue
        return t, a, x, c, p, ns, y, n0, y0, ew, et


def _stem_svg(x):
    """Stem plot of x[0..L-1] with the value at each stem end and the index row BELOW the whole plot
    (svg.stem puts index labels under the axis, where negative stems would cross them)."""
    L = len(x)
    vmax, vmin = max(max(x), 0), min(min(x), 0)
    unit = min(24.0, 130.0 / max(1, vmax - vmin))
    left, right, top = 46.0, 30.0, 30.0
    width = 560
    dx = (width - left - right) / (L + 1)
    y0 = top + vmax * unit
    ybot = y0 - vmin * unit
    height = int(ybot + 58)
    X = lambda n: left + (n + 1) * dx                      # noqa: E731
    Y = lambda v: y0 - v * unit                            # noqa: E731
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" '
           f'style="max-width:{width}px;display:block;margin:0.5em auto;color:inherit" font-family="sans-serif" '
           'font-size="13" role="img" aria-label="stem plot of the input x[n]">']
    out.append(f'<line x1="{left - 6:.1f}" y1="{y0:.1f}" x2="{width - right + 10:.1f}" y2="{y0:.1f}" '
               'stroke="currentColor" stroke-width="1.2"/>')
    out.append(f'<text x="{left - 10:.1f}" y="{y0 + 4:.1f}" text-anchor="end" font-style="italic">x[n]</text>')
    for n, v in enumerate(x):
        cx = X(n)
        if v == 0:
            out.append(f'<circle cx="{cx:.1f}" cy="{y0:.1f}" r="3.2" fill="{svg.ACCENT}"/>')
        else:
            out.append(f'<line x1="{cx:.1f}" y1="{y0:.1f}" x2="{cx:.1f}" y2="{Y(v):.1f}" stroke="{svg.ACCENT}" '
                       'stroke-width="2"/>')
            out.append(f'<circle cx="{cx:.1f}" cy="{Y(v):.1f}" r="4" fill="{svg.ACCENT}"/>')
        dy = -9 if v > 0 else (17 if v < 0 else -9)
        out.append(f'<text x="{cx:.1f}" y="{Y(v) + dy:.1f}" text-anchor="middle" font-size="12">'
                   f'{svg.fmt_num(v)}</text>')
        out.append(f'<text x="{cx:.1f}" y="{ybot + 42:.1f}" text-anchor="middle" fill="currentColor" '
                   f'fill-opacity="0.8">{n}</text>')
    out.append(f'<text x="{left - 10:.1f}" y="{ybot + 42:.1f}" text-anchor="end" font-style="italic">n</text>')
    out.append("</svg>")
    return "".join(out)


def generate(data):
    t, a, x, c, p, ns, y, n0, y0, ew, et = _draw()
    L = len(x)
    rho = y0 / math.sqrt(et * ew)
    # one other index, not the peak, with a nonzero window
    cand = [n for n in ns if n != n0 and abs(n - n0) <= 4]
    n1 = random.choice(cand)
    y1 = y[ns.index(n1)]

    hvals, hstart = t[::-1], -(a + 2)
    p_ = data["params"]
    p_["t"], p_["a"], p_["x"], p_["n0_peak"] = t, a, x, n0
    p_["t_tex"] = fmt.tex_seq(t, a)
    p_["x_tex"] = fmt.tex_seq(x, 0)
    p_["L"], p_["Lm1"] = L, L - 1
    p_["n1"] = n1
    p_["fig_html"] = _stem_svg(x)
    lo_m, hi_m = a, a + 2
    p_["win_lo_tex"] = "n" if a == 0 else (f"n + {a}" if a > 0 else f"n - {-a}")
    p_["win_hi_tex"] = f"n + {a + 2}" if a + 2 != 0 else "n"
    p_["h_tex"] = fmt.tex_seq(hvals, hstart)
    p_["hstart"] = hstart
    p_["et"], p_["ew"], p_["y0"], p_["y1"] = et, ew, y0, y1
    p_["et_tex"] = " + ".join(f"{fmt.tex_num(v, paren=True)}^2" for v in t) + f" = {et}"
    win = [x[p + m] for m in range(3)]
    p_["win_tex"] = r"\{" + r",\ ".join(str(v) for v in win) + r"\}"
    p_["ew_tex"] = " + ".join(f"{fmt.tex_num(v, paren=True)}^2" for v in win) + f" = {ew}"
    p_["k_lo"], p_["k_hi"] = p, p + 2
    p_["y0_terms_tex"] = " + ".join(f"{fmt.tex_num(x[p + m], paren=True)}\\cdot{fmt.tex_num(t[m], paren=True)}"
                                    for m in range(3))
    k1 = [k for k in range(n1 + lo_m, n1 + hi_m + 1)]
    p_["y1_terms_tex"] = " + ".join(
        f"{fmt.tex_num(x[k] if 0 <= k < L else 0, paren=True)}\\cdot{fmt.tex_num(t[k - n1 - a], paren=True)}"
        for k in k1)
    p_["y1_klo"], p_["y1_khi"] = k1[0], k1[-1]
    p_["k_minus_n1_tex"] = "k" if n1 == 0 else (f"k - {n1}" if n1 > 0 else f"k + {-n1}")
    root = math.isqrt(et * ew)
    if root * root == et * ew:
        tail = rf" = \dfrac{{{y0}}}{{{root}}}" + ("" if y0 % root else rf" = {y0 // root}")
        if y0 % root:
            tail += rf" \approx {rho:.4f}"
    else:
        tail = rf" \approx {rho:.4f}"
    p_["rho_tex"] = rf"\dfrac{{{y0}}}{{\sqrt{{{et}}}\,\sqrt{{{ew}}}}} = \dfrac{{{y0}}}{{\sqrt{{{et * ew}}}}}" + tail
    p_["exact_copy"] = (ew * et == y0 * y0)
    p_["noisy_copy"] = not p_["exact_copy"]
    p_["c"] = c
    p_["copy_tex"] = (f"{c}\\,t[n - {p - a}]" if c != 1 else f"t[n - {p - a}]") if p - a > 0 else \
        ((f"{c}\\,t[n + {a - p}]" if c != 1 else f"t[n + {a - p}]") if p - a < 0 else (f"{c}\\,t[n]" if c != 1 else "t[n]"))
    head = "".join(f"<th>${n}$</th>" for n in ns)
    cells = "".join((f"<td><b>${v}$</b></td>" if n == n0 else f"<td>${v}$</td>") for n, v in zip(ns, y))
    p_["ytable_html"] = ('<table class="table table-sm table-bordered text-center" style="width:auto">'
                         f"<tr><th>$n$</th>{head}</tr><tr><th>$y[n]$</th>{cells}</tr></table>")

    cr = data["correct_answers"]
    cr["hfirst"] = hstart
    cr["hrow"] = pl.to_json(np.array([[float(v) for v in hvals]]))
    cr["npeak"] = n0
    cr["ypeak"] = y0
    cr["yn1"] = y1
    cr["rho"] = float(rho)
