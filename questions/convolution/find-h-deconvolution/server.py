"""Find h[n] from one input-output pair of an LTI system (exam family "finding h from input-output
pairs", on all 7 past Midterm 1 exams).

Template "fir"  (SP2025 #3, FA2024 #3, FA2019 #3): finite x (2-3 samples) and finite y = x * h for
                a hidden finite h (2-4 samples, nonzero ends); deconvolve (long division / peeling).
Template "geo"  (FA2023 #4; HW2 #4 / SP2023 #4 "build a delta"): x = a^n u[n] (a = 1: the unit step)
                and a finite y; x[n] - a x[n-1] = delta[n], so h[n] = y[n] - a y[n-1].
Student enters the first index of h (integer) and its samples from first to last nonzero (row)."""

import random

import numpy as np
import prairielearn as pl
from ece310 import convlib as cl
from ece310 import fmt, seq


def _rand_seq(length, lo, hi, ends):
    v = [random.randint(lo, hi) for _ in range(length)]
    v[0], v[-1] = random.choice(ends), random.choice(ends)
    return v


def _peel_tex(x, nx, y, ny, h, nh):
    """Lines of the long division: h[m] = (y[ny+j] - sum_{i>=1} x[nx+i] h[m-i]) / x[nx]."""
    lines = []
    for j, hv in enumerate(h):
        m = nh + j
        terms = [(-x[i], f"h[{m - i}]", h[j - i]) for i in range(1, len(x)) if 0 <= j - i < len(h) and x[i] != 0]
        yv = seq.value(y, ny, ny + j)
        sym = f"y[{ny + j}]" + "".join((" - " if c < 0 else " + ") + cl.tex_coef(abs(c), lab) for c, lab, _ in terms)
        num = fmt.tex_num(yv) + "".join(
            (" - " if c < 0 else " + ") + ((fmt.tex_num(abs(c)) + r"\cdot ") if abs(c) != 1 else "")
            + fmt.tex_num(val, paren=True) for c, lab, val in terms)
        if x[0] == 1:
            lines.append(rf"h[{m}] = {sym}" + (f" = {num}" if terms else "") + f" = {fmt.tex_num(hv)}")
        else:
            lines.append(rf"h[{m}] = \dfrac{{{sym}}}{{{fmt.tex_num(x[0])}}} = \dfrac{{{num}}}{{{fmt.tex_num(x[0])}}}"
                         f" = {fmt.tex_num(hv)}")
    return r"\begin{gathered}" + r"\\ ".join(lines) + r"\end{gathered}"


def _gen_fir(p):
    while True:
        lx, lh = random.choice([(2, 2), (2, 3), (2, 4), (3, 2), (3, 3)])
        x = _rand_seq(lx, -2, 2, [-2, -1, 1, 1, 2])
        h = _rand_seq(lh, -3, 3, [-3, -2, -1, 1, 2, 3])
        if h == h[::-1] and random.random() < 0.7:      # mostly non-palindromic h (the reversal probe bites)
            continue
        nx, nh = random.randint(-2, 2), random.randint(-2, 2)
        y, ny = seq.conv(x, nx, h, nh)
        if max(abs(v) for v in y) <= 12:
            break
    y = [int(v) for v in y]
    p.update(template="fir", x=x, nx=nx, y=y, ny=ny)
    p["x_tex"] = fmt.tex_seq(x, nx)
    p["y_tex"] = fmt.tex_seq(y, ny)
    p["nx_last"], p["ny_last"] = nx + lx - 1, ny + len(y) - 1
    p["lx"], p["ly"], p["lh"] = lx, len(y), lh
    p["nh_tex"] = f"{ny} - ({nx})"
    p["peel_tex"] = _peel_tex(x, nx, y, ny, h, nh)
    rows = [(rf"h[{nh + j}]\,x[{cl.tex_arg(nh + j)}]", [hv * xv for xv in x], nx + nh + j)
            for j, hv in enumerate(h) if hv != 0]
    p["table_html"] = cl.shift_add_table_html(rows, r"x * h", y, ny)
    p["h_tex"] = fmt.tex_seq(h, nh)
    return h, nh


def _gen_geo(p):
    a = random.choice([1, 1, 1, -1, 2, -2, 3])
    while True:
        ly = random.choice([2, 2, 3])
        y = _rand_seq(ly, -3, 3, [-3, -2, -1, 1, 2, 3])
        ny = random.randint(-1, 2)
        h = [y[0]] + [y[j] - a * y[j - 1] for j in range(1, ly)] + [-a * y[-1]]
        if max(abs(v) for v in h) <= 12 and any(v != y[0] for v in y[1:] + [0]):
            break
    nh = ny
    p.update(template="geo", a=a, y=y, ny=ny)
    p["x_tex"] = "u[n]" if a == 1 else cl.tex_coef(1, fmt.tex_pow(a) + r"\,u[n]")
    p["y_tex"] = fmt.tex_seq(y, ny)
    p["a_tex"] = fmt.tex_num(a, paren=True)
    p["xa_tex"] = "u[n]" if a == 1 else fmt.tex_pow(a) + r"\,u[n]"
    p["xa1_tex"] = "u[n-1]" if a == 1 else fmt.tex_pow(a, "n-1") + r"\,u[n-1]"
    p["comb_tex"] = cl.tex_signed_terms([(1, "x[n]"), (-a, "x[n-1]")])
    p["hcomb_tex"] = cl.tex_signed_terms([(1, "y[n]"), (-a, "y[n-1]")])
    p["ya_coef_tex"] = cl.tex_coef(-a, "y[n-1]")
    rows = [("y[n]", y, ny), (cl.tex_coef(-a, "y[n-1]"), [-a * v for v in y], ny + 1)]
    p["table_html"] = cl.shift_add_table_html(rows, "h[n]", h, nh)
    p["h_tex"] = fmt.tex_seq(h, nh)
    p["is_step"] = a == 1
    return h, nh


_DEFAULTS = {"x": [], "nx": 0, "nx_last": 0, "ny_last": 0, "lx": 0, "ly": 0, "lh": 0, "nh_tex": "", "peel_tex": "",
             "a": 0, "a_tex": "", "xa_tex": "", "xa1_tex": "", "comb_tex": "", "hcomb_tex": "", "ya_coef_tex": "",
             "is_step": False}


def generate(data):
    p = data["params"]
    p.update(_DEFAULTS)
    h, nh = _gen_fir(p) if random.random() < 0.6 else _gen_geo(p)
    p["is_fir"] = p["template"] == "fir"
    p["is_geo"] = p["template"] == "geo"
    p["h_len"] = len(h)
    data["correct_answers"]["nh"] = nh
    data["correct_answers"]["h"] = pl.to_json(np.array([[float(v) for v in h]]))
