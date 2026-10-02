"""Finite-length convolution y = x * h of two short integer sequences with arbitrary start
indices (exam family "finite-length convolution", 7 of 7 past Midterm 1 exams; HW2 #5a).

Student enters the start index of y and the samples of y from first to last nonzero.
The first and last samples of x and h are nonzero, so y has exactly len(x)+len(h)-1 samples with
nonzero ends (no ambiguity about leading/trailing zeros)."""

import random

import numpy as np
import prairielearn as pl
from ece310 import fmt, seq


def _rand_seq(length):
    vals = [random.randint(-3, 3) for _ in range(length)]
    vals[0] = random.choice([-3, -2, -1, 1, 2, 3])
    vals[-1] = random.choice([-3, -2, -1, 1, 2, 3])
    return vals


def _table_html(x, nx, h, nh, y, ny):
    """Shift-and-add table: one row x[k] h[n-k] per input sample, then the column sums."""
    ns = list(range(ny, ny + len(y)))
    head = "".join(f"<th>$n={n}$</th>" for n in ns)
    rows = []
    for i, xv in enumerate(x):
        k = nx + i
        if xv == 0:
            continue
        cells = []
        for n in ns:
            j = n - k - nh
            cells.append(f"<td>${fmt.tex_num(xv * h[j])}$</td>" if 0 <= j < len(h) else "<td></td>")
        rows.append(f"<tr><th>$x[{k}]\\,h[n{'-' if k > 0 else '+'}{abs(k)}]$</th>{''.join(cells)}</tr>"
                    if k != 0 else f"<tr><th>$x[0]\\,h[n]$</th>{''.join(cells)}</tr>")
    total = "".join(f"<td><b>${fmt.tex_num(v)}$</b></td>" for v in y)
    return ('<table class="table table-sm table-bordered text-center" style="width:auto">'
            f"<tr><th></th>{head}</tr>{''.join(rows)}<tr><th>$y[n]$</th>{total}</tr></table>")


def generate(data):
    lx, lh = random.choice([(3, 2), (3, 3), (4, 2), (4, 3), (5, 2), (5, 3)])
    x, h = _rand_seq(lx), _rand_seq(lh)
    nx, nh = random.randint(-2, 2), random.randint(-2, 1)
    y, ny = seq.conv(x, nx, h, nh)

    p = data["params"]
    p["x"], p["nx"], p["h"], p["nh"] = x, nx, h, nh
    p["x_tex"] = fmt.tex_seq(x, nx)
    p["h_tex"] = fmt.tex_seq(h, nh)
    p["y_tex"] = fmt.tex_seq(y, ny)
    p["nx_last"], p["nh_last"] = nx + lx - 1, nh + lh - 1
    p["ny_last"] = ny + len(y) - 1
    p["ly"] = len(y)
    p["start_tex"] = f"{nx} {'-' if nh < 0 else '+'} {abs(nh)}"
    p["sumcheck_tex"] = fr"{fmt.tex_num(sum(x), paren=True)}\cdot{fmt.tex_num(sum(h), paren=True)} = {fmt.tex_num(sum(y))}"
    p["table_html"] = _table_html(x, nx, h, nh, y, ny)

    data["correct_answers"]["ny"] = ny
    data["correct_answers"]["y"] = pl.to_json(np.array([[float(v) for v in y]]))
