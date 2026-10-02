"""Response of an LTI system to a combination of shifted, scaled copies of a known input
(linearity + time invariance; HW2 #4 and SP2023 #4 use the same bookkeeping to build h from
input-output pairs; Lecture 3 superposition, Lecture 4 section 1.2).

The system maps x1 -> y1 (finite, shown). Input x2[n] = alpha x1[n-k] + beta x1[n-m]:
  template "formula": x2 is given by that formula (the copies may overlap);
  template "seq":     x2 is given only as a sequence; the copies do not overlap, so the student
                      first decomposes x2 into scaled, shifted copies of x1.
Student enters the first index of y2 (integer) and its samples from first to last nonzero (row)."""

import random

import numpy as np
import prairielearn as pl
from ece310 import convlib as cl
from ece310 import fmt, seq


def _rand_seq(length, lo, hi, ends):
    v = [random.randint(lo, hi) for _ in range(length)]
    v[0], v[-1] = random.choice(ends), random.choice(ends)
    return v


def _combo(vals, start, terms):
    """sum_j c_j vals[n - d_j] as (list, first index), exact; terms = [(c, d), ...]."""
    lo = min(start + d for _, d in terms)
    hi = max(start + d + len(vals) - 1 for _, d in terms)
    out = [0] * (hi - lo + 1)
    for c, d in terms:
        for i, v in enumerate(vals):
            out[start + d + i - lo] += c * v
    return out, lo


def generate(data):
    p = data["params"]
    tmpl = "formula" if random.random() < 0.5 else "seq"
    while True:
        l1 = random.choice([2, 3, 3])
        x1 = _rand_seq(l1, -2, 2, [-2, -1, 1, 1, 2])
        h = _rand_seq(random.choice([2, 2, 3]), -2, 2, [-2, -1, 1, 1, 2])
        n1 = random.randint(-1, 1)
        y1, ny1 = seq.conv(x1, n1, h, 0)
        y1 = [int(v) for v in y1]
        alpha, beta = random.choice([1, 2, -1, -2, 3]), random.choice([1, 2, -1, -2, 3])
        k = random.randint(-2, 2)
        if tmpl == "seq":
            m = k + random.choice([l1, l1 + 1, l1 + 2])        # non-overlapping copies, visible by inspection
        else:
            m = k + random.choice([1, 2, 3, -1, -2])
        if alpha == 1 and beta == 1:
            continue
        x2, nx2 = _combo(x1, n1, [(alpha, k), (beta, m)])
        y2, ny2 = _combo(y1, ny1, [(alpha, k), (beta, m)])
        if max(abs(v) for v in y2) <= 20 and max(abs(v) for v in y1) <= 12 and abs(nx2) <= 4 and len(y2) <= 9:
            break
    terms = sorted([(alpha, k), (beta, m)], key=lambda t: t[1])
    p.update(template=tmpl, x1=x1, n1=n1, y1=y1, ny1=ny1, x2=x2, nx2=nx2,
             alpha=alpha, beta=beta, k=k, m=m)
    p["is_formula"], p["is_seq"] = tmpl == "formula", tmpl == "seq"
    p["x1_tex"] = fmt.tex_seq(x1, n1)
    p["y1_tex"] = fmt.tex_seq(y1, ny1)
    p["x2_tex"] = fmt.tex_seq(x2, nx2)
    p["x2f_tex"] = cl.tex_signed_terms([(c, f"x_1[{cl.tex_arg(d)}]") for c, d in terms])
    p["y2f_tex"] = cl.tex_signed_terms([(c, f"y_1[{cl.tex_arg(d)}]") for c, d in terms])
    # decomposition table for x2 (seq template) and the superposition table for y2
    xrows = [(cl.tex_coef(c, f"x_1[{cl.tex_arg(d)}]"), [c * v for v in x1], n1 + d) for c, d in terms]
    p["xtable_html"] = cl.shift_add_table_html(xrows, "x_2[n]", x2, nx2)
    yrows = [(cl.tex_coef(c, f"y_1[{cl.tex_arg(d)}]"), [c * v for v in y1], ny1 + d) for c, d in terms]
    p["ytable_html"] = cl.shift_add_table_html(yrows, "y_2[n]", y2, ny2)
    p["y2_tex"] = fmt.tex_seq(y2, ny2)
    p["start_tex"] = f"{ny1} + ({terms[0][1]})"
    p["y2_last"] = ny2 + len(y2) - 1
    data["correct_answers"]["ny2"] = ny2
    data["correct_answers"]["y2"] = pl.to_json(np.array([[float(v) for v in y2]]))
