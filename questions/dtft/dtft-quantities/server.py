"""Four DTFT quantities of a short finite sequence, read off without computing X_d(w)
(HW5 Problem 2, Fall 2026; Lecture 13 DTFT pair, Lecture 14 Parseval's relation).

    X_d(0)               = sum_n x[n]
    X_d(pi)              = sum_n (-1)^n x[n]
    int_{-pi}^{pi} X_d   = 2 pi x[0]               (inverse DTFT at n = 0)
    int_{-pi}^{pi} |X_d|^2 = 2 pi sum_n |x[n]|^2   (Parseval)

The integrals are entered as multiples of pi. Design: 3-6 integer samples, the support contains
n = 0 with x[0] != 0 (the integral is not trivially 0), nonzero end samples, and the odd-indexed
samples do not sum to zero, so X_d(pi) != X_d(0) (forgetting (-1)^n is always caught)."""

import random

from ece310 import fmt


def _signed_terms_tex(vals):
    """'2 - 1 + 1' for the list [2, -1, 1] (zeros skipped); '0' if all are zero."""
    return fmt.tex_sum((v, "") for v in vals)


def _table_html(x, start):
    ns = list(range(start, start + len(x)))
    head = "".join(f"<th>${n}$</th>" for n in ns)
    r1 = "".join(f"<td>${fmt.tex_num(v)}$</td>" for v in x)
    r2 = "".join(f"<td>${fmt.tex_num((-1) ** (n % 2) * v)}$</td>" for n, v in zip(ns, x))
    r3 = "".join(f"<td>${fmt.tex_num(v * v)}$</td>" for v in x)
    return ('<table class="table table-sm table-bordered text-center" style="width:auto">'
            f"<tr><th>$n$</th>{head}</tr>"
            f"<tr><th>$x[n]$</th>{r1}</tr>"
            f"<tr><th>$(-1)^n\\,x[n]$</th>{r2}</tr>"
            f"<tr><th>$|x[n]|^2$</th>{r3}</tr></table>")


def generate(data):
    while True:
        L = random.randint(3, 6)
        start = random.randint(-(L - 1), 0)       # the support always contains n = 0
        x = [random.randint(-3, 3) for _ in range(L)]
        x[0] = random.choice([-3, -2, -1, 1, 2, 3])
        x[-1] = random.choice([-3, -2, -1, 1, 2, 3])
        ns = list(range(start, start + L))
        x0 = x[ns.index(0)]
        odd_sum = sum(v for n, v in zip(ns, x) if n % 2)
        if x0 == 0 or odd_sum == 0:
            continue
        break

    X0 = sum(x)
    Xpi = sum((-1) ** (n % 2) * v for n, v in zip(ns, x))
    energy = sum(v * v for v in x)

    p = data["params"]
    p["x"], p["start"] = x, start
    p["x_tex"] = fmt.tex_seq(x, start)
    p["table_html"] = _table_html(x, start)
    p["X0_terms_tex"] = _signed_terms_tex(x)
    p["Xpi_terms_tex"] = _signed_terms_tex([(-1) ** (n % 2) * v for n, v in zip(ns, x)])
    p["energy_terms_tex"] = _signed_terms_tex([v * v for v in x])
    p["X0_tex"] = fmt.tex_num(X0)
    p["Xpi_tex"] = fmt.tex_num(Xpi)
    p["x0_tex"] = fmt.tex_num(x0, paren=True)
    p["energy_tex"] = fmt.tex_num(energy)
    p["intX_tex"] = fmt.tex_num(2 * x0)
    p["intX2_tex"] = fmt.tex_num(2 * energy)

    c = data["correct_answers"]
    c["X0"] = float(X0)
    c["Xpi"] = float(Xpi)
    c["intX"] = float(2 * x0)
    c["intX2"] = float(2 * energy)
