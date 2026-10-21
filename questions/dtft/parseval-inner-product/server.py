"""Generalized Parseval and the inverse DTFT at one n, for two short real sequences (Quiz 2, Oct 2026:
"Parseval without squares"; Lecture 13 inverse DTFT, Lecture 14 Parseval).

    (a) (1/2pi) int X_d Y_d^* dw   = sum_n x[n] y[n]            (y real)
    (b) int X_d e^{j w n0} dw      = 2 pi x[n0]                  (entered as a multiple of pi)
    (c) (1/2pi) int X_d Y_d dw     = sum_n x[n] y[-n] = (x*y)[0] (no conjugate: the trap)
    (d) int |X_d|^2 dw             = 2 pi sum_n x[n]^2           (entered as a multiple of pi)

Design: lengths 3-4, integers in [-3, 3], nonzero end samples, random start indices; (a) and (c)
are nonzero and differ; n0 != 0 lies in the support with x[n0] != 0, x[-n0] != x[n0] and
x[0] != x[n0] (the sign of n0 and the plain n = 0 shortcut are both caught)."""

import random

from ece310 import fmt

NZ = [-3, -2, -1, 1, 2, 3]


def _draw():
    L = random.randint(3, 4)
    s = random.randint(-2, 1)
    v = [random.randint(-3, 3) for _ in range(L)]
    v[0], v[-1] = random.choice(NZ), random.choice(NZ)
    return v, s


def _ejn(n0):
    """e^{j n0 w} in LaTeX."""
    k = "" if abs(n0) == 1 else str(abs(n0))
    return "e^{" + ("-" if n0 < 0 else "") + "j" + k + r"\omega}"


def _prod_terms(pairs):
    """'2\\cdot(-1) + 3\\cdot 1' over the pairs with both factors nonzero; '0' if none."""
    ts = [f"{fmt.tex_num(a, paren=True)}\\cdot {fmt.tex_num(b, paren=True)}" for a, b in pairs if a and b]
    return " + ".join(ts) if ts else "0"


def _table_html(dx, dy, ns):
    def row(head, vals):
        return f"<tr><th>{head}</th>" + "".join(f"<td>${fmt.tex_num(v)}$</td>" for v in vals) + "</tr>"
    hdr = "<tr><th>$n$</th>" + "".join(f"<th>${n}$</th>" for n in ns) + "</tr>"
    return ('<table class="table table-sm table-bordered text-center" style="width:auto">' + hdr
            + row("$x[n]$", [dx.get(n, 0) for n in ns])
            + row("$y[n]$", [dy.get(n, 0) for n in ns])
            + row("$y[-n]$", [dy.get(-n, 0) for n in ns])
            + row("$x[n]\\,y[n]$", [dx.get(n, 0) * dy.get(n, 0) for n in ns])
            + row("$x[n]\\,y[-n]$", [dx.get(n, 0) * dy.get(-n, 0) for n in ns])
            + "</table>")


def generate(data):
    while True:
        x, sx = _draw()
        y, sy = _draw()
        if x == y and sx == sy:
            continue
        dx = {sx + i: v for i, v in enumerate(x)}
        dy = {sy + i: v for i, v in enumerate(y)}
        A = sum(v * dy.get(n, 0) for n, v in dx.items())
        C = sum(v * dy.get(-n, 0) for n, v in dx.items())
        if A == 0 or C == 0 or A == C:
            continue
        nA = sum(1 for n, v in dx.items() if v and dy.get(n, 0))
        nC = sum(1 for n, v in dx.items() if v and dy.get(-n, 0))
        if nA < 2 or nA + nC < 4:                       # real overlaps, not a single product
            continue
        cand = [n for n, v in dx.items() if n != 0 and v != 0 and dx.get(-n, 0) != v and dx.get(0, 0) != v]
        if not cand:
            continue
        n0 = random.choice(cand)
        break

    E = sum(v * v for v in x)
    lo = min(sx, sy, -(sy + len(y) - 1), 0)
    hi = max(sx + len(x) - 1, sy + len(y) - 1, -sy, 0)
    ns = [n for n in range(lo, hi + 1)
          if dx.get(n, 0) or dy.get(n, 0) or dy.get(-n, 0) or n == 0]

    p = data["params"]
    p["x"], p["sx"], p["y"], p["sy"], p["n0"] = x, sx, y, sy, n0
    p["x_tex"] = fmt.tex_seq(x, sx)
    p["y_tex"] = fmt.tex_seq(y, sy)
    p["ejn0_tex"] = _ejn(n0)
    p["n0_tex"] = str(n0)
    p["table_html"] = _table_html(dx, dy, ns)
    p["A_terms_tex"] = _prod_terms([(dx[n], dy.get(n, 0)) for n in sorted(dx)])
    p["C_terms_tex"] = _prod_terms([(dx[n], dy.get(-n, 0)) for n in sorted(dx)])
    p["E_terms_tex"] = fmt.tex_sum((v * v, "") for v in x)
    p["A_tex"], p["C_tex"], p["E_tex"] = fmt.tex_num(A), fmt.tex_num(C), fmt.tex_num(E)
    p["xn0_tex"] = fmt.tex_num(dx[n0], paren=True)
    p["B_tex"] = fmt.tex_num(2 * dx[n0])
    p["D_tex"] = fmt.tex_num(2 * E)

    c = data["correct_answers"]
    c["A"] = float(A)
    c["B"] = float(2 * dx[n0])
    c["C"] = float(C)
    c["D"] = float(2 * E)
