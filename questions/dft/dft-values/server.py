"""DFT values from the definition (Midterm 2 family "DFT computation and properties", 7 of 7 past
exams): SP2023 MT2 #9(a) (X_0 and X_3 of a length-6 signal), SP2025 #7(a) (X_0), FA2021 #3(a)
(X[k] for k = 0, 2, 4, 6), FA2019 #1(i) / SP2021 #1(e) (X[0] = 0?), SP2021 #5 (x[0] from the DFT
values) and SP2023 #1(d) / SP2025 #1(b) (the inverse DFT). Course notes 2019, Sec. 2.5 (definition)
and 2.5.1.9 (Parseval).

One variant: an integer x of length N in {4, 6, 8}; asked
  (a) X[0] = sum x[n] and X[N/2] = sum (-1)^n x[n];
  (b) Re and Im of X[k0] for one other k0 with Im X[k0] != 0 (so the sign of the exponent matters;
      for N = 6 and odd k0 with N = 8 the values are irrational, entered as decimals);
  (c) sum_k |X[k]|^2 = N sum_n |x[n]|^2 (Parseval);
  (d) the inverse direction, one of two templates:
      x0  : v[0] = (1/N') sum_k V[k] from a given Hermitian list of Gaussian integers (N' = 5..8);
      vec : the whole length-4 v[n] (row vector) from its 4-point DFT V[k] (Hermitian, so v is real).
All values are exact (ece310.w8_dft: arithmetic in Q(sqrt2), Q(sqrt3))."""

import random
from fractions import Fraction as Fr

import numpy as np
import prairielearn as pl
from ece310 import fmt, w8_dft as W

def _signed_terms(vals):
    """'1 - 2 + 3' from a list of ints (each as a signed term, negatives in parentheses after +/-)."""
    out = ""
    for v in vals:
        if not out:
            out = fmt.tex_num(v)
        else:
            out += (" - " if v < 0 else " + ") + fmt.tex_num(abs(v))
    return out or "0"


def _alt_terms(x):
    """x[0] - x[1] + x[2] - ... written with the sign of (-1)^n kept visible."""
    out = []
    for n, v in enumerate(x):
        t = fmt.tex_num(v, paren=True)
        out.append(t if n == 0 else ((" - " if n % 2 else " + ") + t))
    return "".join(out)


def _w_tex(N, k0):
    """e^{-j 2 pi k0 n / N} without a written-out factor 1."""
    num = r"2\pi" if k0 == 1 else rf"2\pi\cdot {k0}"
    return rf"e^{{-j\frac{{{num}}}{{{N}}}n}}"


def _bin_table(x, k0):
    N = len(x)
    rows = []
    for n, v in enumerate(x):
        c, s = W.unit_exact(N, k0 * n)
        e_tex = rf"e^{{-j\frac{{2\pi\cdot{k0 * n}}}{{{N}}}}}"
        if v == 0:
            prod = "0"
        else:
            prod = W.tex_cplx(c.scale(v), s.scale(v))
        rows.append(f"<tr><td>${n}$</td><td>${fmt.tex_num(v)}$</td><td>${e_tex}$</td>"
                    f"<td>$= {W.tex_cplx(c, s)}$</td><td>${prod}$</td></tr>")
    w = _w_tex(N, k0)
    head = (f"<tr><th>$n$</th><th>$x[n]$</th><th colspan='2'>${w}$</th><th>$x[n]\\,{w}$</th></tr>")
    return ('<table class="table table-sm table-bordered text-center" style="width:auto">'
            + head + "".join(rows) + "</table>")


def _pick_x(N):
    while True:
        x = [random.randint(-3, 4) for _ in range(N)]
        if sum(1 for v in x if v == 0) > 2:
            continue
        s0 = sum(x)
        sa = sum(v * (-1) ** n for n, v in enumerate(x))
        if sum(abs(v) for v in x) == s0 or s0 == sa:
            continue
        odd = random.random() < 0.5
        cands = []
        for k in range(1, N):
            if 2 * k == N:
                continue
            re, im = W.dft_exact(x, k)
            if im.is_zero() or abs(re.value() - im.value()) < 1e-9:
                continue
            if any(part.b != 0 and abs(part.value()) < 0.3 * abs(float(part.b)) for part in (re, im)):
                continue                                    # e.g. -7 + 5 sqrt2 = 0.071: decimals would be unfair
            cands.append(k)
        if N == 8:
            pref = [k for k in cands if (k % 2 == 1) == odd]
            cands = pref or cands
        if cands:
            return x, s0, sa, random.choice(cands)


def _inverse_vec(p, c):
    while True:
        v = [random.randint(-3, 4) for _ in range(4)]
        if v[1] != v[3] and sum(1 for t in v if t == 0) <= 1 and len(set(v)) >= 3:
            break
    a, b = v[0] - v[2], v[1] - v[3]
    V = [(sum(v), 0), (a, -b), (v[0] - v[1] + v[2] - v[3], 0), (a, b)]
    p["V"] = [list(t) for t in V]
    p["Nv"] = 4
    p["V_tex"] = W.tex_list([W.tex_cplx(*t) for t in V])
    p["V0_tex"], p["V2_tex"] = fmt.tex_num(V[0][0]), fmt.tex_num(V[2][0])
    p["V1_tex"] = W.tex_cplx(*V[1])
    p["V3_tex"] = W.tex_cplx(*V[3])
    p["reV1_tex"] = fmt.tex_num(2 * a)
    p["imV1_tex"] = fmt.tex_num(-2 * b)
    s02, d02 = V[0][0] + V[2][0], V[0][0] - V[2][0]
    p["inv_lines"] = [
        {"line_tex": rf"v[0] = \tfrac14\left[(V_0 + V_2) + (V_1 + V_3)\right] = \tfrac14\left[{fmt.tex_num(s02)} + {fmt.tex_num(2 * a, paren=True)}\right] = {v[0]}"},
        {"line_tex": rf"v[1] = \tfrac14\left[(V_0 - V_2) + j(V_1 - V_3)\right] = \tfrac14\left[{fmt.tex_num(d02)} + {fmt.tex_num(2 * b, paren=True)}\right] = {v[1]}"},
        {"line_tex": rf"v[2] = \tfrac14\left[(V_0 + V_2) - (V_1 + V_3)\right] = \tfrac14\left[{fmt.tex_num(s02)} - {fmt.tex_num(2 * a, paren=True)}\right] = {v[2]}"},
        {"line_tex": rf"v[3] = \tfrac14\left[(V_0 - V_2) - j(V_1 - V_3)\right] = \tfrac14\left[{fmt.tex_num(d02)} - {fmt.tex_num(2 * b, paren=True)}\right] = {v[3]}"},
    ]
    p["jdiff_tex"] = fmt.tex_num(2 * b)       # j (V1 - V3) = j * (2j * Im V1) = -2 Im V1 = 2b
    p["v_tex"] = W.tex_list([str(t) for t in v])
    c["v"] = pl.to_json(np.array([[float(t) for t in v]]))


def _inverse_x0(p, c):
    while True:
        N2 = random.choice([5, 6, 7, 8])
        V = [(random.randint(-4, 6), 0)]
        half = []
        for _ in range(1, (N2 - 1) // 2 + 1):
            half.append((random.randint(-3, 4), random.randint(-4, 4)))
        if not any(b for _, b in half) or all(a == 0 for a, _ in half):
            continue
        V += half
        if N2 % 2 == 0:
            V.append((random.randint(-4, 4), 0))
        V += [(a, -b) for a, b in reversed(half)]
        total = sum(a for a, _ in V)
        if total == 0 or total == V[0][0] or any(a == 0 and b == 0 for a, b in V):
            continue
        break
    v0 = Fr(total, N2)
    p["V"] = [list(t) for t in V]
    p["Nv"] = N2
    p["Nvm1"] = N2 - 1
    p["V_tex"] = W.tex_list([W.tex_cplx(*t) for t in V])
    p["Vsum_terms_tex"] = " + ".join(W.tex_cplx(*t) if i == 0 else r"\left(" + W.tex_cplx(*t) + r"\right)"
                                     for i, t in enumerate(V))
    p["Vsum_tex"] = fmt.tex_num(total)
    p["v0_tex"] = fmt.tex_num(v0)
    p["v0_eq_tex"] = "" if (v0.numerator, v0.denominator) == (total, N2) else " = " + fmt.tex_num(v0)
    p["v0_dec"] = f"{float(v0):.6g}"
    c["v0"] = float(v0)


def generate(data):
    p, c = data["params"], data["correct_answers"]
    N = random.choice([4, 6, 8])
    x, s0, sa, k0 = _pick_x(N)
    re, im = W.dft_exact(x, k0)
    E = sum(v * v for v in x)

    p["N"], p["Nm1"], p["half"], p["k0"] = N, N - 1, N // 2, k0
    p["x"] = x
    p["x_tex"] = W.tex_list([fmt.tex_num(v) for v in x])
    p["X0_terms_tex"] = _signed_terms(x)
    p["XH_terms_tex"] = _alt_terms(x)
    p["X0_tex"], p["XH_tex"] = fmt.tex_num(s0), fmt.tex_num(sa)
    p["table_html"] = _bin_table(x, k0)
    p["re_tex"], p["im_tex"] = re.tex(), im.tex()
    p["Xk_tex"] = W.tex_cplx(re, im)
    p["re_irr"], p["im_irr"] = not re.is_rational(), not im.is_rational()
    p["re_dec"], p["im_dec"] = f"{re.value():.6g}", f"{im.value():.6g}"
    p["W_tex"] = _w_tex(N, k0)
    p["sq_terms_tex"] = " + ".join(f"{fmt.tex_num(v, paren=True)}^2" for v in x if v != 0)
    p["E_tex"], p["P_tex"] = str(E), str(N * E)
    kind = random.choice(["x0", "vec"])
    p["inv_x0"], p["inv_vec"] = kind == "x0", kind == "vec"
    p["inv_kind"] = kind

    c["X0"] = float(s0)
    c["XH"] = float(sa)
    c["re"] = re.value()
    c["im"] = im.value()
    c["P"] = float(N * E)
    if kind == "vec":
        _inverse_vec(p, c)
    else:
        _inverse_x0(p, c)
