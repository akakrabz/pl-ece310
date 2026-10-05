"""Improper causal H(z): numerator degree (in z^-1) >= denominator degree, 1 or 2 real poles.
Find the direct (polynomial) part C0 (, C1) by long division, the partial-fraction residues of
the remainder, and h[0..3].

Model: Lecture 10 Exercise 1 (long division in z^-1, then PFE of the remainder), HW4 #1(b)
(an improper transform), Lecture 8 (PFE).

Built backwards from nice parts: H = sum_k C_k z^-k + sum_k A_k / (1 - p_k z^-1), so
numerator = C(z^-1) D(z^-1) + R(z^-1) with D = prod (1 - p_k z^-1)."""

import random
from fractions import Fraction as F

from ece310 import fmt, poly, zinv

CS = [1, 2, 3, -1, -2, -3, 4, -4, F(1, 2), F(-1, 2), F(3, 2), F(-3, 2)]
RES = [1, 2, 3, -1, -2, -3, F(1, 2), F(-1, 2), F(3, 2), F(-3, 2)]


def _remainder(poles, A):
    rem = [F(0)]
    for k, a in enumerate(A):
        t = [F(1)]
        for j, q in enumerate(poles):
            if j != k:
                t = poly.pmul(t, [F(1), -q])
        rem = poly.padd(rem, poly.pscale(t, a))
    return poly.trim(rem)


def _pick():
    while True:
        m = random.choice([1, 2, 2])
        L = random.choice([1, 2])
        poles = sorted(random.sample(zinv.NICE_POLES, m), key=lambda q: (abs(q), q))
        C = [F(random.choice(CS)) for _ in range(L)]
        A = [F(random.choice(RES)) for _ in range(m)]
        if m == 2 and A[0] == A[1]:
            continue
        if sum(A) == 0:
            continue          # then h[0] = C0 and the "front division" mistake is invisible
        D = poly.from_roots(poles)
        R = _remainder(poles, A)
        num = poly.padd(poly.pmul(C, D), R)
        if len(num) != m + L:
            continue
        if not all(zinv.nice(x, maxden=12, maxnum=30, bad_dens=(5, 7, 10, 11)) for x in num + D):
            continue
        h = [(C[n] if n < L else 0) + sum(a * q**n for a, q in zip(A, poles)) for n in range(4)]
        if not all(zinv.nice(v, maxden=36, maxnum=150, bad_dens=()) for v in h):
            continue
        return poles, C, A, D, R, num, h


def _division_steps(num, den):
    """LaTeX lines of the long division in w = z^-1, highest power first (Lecture 10)."""
    r = list(num)
    lines = []
    dm = len(den) - 1
    for i in range(len(num) - len(den), -1, -1):
        lead = r[i + dm]
        q = lead / den[-1]
        sub = [F(0)] * i + [q * d for d in den]
        newr = poly.padd(r, poly.pscale(sub, -1))
        newr = (newr + [F(0)] * len(r))[: i + dm]
        lines.append(
            r"\dfrac{" + fmt.tex_sum([(lead, fmt.zpow_tex(i + dm))]) + "}{" + fmt.tex_sum([(den[-1], fmt.zpow_tex(dm))])
            + "} = " + fmt.tex_sum([(q, fmt.zpow_tex(i))]) + r":\quad "
            + r"\text{subtract } " + fmt.tex_sum([(q, fmt.zpow_tex(i))]) + r"\cdot D \;\Rightarrow\; \text{remainder } "
            + (fmt.tex_poly_zinv(newr) if any(x != 0 for x in newr) else "0"))
        r = newr
    return lines


def generate(data):
    poles, C, A, D, R, num, h = _pick()
    m, L = len(poles), len(C)
    p = data["params"]
    p["num"] = [str(x) for x in num]
    p["poles"] = [str(q) for q in poles]
    p["L"] = L
    p["has_C1"] = L == 2
    p["two_poles"] = m == 2
    p["one_pole"] = m == 1
    # one pole: the lone factor without parentheses; two poles: the product of the factors
    den_disp = fmt.tex_poly_zinv(D) if m == 1 else "".join(fmt.tex_factor(q) for q in poles)
    p["H_tex"] = fmt.tex_frac(fmt.tex_poly_zinv(num), den_disp)
    terms = ["C_0"] + (["C_1 z^{-1}"] if L == 2 else [])
    for k in range(m):
        terms.append(rf"\dfrac{{A_{k + 1}}}{{1 - p_{k + 1} z^{{-1}}}}")
    p["form_tex"] = " + ".join(terms)
    p["poles_named_tex"] = r",\ ".join(f"p_{k + 1} = {fmt.tex_num(q)}" for k, q in enumerate(poles))
    p["p1_tex"] = fmt.tex_num(poles[0])
    p["p2_tex"] = fmt.tex_num(poles[1]) if m == 2 else ""
    p["roc_tex"] = r"\lvert z\rvert \gt " + fmt.tex_num(max(abs(q) for q in poles))
    p["deg_num"], p["deg_den"] = len(num) - 1, m
    p["D_tex"] = fmt.tex_poly_zinv(D)
    p["Dfact_tex"] = den_disp
    p["D_line_tex"] = fmt.tex_poly_zinv(D) if m == 1 else den_disp + " = " + fmt.tex_poly_zinv(D)
    p["div_steps"] = [{"tex": s} for s in _division_steps(num, D)]
    p["Q_tex"] = fmt.tex_poly_zinv(C)
    p["R_tex"] = fmt.tex_poly_zinv(R)
    p["split_tex"] = fmt.tex_poly_zinv(C) + " + " + fmt.tex_frac(fmt.tex_poly_zinv(R), p["Dfact_tex"])
    p["cover_tex"] = r" \\[2pt] ".join(zinv.coverup_tex(R, poles, k, f"A_{k + 1}") for k in range(m)) if m == 2 else ""
    p["R0_tex"] = fmt.tex_num(R[0])
    p["R0_frac_tex"] = fmt.tex_frac(fmt.tex_num(R[0]), fmt.tex_poly_zinv(D))
    full = [(c, fmt.zpow_tex(k)) for k, c in enumerate(C)]
    full += [(a, fmt.tex_frac("1", fmt.tex_sum([(1, ""), (-q, "z^{-1}")]))) for a, q in zip(A, poles)]
    p["full_tex"] = fmt.tex_sum(full)
    hterms = [(c, r"\delta[n]" if k == 0 else rf"\delta[n-{k}]") for k, c in enumerate(C)]
    hterms += [(a, fmt.tex_pow(q, "n") + r"\,u[n]" if q != 1 else "u[n]") for a, q in zip(A, poles)]
    p["h_tex"] = fmt.tex_sum(hterms)
    rows = []
    for n in range(4):
        dpart = C[n] if n < L else F(0)
        epart = sum(a * q**n for a, q in zip(A, poles))
        rows.append(f"<tr><td>${n}$</td><td>${fmt.tex_num(dpart)}$</td><td>${fmt.tex_num(epart)}$</td>"
                    f"<td><b>${fmt.tex_num(h[n])}$</b></td></tr>")
    p["table_html"] = ('<table class="table table-sm table-bordered text-center" style="width:auto">'
                       r"<tr><th>$n$</th><th>$\delta$ terms</th><th>$\sum_k A_k p_k^n$</th><th>$h[n]$</th></tr>"
                       + "".join(rows) + "</table>")
    p["h0_check_tex"] = fr"h[0] = \lim_{{z\to\infty}} H(z) = {fmt.tex_num(num[0])}"

    c = data["correct_answers"]
    c["C0"] = float(C[0])
    if L == 2:
        c["C1"] = float(C[1])
    c["A1"] = float(A[0])
    if m == 2:
        c["A2"] = float(A[1])
    for n in range(4):
        c[f"h{n}"] = float(h[n])
