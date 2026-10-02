"""Output of a causal LTI system to x[n] = a^n u[n] (u[n] when a = 1) via Y(z) = H(z) X(z) and a
partial-fraction expansion (Lecture 9 §1, Lecture 10; HW4 #6(c); past exams FA2023 #7(b),
FA2025 #6(b), SP2025 #8(b)).

Kinds (the student is not told which):
  plain1  : H = K(1 - q z^-1)/(1 - p z^-1)                         -> 2-term PFE
  plain2  : H = K(1 - q z^-1)/((1 - p1 z^-1)(1 - p2 z^-1))         -> 3-term PFE
  xcancel : a zero of H sits on the input pole (FA2023 #7: zero at z = 1, input u[n])
  hcancel : H has a common factor (HW4 #6: the LCCDE's H reduces to first order)
In every kind Y(z) is proper after cancelling, so y[n] = sum_k A_k p_k^n holds for all n >= 0.
The system is given as an LCCDE (recursion form) or as H(z) with multiplied-out polynomials, so a
cancellation has to be found by factoring."""

import random
from fractions import Fraction as F

from ece310 import fmt, lccde, poly, zt

POLES = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(3, 4), F(-3, 4), F(2, 3), F(-2, 3),
         F(2), F(-2)]
INPUTS = [F(1), F(1), F(-1), F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(1, 4)]
ZEROS = [F(0), F(0), F(1), F(-1), F(2), F(-2), F(1, 2), F(-1, 2), F(3), F(-3)]


NICE_DEN = {1, 2, 3, 4, 5, 6, 8, 9, 10, 12}


def _nice(vals):
    return all(v != 0 and v.denominator in NICE_DEN and abs(v.numerator) <= 40 for v in vals)


def _subst_num(num, w):
    """LaTeX of num(z^-1) with z^-1 = w substituted: c0 + c1 (w) (num has degree <= 1 here)."""
    out = fmt.tex_num(num[0])
    for i, c in enumerate(num[1:], start=1):
        if c == 0:
            continue
        wt = fmt.tex_num(w, paren=True) + ("" if i == 1 else f"^{{{i}}}")
        mag = abs(c)
        body = wt if mag == 1 else fmt.tex_num(mag) + r"\cdot " + wt
        out += (" - " if c < 0 else " + ") + body
    return out


def _pick():
    while True:
        kind = random.choice(["plain1", "plain2", "xcancel", "xcancel", "hcancel", "hcancel"])
        a_in = random.choice(INPUTS)
        K = F(random.choice([1, 1, 2, 3, -1, -2]))
        q = random.choice(ZEROS)
        if kind == "plain1":
            hz, hp = [q], [random.choice(POLES)]
            yz, yp = [q], hp + [a_in]
        elif kind == "plain2":
            hp = random.sample(POLES, 2)
            hz = [q]
            yz, yp = [q], hp + [a_in]
        elif kind == "xcancel":
            hp = random.sample(POLES, 2)
            hz = [a_in, q]
            yz, yp = [q], list(hp)
        else:
            hp = random.sample(POLES, 2)
            hz = [hp[1], q]
            yz, yp = [q], [hp[0], a_in]
        if len(set(yp)) != len(yp) or q in yp:
            continue
        if kind != "hcancel" and any(z in hp for z in hz if z != 0):
            continue
        if kind == "hcancel" and (q in hp or a_in == hp[1]):
            continue
        if kind in ("plain1", "plain2", "hcancel") and q == a_in:
            continue
        if kind == "xcancel" and a_in in hp:
            continue
        num_y = poly.pscale(poly.from_roots([z for z in yz if z != 0]), K)
        if q == 0 and kind in ("xcancel", "hcancel"):
            num_y = [K]
        C, A = zt.pfe(num_y, yp)
        if any(c != 0 for c in C) or not _nice(A):
            continue
        b = poly.pscale(poly.from_roots([z for z in hz if z != 0]), K)
        a = poly.from_roots(hp)
        y = lccde.run(b, a, lambda n: a_in ** n, 4)
        if not all(v.denominator <= 64 and abs(v.numerator) <= 200 for v in y[:2]):
            continue
        return kind, a_in, K, hz, hp, yz, yp, num_y, A, b, a


def _H_factored(K, zeros, poles):
    nz = [z for z in zeros if z != 0]
    num = "".join(fmt.tex_factor(z) for z in nz)
    if K != 1 or not nz:
        num = (fmt.tex_num(K) if not (K == -1 and nz) else "-") + num
    return fmt.tex_frac(num, "".join(fmt.tex_factor(p) for p in poles))


def generate(data):
    kind, a_in, K, hz, hp, yz, yp, num_y, A, b, a = _pick()
    form = random.choice(["lccde", "H"])
    p = data["params"]
    p["kind"] = kind
    p["b"], p["a"] = [str(v) for v in b], [str(v) for v in a]
    p["a_in"] = str(a_in)
    p["poles_y"] = [str(v) for v in yp]
    p["res_y"] = [str(v) for v in A]
    p["given_lccde"] = form == "lccde"
    p["given_H"] = form == "H"
    p["eq_tex"] = lccde.tex_lccde(b, a, "rec")
    p["H_tex"] = lccde.tex_H(b, a)
    p["Hfact_tex"] = _H_factored(K, hz, hp)
    p["x_tex"] = r"u[n]" if a_in == 1 else fmt.tex_pow(a_in, "n") + r"\,u[n]"
    p["X_tex"] = fmt.tex_frac("1", fmt.tex_sum([(1, ""), (-a_in, "z^{-1}")]))
    p["Xroc_tex"] = r"\lvert z\rvert \gt " + fmt.tex_num(abs(a_in))
    p["Y_tex"] = fmt.tex_frac(fmt.tex_poly_zinv(num_y), "".join(fmt.tex_factor(v) for v in yp))
    yroc = max(abs(v) for v in yp)
    p["Yroc_tex"] = r"\lvert z\rvert \gt " + fmt.tex_num(yroc)
    p["xcancel"] = kind == "xcancel"
    p["hcancel"] = kind == "hcancel"
    p["nocancel"] = kind in ("plain1", "plain2")
    p["cancel_tex"] = ""
    if kind == "xcancel":
        p["cancel_tex"] = fmt.tex_factor(a_in)
    elif kind == "hcancel":
        p["cancel_tex"] = fmt.tex_factor(hp[1])
    # cover-up for each residue
    lines = []
    for k, (pk, Ak) in enumerate(zip(yp, A)):
        w = 1 / pk
        numv = sum(c * w ** i for i, c in enumerate(num_y))
        others = [pj for j, pj in enumerate(yp) if j != k]
        den_tex = "".join(r"\left(1 - " + fmt.tex_num(pj, paren=True) + r"\cdot " + fmt.tex_num(w, paren=True) + r"\right)"
                          for pj in others)
        denv = 1
        for pj in others:
            denv *= 1 - pj * w
        num_sub = _subst_num(num_y, w)
        line = (f"A_{k + 1} &= \\left. Y(z)\\,{fmt.tex_factor(pk)} \\right|_{{z = {fmt.tex_num(pk)}}}"
                f" = \\frac{{{num_sub}}}{{{den_tex}}}")
        line += f" = \\frac{{{fmt.tex_num(numv)}}}{{{fmt.tex_num(denv)}}}"
        lines.append(line + f" = {fmt.tex_num(Ak)}")
    p["cover_tex"] = r"\begin{aligned}" + r"\\".join(lines) + r"\end{aligned}"
    p["pfe_tex"] = " + ".join(fmt.tex_frac(fmt.tex_num(Ak), fmt.tex_sum([(1, ""), (-pk, "z^{-1}")])) for pk, Ak in zip(yp, A))
    p["y_tex"] = fmt.tex_sum((Ak, zt.exp_body_tex(pk, "n")) for pk, Ak in zip(yp, A))
    y0 = sum(A)
    y1 = sum(Ak * pk for pk, Ak in zip(yp, A))
    rec = lccde.run(b, a, lambda n: a_in ** n, 2)
    assert rec[0] == y0 and rec[1] == y1
    p["y0_tex"], p["y1_tex"] = fmt.tex_num(y0), fmt.tex_num(y1)
    p["y1sum_tex"] = fmt.tex_sum((Ak * pk, "") for pk, Ak in zip(yp, A))
    p["y0sum_tex"] = fmt.tex_sum((Ak, "") for Ak in A)

    ca = data["correct_answers"]
    ca["y"] = fmt.sym_sum((Ak, fmt.sym_pow(pk, "n")) for pk, Ak in zip(yp, A))
    ca["y0"] = float(y0)
    ca["y1"] = float(y1)


def grade(data):
    """y[n] is a sequence: accept any expression that agrees with the answer at n = 0, 1, ..., 15.
    (SymPy's equality test treats n as complex, so it cannot prove e.g. (-2)**(-n) == (-1/2)**n,
    and sometimes gives up on sums of shifted exponentials; both are equal on the integers.)"""
    ps = data["partial_scores"].get("y")
    sub = data["submitted_answers"].get("y")
    if not ps or (ps.get("score") or 0) >= 1 or "y" in data["format_errors"] or not isinstance(sub, dict):
        return
    try:
        import prairielearn as pl
        from prairielearn import sympy_utils as psu
        expr = psu.json_to_sympy(sub, allow_complex=False)
        syms = list(expr.free_symbols)
        poles = [F(v) for v in data["params"]["poles_y"]]
        res = [F(v) for v in data["params"]["res_y"]]
        for k in range(16):
            want = float(sum(A * p ** k for A, p in zip(res, poles)))
            got = complex(expr.subs({s: k for s in syms}).evalf(30))
            # written so that NaN / inf (0/0, zoo, log(0), ...) fail: every comparison with NaN is False
            if not (abs(got.imag) <= 1e-9 and abs(got.real - want) <= 1e-9 * max(1.0, abs(want))):
                return
    except Exception:
        return
    ps["score"] = 1.0
    pl.set_weighted_score_data(data)
