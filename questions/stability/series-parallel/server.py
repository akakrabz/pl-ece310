"""Two causal first-order systems in series or in parallel (HW4 #5; Lecture 10 system algebra,
Lecture 11 BIBO stability): the overall H(z) and whether each system and the overall system are
BIBO stable.

Subsystem i: H_i(z) = K_i (1 - q_i z^-1) / (1 - p_i z^-1) = C_i + A_i/(1 - p_i z^-1), shown as an
impulse response C_i delta[n] + A_i p_i^n u[n], as a difference equation, or as H_i(z).
Series variants include pole-zero cancellations between the two systems (HW4 #5: the zero of h2
at z = 1 cancels the unstable pole of h1, and the cascade is stable); |p| = 1 counts as unstable."""

import random
from fractions import Fraction as F

from ece310 import fmt, lccde, mc, poly

STABLE = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(3, 4), F(-3, 4), F(2, 3), F(-2, 3)]
UNSTABLE = [F(1), F(-1), F(2), F(-2), F(3, 2), F(-3, 2)]
ZEROS = [F(0), F(0), F(1), F(-1), F(2), F(-2), F(1, 2), F(-1, 2), F(3), F(-3)]
GAINS = [F(1), F(1), F(2), F(3), F(-1), F(-2)]


def _sub_tex(i, K, q, p, form):
    b = [K] if q == 0 else [K, -K * q]
    a = [F(1), -p]
    if form == "h":
        C = -b[1] / p if len(b) > 1 else F(0)
        A = b[0] - C
        body = fmt.tex_sum([(C, r"\delta[n]"), (A, (fmt.tex_pow(p, "n") + r"\,u[n]") if p != 1 else "u[n]")])
        return f"h_{i}[n] = " + body
    if form == "lccde":
        return lccde.tex_lccde(b, a, "rec").replace("y[", f"y_{i}[").replace("x[", f"x_{i}[")
    return f"H_{i}(z) = " + lccde.tex_H(b, a)


def _num_fact(K, q):
    r"""K (1 - q z^-1) as LaTeX: '-\left(1 - 2z^{-1}\right)' for K = -1, never '-1\left(...\right)'."""
    if q == 0:
        return fmt.tex_num(K)
    f = fmt.tex_factor(q)
    return f if K == 1 else ("-" + f if K == -1 else fmt.tex_num(K) + f)


def _derive(i, K, q, p, form):
    """One display line: how H_i(z) follows from the given form of system i."""
    b = [K] if q == 0 else [K, -K * q]
    a = [F(1), -p]
    den = fmt.tex_sum([(1, ""), (-p, "z^{-1}")])
    H = lccde.tex_H(b, a)
    if form == "h":
        C = -b[1] / p if len(b) > 1 else F(0)
        A = b[0] - C
        if C == 0:
            line = f"H_{i}(z) = {H}"
        else:
            cden = (r"\left(" + den + r"\right)") if C == 1 else (("-" if C == -1 else fmt.tex_num(C)) + r"\left(" + den + r"\right)")
            mid = fmt.tex_frac(cden + (" - " if A < 0 else " + ") + fmt.tex_num(abs(A)), den)
            line = f"H_{i}(z) = {fmt.tex_num(C)} + " + fmt.tex_frac(fmt.tex_num(A), den) + f" = {mid} = {H}"
    elif form == "lccde":
        line = (f"Y_{i}(z)\\left({den}\\right) = X_{i}(z)\\left({lccde.tex_poly_zinv(b)}\\right)"
                f"\\ \\Rightarrow\\ H_{i}(z) = {H}")
    else:
        line = f"H_{i}(z) = {H}"
    if K != 1 and q != 0:
        line += " = " + fmt.tex_frac(_num_fact(K, q), den)
    return line


def _pick():
    while True:
        conn = random.choice(["series", "series", "series", "parallel", "parallel"])
        kind = "plain"
        if conn == "series":
            kind = random.choice(["plain", "cancel_unstable", "cancel_unstable", "cancel_stable"])
        p1 = random.choice(STABLE + UNSTABLE)
        p2 = random.choice(STABLE + UNSTABLE)
        q1, q2 = random.choice(ZEROS), random.choice(ZEROS)
        if kind == "cancel_unstable":
            p1, p2 = random.choice(UNSTABLE), random.choice(STABLE)
            q2 = p1
        elif kind == "cancel_stable":
            p1, p2 = random.choice(STABLE), random.choice(UNSTABLE)
            q2 = p1
        if p1 == p2 or q1 == p1 or q2 == p2:
            continue
        if kind == "plain" and (q1 == p2 or q2 == p1):
            continue
        if kind != "plain" and q1 == p2:
            continue
        K1, K2 = random.choice(GAINS), random.choice(GAINS)
        if random.random() < 0.5:          # randomly swap the roles of the two systems
            p1, p2, q1, q2, K1, K2 = p2, p1, q2, q1, K2, K1
        return conn, kind, (K1, q1, p1), (K2, q2, p2)


def generate(data):
    conn, kind, s1, s2 = _pick()
    (K1, q1, p1), (K2, q2, p2) = s1, s2
    b1 = [K1] if q1 == 0 else [K1, -K1 * q1]
    b2 = [K2] if q2 == 0 else [K2, -K2 * q2]
    a1, a2 = [F(1), -p1], [F(1), -p2]
    if conn == "series":
        B, A = poly.pmul(b1, b2), poly.pmul(a1, a2)
    else:
        B, A = poly.padd(poly.pmul(b1, a2), poly.pmul(b2, a1)), poly.pmul(a1, a2)
    # cancellation between the systems (series only by construction)
    cancelled = [p for p in (p1, p2) if p in [z for z in (q1, q2) if z != 0]] if conn == "series" else []
    remaining = [p for p in (p1, p2) if p not in cancelled]
    st1, st2 = abs(p1) < 1, abs(p2) < 1
    st = all(abs(p) < 1 for p in remaining)

    p = data["params"]
    forms = [random.choice(["h", "lccde", "H"]) for _ in range(2)]
    p["sys1_tex"] = _sub_tex(1, K1, q1, p1, forms[0])
    p["sys2_tex"] = _sub_tex(2, K2, q2, p2, forms[1])
    p["series"] = conn == "series"
    p["parallel"] = conn == "parallel"
    p["conn_word"] = "in series (cascade)" if conn == "series" else "in parallel (outputs added)"
    p["spec"] = {"conn": conn, "s1": [str(K1), str(q1), str(p1)], "s2": [str(K2), str(q2), str(p2)]}
    p["H1_tex"], p["H2_tex"] = lccde.tex_H(b1, a1), lccde.tex_H(b2, a2)
    p["H1_deriv_tex"] = _derive(1, K1, q1, p1, forms[0])
    p["H2_deriv_tex"] = _derive(2, K2, q2, p2, forms[1])
    p["H_tex"] = lccde.tex_H(B, A)
    p["p1_tex"], p["p2_tex"] = fmt.tex_num(p1), fmt.tex_num(p2)
    p["absp1_tex"], p["absp2_tex"] = fmt.tex_num(abs(p1)), fmt.tex_num(abs(p2))
    p["st1_word"] = "stable" if st1 else "not stable"
    p["st2_word"] = "stable" if st2 else "not stable"
    p["cancel"] = bool(cancelled)
    p["cancel_tex"] = fmt.tex_factor(cancelled[0]) if cancelled else ""
    if cancelled:
        Br = poly.pmul([K1 * K2], poly.from_roots([z for z in (q1, q2) if z != 0 and z not in cancelled]))
        Ar = poly.from_roots(remaining)
        p["Hred_tex"] = lccde.tex_H(Br, Ar)
    else:
        p["Hred_tex"] = p["H_tex"]
    p["rem_tex"] = r",\ ".join(fmt.tex_num(v) for v in remaining)
    p["rem_one"] = len(remaining) == 1
    p["roc_tex"] = r"\lvert z\rvert \gt " + fmt.tex_num(max(abs(v) for v in remaining))
    p["st_word"] = "is" if st else "is not"
    p["st_reason"] = ("every remaining pole lies inside the unit circle" if st
                      else "a remaining pole lies on or outside the unit circle")
    p["st1_choices"] = mc.yes_no(st1)
    p["st2_choices"] = mc.yes_no(st2)
    p["st_choices"] = mc.yes_no(st)
    data["correct_answers"]["H"] = lccde.sym_H(B, A)
