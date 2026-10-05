"""Inverse z-transform of a second-order H(z) by partial fractions, with the ROC fixed by a stated
property (causal / BIBO stable / left-sided).

Model: Lecture 8 Exercise 1 (PFE, right- and left-sided answers), HW3 #3(c)-(d) and #4, and the
"all possible ROCs" exam family (FA2025 #7a stable h, FA2024 #7b causal h, SP2025 #8c).

H(z) = (b0 + b1 z^-1) / ((1 - p1 z^-1)(1 - p2 z^-1)), real poles with |p1| < |p2| (never on the
unit circle, so exactly one ROC is stable), shown factored or multiplied out (then the student
also finds the poles). Asked: (poles), A1, A2, the ROC (multiple choice over the three possible
ROCs), and h[n] for n >= 0 and for n <= -1 (pl-symbolic-input in an integer n)."""

import random
from fractions import Fraction as F

from ece310 import fmt, poly, zinv, zt

B0 = [1, 2, 3, -1, -2, F(1, 2), F(3, 2), F(-1, 2), F(-3, 2)]
B1 = [0, 1, -1, 2, -2, 3, -3, F(1, 2), F(-1, 2), F(3, 2), F(-3, 2)]
PROPS = {
    "causal": "causal",
    "stable": "BIBO stable",
    "left": "left-sided",
}


def _pick():
    while True:
        pa, pb = random.sample(zinv.NICE_POLES, 2)
        if abs(pa) == abs(pb):
            continue
        poles = sorted([pa, pb], key=abs)
        num = [F(random.choice(B0)), F(random.choice(B1))]
        A = [zinv.coverup(num, poles, k)[2] for k in range(2)]
        if A[0] == 0 or A[1] == 0 or A[0] == A[1]:
            continue          # A = 0 means a zero cancels that pole
        if not all(zinv.nice(a) for a in A):
            continue
        return num, poles, A


def _tex_terms(terms):
    return fmt.tex_sum((a, fmt.tex_pow(p, "n")) for a, p in terms)


def generate(data):
    num, poles, A = _pick()
    p1, p2 = poles
    m1, m2 = abs(p1), abs(p2)
    den = poly.from_roots(poles)
    multiplied = random.random() < 0.5 and all(zinv.nice(c, maxden=12, maxnum=20, bad_dens=()) for c in den)
    prop = random.choice(["stable", "stable", "causal", "left"])

    rocs = zt.rocs_for_poles(poles)            # |z| < m1,  m1 < |z| < m2,  |z| > m2
    if prop == "causal":
        roc = rocs[2]
    elif prop == "left":
        roc = rocs[0]
    else:
        roc = next(r for r in rocs if r.contains_unit_circle())
    inv = zt.inverse_pfe([F(0)], poles, A, roc)

    p = data["params"]
    p["num"] = [str(c) for c in num]
    p["poles"] = [str(q) for q in poles]
    p["prop"] = prop
    p["prop_text"] = PROPS[prop]
    p["multiplied"] = multiplied
    p["factored"] = not multiplied
    factors = [fmt.tex_factor(q) for q in poles]
    random.shuffle(factors)
    p["H_tex"] = fmt.tex_frac(fmt.tex_poly_zinv(num), fmt.tex_poly_zinv(den) if multiplied else "".join(factors))
    p["num_tex"] = fmt.tex_poly_zinv(num)
    p["den_tex"] = fmt.tex_poly_zinv(den)
    p["den_fact_tex"] = fmt.tex_factor(p1) + fmt.tex_factor(p2)
    p["den_z_tex"] = fmt.tex_sum([(1, "z^{2}"), (den[1], "z"), (den[2], "")])
    p["den_zfact_tex"] = (r"\left(" + fmt.tex_sum([(1, "z"), (-p1, "")]) + r"\right)\left("
                          + fmt.tex_sum([(1, "z"), (-p2, "")]) + r"\right)")
    p["p1_tex"], p["p2_tex"] = fmt.tex_num(p1), fmt.tex_num(p2)
    p["m1_tex"], p["m2_tex"] = fmt.tex_num(m1), fmt.tex_num(m2)
    p["letters"] = {"pfe": "b" if multiplied else "a", "roc": "c" if multiplied else "b",
                    "h": "d" if multiplied else "c"}
    p["A1_cover_tex"] = zinv.coverup_tex(num, poles, 0, "A_1")
    p["A2_cover_tex"] = zinv.coverup_tex(num, poles, 1, "A_2")
    p["Asum_tex"] = fr"A_1 + A_2 = {fmt.tex_num(A[0])} {'-' if A[1] < 0 else '+'} {fmt.tex_num(abs(A[1]))} = {fmt.tex_num(A[0] + A[1])} = b_0"
    lin = lambda q: fmt.tex_sum([(1, ""), (-q, "z^{-1}")])  # noqa: E731  1 - q z^{-1}
    p["pfe_tex"] = fmt.tex_frac(fmt.tex_num(A[0]), lin(p1)) + (" - " if A[1] < 0 else " + ") \
        + fmt.tex_frac(fmt.tex_num(abs(A[1])), lin(p2))

    kinds = ["left-sided (anti-causal)", "two-sided", "right-sided (causal)"]
    # multiplied-out H: the options must not reveal the pole magnitudes, so they name |p1|, |p2|
    z, P1, P2 = r"\lvert z\rvert", r"\lvert p_1\rvert", r"\lvert p_2\rvert"
    sym_rocs = [z + r" \lt " + P1, P1 + r" \lt " + z + r" \lt " + P2, z + r" \gt " + P2]
    p["rocs_list_tex"] = [{"tex": (sr + r",\ \text{i.e. } " + r.tex()) if multiplied else r.tex(), "kind": k}
                          for r, sr, k in zip(rocs, sym_rocs, kinds)]
    p["roc_tex"] = roc.tex()
    if prop == "causal":
        why = (r"A causal system has its ROC outside the outermost pole (and including $z = \infty$): "
               rf"$\lvert z\rvert \gt {fmt.tex_num(m2)}$.")
    elif prop == "left":
        why = rf"A left-sided $h[n]$ has its ROC inside the innermost pole: $\lvert z\rvert \lt {fmt.tex_num(m1)}$."
    else:
        where = ("between the two pole circles" if m1 < 1 < m2 else
                 "outside both poles (both are inside the unit circle)" if m2 < 1 else
                 "inside both poles (both are outside the unit circle)")
        why = (r"BIBO stable $\Leftrightarrow$ the ROC contains the unit circle $\lvert z\rvert = 1$, "
               f"which lies {where}: ${roc.tex()}$.")
    p["why_html"] = why
    side_lines = []
    for k, (q, a, s) in enumerate(zip(poles, A, inv.sides), start=1):
        if s == "R":
            side_lines.append(rf"$p_{k} = {fmt.tex_num(q)}$ lies inside the ROC's inner circle: right-sided, "
                              rf"$\dfrac{{A_{k}}}{{1 - p_{k}z^{{-1}}}} \to {fmt.tex_sum([(a, zt.exp_body_tex(q, 'n'))])}$")
        else:
            side_lines.append(rf"$p_{k} = {fmt.tex_num(q)}$ lies outside the ROC's outer circle: left-sided, "
                              rf"$\dfrac{{A_{k}}}{{1 - p_{k}z^{{-1}}}} \to {fmt.tex_sum([(-a, zt.exp_body_tex(q, '-n-1'))])}$")
    p["side_lines"] = [{"line_html": s} for s in side_lines]
    p["h_tex"] = inv.tex()
    p["hpos_tex"] = _tex_terms(inv.right_terms())
    p["hneg_tex"] = _tex_terms([(-a, q) for a, q in inv.left_terms()])
    p["roc_choices"] = [{"text": f"${sr if multiplied else r.tex()}$", "correct": r == roc}
                        for r, sr in zip(rocs, sym_rocs)]

    c = data["correct_answers"]
    if multiplied:
        c["p1"], c["p2"] = float(p1), float(p2)
    c["A1"], c["A2"] = float(A[0]), float(A[1])
    c["hpos"] = zinv.sym_json(inv.sym_nonneg())
    c["hneg"] = zinv.sym_json(inv.sym_neg())
