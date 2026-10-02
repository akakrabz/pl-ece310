"""All possible ROCs of a rational H(z) with 2 or 3 real poles of distinct magnitudes (none on the
unit circle): how many ROCs, which one is BIBO stable, and h[-2], h[0], h[2] for one given ROC.

Model: the "all possible ROCs" exam family (FA2023 Midterm 1 #6a, FA2019 #7, SP2025 #8b-c),
HW4 #1 (all possible ROCs, including the cancelled pole in #1b) and HW3 #3(c) (two-sided h for a
given annulus); Lecture 8 (PFE) and Lecture 11 (ROC <-> stability).

In about a third of the variants the numerator (shown multiplied out) contains a factor that
cancels one of the displayed denominator factors: that pole does not count (HW4 #1b trap).
The ROC asked in part (c) is never the stable one, so it does not give away part (b)."""

import random
from fractions import Fraction as F

from ece310 import fmt, poly, zinv, zt

RES = [1, 2, 3, -1, -2, -3, F(1, 2), F(-1, 2), F(3, 2), F(-3, 2)]


def _num_from_residues(poles, A):
    """sum_k A_k prod_{j != k} (1 - p_j z^-1)."""
    num = [F(0)]
    for k, a in enumerate(A):
        term = [F(1)]
        for j, q in enumerate(poles):
            if j != k:
                term = poly.pmul(term, [F(1), -q])
        num = poly.padd(num, poly.pscale(term, a))
    return poly.trim(num)


def _pick():
    while True:
        m = random.choice([2, 2, 3])
        poles = random.sample(zinv.NICE_POLES, m)
        mags = {abs(q) for q in poles}
        if len(mags) != m:
            continue
        c = None
        if random.random() < 0.35:
            c = random.choice([q for q in zinv.NICE_POLES if abs(q) not in mags])
        A = [F(random.choice(RES)) for _ in range(m)]
        num = _num_from_residues(poles, A)
        shown = poly.pmul(num, [F(1), -c]) if c is not None else num
        if len(shown) < 2 or not all(zinv.nice(x, maxden=12, maxnum=24, bad_dens=(5, 7, 10, 11)) for x in shown):
            continue
        return poles, A, num, c, shown


def _roc_sort_key(r):
    return (r.inner, F(10**6) if r.outer is None else r.outer)


def _h_tex(n, contrib, value):
    """h[n] = c1 (p1)^n + ... = value (contrib: list of (coef, pole))."""
    if not contrib:
        return f"h[{n}] = 0"
    parts = []
    for i, (cf, q) in enumerate(contrib):
        t = fmt.tex_num(abs(cf)) + r"\cdot" + fmt.tex_num(q, paren=True) + "^{" + str(n) + "}"
        parts.append(("-" if cf < 0 else "") + t if i == 0 else (" - " if cf < 0 else " + ") + t)
    return f"h[{n}] = " + "".join(parts) + " = " + fmt.tex_num(value)


def _choose():
    """Pick H and the part-(c) ROC until h[-2], h[0], h[2] are exam-nice."""
    while True:
        poles, A, num, c, shown = _pick()
        rocs = zt.rocs_for_poles(poles)
        stable = next(r for r in rocs if r.contains_unit_circle())
        others = [r for r in rocs if r != stable]
        weights = [2 if r.kind() == "two-sided" else 1 for r in others]
        spec = random.choices(others, weights=weights)[0]
        if c is not None and spec.inner < abs(c) and (spec.outer is None or abs(c) < spec.outer):
            continue          # an ROC across the cancelled pole's circle would reveal the cancellation
        inv = zt.inverse_pfe([F(0)], poles, A, spec)
        hv = {n: inv.value(n) for n in (-2, 0, 2)}
        if all(zinv.nice(v, maxden=18, maxnum=80, bad_dens=()) for v in hv.values()):
            return poles, A, num, c, shown, rocs, stable, spec, inv, hv


def generate(data):
    poles, A, num, c, shown, rocs, stable, spec, inv, hv = _choose()
    m = len(poles)
    order = poles + ([c] if c is not None else [])
    random.shuffle(order)

    opts = {r.tex(): r for r in rocs}
    if c is not None:
        for r in zt.rocs_for_poles(poles + [c]):
            opts.setdefault(r.tex(), r)
    # 1-2 regions that contain a pole (never ROCs), so the number of options does not answer (a)
    mags = sorted(abs(q) for q in poles)
    bad = [zt.ROC(inner=F(0), outer=mags[k], has0=True, hasinf=False) for k in range(1, len(mags))]
    bad += [zt.ROC(inner=mags[k], outer=None, has0=False, hasinf=True) for k in range(len(mags) - 1)]
    bad += [zt.ROC(inner=mags[i], outer=mags[j], has0=False, hasinf=False)
            for i in range(len(mags)) for j in range(i + 2, len(mags))]
    bad = [r for r in bad if r.tex() not in opts]
    for r in random.sample(bad, min(len(bad), random.choice([1, 2]))):
        opts[r.tex()] = r
    ordered = sorted(opts.values(), key=_roc_sort_key)

    p = data["params"]
    p["shown_num"] = [str(x) for x in shown]
    p["den_poles"] = [str(q) for q in order]
    p["spec_inner"] = str(spec.inner)
    p["spec_outer"] = None if spec.outer is None else str(spec.outer)
    p["H_tex"] = fmt.tex_frac(fmt.tex_poly_zinv(shown), "".join(fmt.tex_factor(q) for q in order))
    p["spec_tex"] = spec.tex()
    p["stab_choices"] = [{"text": f"${r.tex()}$", "correct": r == stable} for r in ordered]

    # ---------------------------------------------------------------- worked solution
    p["cancel"] = c is not None
    p["nocancel"] = c is None
    p["c_tex"] = fmt.tex_num(c) if c is not None else ""
    p["cmag_tex"] = fmt.tex_num(abs(c)) if c is not None else ""
    p["cinv_tex"] = fmt.tex_num(1 / c) if c is not None else ""
    p["cfac_tex"] = fmt.tex_factor(c) if c is not None else ""
    p["shown_tex"] = fmt.tex_poly_zinv(shown)
    p["red_num_tex"] = fmt.tex_poly_zinv(num)
    p["Hred_tex"] = fmt.tex_frac(fmt.tex_poly_zinv(num), "".join(fmt.tex_factor(q) for q in sorted(poles, key=abs)))
    sp = sorted(poles, key=abs)
    p["poles_list_tex"] = r",\ ".join(fmt.tex_num(q) for q in sp)
    p["mags_tex"] = r" \lt ".join(fmt.tex_num(abs(q)) for q in sp)
    p["m"] = m
    p["nroc"] = m + 1
    p["rocs_tex"] = [{"tex": r.tex(), "stable": r == stable} for r in rocs]
    p["stable_tex"] = stable.tex()
    p["cover_tex"] = r" \\[2pt] ".join(zinv.coverup_tex(num, poles, k, f"A_{{{k + 1}}}") for k in range(m))
    p["two_m_tex"] = f"2^{{{m}}}"
    p["pfe_tex"] = fmt.tex_sum((a, fmt.tex_frac("1", fmt.tex_sum([(1, ""), (-q, "z^{-1}")]))) for a, q in zip(A, poles))
    p["poles_named_tex"] = r",\ ".join(f"p_{k + 1} = {fmt.tex_num(q)}" for k, q in enumerate(poles))
    lines = []
    for q, a, s in zip(poles, A, inv.sides):
        if s == "R":
            lines.append(rf"pole ${fmt.tex_num(q)}$: $\lvert p\rvert = {fmt.tex_num(abs(q))}$ is inside the ROC's inner circle, "
                         rf"right-sided term ${fmt.tex_sum([(a, zt.exp_body_tex(q, 'n'))])}$")
        else:
            lines.append(rf"pole ${fmt.tex_num(q)}$: $\lvert p\rvert = {fmt.tex_num(abs(q))}$ is outside the ROC's outer circle, "
                         rf"left-sided term ${fmt.tex_sum([(-a, zt.exp_body_tex(q, '-n-1'))])}$")
    p["side_lines"] = [{"line_html": s} for s in lines]
    p["h_tex"] = inv.tex()
    right = [(a, q) for a, q in inv.right_terms()]
    left = [(-a, q) for a, q in inv.left_terms()]
    p["hm2_tex"] = _h_tex(-2, left, hv[-2])
    p["h0_tex"] = _h_tex(0, right, hv[0])
    p["h2_tex"] = _h_tex(2, right, hv[2])

    cr = data["correct_answers"]
    cr["nroc"] = m + 1
    cr["hm2"], cr["h0"], cr["h2"] = float(hv[-2]), float(hv[0]), float(hv[2])
