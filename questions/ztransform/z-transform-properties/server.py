"""z-transform properties (Lecture 7 table; HW3 #2; SP2021 #4): given x[n] <-> X(z) = A/(1 - a z^-1)
with a right-sided ROC |z| > |a| or a left-sided ROC |z| < |a|, find the transform and ROC of one of
x[n-k], n x[n], c^n x[n], x[-n], x[n]*x[n], x[n] - x[n-1].

The student types Y(z) (pl-symbolic-input) and picks the ROC (multiple choice). ROC effects that are
tested: shifts can drop z = 0 or z = infinity, scaling multiplies the radius by |c|, reversal inverts it
and swaps the side, the others keep R_x."""

import random
from fractions import Fraction as F

from ece310 import fmt, mc, zt

MAGS = [F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(4, 3), F(3, 2), F(2), F(3)]
SCALES = [F(2), F(3), F(1, 2), F(1, 3), F(-2), F(-1, 2), F(3, 2), F(2, 3), F(-3)]
LZ = r"\lvert z\rvert"


def _frac_tex(num, den):
    return fmt.tex_frac(num, den)


def _sfrac(coef, body, den):
    """coef*body/den with the sign in front of the fraction (library convention)."""
    frac = fmt.tex_frac(fmt.tex_sum([(abs(coef), body)]), den)
    return "-" + frac if coef < 0 else frac


def _y_terms(op, side, A, a, k, c):
    """y[n] in the time domain as terms c * P(n) * b^(n-s) * step; step ('R', e) = u[n-e], ('L', e) = u[-n+e]."""
    T = lambda co, poly, b, sh, sd, e: {"c": str(co), "poly": poly, "b": str(b), "s": sh, "side": sd, "e": e}
    R = side == "R"
    if op == "shift":
        return [T(A, "1", a, k, "R", k)] if R else [T(-A, "1", a, k, "L", k - 1)]
    if op == "mult_n":
        return [T(A, "n", a, 0, "R", 0)] if R else [T(-A, "n", a, 0, "L", -1)]
    if op == "scale":
        return [T(A, "1", a * c, 0, "R", 0)] if R else [T(-A, "1", a * c, 0, "L", -1)]
    if op == "reverse":
        return [T(A, "1", 1 / a, 0, "L", 0)] if R else [T(-A, "1", 1 / a, 0, "R", 1)]
    if op == "conv":
        return [T(A * A, "n+1", a, 0, "R", 0)] if R else [T(-A * A, "n+1", a, 0, "L", -2)]
    return ([T(A, "1", a, 0, "R", 0), T(-A, "1", a, 1, "R", 1)] if R
            else [T(-A, "1", a, 0, "L", -1), T(A, "1", a, 1, "L", 0)])


def _y_terms_tex(terms):
    out = []
    for t in terms:
        b, sh = F(t["b"]), t["s"]
        poly = {"1": "", "n": r"n\,", "n+1": r"(n+1)\,"}[t["poly"]]
        expo = "n" if sh == 0 else (f"n-{sh}" if sh > 0 else f"n+{-sh}")
        power = "" if b == 1 else fmt.tex_pow(b, expo) + r"\,"
        if poly and power[:1].isdigit():
            poly = poly[:-2] + r" \cdot "          # n \cdot 3^{n}, not n\,3^{n}
        step = "u[" + (zt.u_arg_right(t["e"]) if t["side"] == "R" else zt.u_arg_left(t["e"])) + "]"
        out.append((F(t["c"]), poly + power + step))
    return fmt.tex_sum(out)


def _pole_den(p, power=1):
    d = fmt.tex_sum([(1, ""), (-p, "z^{-1}")])
    return d if power == 1 else r"\left(" + d + r"\right)^{" + str(power) + "}"


def roc_choices(correct, extra, radii):
    """correct ROC + property-specific wrong ROCs (extra, in priority order) + generic ones."""
    R = zt.ROC
    generic = []
    k = correct.kind()
    if k == "right":
        r = correct.inner
        generic += [R(inner=r, outer=None, has0=False, hasinf=not correct.hasinf), R(inner=F(0), outer=r, has0=True, hasinf=False)]
        generic += [R(inner=w, outer=None, has0=False, hasinf=correct.hasinf) for w in radii if w != r]
    else:
        r = correct.outer
        generic += [R(inner=F(0), outer=r, has0=not correct.has0, hasinf=False), R(inner=r, outer=None, has0=False, hasinf=True)]
        generic += [R(inner=F(0), outer=w, has0=correct.has0, hasinf=False) for w in radii if w != r]
    pool = list(extra)
    rest = generic[2:]
    random.shuffle(rest)
    pool += generic[:2] + rest
    by_text = {"$" + o.tex() + "$": o for o in reversed([correct] + pool)}
    ch = mc.choices("$" + correct.tex() + "$", ["$" + o.tex() + "$" for o in pool], n=5)
    for c in ch:
        c["roc"] = by_text[c["text"]].to_json()
    return ch


def generate(data):
    a = random.choice(MAGS) * (-1 if random.random() < 0.3 else 1)
    A = random.choice([F(1), F(1), F(1), F(2), F(3)])
    side = "R" if random.random() < 0.7 else "L"
    op = random.choice(["shift", "shift", "mult_n", "mult_n", "scale", "reverse", "conv", "diff"])
    ra = abs(a)
    Rx = zt.ROC.right(ra, 0) if side == "R" else zt.ROC.left(ra, -1)
    X_tex = _frac_tex(fmt.tex_num(A), _pole_den(a))
    p = data["params"]
    radii = {ra, 1 / ra, F(1)}
    k, c = 0, F(1)

    if op == "shift":
        k = random.choice([-3, -2, -1, 1, 2, 3])
        y_tex = rf"x[{zt.u_arg_right(k)}]"
        Y_sym = f"{fmt.sym_num(A)}*z**({-k})/(1 - {fmt.sym_num(a)}*z**(-1))"
        Y_tex = _frac_tex(fmt.tex_sum([(A, fmt.zpow_tex(k))]), _pole_den(a))
        if side == "R":
            roc = zt.ROC.right(ra, k)
            note = (rf"$x[n]$ is right-sided and starts at $n = 0$, so $y[n]$ starts at "
                    + (rf"$n = {k} \lt 0$: the factor $z^{{{-k}}}$ blows up at $z = \infty$, which is excluded." if k < 0
                       else rf"$n = {k} \gt 0$: the factor $z^{{{-k}}}$ only matters at $z = 0$, which was not in the ROC anyway."))
        else:
            roc = zt.ROC.left(ra, k - 1)
            note = (rf"$x[n]$ is left-sided and ends at $n = -1$, so $y[n]$ ends at "
                    + (rf"$n = {k - 1} \gt 0$: negative powers of $z$ appear, so $z = 0$ is excluded." if k - 1 > 0
                       else rf"$n = {k - 1} \le 0$: no negative powers of $z$, so $z = 0$ stays in the ROC."))
        extra = [Rx] if Rx.tex() != roc.tex() else []
        prop = rf"Time shift: $x[n-k] \leftrightarrow z^{{-k}}X(z)$ with $k = {k}$."
        steps = rf"Y(z) = {fmt.zpow_tex(k)}\cdot {X_tex} = {Y_tex}"
    elif op == "mult_n":
        y_tex = r"n\,x[n]"
        Y_sym = f"{fmt.sym_num(A * a)}*z**(-1)/(1 - {fmt.sym_num(a)}*z**(-1))**2"
        Y_tex = _sfrac(A * a, "z^{-1}", _pole_den(a, 2))
        roc = Rx
        note = r"Multiplying by $n$ does not change the ROC: $R_y = R_x$."
        extra = []
        prop = r"Differentiation: $n\,x[n] \leftrightarrow -z\,\dfrac{dX(z)}{dz}$."
        dX = _sfrac(-A * a, "z^{-2}", _pole_den(a, 2))
        steps = rf"\frac{{dX}}{{dz}} = {dX}, \qquad Y(z) = -z\,\frac{{dX}}{{dz}} = {Y_tex}"
    elif op == "scale":
        c = random.choice([s for s in SCALES if abs((s * a).numerator) <= 9 and (s * a).denominator <= 9 and abs(s * a) != 1])
        y_tex = fmt.tex_pow(c, "n") + r"\,x[n]"
        Y_sym = f"{fmt.sym_num(A)}/(1 - {fmt.sym_num(a * c)}*z**(-1))"
        Y_tex = _frac_tex(fmt.tex_num(A), _pole_den(a * c))
        roc = zt.ROC.right(abs(a * c), 0) if side == "R" else zt.ROC.left(abs(a * c), -1)
        note = (rf"The ROC scales by $\lvert c\rvert = {fmt.tex_num(abs(c))}$: the boundary moves from "
                rf"${fmt.tex_num(ra)}$ to ${fmt.tex_num(abs(a * c))}$ (the pole moves from ${fmt.tex_num(a)}$ to ${fmt.tex_num(a * c)}$).")
        wrong = abs(a) / abs(c)
        radii |= {abs(a * c), wrong}
        extra = []
        if wrong != abs(a * c):
            extra.append(zt.ROC.right(wrong, 0) if side == "R" else zt.ROC.left(wrong, -1))
        if abs(c) != 1:
            extra.append(Rx)
        prop = rf"Scaling: $c^{{n}}x[n] \leftrightarrow X(z/c)$, ROC $\lvert c\rvert R_x$, with $c = {fmt.tex_num(c)}$."
        steps = (rf"Y(z) = X\!\left(\frac{{z}}{{{fmt.tex_num(c)}}}\right) = \frac{{{fmt.tex_num(A)}}}{{1 - {fmt.tex_num(a, paren=True)}\left(\frac{{z}}{{{fmt.tex_num(c)}}}\right)^{{-1}}}} = {Y_tex}"
                 if c > 0 else
                 rf"Y(z) = X\!\left(\frac{{z}}{{{fmt.tex_num(c, paren=True)}}}\right) = {Y_tex}")
    elif op == "reverse":
        y_tex = r"x[-n]"
        Y_sym = f"{fmt.sym_num(A)}/(1 - {fmt.sym_num(a)}*z)"
        Y_tex = _frac_tex(fmt.tex_num(A), fmt.tex_sum([(1, ""), (-a, "z")]))
        alt = _sfrac(-A / a, "z^{-1}", _pole_den(1 / a))
        if side == "R":
            roc = zt.ROC.left(1 / ra, 0)
            note = (rf"The ROC inverts: ${LZ} \gt {fmt.tex_num(ra)}$ becomes ${LZ} \lt {fmt.tex_num(1 / ra)}$; "
                    r"$y[n]$ is left-sided and ends at $n = 0$, so $z = 0$ is included.")
            extra = [zt.ROC.right(1 / ra, 0), Rx]
        else:
            roc = zt.ROC.right(1 / ra, 1)
            note = (rf"The ROC inverts: ${LZ} \lt {fmt.tex_num(ra)}$ becomes ${LZ} \gt {fmt.tex_num(1 / ra)}$; "
                    r"$y[n]$ is right-sided and starts at $n = 1$, so $z = \infty$ is included.")
            extra = [zt.ROC.left(1 / ra, -1), Rx]
        prop = r"Time reversal: $x[-n] \leftrightarrow X(z^{-1})$, ROC $1/R_x$."
        steps = rf"Y(z) = X(z^{{-1}}) = {Y_tex} = {alt}"
    elif op == "conv":
        y_tex = r"x[n] * x[n]"
        Y_sym = f"{fmt.sym_num(A * A)}/(1 - {fmt.sym_num(a)}*z**(-1))**2"
        Y_tex = _frac_tex(fmt.tex_num(A * A), _pole_den(a, 2))
        roc = Rx
        note = r"At least $R_x \cap R_x = R_x$, and nothing cancels the (now double) pole, so $R_y = R_x$."
        extra = []
        prop = r"Convolution: $x_1[n] * x_2[n] \leftrightarrow X_1(z)X_2(z)$."
        steps = rf"Y(z) = X(z)^{{2}} = {Y_tex}"
    else:  # diff
        y_tex = r"x[n] - x[n-1]"
        Y_sym = f"{fmt.sym_num(A)}*(1 - z**(-1))/(1 - {fmt.sym_num(a)}*z**(-1))"
        Y_tex = _frac_tex(fmt.tex_num(A) + r"\left(1 - z^{-1}\right)" if A != 1 else r"1 - z^{-1}", _pole_den(a))
        roc = Rx
        if side == "R":
            note = r"Both terms converge on $R_x$ (the delay only affects $z = 0$, already outside), and the zero at $z = 1$ cancels no pole: $R_y = R_x$."
        else:
            note = r"$x[n-1]$ ends at $n = 0$, so $y[n]$ has no samples at $n \gt 0$ and $z = 0$ stays in the ROC: $R_y = R_x$."
        extra = [zt.ROC.left(ra, 1)] if side == "L" else [zt.ROC.right(F(1), 0)] if ra < 1 else []
        prop = r"Linearity and time shift: $x[n] - x[n-1] \leftrightarrow X(z) - z^{-1}X(z)$."
        steps = rf"Y(z) = \left(1 - z^{{-1}}\right)X(z) = {Y_tex}"

    choices = roc_choices(roc, extra, sorted(radii - {F(0)}))
    y_terms = _y_terms(op, side, A, a, k, c)
    p.update({
        "A": str(A), "a": str(a), "side": side, "op": op, "k": k, "c": str(c),
        "X_tex": X_tex, "rocx_tex": Rx.tex(), "y_tex": y_tex, "Y_tex": Y_tex, "roc_tex": roc.tex(),
        "prop": prop, "steps_tex": steps, "note": note, "roc_choices": choices,
        "y_terms": y_terms, "y_time_tex": _y_terms_tex(y_terms),
        "x_tex": (fmt.tex_sum([(A, zt.exp_body_tex(a, "n"))]) if side == "R"
                  else fmt.tex_sum([(-A, zt.exp_body_tex(a, "-n-1"))])),
    })
    data["correct_answers"]["Y"] = Y_sym
