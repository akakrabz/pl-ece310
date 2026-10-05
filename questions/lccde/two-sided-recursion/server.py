"""A BIBO-stable LTI system with one pole inside and one outside the unit circle (ROC an annulus),
implemented as two recursions whose outputs are added (exam family "two-sided systems as recursions":
FA2025 #7b; FA2024 #8a and FA2023 #8b for "non-causal / two-sided LCCDE"; Lecture 11 h = h_l + h_r,
Lecture 10 parallel connection, Lecture 8 inverse z-transform with an annulus ROC).

    H(z) = K(1 - q z^-1) / ((1 - p_in z^-1)(1 - p_out z^-1)),   |p_in| < 1 < |p_out|
         = A_in/(1 - p_in z^-1) + A_out/(1 - p_out z^-1)

shown as H(z) (FA2025 #7 framing) or as the LCCDE of "the BIBO-stable (two-sided) system" (FA2023 #8b).

Parts:
  (a) MC  which PFE term runs forward (causal) and which backward (anti-causal);
  (b) numbers  alpha, beta of the backward recursion  y_a[n-1] = alpha y_a[n] + beta x[n]
      (the form of the FA2025 key: alpha = 1/p_out, beta = -A_out/p_out);
  (c) numbers  y[-1], y[0], y[1] for the two-sample input x[n] = x0 delta[n] + x1 delta[n-1].
Part (c) uses an input with two samples (not delta[n]) on purpose: with x = delta, h[-1] = beta and
h[-2] = alpha*beta would repeat part (b); here every requested output sample needs both recursions
(y[0] = A_in x0 + beta x1 mixes them) and none equals a coefficient."""

import random
from fractions import Fraction as F

from ece310 import fmt, lccde, mc, poly, zt

P_IN = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(2, 3), F(-2, 3), F(3, 4), F(-3, 4)]
P_OUT = [F(2), F(-2), F(3), F(-3), F(3, 2), F(-3, 2), F(4), F(-4), F(4, 3), F(-4, 3)]
ZEROS = [F(0), F(0), F(1), F(-1), F(2), F(-2), F(3), F(-3), F(1, 2), F(-1, 2)]
GAINS = [F(1), F(1), F(2), F(3), F(-1), F(-2)]
X0 = [F(1), F(1), F(2), F(-1), F(-2), F(3)]
X1 = [F(1), F(-1), F(2), F(-2), F(3), F(-3), F(1, 2), F(-1, 2)]
NICE_DEN = {1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 18, 20, 24}


def _nice(v, den=NICE_DEN, top=60):
    return v.denominator in den and abs(v.numerator) <= top


def _pick():
    while True:
        p_in, p_out = random.choice(P_IN), random.choice(P_OUT)
        K, q = random.choice(GAINS), random.choice(ZEROS)
        if q in (p_in, p_out):
            continue
        num = [K] if q == 0 else [K, -K * q]
        C, (A_in, A_out) = zt.pfe(num, [p_in, p_out])
        if any(c != 0 for c in C) or A_in == 0 or A_out == 0:
            continue
        if not (_nice(A_in, {1, 2, 3, 4, 5, 6, 8, 9, 10, 12}, 30) and _nice(A_out, {1, 2, 3, 4, 5, 6, 8, 9, 10, 12}, 30)):
            continue
        alpha, beta = 1 / p_out, -A_out / p_out
        x0, x1 = random.choice(X0), random.choice(X1)
        # outputs: forward (causal) part and backward (anti-causal) part, from rest
        yc = {-1: F(0), 0: A_in * x0}
        yc[1] = p_in * yc[0] + A_in * x1
        ya = {1: F(0), 0: beta * x1}
        ya[-1] = alpha * ya[0] + beta * x0
        y = {n: yc[n] + ya[n] for n in (-1, 0, 1)}
        if not all(_nice(v) and v != 0 for v in list(y.values()) + [beta]):
            continue
        if len({y[-1], y[0], y[1], alpha, beta}) < 5:      # no answer repeats another
            continue
        return p_in, p_out, K, q, num, A_in, A_out, alpha, beta, x0, x1, yc, ya, y


def _cover(name, p, other, num):
    w = 1 / p
    numv = sum(c * w ** i for i, c in enumerate(num))
    denv = 1 - other * w
    sub = fmt.tex_num(num[0])
    if len(num) > 1:
        sub += (" - " if num[1] < 0 else " + ") + (fmt.tex_num(abs(num[1])) + r"\cdot " if abs(num[1]) != 1 else "") \
            + fmt.tex_num(w, paren=True)
    den = r"1 - " + fmt.tex_num(other, paren=True) + r"\cdot " + fmt.tex_num(w, paren=True)
    return (f"{name} &= \\left. H(z)\\,{fmt.tex_factor(p)}\\right|_{{z = {fmt.tex_num(p)}}}"
            f" = \\frac{{{sub}}}{{{den}}} = \\frac{{{fmt.tex_num(numv)}}}{{{fmt.tex_num(denv)}}} = {fmt.tex_num(numv / denv)}")


def _subst(terms):
    """c1 \\cdot (v1) + c2 \\cdot (v2) ... (values in parentheses), '0' if empty."""
    out = ""
    for c, v in terms:
        if c == 0 or v == 0:
            continue
        body = (fmt.tex_num(abs(c)) + r"\cdot " if abs(c) != 1 else "") + fmt.tex_num(v, paren=True)
        out += (("-" if c < 0 else "") if not out else (" - " if c < 0 else " + ")) + body
    return out or "0"


def generate(data):
    p_in, p_out, K, q, num, A_in, A_out, alpha, beta, x0, x1, yc, ya, y = _pick()
    den = poly.from_roots([p_in, p_out])
    p = data["params"]
    order = [p_in, p_out] if random.random() < 0.5 else [p_out, p_in]
    num_tex = fmt.tex_num(K) if q == 0 else ((("-" if K == -1 else fmt.tex_num(K)) if K != 1 else "") + fmt.tex_factor(q))
    if q != 0 and K == 1:
        num_tex = fmt.tex_sum([(1, ""), (-q, "z^{-1}")])
    p["Hgiven_tex"] = fmt.tex_frac(num_tex, "".join(fmt.tex_factor(v) for v in order))
    p["given_H"] = random.random() < 0.6
    p["given_lccde"] = not p["given_H"]
    p["eq_tex"] = lccde.tex_lccde(num, den, "std")
    p["H_tex"] = lccde.tex_H(num, den)
    p["b"], p["a"] = [str(v) for v in num], [str(v) for v in den]
    p["x"] = [str(x0), str(x1)]
    p["x_tex"] = fmt.tex_sum([(x0, r"\delta[n]"), (x1, r"\delta[n-1]")])
    p["pin_tex"], p["pout_tex"] = fmt.tex_num(p_in), fmt.tex_num(p_out)
    p["mpout_tex"] = fmt.tex_num(-p_out)
    p["roc_tex"] = fmt.tex_num(abs(p_in)) + r" \lt \lvert z\rvert \lt " + fmt.tex_num(abs(p_out))
    first, second = order
    p["pole1_tex"], p["pole2_tex"] = fmt.tex_num(first), fmt.tex_num(second)

    # (a) which term runs in which direction
    def opt(fwd, bwd):
        return (f"The term with the pole at $z = {fmt.tex_num(fwd)}$ forward in time (causal), the term with the pole at "
                f"$z = {fmt.tex_num(bwd)}$ backward in time (anti-causal)")
    correct = opt(p_in, p_out)
    swapped = opt(p_out, p_in)
    both_f = "Both terms forward in time (causal): the difference equation solved for the newest output"
    both_b = "Both terms backward in time (anti-causal): the difference equation solved for the oldest output"
    p["dir_choices"] = mc.choices(correct, [swapped, both_f, both_b], n=4)
    p["dir_meta"] = {correct: ["fwd_in"], swapped: ["fwd_out"], both_f: ["both_f"], both_b: ["both_b"]}

    # worked solution
    p["cover_tex"] = (r"\begin{aligned}" + _cover("A_{\\text{in}}", p_in, p_out, num) + r"\\"
                      + _cover("A_{\\text{out}}", p_out, p_in, num) + r"\end{aligned}")
    p["pfe_tex"] = (fmt.tex_frac(fmt.tex_num(A_in), fmt.tex_sum([(1, ""), (-p_in, "z^{-1}")])) + " + "
                    + fmt.tex_frac(fmt.tex_num(A_out), fmt.tex_sum([(1, ""), (-p_out, "z^{-1}")])))
    p["h_tex"] = fmt.tex_sum([(A_in, zt.exp_body_tex(p_in, "n")), (-A_out, zt.exp_body_tex(p_out, "-n-1"))])
    p["Ain_tex"], p["Aout_tex"] = fmt.tex_num(A_in), fmt.tex_num(A_out)
    p["fwd_eq_tex"] = "y_c[n] = " + fmt.tex_sum([(p_in, "y_c[n-1]"), (A_in, "x[n]")])
    p["aout_eq_tex"] = fmt.tex_sum([(1, "y_a[n]"), (-p_out, "y_a[n-1]")]) + " = " + fmt.tex_sum([(A_out, "x[n]")])
    p["bwd_eq_tex"] = "y_a[n-1] = " + fmt.tex_sum([(alpha, "y_a[n]"), (beta, "x[n]")])
    p["alpha_tex"], p["beta_tex"] = fmt.tex_num(alpha), fmt.tex_num(beta)
    p["beta_expr_tex"] = r"-\frac{A_{\text{out}}}{p_{\text{out}}} = -\frac{" + fmt.tex_num(A_out) + "}{" + fmt.tex_num(p_out) + "}"
    rows = [
        r"y_c[0] &= " + fmt.tex_num(A_in) + r"\,x[0] = " + fmt.tex_num(yc[0]),
        r"y_c[1] &= " + fmt.tex_sum([(p_in, "y_c[0]"), (A_in, "x[1]")]) + " = "
        + _subst([(p_in, yc[0]), (A_in, x1)]) + " = " + fmt.tex_num(yc[1]),
        r"y_a[0] &= " + fmt.tex_sum([(alpha, "y_a[1]"), (beta, "x[1]")]) + " = "
        + _subst([(beta, x1)]) + " = " + fmt.tex_num(ya[0]),
        r"y_a[-1] &= " + fmt.tex_sum([(alpha, "y_a[0]"), (beta, "x[0]")]) + " = "
        + _subst([(alpha, ya[0]), (beta, x0)]) + " = " + fmt.tex_num(ya[-1]),
    ]
    p["run_tex"] = r"\begin{aligned}" + r"\\".join(rows) + r"\end{aligned}"
    p["sum_tex"] = (r"\begin{aligned}"
                    + f"y[-1] &= y_c[-1] + y_a[-1] = 0 + {fmt.tex_num(ya[-1], paren=True)} = {fmt.tex_num(y[-1])}\\\\"
                    + f"y[0] &= y_c[0] + y_a[0] = {fmt.tex_num(yc[0])} + {fmt.tex_num(ya[0], paren=True)} = {fmt.tex_num(y[0])}\\\\"
                    + f"y[1] &= y_c[1] + y_a[1] = {fmt.tex_num(yc[1])} + 0 = {fmt.tex_num(y[1])}"
                    + r"\end{aligned}")
    h = {n: (A_in * p_in ** n if n >= 0 else -A_out * p_out ** n) for n in (-2, -1, 0, 1)}
    p["check_tex"] = ("y[n] = " + fmt.tex_sum([(x0, "h[n]"), (x1, "h[n-1]")])
                      + f" \\text{{ with }} h[-2] = {fmt.tex_num(h[-2])},\\ h[-1] = {fmt.tex_num(h[-1])},\\ "
                      f"h[0] = {fmt.tex_num(h[0])},\\ h[1] = {fmt.tex_num(h[1])}")
    assert all(x0 * h[n] + x1 * h[n - 1] == y[n] for n in (-1, 0, 1))

    ca = data["correct_answers"]
    ca["alpha"], ca["beta"] = float(alpha), float(beta)
    ca["ym1"], ca["y0"], ca["y1"] = float(y[-1]), float(y[0]), float(y[1])
