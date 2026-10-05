"""First samples h[0..3] of the impulse response of a causal LCCDE at initial rest, computed by
running the recursion with x[n] = delta[n] (Lecture 5 §1.2 and Exercise 1; SP2021 #6; HW4 #6(b)),
plus FIR or IIR, i.e. whether h[n] is finite.

Three kinds of variants (shown in the Lecture 9 form or as a recursion "y[n] = ..."):
  fir_plain : no feedback terms, y[n] = b0 x[n] + b1 x[n-1] (+ b2 x[n-2])                 -> FIR
  fir_rec   : feedback that cancels, y[n] = p y[n-1] + c (x[n] - p^L x[n-L]) (Lecture 5 Eq. (4),
              the recursive running sum / moving average, for p = 1): H(z) = c(1 - p^L z^-L)/(1 - p z^-1)
              = c(1 + p z^-1 + ... + p^(L-1) z^-(L-1)), so h[n] = c p^n for 0 <= n < L, then 0  -> FIR
  iir       : feedback terms with B, A coprime (no pole-zero cancellation)                 -> IIR
The question asks whether the impulse response is finite, so the answer follows Lecture 9 §1.2
("IIR iff a finite pole is not cancelled by a zero"); Lecture 5's shorthand "K > 0 feedback terms
=> IIR" holds only when nothing cancels (Lecture 5 itself computes its FIR moving average with a
feedback term in Eq. (4))."""

import random
from fractions import Fraction as F

from ece310 import fmt, lccde, mc

A1_CHOICES = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4), F(1), F(-1), F(2), F(-2), F(1, 3), F(-1, 3), F(3, 2), F(-3, 2)]
A2_PAIRS_A1 = [F(0), F(1, 2), F(-1, 2), F(1), F(-1)]
A2_CHOICES = [F(1, 4), F(-1, 4), F(1, 2), F(-1, 2), F(1), F(-1)]
B0_CHOICES = [1, 1, 2, 3, 4, -1, -2, -3]


def _nice(vals):
    return all(v.denominator <= 32 and abs(v.numerator) <= 99 for v in vals)


FIRREC_P = [F(1), F(1), F(-1), F(1, 2), F(-1, 2), F(2), F(-2)]
FIRREC_C = [F(1), F(1), F(2), F(3), F(-1), F(-2), F(-3), F(1, 2), F(-1, 2), F(1, 3), F(1, 4), F(3, 2)]


def _pick():
    while True:
        r = random.random()
        kind = "fir_plain" if r < 0.2 else ("fir_rec" if r < 0.45 else "iir")
        extra = {}
        if kind == "fir_plain":
            a = [F(1)]
            nb = random.choice([2, 3, 3, 3, 4, 4])
            b = [F(random.choice(B0_CHOICES))] + [F(random.randint(-4, 4)) for _ in range(nb - 1)]
        elif kind == "fir_rec":
            pp = random.choice(FIRREC_P)
            L = random.choice([2, 3, 3, 4])
            c = random.choice(FIRREC_C)
            if pp == 1 and random.random() < 0.5:
                c = F(1, L)                  # Lecture 5's moving average
            d = random.choice([0, 0, 1])     # optional extra delay: c(x[n-d] - p^L x[n-d-L])
            a = [F(1), -pp]
            b = [F(0)] * d + [c] + [F(0)] * (L - 1) + [-c * pp ** L]
            extra = {"p": pp, "L": L, "c": c, "d": d}
            if (c * pp ** L).denominator > 16:
                continue
        else:
            order = random.choice([1, 1, 2])
            if order == 1:
                a = [F(1), random.choice(A1_CHOICES)]
            else:
                a = [F(1), random.choice(A2_PAIRS_A1), random.choice(A2_CHOICES)]
            nb = random.choice([1, 2, 2, 2, 3, 3, 3])
            b = [F(random.choice(B0_CHOICES))] + [F(random.randint(-3, 3)) for _ in range(nb - 1)]
            if not lccde.coprime(b, a):
                continue
        if len(b) > 1 and b[-1] == 0:
            continue
        h = lccde.impulse(b, a, 4)
        if not _nice(h) or all(v == 0 for v in h[1:]):
            continue
        return kind, b, a, h, extra


def _subst(terms):
    """LaTeX of sum(coef * value) with the values substituted: c \\cdot (v)."""
    out = ""
    for c, v in terms:
        c = F(c)
        if c == 0:
            continue
        mag = abs(c)
        if v is None:      # a constant term (input coefficient b_n)
            body = fmt.tex_num(mag)
        else:
            v = F(v)
            body = fmt.tex_num(v, paren=True) if mag == 1 else fmt.tex_num(mag) + r"\cdot " + fmt.tex_num(v, paren=True)
        if not out:
            out = ("-" if c < 0 else "") + body
        else:
            out += (" - " if c < 0 else " + ") + body
    return out or "0"


def _firrec_rec_tex(pp, L, c, d):
    """y[n] = p y[n-1] + c(x[n-d] - p^L x[n-d-L]) with c factored out, as Lecture 5 Eq. (4) writes it."""
    inner = fmt.tex_sum([(1, lccde.sig("x", d)), (-pp ** L, lccde.sig("x", d + L))])
    if c == 1:
        tail = " + " + inner
    else:
        tail = (" - " if c < 0 else " + ") + ("" if abs(c) == 1 else fmt.tex_num(abs(c))) + r"\left(" + inner + r"\right)"
    return "y[n] = " + fmt.tex_sum([(pp, "y[n-1]")]) + tail


def generate(data):
    kind, b, a, h, extra = _pick()
    fir = kind != "iir"
    form = random.choice(["rec", "std"])
    p = data["params"]
    p["kind"] = kind
    p["b"] = [str(v) for v in b]
    p["a"] = [str(v) for v in a]
    p["fir"] = fir
    p["fir_plain"] = kind == "fir_plain"
    p["fir_rec"] = kind == "fir_rec"
    p["iir"] = kind == "iir"
    p["has_feedback"] = len(a) > 1
    p["eq_tex"] = lccde.tex_lccde(b, a, form)
    if kind == "fir_rec" and form == "rec":
        p["eq_tex"] = _firrec_rec_tex(extra["p"], extra["L"], extra["c"], extra["d"])
    p["std_tex"] = lccde.tex_lccde(b, a, "std")
    p["rec_tex"] = lccde.tex_lccde(b, a, "rec")

    # worked recursion with x[n] = delta[n]: h[n] = -a1 h[n-1] - a2 h[n-2] + b_n
    steps = []
    for n in range(4):
        sym = [(-a[k], f"h[{n - k}]") for k in range(1, len(a)) if n - k >= 0]
        num = [(-a[k], h[n - k]) for k in range(1, len(a)) if n - k >= 0]
        bn = b[n] if n < len(b) else F(0)
        sym.append((bn, ""))
        num.append((bn, None))
        parts = [fmt.tex_sum(sym)]
        if len([t for t in num if t[0] != 0]) > 1:
            parts.append(_subst(num))
        parts.append(fmt.tex_num(h[n]))
        dedup = [s for i, s in enumerate(parts) if i == 0 or s != parts[i - 1]]
        steps.append(f"h[{n}] &= " + " = ".join(dedup))
    p["steps_tex"] = r"\begin{aligned}" + r"\\".join(steps) + r"\end{aligned}"
    p["hrec_tex"] = ("h[n] = " + fmt.tex_sum([(-a[k], f"h[n-{k}]") for k in range(1, len(a))]
                                              + [(c, r"\delta[n]" if k == 0 else rf"\delta[n-{k}]") for k, c in enumerate(b)]))
    p["H_tex"] = lccde.tex_H(b, a)
    p["h_seq_tex"] = fmt.tex_seq(h, 0) + (r"\ \text{(then zeros)}" if fir else r",\ \dots")
    # fir_rec: the cancellation and the first zero sample h[L]
    p["L"], p["L_minus_1"], p["Hpoly_tex"], p["pole_tex"], p["hL_tex"], p["factor_tex"] = 0, 0, "", "", "", ""
    p["nz"], p["Hnumfact_tex"] = 0, ""
    if kind == "fir_rec":
        pp, L, c, d = extra["p"], extra["L"], extra["c"], extra["d"]
        nz = d + L                      # first sample that is zero for good
        hfull = lccde.impulse(b, a, nz + 1)
        p["L"], p["L_minus_1"], p["nz"] = L, L - 1, nz
        cpart = "" if c == 1 else ("-" if c == -1 else fmt.tex_num(c))
        p["Hnumfact_tex"] = fmt.tex_frac(
            cpart + (f"z^{{-{d}}}" if d else "") + "\\left(" + fmt.tex_sum([(1, ""), (-pp ** L, f"z^{{-{L}}}")]) + "\\right)",
            fmt.tex_sum([(1, ""), (-pp, "z^{-1}")]))
        p["Hpoly_tex"] = lccde.tex_poly_zinv([F(0)] * d + [c * pp ** k for k in range(L)])
        p["pole_tex"] = fmt.tex_num(pp)
        p["factor_tex"] = (fmt.tex_sum([(1, ""), (-pp ** L, f"z^{{-{L}}}")]) + r" = " + fmt.tex_factor(pp)
                           + r"\left(" + lccde.tex_poly_zinv([pp ** k for k in range(L)]) + r"\right)")
        p["hL_tex"] = (f"h[{nz}] = " + fmt.tex_sum([(pp, f"h[{nz - 1}]"), (-c * pp ** L, "")]) + " = "
                       + _subst([(pp, hfull[nz - 1]), (-c * pp ** L, None)]) + " = " + fmt.tex_num(hfull[nz]))
        assert hfull[nz] == 0 and all(v == 0 for v in lccde.impulse(b, a, nz + 6)[nz:])
    p["nb"] = len(b)
    p["nb_last"] = len(b) - 1
    p["firiir_choices"] = mc.yes_no(fir, "FIR", "IIR")

    for n in range(4):
        data["correct_answers"][f"h{n}"] = float(h[n])
