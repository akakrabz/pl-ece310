"""Read a difference equation off a direct form I block diagram (Lecture 5 §2, Figs. 1-4: delay,
gain and adder blocks; the direct form I diagram of eq. 12).

The diagram (inline SVG from ece310.w2_diagrams) draws every gain as it acts: each adder input
enters with a plus sign, so with the drawn feedback gains g_i and feed-forward gains b_j the
system is y[n] = sum_i g_i y[n-i] + sum_j b_j x[n-j], and in the Lecture 9 form
y[n] + sum_k a_k y[n-k] = sum_k b_k x[n-k] the feedback coefficients are a_i = -g_i. (g is a
neutral letter: Lecture 5 itself writes sum b_i y[n-i] + sum c_j x[n-j], with b for the feedback.)
The question states the convention explicitly. Asked: b and a of the equation the diagram
implements term by term (no cancellation), FIR or IIR, and h[0], h[1], h[2].

Families: fir (no feedback), iir1 / iir2 (first / second-order feedback, B and A coprime), cancel
(a feedback loop whose pole is cancelled by a zero, e.g. y[n] = y[n-1] + x[n] - x[n-2]: FIR although
the diagram has a loop - the "recursive moving average" of Lecture 5 §1.1). Zero taps are allowed
(a delay without a gain branch)."""

import random
from fractions import Fraction as F

import numpy as np
import prairielearn as pl
from ece310 import fmt, lccde, mc, poly, w2_diagrams

FF = [1, 2, 3, -1, -2, -3, F(1, 2), F(-1, 2), F(3, 2), F(-3, 2)]
FB1 = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4), F(3, 4), F(-3, 4), 1, -1, F(3, 2), F(-3, 2)]
FB2 = [F(1, 2), F(-1, 2), F(1, 4), F(-1, 4)]


def dec(q):
    """Plain decimal string of a terminating fraction: 2, -0.5, 0.25."""
    q = F(q)
    return str(q.numerator) if q.denominator == 1 else f"{float(q):g}"


def _draw():
    fam = random.choice(["fir", "iir1", "iir1", "iir2", "iir2", "cancel"])
    if fam == "fir":
        nff = random.choice([1, 2])
        ff = [random.choice(FF) for _ in range(nff + 1)]
        if nff == 2 and random.random() < 0.3:
            ff[1] = 0
        return fam, ff, []
    if fam == "iir1":
        nff = random.choice([0, 1, 1, 2, 2])
        ff = [random.choice(FF) for _ in range(nff + 1)]
        if nff >= 1 and random.random() < 0.2:
            ff[0] = 0
        return fam, ff, [random.choice(FB1)]
    if fam == "iir2":
        nff = random.choice([0, 1, 2])
        ff = [random.choice(FF) for _ in range(nff + 1)]
        fb = [random.choice(FB2 + [0]), random.choice(FB2)]
        return fam, ff, fb
    g = random.choice([1, -1, F(1, 2), F(-1, 2)])             # cancel: B = (1 - g z^-1)(beta0 + beta1 z^-1)
    beta = [random.choice([1, 2, -1, -2, 3]), random.choice([1, -1, 2, -2, 3, -3])]
    ff = poly.pmul([1, -g], beta)
    return fam, ff, [g]


def _h_work(n, b, g, h):
    """h[n] = g_1 h[n-1] + ... + b_n with numbers (recursion from rest, x = delta)."""
    sym, val = [], []
    for i in range(1, len(g) + 1):
        if n - i >= 0:
            sym.append(rf"g_{i}\,h[{n - i}]")
            val.append(rf"{fmt.tex_num(g[i - 1], paren=True)}\cdot{fmt.tex_num(h[n - i], paren=True)}")
    bn = b[n] if n < len(b) else F(0)
    sym.append(f"b_{n}")
    val.append(fmt.tex_num(bn, paren=True))
    if len(sym) == 1:
        return rf"h[{n}] = b_{n} = {fmt.tex_num(h[n])}"
    return rf"h[{n}] = {' + '.join(sym)} = {' + '.join(val)} = {fmt.tex_num(h[n])}"


def generate(data):
    while True:
        fam, ff, fb = _draw()
        ff = [F(v) for v in ff]
        fb = [F(v) for v in fb]
        if ff[-1] == 0 or all(v == 0 for v in ff):
            continue
        b = ff
        a = [F(1)] + [-v for v in fb]
        _, r = poly.pdivmod(b, a)
        fir = all(v == 0 for v in r)
        if fam in ("iir1", "iir2") and (not lccde.coprime(b, a) or fir):
            continue
        if fam == "cancel" and not fir:
            continue
        if any(abs(v) > 4 for v in b):
            continue
        break
    K, Mm1 = len(fb), len(ff) - 1

    p = data["params"]
    p["family"], p["form"] = fam, "I"
    p["ff"] = [dec(v) for v in ff]
    p["fb"] = [dec(v) for v in fb]
    p["fig_html"] = w2_diagrams.direct_form(ff, fb, form="I")
    p["rec_tex"] = lccde.tex_lccde(b, a, form="rec")
    p["std_tex"] = lccde.tex_lccde(b, a, form="std")
    p["b_plain"] = "[" + ", ".join(dec(v) for v in b) + "]"
    p["a_plain"] = "[" + ", ".join(dec(v) for v in a) + "]"
    p["K"], p["Mm1"] = K, Mm1
    p["has_fb"] = K > 0
    p["no_fb"] = K == 0
    p["fir"], p["iir"] = fir, not fir
    p["fir_choices"] = mc.choices("FIR (finite-length impulse response)", ["IIR (infinite-length impulse response)"],
                                  n=2, shuffle=False)
    for o in p["fir_choices"]:
        o["correct"] = (o["text"].startswith("FIR") == fir)
    h = lccde.impulse(b, a, 6)
    p["hwork"] = [{"line_tex": _h_work(n, b, fb, h)} for n in range(3)]
    p["h_list_tex"] = r",\ ".join(fmt.tex_num(v) for v in h)
    p["is_cancel"] = fam == "cancel"
    if fir:
        q, _ = poly.pdivmod(b, a)
        q = poly.trim(q)
        p["Hred_tex"] = fmt.tex_poly_zinv(q)
        p["fir_why_html"] = (
            (r"There is no feedback branch, so $h[n]$ is just the list of feed-forward gains: FIR." if K == 0 else
             rf"The diagram has a feedback loop, but its pole is cancelled: $H(z) = {lccde.tex_H(b, a)} = {fmt.tex_poly_zinv(q)}$, "
             r"a polynomial in $z^{-1}$. So $h[n]$ has finitely many nonzero samples: FIR (like the recursive moving "
             r"average of Lecture 5 &sect;1.1)."))
        if fam == "cancel":
            p["cancel_note_html"] = (
                rf"(Term by term as drawn: $B(z)$ and $A(z)$ share a factor here, which cancels in $H(z) = {fmt.tex_poly_zinv(q)}$, "
                r"but the diagram computes the uncancelled equation; part (b) is where the cancellation matters.)")
    else:
        p["fir_why_html"] = (rf"$H(z) = {lccde.tex_H(b, a)}$ keeps its pole(s) (no common factor with the numerator), "
                             r"so the feedback keeps the impulse response going forever: IIR.")

    c = data["correct_answers"]
    c["bvec"] = pl.to_json(np.array([[float(v) for v in b]]))
    c["avec"] = pl.to_json(np.array([[float(v) for v in a]]))
    c["h0"], c["h1"], c["h2"] = float(h[0]), float(h[1]), float(h[2])
