"""Causal system with a real parameter K: for which K is it BIBO stable? (SP2025 #7, FA2024 #8b,
FA2023 #8a; study-site "parameters for stability", Practice 2.)

Three answer shapes:
  interval  K multiplies a pole: y[n] = ((K - m)/c) y[n-1] + ..., a pole pair +-sqrt(K/c), a pole (K - m)/c
            next to a fixed stable pole. The student enters both endpoints and says which endpoints belong
            to the set: an endpoint belongs exactly when the unit-circle pole there is cancelled by a zero of
            the numerator (zeros at 1, at -1, or at both, in some difference-equation variants);
  single    K sits in the numerator (FA2023 #8a / FA2024 #8b): only the cancelling value works;
  set       pole K/c and a zero q (Practice 2), multiple choice over set expressions. The shape of the
            answer depends on q: |q| > 1 -> interval or K = cq; |q| < 1 -> the open interval only (the
            cancelling value is already inside); q = +-1 -> one endpoint included.
params["model"] gives every coefficient of b(K), a(K) (powers of z^-1) as [const, slope] for the checker."""

import random
from fractions import Fraction as F

from ece310 import fmt

STAB = [F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(1, 4), F(-1, 4), F(2, 3), F(3, 4), F(-2, 3), F(-3, 4)]
BAD = [F(2), F(-2), F(3), F(-3), F(3, 2), F(-3, 2), F(4, 3), F(-4, 3), F(5, 2)]
ENDS = [("Neither endpoint", False, False), (r"Only $K_{\text{lo}}$", True, False),
        (r"Only $K_{\text{hi}}$", False, True), ("Both endpoints", True, True)]


def lin(c0, c1=0):
    return [str(F(c0)), str(F(c1))]


def kfrac(m, c):
    """LaTeX of (K - m)/c."""
    top = "K" if m == 0 else (f"K - {m}" if m > 0 else f"K + {-m}")
    return r"\dfrac{" + top + "}{" + str(c) + "}" if c != 1 else (top if m == 0 else "(" + top + ")")


def abs_km(m):
    return r"\lvert K\rvert" if m == 0 else (r"\lvert K - " + str(m) + r"\rvert" if m > 0 else r"\lvert K + " + str(-m) + r"\rvert")


def set_tex(lo, hi, lo_in, hi_in):
    le, lt = r"\le", r"\lt"
    return f"{lo} {le if lo_in else lt} K {le if hi_in else lt} {hi}"


def gen_interval(p):
    kind = random.choice(["lccde", "lccde", "lccde", "pair", "cascade"])
    c = random.choice([2, 3, 4, 5, 6])
    lo_in = hi_in = False
    if kind == "lccde":
        m = random.choice([0, 0, 1, -1, 2, -2, 3, -3])
        ends = random.choices(["none", "inside", "upper", "lower", "both"], weights=[0.7, 1.3, 1.4, 1.4, 1.4])[0]
        b0 = random.choice([F(1), F(1), F(2), F(3)])
        b = {"none": [F(1)], "upper": [F(1), F(-1)], "lower": [F(1), F(1)], "both": [F(1), F(0), F(-1)]}.get(ends)
        if ends == "inside":
            b = [F(1), random.choice([F(1, 2), F(-1, 2), F(1, 3), F(-1, 4), F(1, 4), F(-1, 3), F(2, 3)])]
        b = [b0 * bk for bk in b]
        g = fmt.tex_num(b0) if b0 != 1 else ""                     # gain left after a cancellation
        lo, hi = m - c, m + c
        lo_in, hi_in = ends in ("lower", "both"), ends in ("upper", "both")
        rhs = fmt.tex_sum([(bk, "x[n]" if k == 0 else f"x[n-{k}]") for k, bk in enumerate(b)])
        p["system_html"] = fr"the causal LTI system described by the difference equation $$y[n] = {kfrac(m, c)}\,y[n-1] + {rhs}.$$"
        p["model"] = {"b": [lin(bk) for bk in b], "a": [lin(1), lin(F(m, c), F(-1, c))]}
        num = fmt.tex_poly_zinv(b)
        endtxt = {
            "none": f"The numerator ${num}$ has no zero: nothing can cancel the pole on the unit circle, so neither endpoint belongs to the set "
                    "(a pole on $\\lvert z\\rvert = 1$ means marginally stable, which is <b>not</b> BIBO stable).",
            "inside": (f"The only zero, $z = {fmt.tex_num(-b[-1] / b[0])}$, is inside the unit circle, so it cannot cancel a pole on the circle: "
                       "neither endpoint belongs to the set. (It cancels the pole for one $K$ inside the interval, which changes nothing.)"),
            "upper": (f"The numerator ${num}$ vanishes at $z = 1$: at $K = {hi}$ the pole cancels, $H(z) = {g or '1'}$, $h[n] = {g}\\delta[n]$, "
                      f"stable, so $K = {hi}$ belongs to the set. At $K = {lo}$ the pole $z = -1$ survives: not stable."),
            "lower": (f"The numerator ${num}$ vanishes at $z = -1$: at $K = {lo}$ the pole cancels, $H(z) = {g or '1'}$, $h[n] = {g}\\delta[n]$, "
                      f"stable, so $K = {lo}$ belongs to the set. At $K = {hi}$ the pole $z = 1$ survives: not stable."),
            "both": (f"The numerator ${num} = {g}(1 - z^{{-1}})(1 + z^{{-1}})$ vanishes at $z = \\pm 1$: at $K = {hi}$ the pole $1$ cancels and "
                     f"leaves the FIR filter $H(z) = {g}(1 + z^{{-1}})$; at $K = {lo}$ the pole $-1$ cancels and leaves $H(z) = {g}(1 - z^{{-1}})$. "
                     "Both are stable, so both endpoints belong to the set."),
        }[ends]
        p["sol_html"] = (fr"$H(z) = \dfrac{{{num}}}{{1 - {kfrac(m, c)}\,z^{{-1}}}}$ has one pole, $p = {kfrac(m, c)}$. "
                         r"Causal $\Rightarrow$ ROC $\lvert z\rvert \gt \lvert p\rvert$, stable when $\lvert p\rvert \lt 1$:"
                         fr"$$\left\lvert {kfrac(m, c)}\right\rvert \lt 1 \iff {abs_km(m)} \lt {c} \iff {lo} \lt K \lt {hi}.$$"
                         fr"At the endpoints the pole is on the unit circle ($K = {hi}$: $p = 1$; $K = {lo}$: $p = -1$). " + endtxt)
    elif kind == "pair":
        b0 = random.choice([1, 1, 2, 3])
        b1 = random.choice([F(1, 2), F(-1, 2), F(1, 3), F(2), F(-1), F(1), F(-2), F(3)])
        num = fmt.tex_sum([(b0, ""), (b1, "z^{-1}")])
        p["system_html"] = fr"the causal LTI system with transfer function $$H(z) = \dfrac{{{num}}}{{1 - \tfrac{{K}}{{{c}}}\,z^{{-2}}}}.$$"
        lo, hi = -c, c
        p["model"] = {"b": [lin(b0), lin(b1)], "a": [lin(1), lin(0), lin(0, F(-1, c))]}
        p["sol_html"] = (fr"The poles solve $z^2 = \tfrac{{K}}{{{c}}}$: for $K \ge 0$ they are $z = \pm\sqrt{{K/{c}}}$ (real), for $K \lt 0$ they are "
                         fr"$z = \pm j\sqrt{{\lvert K\rvert/{c}}}$ (imaginary). Either way $\lvert p\rvert = \sqrt{{\lvert K\rvert/{c}}}$, and causal "
                         fr"stability needs $\lvert p\rvert \lt 1$: $$\sqrt{{\lvert K\rvert/{c}}} \lt 1 \iff \lvert K\rvert \lt {c} \iff -{c} \lt K \lt {c}.$$"
                         fr"Negative $K$ is allowed: imaginary poles are fine as long as they are inside the unit circle. At $K = {c}$ the poles "
                         fr"are $\pm 1$ and at $K = -{c}$ they are $\pm j$, all on the unit circle; the single zero $z = {fmt.tex_num(-b1 / b0)}$ can "
                         r"cancel at most one of them, so neither endpoint belongs to the set.")
    else:
        p0 = random.choice(STAB)
        m = random.choice([0, 0, 1, -1, 2])
        p["system_html"] = (fr"a series connection of two causal LTI systems with impulse responses "
                            fr"$$h_1[n] = \left({kfrac(m, c)}\right)^{{n}}u[n] \quad\text{{and}}\quad h_2[n] = {fmt.tex_pow(p0)}u[n]$$"
                            r"(the overall system is $h = h_1 * h_2$).")
        lo, hi = m - c, m + c
        # a(K) = (1 - ((K - m)/c) z^-1)(1 - p0 z^-1)
        p["model"] = {"b": [lin(1)], "a": [lin(1), lin(-p0 + F(m, c), F(-1, c)), lin(-p0 * F(m, c), p0 / c)]}
        p["sol_html"] = (fr"$H(z) = H_1(z)H_2(z) = \dfrac{{1}}{{\left(1 - {kfrac(m, c)}\,z^{{-1}}\right){fmt.tex_factor(p0)}}}$ with poles "
                         fr"${kfrac(m, c)}$ and ${fmt.tex_num(p0)}$, and no zeros that could cancel anything. The pole ${fmt.tex_num(p0)}$ is inside "
                         fr"the unit circle; the causal cascade is stable iff also $\left\lvert {kfrac(m, c)}\right\rvert \lt 1$, i.e. "
                         fr"$${abs_km(m)} \lt {c} \iff {lo} \lt K \lt {hi}.$$"
                         r"At either endpoint the pole is on the unit circle and nothing cancels it, so neither endpoint belongs to the set.")
    p["lo"], p["hi"], p["lo_in"], p["hi_in"] = lo, hi, lo_in, hi_in
    p["set_tex"] = set_tex(lo, hi, lo_in, hi_in)
    p["ends_choices"] = [{"text": t, "correct": (li == lo_in and hj == hi_in), "lo_in": li, "hi_in": hj} for t, li, hj in ENDS]
    return {"K_lo": float(lo), "K_hi": float(hi)}


def gen_single(p):
    pb, pg = random.choice(BAD), random.choice(STAB)
    b0 = random.choice([1, 2, 3, 4, 5])
    a = [F(1), -(pb + pg), pb * pg]
    lhs = fmt.tex_sum([(1, "y[n]"), (a[1], "y[n-1]"), (a[2], "y[n-2]")])
    p["system_html"] = fr"the causal LTI system described by the difference equation $$ {lhs} = {fmt.tex_sum([(b0, 'x[n]')])} + K\,x[n-1].$$"
    Kv = -b0 * pb
    p["model"] = {"b": [lin(b0), lin(0, 1)], "a": [lin(v) for v in a]}
    p["Kval"] = str(Kv)
    p["Kwrong"] = [str(b0 * pb), str(-b0 / pb), str(-b0 * pg)]
    num_after = fmt.tex_num(b0) if b0 != 1 else ""
    p["sol_html"] = (fr"$$H(z) = \dfrac{{{b0} + K z^{{-1}}}}{{{fmt.tex_sum([(1, ''), (a[1], 'z^{-1}'), (a[2], 'z^{-2}')])}}}"
                     fr" = \dfrac{{{b0} + K z^{{-1}}}}{{{fmt.tex_factor(pb)}{fmt.tex_factor(pg)}}}.$$"
                     fr"Poles ${fmt.tex_num(pb)}$ and ${fmt.tex_num(pg)}$; the causal ROC $\lvert z\rvert \gt {fmt.tex_num(abs(pb))}$ misses the unit "
                     fr"circle unless the bad pole $z = {fmt.tex_num(pb)}$ is cancelled. $K$ only moves the zero, so the zero must sit on it: "
                     fr"substitute $z^{{-1}} = 1/({fmt.tex_num(pb)}) = {fmt.tex_num(1 / pb)}$ into the numerator,"
                     fr"$${b0} + K\cdot{fmt.tex_num(1 / pb, paren=True)} = 0 \iff K = {fmt.tex_num(Kv)}.$$"
                     fr"Check: ${b0} {'+' if Kv >= 0 else '-'} {fmt.tex_num(abs(Kv))}z^{{-1}} = {num_after}{fmt.tex_factor(pb)}$, leaving "
                     fr"$H(z) = \dfrac{{{num_after or '1'}}}{{{fmt.tex_factor(pg)}}}$ with ROC $\lvert z\rvert \gt {fmt.tex_num(abs(pg))}$, which contains the unit "
                     fr"circle. Every other $K$ leaves the pole ${fmt.tex_num(pb)}$ in place: unstable.")
    return {"K_val": float(Kv)}


def gen_set(p):
    c = random.choice([2, 3, 4, 5])
    shape = random.choice(["point", "inside", "endpoint"])
    if shape == "point":
        q = random.choice([F(2), F(-2), F(3), F(3, 2), F(-3, 2), F(5, 2), F(-3)])
    elif shape == "inside":
        q = random.choice([F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(1, 4), F(2, 3), F(-2, 3)])
    else:
        q = random.choice([F(1), F(-1)])
    rhs = fmt.tex_sum([(1, "x[n]"), (-q, "x[n-1]")])
    p["system_html"] = (fr"the causal LTI system described by the difference equation "
                        fr"$$y[n] = \dfrac{{K}}{{{c}}}\,y[n-1] + {rhs}.$$")
    p["model"] = {"b": [lin(1), lin(-q)], "a": [lin(1), lin(0, F(-1, c))]}
    Kc = c * q
    K_tex, ct = fmt.tex_num(Kc), str(c)

    def iv(li, hi_, txt=None):
        return (txt or fr"${set_tex(-c, c, li, hi_)}$", {"lo": -c, "hi": c, "lo_in": li, "hi_in": hi_, "point": None})

    A, B, C, D = iv(False, False), iv(True, True), iv(False, True), iv(True, False)
    E = (fr"$-{ct} \lt K \lt {ct}$ or $K = {K_tex}$", {"lo": -c, "hi": c, "lo_in": False, "hi_in": False, "point": str(Kc)})
    Fo = (fr"$K = {K_tex}$ only", {"lo": None, "hi": None, "lo_in": False, "hi_in": False, "point": str(Kc)})
    G = (fr"$\lvert K\rvert \gt {ct}$", {"outside": True, "c": c})
    if shape == "point":
        correct, others = E, [A, Fo, B, random.choice([C, D, G])]
    elif shape == "inside":
        correct, others = A, [Fo, B, G, random.choice([C, D])]
    else:
        correct = C if q == 1 else D
        others = [A, B, D if q == 1 else C, Fo]
    opts = [(correct, True)] + [(o, False) for o in others]
    p["K_choices"] = [{"text": t, "correct": ok, "set": s} for (t, s), ok in opts]
    random.shuffle(p["K_choices"])
    p["Kc"] = str(Kc)
    p["shape"] = shape
    head = (fr"$H(z) = \dfrac{{{fmt.tex_sum([(1, ''), (-q, 'z^{-1}')])}}}{{1 - \frac{{K}}{{{c}}}z^{{-1}}}}$: one pole at "
            fr"$\tfrac{{K}}{{{c}}}$ and a zero at $z = {fmt.tex_num(q)}$. If the pole does not land on the zero it survives; the causal ROC "
            fr"$\lvert z\rvert \gt \lvert K/{c}\rvert$ contains the unit circle iff $\lvert K\rvert \lt {c}$ (at $\lvert K\rvert = {c}$ the pole is "
            r"on the unit circle: marginally stable, <b>not</b> BIBO stable). The pole lands on the zero, and cancels, for "
            fr"$\tfrac{{K}}{{{c}}} = {fmt.tex_num(q)}$, i.e. $K = {K_tex}$: then $H(z) = 1$, $h[n] = \delta[n]$, stable. ")
    if shape == "point":
        tail = (fr"That value lies outside $-{ct} \lt K \lt {ct}$, so it is an extra isolated stable value: the system is BIBO stable exactly for "
                fr"$-{ct} \lt K \lt {ct}$ or $K = {K_tex}$.")
    elif shape == "inside":
        tail = (fr"That value already lies inside $-{ct} \lt K \lt {ct}$ (the zero is inside the unit circle), so it adds nothing: the system is "
                fr"BIBO stable exactly for $-{ct} \lt K \lt {ct}$.")
    else:
        tail = (fr"That value is the endpoint $K = {K_tex}$, where the pole sits on the unit circle at $z = {fmt.tex_num(q)}$ but is cancelled; at the "
                fr"other endpoint $K = {-Kc}$ the pole $z = {fmt.tex_num(-q)}$ survives. So the system is BIBO stable exactly for "
                fr"${set_tex(-c, c, q == -1, q == 1)}$.")
    p["sol_html"] = head + tail + (r" (Practice 2 of the study site; the cancellation logic is the one of FA2024 #8b and FA2023 #8a.)")
    return {}


def generate(data):
    p = data["params"]
    fam = random.choices(["interval", "single", "set"], weights=[5, 3, 2])[0]
    p["fam_interval"], p["fam_single"], p["fam_set"] = fam == "interval", fam == "single", fam == "set"
    p["family"] = fam
    ans = {"interval": gen_interval, "single": gen_single, "set": gen_set}[fam](p)
    data["correct_answers"].update(ans)
