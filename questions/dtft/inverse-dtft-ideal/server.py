"""Inverse DTFT of a sketched piecewise-constant spectrum (Lecture 14 Table 1 rect <-> sinc and the time
shift; study-site Lecture 13 page §8; site family "DTFTs and inverse DTFTs"; FA2024 MT2 #2, SP2021 MT2 #7, FA2019 MT2 #5).

X_d(w) = G(w) e^{-j w n0} on [-pi, pi] with G drawn as a figure (ticks at the cutoffs):
  lp     A on |w| <= wc
  hp     A on wc < |w| <= pi                  (all-pass minus low-pass, FA2024 #2)
  bp     A on w1 < |w| < w2                   (difference of two low-passes = modulated low-pass)
  stair  A1 on |w| <= wa, A2 on wa < |w| <= wb (two stacked low-passes, SP2021 #7)
Every G is a sum of centred rectangles (h_i, W_i): height h_i on |w| <= W_i pi (W = 1 is the all-pass,
whose inverse is delta[n]), so g[n] = sum_i h_i sin(W_i pi n)/(pi n) for n != 0 and
g[0] = sum_i h_i W_i, and x[n] = g[n - n0].
Asked: x[n] for n != n0 (pl-symbolic-input in n, graded NUMERICALLY at n0 +- 1..40 in grade(), since
SymPy cannot prove trigonometric identities such as the modulated form), x[n0] (exact) and one more
sample x[n0 + m] (decimal)."""

import math
import random
from fractions import Fraction as F

import prairielearn as pl
import prairielearn.sympy_utils as psu
import sympy
from ece310 import fmt
from ece310.w4_dtft import step_svg, sym_pi, tex_pi

CUTS = [F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4)]
# two cutoffs at least pi/6 apart, so their tick labels (and the step between them) stay readable
CUT_PAIRS = [(a, b) for i, a in enumerate(CUTS) for b in CUTS[i + 1:] if b - a >= F(1, 6)]
K_TEST = [k for k in range(-40, 41) if k != 0]


# ----------------------------------------------------------------------------- formatting
def shifted(n0):
    """'n', '(n - 2)', '(n + 3)' in LaTeX and SymPy."""
    if n0 == 0:
        return "n", "n"
    s = f"n {'-' if n0 > 0 else '+'} {abs(n0)}"
    return rf"\left({s}\right)", f"({s})"


def sin_tex(W, n0):
    """sin(W pi (n - n0)) in LaTeX."""
    m_tex, _ = shifted(n0)
    if n0 == 0:
        arg = tex_pi(W) + ("" if tex_pi(W).endswith("}") else " ") + "n"
    else:
        arg = tex_pi(W) + m_tex
    return rf"\sin\left({arg}\right)"


def den_tex(n0):
    m_tex, _ = shifted(n0)
    return r"\pi n" if n0 == 0 else r"\pi" + m_tex


def sym_sinc_sum(rects, n0):
    """SymPy string of sum_i h_i sin(W_i pi (n - n0)) / (pi (n - n0)); the all-pass (W = 1) is
    dropped (it is delta[n - n0], zero for n != n0)."""
    _, m = shifted(n0)
    parts = [f"{fmt.sym_num(h)}*sin({sym_pi(W)}*{m})" for h, W in rects if W != 1]
    return "(" + " + ".join(parts) + f")/(pi*{m})"


def g_value(rects, k):
    """g[k] = sum_i h_i sin(W_i pi k)/(pi k), k != 0 (float)."""
    return sum(float(h) * math.sin(float(W) * math.pi * k) for h, W in rects) / (math.pi * k)


def g0(rects):
    return sum((F(h) * F(W) for h, W in rects), F(0))


# ----------------------------------------------------------------------------- variants
def _draw():
    fam = random.choices(["lp", "hp", "bp", "stair"], weights=[25, 25, 25, 25])[0]
    if fam == "lp":
        wc, A = random.choice(CUTS), random.choice([1, 2, 3, 4])
        rects = [(F(A), wc)]
        levels = [(-wc, wc, F(A))]
        cuts, heights = [wc], [A]
    elif fam == "hp":
        wc, A = random.choice(CUTS), random.choice([1, 2, 3, 4])
        rects = [(F(A), F(1)), (F(-A), wc)]
        levels = [(F(-1), -wc, F(A)), (wc, F(1), F(A))]
        cuts, heights = [wc], [A]
    elif fam == "bp":
        w1, w2 = sorted(random.choice(CUT_PAIRS))
        A = random.choice([1, 2, 3, 4])
        rects = [(F(A), w2), (F(-A), w1)]
        levels = [(-w2, -w1, F(A)), (w1, w2, F(A))]
        cuts, heights = [w1, w2], [A]
    else:
        wa, wb = sorted(random.choice(CUT_PAIRS))
        A1, A2 = random.sample([1, 2, 3, 4], 2)
        rects = [(F(A1 - A2), wa), (F(A2), wb)]
        levels = [(-wb, -wa, F(A2)), (-wa, wa, F(A1)), (wa, wb, F(A2))]
        cuts, heights = [wa, wb], sorted({A1, A2})
    return fam, rects, levels, cuts, heights


def _curve(levels):
    """Polyline of the step function on [-1, 1] (units of pi)."""
    bps = sorted({F(-1), F(1)} | {a for a, _, _ in levels} | {b for _, b, _ in levels})
    pts = []
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        h = next((hh for a, b, hh in levels if a < mid < b), F(0))
        pts += [(x0, h), (x1, h)]
    return [(float(a), float(b)) for a, b in pts]


def generate(data):
    while True:
        fam, rects, levels, cuts, heights = _draw()
        n0 = random.choice([0, 0, 1, 2, 3, -1, -2, -3])
        ms = [m for m in (1, 2, 3) if abs(g_value(rects, m)) > 1e-3]
        if ms:
            m = random.choice(ms)
            break
    p = data["params"]
    p["family"] = fam
    p["rects"] = [[fmt.plain_num(h), fmt.plain_num(W)] for h, W in rects]
    p["levels"] = [[fmt.plain_num(a), fmt.plain_num(b), fmt.plain_num(h)] for a, b, h in levels]
    p["n0"] = n0
    p["m"] = m
    p["n1"] = n0 + m
    p["gshift_tex"] = "n" if n0 == 0 else f"n {'-' if n0 > 0 else '+'} {abs(n0)}"
    p["has_delay"] = n0 != 0
    p["phase_tex"] = (rf"e^{{-j{n0 if n0 != 1 else ''}\omega}}" if n0 > 0 else rf"e^{{j{-n0 if n0 != -1 else ''}\omega}}") if n0 else ""
    ymax = max(heights)
    xt = sorted({F(-1), F(1), F(0)} | {c for c in cuts} | {-c for c in cuts})
    p["fig_html"] = step_svg([_curve(levels)], ymax, xt, heights, ylabel="G(ω)")

    # closed forms
    _, msym = shifted(n0)
    p["x_sym"] = sym_sinc_sum(rects, n0)
    one_rect_tex = []
    for h, W in rects:
        if W == 1:
            one_rect_tex.append((h, r"\delta" + ("[n]" if n0 == 0 else f"[n {'-' if n0 > 0 else '+'} {abs(n0)}]")))
        else:
            one_rect_tex.append((h, r"\dfrac{" + sin_tex(W, n0) + "}{" + den_tex(n0) + "}"))
    p["x_full_tex"] = fmt.tex_sum(one_rect_tex)
    p["x_tex"] = fmt.tex_sum([(h, r"\dfrac{" + sin_tex(W, n0) + "}{" + den_tex(n0) + "}") for h, W in rects if W != 1])
    dec = []
    for h, W in rects:
        if W == 1:
            dec.append({"line_html": rf"height ${fmt.tex_num(h)}$ on the whole period $\leftrightarrow$ ${fmt.tex_sum([(h, r'\delta[n]')])}$"})
        else:
            dec.append({"line_html": rf"a rectangle of height ${fmt.tex_num(h)}$ on $|\omega| \le {tex_pi(W)}$ $\leftrightarrow$ "
                                     rf"${fmt.tex_sum([(h, r'\dfrac{' + sin_tex(W, 0) + r'}{\pi n}')])}$"})
    p["dec_lines"] = dec
    if fam == "bp":
        w1, w2 = cuts
        wm, dw = (w1 + w2) / 2, (w2 - w1) / 2
        A = rects[0][0]
        m_tex, _ = shifted(n0)
        cos_arg = (tex_pi(wm) + ("" if tex_pi(wm).endswith("}") else " ") + "n") if n0 == 0 else tex_pi(wm) + m_tex
        p["mod_tex"] = (rf"{fmt.tex_num(2 * A)}\cos\left({cos_arg}\right)\dfrac{{{sin_tex(dw, n0)}}}{{{den_tex(n0)}}}")
        p["wm_tex"], p["dw_tex"] = tex_pi(wm), tex_pi(dw)
        p["mod_sym"] = f"2*{fmt.sym_num(A)}*cos({sym_pi(wm)}*{msym})*sin({sym_pi(dw)}*{msym})/(pi*{msym})"
    p["is_bp"] = fam == "bp"
    p["is_hp"] = fam == "hp"
    p["is_stair"] = fam == "stair"
    if fam == "hp":
        A = rects[0][0]
        p["alt_tex"] = (rf"(-1)^{{{'n' if n0 == 0 else '(' + ('n-' if n0 > 0 else 'n+') + str(abs(n0)) + ')'}}}"
                        rf"\,{fmt.tex_num(A)}\dfrac{{{sin_tex(1 - rects[1][1], n0)}}}{{{den_tex(n0)}}}")
    # x[n0] and x[n0 + m]
    v0 = g0(rects)
    p["x0_tex"] = fmt.tex_num(v0)
    p["area_tex"] = " + ".join(rf"{fmt.tex_num(h, paren=True)}\cdot {tex_pi(2 * W)}" for h, W in rects)
    vm = g_value(rects, m)
    p["xm_terms_tex"] = fmt.tex_sum([(h, rf"\sin\left({tex_pi(W * m)}\right)") for h, W in rects if W != 1])
    p["xm_den_tex"] = rf"{m}\pi" if m != 1 else r"\pi"
    p["xm_dec"] = f"{vm:.4f}"
    p["delay_word"] = (f"delayed by {n0}" if n0 > 0 else f"advanced by {-n0}") if n0 else ""
    data["correct_answers"]["x0"] = float(v0)
    data["correct_answers"]["xm"] = vm


# ----------------------------------------------------------------------------- numerical grading of x[n]
def _true(params, n):
    rects = [(F(h), F(W)) for h, W in params["rects"]]
    return g_value(rects, n - params["n0"])


def _value(expr, nsym, k):
    v = complex(expr.subs(nsym, k).evalf(30)) if nsym is not None else complex(expr.evalf(30))
    if not (math.isfinite(v.real) and math.isfinite(v.imag)):
        raise ValueError("not finite")
    return v


def matches(expr, params):
    syms = {s.name: s for s in expr.free_symbols}
    if set(syms) - {"n"}:
        return False
    nsym = syms.get("n")
    for k in K_TEST:
        n = params["n0"] + k
        try:
            v = _value(expr, nsym, n)
        except Exception:
            return False
        t = _true(params, n)
        if not abs(v - t) <= 1e-7:
            return False
    return True


def _with_sinc(expr):
    from sympy.core.function import AppliedUndef
    return expr.replace(lambda e: isinstance(e, AppliedUndef) and e.func.__name__ == "sinc",
                        lambda e: sympy.sinc(*e.args))


def grade(data):
    sub = data["submitted_answers"].get("x")
    if sub is None or "x" in data["format_errors"]:
        return
    try:
        expr = (psu.json_to_sympy(sub, allow_complex=False) if isinstance(sub, dict)
                else psu.convert_string_to_sympy(sub, ["n"], custom_functions=["sinc"]))
        ok = matches(_with_sinc(expr), data["params"])
    except Exception:
        ok = False
    data["partial_scores"]["x"] = {"score": 1.0 if ok else 0.0, "weight": 1}
    pl.set_weighted_score_data(data)


def test(data):
    tt = data["test_type"]
    p = data["params"]
    if tt == "correct":
        data["raw_submitted_answers"]["x"] = p["x_sym"]
        data["partial_scores"]["x"] = {"score": 1.0, "weight": 1}
    elif tt == "incorrect":
        data["raw_submitted_answers"]["x"] = p["x_sym"].replace("/(pi*", "/(2*pi*")
        data["partial_scores"]["x"] = {"score": 0.0, "weight": 1}
    if tt in ("correct", "incorrect"):
        pl.set_weighted_score_data(data)
