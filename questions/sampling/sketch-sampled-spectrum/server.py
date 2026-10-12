"""The "sketch X_d(omega)" exam problem answered by drawing (pl-sketch), Lecture 17.

X_a(Omega) is a triangle or a trapezoid (continuous, so X_d is a continuous polyline) of height A and
band edge Omega_B = 2 pi B, shown as a figure; x[n] = x_a(nT). The student draws X_d(omega) on one
period, -pi <= omega <= pi, with a polyline on a canvas whose horizontal axis is omega/pi on [-1, 1]
(neighbouring copies matter only where they reach into that period, i.e. when they overlap) and
whose vertical axis is X_d itself (the factor 1/T included).

Grading (all or nothing; R5 review): the CORNERS of X_d are graded with the point tool, whose
submission is exact (a browser polyline is sent as a smoothed fitCurve(points, 0) Bezier, which bulges
past any tight band at the apex). Stage 1, weight 0: `count` (at least as many points as corners,
the ends x = -1 and x = 1 included); on each quarter of [-1, 1] every marked point must lie neither
above nor below the exact X_d (sum of the copies, written with abs()) by more than 15 px (`less-than`
+ `greater-than`); a `match` at every corner. Stage 2, weight 1: `count >= 0`, which scores only if
all of stage 1 passed. An ungraded helper polyline is offered for sketching the shape.
Variants: b = 2BT and the trapezoid corner are multiples of 1/60 (within 0.017
of a 0.05 grid line), ramps at least 1/3 wide (no near-vertical edges), peak at <= 5/6 of the
canvas height. Models: Lecture 17 Example 2, FA2024 MT2 #5(b), SP2025 MT2 #5(a)."""

import random
from fractions import Fraction as F

from ece310 import fmt
from ece310 import w6_sampling as w

NAMED = [
    ("tri", 1, 50, 150, None, "Lecture 17, Example 2(a)"),
    ("tri", 1, 50, 100, None, "Lecture 17, Example 2(b)"),
    ("tri", 1, 50, 75, None, "Lecture 17, Example 2(c)"),
    ("tri", 1, 500, 1000, None, "FA2024 Midterm 2 #5(b) with $T = T_N$"),
    ("tri", 1, 500, 500, None, "FA2024 Midterm 2 #5(b) with $T = 2T_N$"),
    ("tri", 1, 100, 300, None, "SP2025 Midterm 2 #5(a)"),
]
B_CHOICES = [5, 10, 20, 25, 30, 40, 50, 60, 75, 80, 100, 120, 125, 150, 200, 250, 300, 400, 500, 600, 750, 800,
             1000, 1200, 1500, 2000, 2500, 3000, 4000, 5000]
B_ABOVE = [F(1, 3), F(2, 5), F(1, 2), F(3, 5), F(2, 3), F(3, 4), F(4, 5), F(5, 6)]
B_BELOW = [F(6, 5), F(7, 6), F(5, 4), F(4, 3), F(7, 5), F(3, 2), F(8, 5), F(5, 3), F(7, 4), F(2)]
FN_RANGES = [("-1", "-0.5"), ("-0.5", "0"), ("0", "0.5"), ("0.5", "1")]


def _ok(sp):
    if sp.shape == "trap":
        if sp.b - sp.c < F(1, 3) or (sp.c * 60).denominator != 1:
            return False
    if not ((sp.b * 60).denominator == 1 and sp.b >= F(1, 3)):
        return False
    xs = [x for x, _ in sp.polyline()]
    return min(b - a for a, b in zip(xs, xs[1:])) >= F(1, 20)      # corners at least 17 px apart


def _pick():
    while True:
        if random.random() < 0.03:
            shape, A, B, fs, r, src = random.choice(NAMED)
            sp = w.Spec(shape, A, B, fs, r)
            return sp, src
        shape = random.choice(["tri", "trap"])
        reg = random.choices(["above", "at", "below"], weights=[4, 1.5, 4.5])[0]
        b = random.choice(B_ABOVE) if reg == "above" else (F(1) if reg == "at" else random.choice(B_BELOW))
        B = F(random.choice(B_CHOICES))
        fs = 2 * B / b
        if fs.denominator != 1 or fs < 10:
            continue
        A = F(random.choice([1, 1, 2, 3, 4, 5]))
        r = random.choice([F(1, 4), F(1, 3), F(1, 2), F(1, 5), F(2, 5)]) if shape == "trap" else None
        if r is not None and (2 * B * r).denominator != 1:
            continue
        sp = w.Spec(shape, A, B, fs, r)
        if _ok(sp):
            return sp, None


def _num(q):
    """Coordinate for the solution string (exact enough: 10 significant digits)."""
    return f"{float(q):.10g}"


def generate(data):
    sp, src = _pick()
    poly = sp.polyline()
    peak = max(y for _, y in poly)
    ymax = w.nice_ceil(F(6, 5) * peak)
    p = data["params"]
    p["spec"] = {"shape": sp.shape, "A": str(sp.A), "B": str(sp.B), "fs": str(sp.fs),
                 "r": str(sp.r) if sp.r is not None else None}
    p["named"] = src is not None
    p["source_html"] = src or ""
    p["fig_html"] = sp.figure()
    p["words_html"] = sp.words()
    p["xa_tex"] = sp.tex_xa()
    p["T_tex"] = sp.tex_T()
    p["ymax"] = w.plain_dec(ymax)
    # the axis y = 0 sits ~20 px above the bottom of the range (~30 px above the canvas edge): a zero
    # corner dragged onto the clamped bottom edge is then clearly wrong instead of silently dropped
    # (pl-sketch's pixel -> graph conversion can land one rounding step below y_start there)
    p["ymin"] = _num(-ymax / 19)
    p["ybig"] = w.plain_dec(10 * ymax)
    p["fn"] = sp.sketch_function()
    p["fn_ranges"] = [{"range": f"{lo},{hi}", "lo": lo, "hi": hi} for lo, hi in FN_RANGES]
    pts = {x: y for x, y in poly}
    keys = []
    for x in sorted(pts):
        fb = ("No point at the corner at x = 0, or its height is wrong: remember the factor 1/T and any copies "
              "that overlap there." if x == 0 else
              "A corner of X_d is missing or misplaced: check where each copy starts and ends "
              "(the band edge in omega is Omega_B T), add the copies where they overlap, and mark both ends.")
        keys.append({"x": _num(x), "y": _num(pts[x]), "fb": fb})
    p["key_points"] = keys
    p["sol_pts"] = [{"c": f"({_num(x)},{_num(y)})"} for x, y in poly]
    p["ncorners"] = len(poly)
    p["vertices"] = [[str(x), str(y)] for x, y in poly]

    # worked solution
    reg = sp.regime()
    p["B_tex"] = w.tex_dec(sp.B)
    p["twoB_tex"] = w.tex_dec(2 * sp.B)
    p["fs_tex"] = w.tex_dec(sp.fs)
    p["cmp_tex"] = {"above": r"\gt", "at": "=", "below": r"\lt"}[reg]
    p["WB_tex"] = sp.tex_omega()
    p["b_work_tex"] = rf"{sp.tex_omega()}\cdot\tfrac{{1}}{{{w.tex_dec(sp.fs)}}} = {w.tex_pi(sp.b)}"
    p["H_work_tex"] = rf"\tfrac{{1}}{{T}}\cdot{w.tex_dec(sp.A)} = {w.tex_dec(sp.fs)}\cdot{w.tex_dec(sp.A)} = {w.tex_dec(sp.H)}"
    p["central_tex"] = w.central_tex(sp)
    p["trap_note"] = sp.shape == "trap"
    if sp.shape == "trap":
        p["c_work_tex"] = rf"{sp.tex_omega(sp.r)}\cdot\tfrac{{1}}{{{w.tex_dec(sp.fs)}}} = {w.tex_pi(sp.c)}"
    if reg == "above":
        p["regime_html"] = (rf"above the Nyquist rate ($\omega_B \lt \pi$): one copy on the period, zero (a guard band) for "
                            rf"${w.tex_pi(sp.b)} \lt \lvert\omega\rvert \le \pi$; the next copy only starts at "
                            rf"$2\pi - \omega_B = {w.tex_pi(2 - sp.b)}$")
    elif reg == "at":
        p["regime_html"] = (r"exactly at the Nyquist rate ($\omega_B = \pi$): the central copy fills the period and touches "
                            r"its neighbours at $\pm\pi$ without overlapping")
    else:
        p["regime_html"] = (rf"below the Nyquist rate ($\omega_B \gt \pi$): the copies centred at $\pm 2\pi$ reach into the "
                            rf"period down to $\lvert\omega\rvert = 2\pi - \omega_B = {w.tex_pi(2 - sp.b)}$, and there "
                            r"the copies must be <b>added</b>")
    p["corner_lines"] = [{"line": w.value_tex(sp, x, rf"X_d\!\left({w.tex_pi(x)}\right)")}
                         for x in sorted(set(pts) | {F(0)}) if x >= 0]
    p["vertex_list_tex"] = r",\ ".join(rf"\left({fmt.tex_num(x)},\ {w.tex_dec(y)}\right)" for x, y in poly)
