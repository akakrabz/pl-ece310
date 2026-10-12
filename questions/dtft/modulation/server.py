"""Modulation: shifted copies of a sketched spectrum, with overlaps and wrap-around (Lecture 14
frequency shift / modulation; SP2025 MT2 #2; site family "DTFTs and inverse DTFTs", Practice 4).

X_d(w) on [-pi, pi]: an ideal low-pass (height A, cutoff W) or a triangle (height A, edge W), drawn
as a figure, 2pi-periodic. y[n] is
    cos   x[n] cos(w0 n)        Y = 1/2 X(w - w0) + 1/2 X(w + w0)
    cos2  x[n] cos^2(w0 n)      Y = 1/2 X(w) + 1/4 X(w - 2w0) + 1/4 X(w + 2w0)   (2 w0 < pi)
    alt   (-1)^n x[n]           Y = X(w - pi)
with the PERIODIC X, so copies that leave [-pi, pi] come back in from the other side.
Asked: Y_d at three frequencies in [0, pi] (never on a jump of a rectangle; an overlap point and a
wrap-around point whenever the variant has them) and the band edges of the copy whose centre lies in
0 < w <= pi (before folding: the upper edge may exceed pi). All frequencies are exact multiples of pi
(Fractions); values are exact rationals."""

import random
from fractions import Fraction as F

from ece310 import fmt
from ece310.w4_dtft import plain, step_svg, tex_pi

WS = [F(1, 6), F(1, 5), F(1, 4), F(1, 3), F(2, 5), F(1, 2)]
W0_COS = [F(1, 6), F(1, 4), F(1, 3), F(5, 12), F(1, 2), F(7, 12), F(2, 3), F(3, 4), F(5, 6), F(11, 12),
          F(1, 5), F(2, 5), F(3, 5), F(4, 5)]
W0_COS2 = [F(1, 5), F(1, 4), F(1, 3), F(2, 5), F(1, 6), F(3, 8), F(5, 12), F(1, 12)]


def base(shape, A, W, w):
    """Non-periodic X(w) (w in units of pi)."""
    if abs(w) > W:
        return F(0)
    return F(A) if shape == "rect" else F(A) * (1 - abs(w) / W)


def terms_of(kind, w0):
    if kind == "cos":
        return [(F(1, 2), w0), (F(1, 2), -w0)]
    if kind == "cos2":
        return [(F(1, 2), F(0)), (F(1, 4), 2 * w0), (F(1, 4), -2 * w0)]
    return [(F(1), F(1))]


def contributions(shape, A, W, terms, w, periodic=True):
    """[(coef, shift, k, value)] of the nonzero copies at w."""
    out = []
    for c, s in terms:
        for k in ((-1, 0, 1) if periodic else (0,)):
            v = base(shape, A, W, w - s - 2 * k)
            if v != 0:
                out.append((c, s, k, c * v))
    return out


def Y(shape, A, W, terms, w, periodic=True):
    return sum((v for *_, v in contributions(shape, A, W, terms, w, periodic)), F(0))


def on_jump(shape, W, terms, w):
    if shape != "rect":
        return False
    return any(abs(w - s - 2 * k) == W for _, s in terms for k in (-1, 0, 1))


def _candidates(shape, W, terms):
    """Midpoints between consecutive band edges in [0, 1] (and the edges/centres for triangles)."""
    marks = {F(0), F(1)}
    for _, s in terms:
        for k in (-1, 0, 1):
            for e in ((s + 2 * k - W, s + 2 * k + W) + ((s + 2 * k,) if shape == "tri" else ())):
                if 0 < e < 1:
                    marks.add(e)
    marks = sorted(marks)
    mids = [a + (b - a) * r for a, b in zip(marks, marks[1:]) for r in (F(1, 4), F(1, 2), F(3, 4))]
    return sorted(set(mids) | (set(marks) if shape == "tri" else {F(0), F(1)}))


def _pick_points(shape, A, W, terms):
    """Three test frequencies in [0, pi] with at least two different values of Y_d (None if impossible)."""
    if terms == [(F(1), F(1))]:                       # (-1)^n x[n]: one copy, centred at pi
        lows = [W / 4, W / 2, 3 * W / 4]               # Y = 0 there, although X_d(w) != 0 (forgot the shift)
        if shape == "rect":                            # values (0, A, A) or (0, 0, A)
            if random.random() < 0.5:
                return [random.choice(lows), random.choice([1 - W / 4, 1 - W / 2, 1 - 3 * W / 4]), F(1)]
            return sorted(random.sample(lows, 2)) + [F(1)]
        return sorted([random.choice(lows + [1 - 3 * W / 4]), random.choice([1 - W / 2, 1 - W / 4]), F(1)])
    cands = [w for w in _candidates(shape, W, terms)
             if not on_jump(shape, W, terms, w) and Y(shape, A, W, terms, w) != 0 and w.denominator <= 24]
    overlap = [w for w in cands if len(contributions(shape, A, W, terms, w)) >= 2]
    wrap = [w for w in cands if any(k != 0 for _, _, k, _ in contributions(shape, A, W, terms, w))]
    for _ in range(30):
        pts = []
        for pool in (overlap, wrap):
            rest = [w for w in pool if w not in pts]
            if rest:
                pts.append(random.choice(rest))
        rest = [w for w in cands if w not in pts]
        random.shuffle(rest)
        for w in rest:
            if len(pts) >= 3:
                break
            if all(abs(w - q) >= F(1, 24) for q in pts):
                pts.append(w)
        if len(pts) == 3 and len({Y(shape, A, W, terms, w) for w in pts}) >= 2:
            return sorted(pts)
    return None


def generate(data):
    kind = random.choices(["cos", "cos2", "alt"], weights=[55, 33, 12])[0]
    while True:
        shape = random.choice(["rect", "rect", "tri"])
        W = random.choice(WS)
        A = random.choice([1, 2, 3, 4])
        w0 = random.choice(W0_COS if kind == "cos" else W0_COS2) if kind != "alt" else F(1)
        if kind == "cos" and w0 == W:
            continue
        terms = terms_of(kind, w0)
        pts = _pick_points(shape, A, W, terms)
        if pts is None:
            continue
        special = any(len(contributions(shape, A, W, terms, w)) >= 2
                      or any(k for _, _, k, _ in contributions(shape, A, W, terms, w)) for w in pts)
        if kind != "alt" and not special and random.random() < 0.85:
            continue                                   # mostly variants with an overlap or a wrap-around
        centre = {"cos": w0, "cos2": 2 * w0, "alt": F(1)}[kind]
        if kind == "cos2" and not 0 < centre < 1:
            continue
        break
    p = data["params"]
    p["shape"], p["kind"], p["A"], p["W"], p["w0"] = shape, kind, A, plain(W), plain(w0)
    p["is_rect"], p["is_tri"] = shape == "rect", shape == "tri"
    p["is_cos"], p["is_cos2"], p["is_alt"] = kind == "cos", kind == "cos2", kind == "alt"
    p["W_tex"], p["w0_tex"], p["A_tex"] = tex_pi(W), tex_pi(w0), str(A)
    p["two_w0_tex"] = tex_pi(2 * w0)
    p["pts"] = [plain(w) for w in pts]
    for i, w in enumerate(pts):
        p[f"wpt{i + 1}_tex"] = tex_pi(w)
    curve = ([(-1, 0), (-W, 0), (-W, A), (W, A), (W, 0), (1, 0)] if shape == "rect"
             else [(-1, 0), (-W, 0), (0, A), (W, 0), (1, 0)])
    p["fig_html"] = step_svg([curve], A, [F(-1), -W, F(0), W, F(1)], [A], ylabel="Xd(ω)")
    if kind == "cos":
        p["y_tex"] = rf"x[n]\cos\left({tex_pi(w0)}{'' if tex_pi(w0).endswith('}') else ' '}n\right)"
        p["Y_tex"] = (rf"\tfrac12 X_d\left(\omega - {tex_pi(w0)}\right) + \tfrac12 X_d\left(\omega + {tex_pi(w0)}\right)")
    elif kind == "cos2":
        p["y_tex"] = rf"x[n]\cos^2\left({tex_pi(w0)}{'' if tex_pi(w0).endswith('}') else ' '}n\right)"
        p["Y_tex"] = (rf"\tfrac12 X_d(\omega) + \tfrac14 X_d\left(\omega - {tex_pi(2 * w0)}\right)"
                      rf" + \tfrac14 X_d\left(\omega + {tex_pi(2 * w0)}\right)")
    else:
        p["y_tex"] = r"(-1)^n x[n]"
        p["Y_tex"] = r"X_d(\omega - \pi)"
    centre = {"cos": w0, "cos2": 2 * w0, "alt": F(1)}[kind]
    wl, wh = centre - W, centre + W
    p["centre_tex"] = tex_pi(centre)
    p["wl_tex"], p["wh_tex"] = tex_pi(wl), tex_pi(wh)
    p["copy_height_tex"] = fmt.tex_num({"cos": F(A, 2), "cos2": F(A, 4), "alt": F(A)}[kind])
    p["wraps"] = wh > 1
    p["wrap_lo_tex"] = tex_pi(wh - 2)

    # copies list for the solution
    copies = []
    for c, s in terms:
        lo, hi = s - W, s + W
        line = rf"height ${fmt.tex_num(c * A)}$, centred at ${tex_pi(s)}$, on $[{tex_pi(lo)}, {tex_pi(hi)}]$"
        if hi > 1:
            line += rf"; the part beyond $\pi$ re-enters on $[-\pi, {tex_pi(hi - 2)}]$"
        if lo < -1:
            line += rf"; the part below $-\pi$ re-enters on $[{tex_pi(lo + 2)}, \pi]$"
        copies.append({"line_html": line})
    p["copies"] = copies

    # evaluation lines
    evals = []
    for i, w in enumerate(pts):
        cs = contributions(shape, A, W, terms, w)
        parts = []
        for c, s, k, v in cs:
            arg = w - s - 2 * k
            where = ("" if c == 1 else rf"{fmt.tex_num(c)}\,") + rf"X_d\!\left({tex_pi(arg)}\right)"
            if k:
                where += r"_{\,\text{wrapped}}"
            parts.append((where, v))
        val = Y(shape, A, W, terms, w)
        expr = " + ".join(t for t, _ in parts)
        nums = " + ".join(fmt.tex_num(v) for _, v in parts)
        note = []
        if len(cs) >= 2:
            note.append("copies overlap here: add them")
        if any(k for _, _, k, _ in cs):
            note.append(r"a copy wraps around $\pm\pi$")
        if not parts:
            line = (rf"$Y_d\left({tex_pi(w)}\right) = 0$: no copy covers this frequency"
                    + (r" (the copy now sits around $\pm\pi$)" if kind == "alt" else ""))
        elif len(parts) == 1:
            line = rf"$Y_d\left({tex_pi(w)}\right) = {expr} = {fmt.tex_num(val)}$"
        else:
            line = rf"$Y_d\left({tex_pi(w)}\right) = {expr} = {nums} = {fmt.tex_num(val)}$"
        evals.append({"line_html": line + (f" ({'; '.join(note)})" if note else "")})
        data["correct_answers"][f"y{i + 1}"] = float(val)
    p["evals"] = evals
    data["correct_answers"]["wl"] = float(wl)
    data["correct_answers"]["wh"] = float(wh)
