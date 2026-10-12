"""Helpers shared by the Lecture 13-14 DTFT questions (dtft/periodicity, dtft/inverse-dtft-ideal,
dtft/symmetry-and-existence, dtft/modulation).

Frequencies are exact multiples of pi held as ``fractions.Fraction`` (q stands for q*pi).

    tex_pi(q)            LaTeX of q*pi: 0, \\pi, -\\pi, 2\\pi, \\tfrac{3\\pi}{4}, -\\tfrac{\\pi}{3}
    tex_pi_n(q, var)     LaTeX of the argument q*pi*n: \\tfrac{3\\pi}{4}n, \\pi n, -\\tfrac{\\pi}{3}n
    sym_pi(q)            SymPy string of q*pi: (3/4)*pi, pi, (-1/3)*pi, 0
    plain(q)             what a student types for q: 3/4, -1/3, 2
    fold(q)              the frequency in [0, 1] (times pi) whose cosine equals cos(q*pi*n) for all n
    same_cos(a, b)       cos(a*pi*n) == cos(b*pi*n) for every integer n  (a == +-b mod 2)
    step_svg(...)        inline SVG of a piecewise-constant / piecewise-linear spectrum on [-pi, pi]
"""

from __future__ import annotations

from fractions import Fraction as F

from . import svg


def Fr(q) -> F:
    return q if isinstance(q, F) else F(q)


def tex_pi(q) -> str:
    """LaTeX of q*pi."""
    q = Fr(q)
    if q == 0:
        return "0"
    s = "-" if q < 0 else ""
    a, b = abs(q.numerator), q.denominator
    num = r"\pi" if a == 1 else f"{a}\\pi"
    return s + (num if b == 1 else r"\tfrac{" + num + "}{" + str(b) + "}")


def tex_pi_n(q, var: str = "n") -> str:
    """LaTeX of q*pi*var, e.g. \\tfrac{3\\pi}{4}n, \\pi n, 2\\pi n, -\\tfrac{\\pi}{3}n."""
    q = Fr(q)
    if q == 0:
        return "0"
    t = tex_pi(q)
    return t + (" " if t.endswith(r"\pi") else "") + var


def sym_pi(q) -> str:
    """SymPy string of q*pi."""
    q = Fr(q)
    if q == 0:
        return "0"
    if q == 1:
        return "pi"
    if q.denominator == 1:
        return f"({q.numerator})*pi"
    return f"({q.numerator}/{q.denominator})*pi"


def plain(q) -> str:
    q = Fr(q)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def mod2(q) -> F:
    """q reduced to [0, 2)."""
    q = Fr(q)
    return q - 2 * (q.numerator // (2 * q.denominator))


def fold(q) -> F:
    """The frequency in [0, 1] (times pi) with the same cosine samples as cos(q*pi*n)."""
    r = mod2(q)
    return r if r <= 1 else 2 - r


def same_cos(a, b) -> bool:
    """cos(a pi n) = cos(b pi n) for all integers n  <=>  a = +-b (mod 2)."""
    return mod2(Fr(a) - Fr(b)) == 0 or mod2(Fr(a) + Fr(b)) == 0


XLIM = (-1.12, 1.12)
PLOT_PX = 560 - 48 - 34          # svg.plot's drawing width for width=560 (left 48, right 34 px)


def tick_px(x) -> float:
    """Horizontal pixel position of x (units of pi) in a step_svg figure."""
    return 48 + (float(x) - XLIM[0]) * PLOT_PX / (XLIM[1] - XLIM[0])


def step_svg(curves, ymax, xticks, yticks, ylabel="Xd(ω)", *, guides=None, labels=None, width=560,
             height=230, title=None) -> str:
    """A spectrum on [-pi, pi] (x in units of pi): curves are lists of (x, y) points; xticks/yticks are
    numbers (multiples of pi / heights) labelled with svg.fmt_pi / svg.fmt_num.

    The height labels are written at the left edge of the frame (left of -pi), each joined to the
    leftmost point of the curve at that height by a dashed guide, so no plateau runs through them;
    the vertical axis keeps unlabelled tick marks."""
    xt = [(float(x), svg.fmt_pi(x)) for x in xticks]
    yt = [(float(y), "") for y in yticks]
    pts = [[(float(x), float(y)) for x, y in c] for c in curves]
    lab = list(labels or [])
    gd = list(guides or [])
    for y in yticks:
        yf = float(y)
        lab.append((-1.075, yf - 0.035 * float(ymax), svg.fmt_num(Fr(y))))
        xs = [x for c in pts for x, yy in c if abs(yy - yf) < 1e-12]
        if xs and min(xs) > -0.98:                     # (a plateau that starts at -pi needs no guide)
            gd.append(((-1.04, yf), (min(xs), yf)))
    return svg.plot(pts, xlim=XLIM, ylim=(0, float(ymax) * 1.18), xticks=xt, yticks=yt,
                    xlabel="ω", ylabel=ylabel, guides=gd, labels=lab, width=width, height=height,
                    title=title)
