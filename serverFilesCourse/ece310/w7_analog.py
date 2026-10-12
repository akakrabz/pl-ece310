"""Helpers for the A/D -> H_d -> D/A and ideal-reconstruction questions (sampling/ideal-reconstruction,
sampling/ad-filter-da): multiples of pi, principal values, frequencies in rad/s, a block diagram.

Digital frequencies are handled as exact multiples of pi (``Fraction`` c means c*pi).

    tex_pi(c)          LaTeX for c*pi: 0, \\pi, -\\tfrac{3\\pi}{4}, \\tfrac{9\\pi}{4}
    principal(c)       c reduced to (-1, 1]  (the principal value of c*pi)
    rad_tex(f_hz)      Omega = 2 pi f in rad/s as LaTeX: 300\\pi (for 150 Hz)
    plain(c)           what a student types: 3/4, -1/2, 2
    chain_svg(...)     inline SVG of x_a(t) -> ideal A/D -> [H_d(w)] -> ideal D/A -> y_a(t)
"""

from __future__ import annotations

from fractions import Fraction

__all__ = ["tex_pi", "principal", "rad_tex", "plain", "chain_svg", "num_tex"]


def tex_pi(c) -> str:
    """c*pi in LaTeX: 0, \\pi, -\\pi, 2\\pi, \\tfrac{\\pi}{4}, -\\tfrac{3\\pi}{4}, \\tfrac{9\\pi}{4}."""
    c = Fraction(c)
    if c == 0:
        return "0"
    s = "-" if c < 0 else ""
    p, q = abs(c.numerator), c.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return s + (num if q == 1 else r"\tfrac{" + num + "}{" + str(q) + "}")


def principal(c) -> Fraction:
    """Reduce c (a multiple of pi) to the principal interval (-1, 1]."""
    c = Fraction(c)
    while c <= -1:
        c += 2
    while c > 1:
        c -= 2
    return c


def rad_tex(f_hz) -> str:
    """Omega = 2*pi*f in rad/s as LaTeX, e.g. 150 Hz -> '300\\pi', 1/2 Hz -> '\\pi'."""
    w = 2 * Fraction(f_hz)
    if w == 1:
        return r"\pi"
    if w.denominator == 1:
        return f"{w.numerator}\\pi"
    return r"\tfrac{" + str(w.numerator) + r"\pi}{" + str(w.denominator) + "}"


def plain(c) -> str:
    """Plain-text form of a rational a student types: 3/4, -1/2, 2."""
    c = Fraction(c)
    return str(c.numerator) if c.denominator == 1 else f"{c.numerator}/{c.denominator}"


def num_tex(c) -> str:
    """A rational in LaTeX: 3, -3, \\tfrac{3}{4}."""
    c = Fraction(c)
    if c.denominator == 1:
        return str(c.numerator)
    return ("-" if c < 0 else "") + r"\tfrac{" + str(abs(c.numerator)) + "}{" + str(c.denominator) + "}"


def _box(x, y, w, h, label, sub=None):
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="none" stroke="#1f6fd1" stroke-width="1.8"/>']
    if sub:
        out.append(f'<text x="{x + w / 2}" y="{y + h / 2 + 5}" text-anchor="middle">{label}'
                   f'<tspan font-size="10" dy="4">{sub}</tspan><tspan dy="-4">(ω)</tspan></text>')
    else:
        out.append(f'<text x="{x + w / 2}" y="{y + h / 2 + 5}" text-anchor="middle">{label}</text>')
    return out


def _arrow(x1, x2, y, label=None):
    out = [f'<line x1="{x1}" y1="{y}" x2="{x2 - 7}" y2="{y}" stroke="currentColor" stroke-width="1.4"/>',
           f'<path d="M{x2 - 8},{y - 4.5} L{x2},{y} L{x2 - 8},{y + 4.5} z" fill="currentColor"/>']
    if label:
        out.append(f'<text x="{(x1 + x2) / 2}" y="{y - 9}" text-anchor="middle" font-style="italic">{label}</text>')
    return out


def chain_svg(with_filter: bool = True, title: str | None = None) -> str:
    """Block diagram: x_a(t) -> ideal A/D (T) -> x[n] -> H_d(w) -> y[n] -> ideal D/A (T) -> y_a(t).
    Without the filter: x_a(t) -> ideal A/D -> x[n] -> ideal D/A -> y_a(t)."""
    W = 640 if with_filter else 500
    H = 110
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="100%" '
           f'style="max-width:{W}px;display:block;margin:0.5em auto;color:inherit" font-family="sans-serif" '
           f'font-size="14" role="img" aria-label="{title or "block diagram of the A/D, digital filter, D/A chain"}">']
    y, bh, bw = 40, 40, 100
    x = 10
    out += _arrow(x, x + 80, y, "xₐ(t)")
    x += 80
    out += _box(x, y - bh / 2, bw, bh, "ideal A/D")
    out.append(f'<line x1="{x + bw / 2}" y1="{y + bh / 2 + 28}" x2="{x + bw / 2}" y2="{y + bh / 2 + 6}" '
               'stroke="currentColor" stroke-width="1.2"/>'
               f'<path d="M{x + bw / 2 - 4},{y + bh / 2 + 8} L{x + bw / 2},{y + bh / 2} L{x + bw / 2 + 4},{y + bh / 2 + 8} z" '
               'fill="currentColor"/>'
               f'<text x="{x + bw / 2}" y="{y + bh / 2 + 44}" text-anchor="middle" font-style="italic">T</text>')
    x += bw
    if with_filter:
        out += _arrow(x, x + 70, y, "x[n]")
        x += 70
        out += _box(x, y - bh / 2, 90, bh, "H", sub="d")
        x += 90
        out += _arrow(x, x + 70, y, "y[n]")
        x += 70
    else:
        out += _arrow(x, x + 90, y, "x[n]")
        x += 90
    out += _box(x, y - bh / 2, bw, bh, "ideal D/A")
    out.append(f'<line x1="{x + bw / 2}" y1="{y + bh / 2 + 28}" x2="{x + bw / 2}" y2="{y + bh / 2 + 6}" '
               'stroke="currentColor" stroke-width="1.2"/>'
               f'<path d="M{x + bw / 2 - 4},{y + bh / 2 + 8} L{x + bw / 2},{y + bh / 2} L{x + bw / 2 + 4},{y + bh / 2 + 8} z" '
               'fill="currentColor"/>'
               f'<text x="{x + bw / 2}" y="{y + bh / 2 + 44}" text-anchor="middle" font-style="italic">T</text>')
    x += bw
    out += _arrow(x, x + 80, y, "yₐ(t)")
    out.append("</svg>")
    return "".join(out)
