"""Small inline-SVG plots for question text: spectra, stem plots, pole-zero diagrams.

Every function returns one ``<svg>`` string to put in ``params["..._html"]`` and insert with triple
braces (``{{{params.fig_html}}}``). Pure string building, no dependencies, deterministic output.

Conventions (keep figures consistent across the course):
- labels are plain text with Unicode (``"−π"``, ``"2π/3"``, ``"ω"``, ``"X_d(ω)"`` is written
  ``"Xₓ"``-free: use ``"Xd(ω)"`` or pass an italic label); math is NOT typeset inside SVG;
- dark strokes (``currentColor``) on a transparent background, so the figure follows the page's
  text colour; accent colour only for the data;
- the caller decides every tick and label: a figure must show what students need to read off and
  nothing that another part of the question asks for.

Functions
    plot(curves, xlim, ylim, xticks, yticks, ...)   piecewise-linear curves, impulses (arrows), dashed guides
    stem(n, x, ...)                                  stem plot of a sequence (or |X[k]|)
    pole_zero(poles, zeros, ...)                     z-plane with unit circle, × poles, ○ zeros, multiplicities
    fmt_pi(fraction)                                 "2π/3", "−π", "0", "π/4" for a multiple of π

Example::

    from ece310 import svg
    fig = svg.plot([[(-1, 0), (-0.5, 0), (-0.5, 2), (0.5, 2), (0.5, 0), (1, 0)]],
                   xlim=(-1, 1), ylim=(0, 2.5),
                   xticks=[(-1, "−π"), (-0.5, "−π/2"), (0, "0"), (0.5, "π/2"), (1, "π")],
                   yticks=[(2, "2")], xlabel="ω", ylabel="X_d(ω)")
"""

from __future__ import annotations

from fractions import Fraction
from html import escape

ACCENT = "#1f6fd1"
ACCENT2 = "#d1501f"

__all__ = ["plot", "stem", "pole_zero", "fmt_pi", "fmt_num"]


def fmt_num(v) -> str:
    """Short label for a number: integers plain, simple fractions as a/b, else 3 significant digits.
    Uses the Unicode minus sign."""
    if isinstance(v, Fraction):
        s = str(v.numerator) if v.denominator == 1 else f"{abs(v.numerator)}/{v.denominator}"
        if v.denominator != 1 and v < 0:
            s = "-" + s
    elif isinstance(v, int) or (isinstance(v, float) and v.is_integer()):
        s = str(int(v))
    else:
        f = Fraction(v).limit_denominator(64)
        if abs(float(f) - v) < 1e-9:
            return fmt_num(f)
        s = f"{v:.3g}"
    return s.replace("-", "−")


def fmt_pi(c) -> str:
    """Label for c·π with c a Fraction/int/float: 0, π, −π, 2π/3, −π/4, 3π."""
    c = Fraction(c).limit_denominator(1000)
    if c == 0:
        return "0"
    sign = "−" if c < 0 else ""
    a, b = abs(c.numerator), c.denominator
    num = "π" if a == 1 else f"{a}π"
    return sign + (num if b == 1 else f"{num}/{b}")


class _Frame:
    def __init__(self, xlim, ylim, width, height, left, right, top, bottom):
        self.x0, self.x1 = map(float, xlim)
        self.y0, self.y1 = map(float, ylim)
        self.w, self.h = width, height
        self.l, self.r, self.t, self.b = left, right, top, bottom

    def X(self, x) -> float:
        return self.l + (float(x) - self.x0) * (self.w - self.l - self.r) / (self.x1 - self.x0)

    def Y(self, y) -> float:
        return self.h - self.b - (float(y) - self.y0) * (self.h - self.t - self.b) / (self.y1 - self.y0)


def _open(width, height, max_width, title) -> list[str]:
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" '
           f'style="max-width:{max_width}px;display:block;margin:0.5em auto;color:inherit" '
           f'font-family="sans-serif" font-size="13" role="img"'
           + (f' aria-label="{escape(title)}"' if title else "") + ">"]
    out.append('<defs><marker id="arrowhead" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="7" '
               'markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/>'
               '</marker><marker id="arrowhead-accent" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="5" '
               f'markerHeight="5" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{ACCENT}"/>'
               '</marker></defs>')
    return out


def _axes(f: _Frame, xticks, yticks, xlabel, ylabel, grid: bool, xlabel_row: float | None = None) -> list[str]:
    """Axes, ticks and labels. xlabel_row: pixel y of the x tick labels (default just under the axis)."""
    out: list[str] = []
    yz = f.Y(0) if f.y0 <= 0 <= f.y1 else f.Y(f.y0)
    xz = f.X(0) if f.x0 <= 0 <= f.x1 else f.X(f.x0)
    # axes with arrows
    out.append(f'<line x1="{f.l - 6:.1f}" y1="{yz:.1f}" x2="{f.w - f.r + 14:.1f}" y2="{yz:.1f}" '
               'stroke="currentColor" stroke-width="1.2" marker-end="url(#arrowhead)"/>')
    out.append(f'<line x1="{xz:.1f}" y1="{f.h - f.b + 4:.1f}" x2="{xz:.1f}" y2="{f.t - 12:.1f}" '
               'stroke="currentColor" stroke-width="1.2" marker-end="url(#arrowhead)"/>')
    if xlabel:
        out.append(f'<text x="{f.w - f.r + 18:.1f}" y="{yz + 4:.1f}" font-style="italic">{escape(xlabel)}</text>')
    if ylabel:
        out.append(f'<text x="{xz + 6:.1f}" y="{f.t - 14:.1f}" font-style="italic">{escape(ylabel)}</text>')
    for v, lab in xticks or []:
        x = f.X(v)
        if grid:
            out.append(f'<line x1="{x:.1f}" y1="{f.t:.1f}" x2="{x:.1f}" y2="{f.h - f.b:.1f}" stroke="currentColor" '
                       'stroke-opacity="0.15" stroke-width="1"/>')
        out.append(f'<line x1="{x:.1f}" y1="{yz - 4:.1f}" x2="{x:.1f}" y2="{yz + 4:.1f}" stroke="currentColor"/>')
        ty = yz + 19 if xlabel_row is None else xlabel_row
        out.append(f'<text x="{x:.1f}" y="{ty:.1f}" text-anchor="middle">{escape(lab)}</text>')
    for v, lab in yticks or []:
        y = f.Y(v)
        if grid:
            out.append(f'<line x1="{f.l:.1f}" y1="{y:.1f}" x2="{f.w - f.r:.1f}" y2="{y:.1f}" stroke="currentColor" '
                       'stroke-opacity="0.15" stroke-width="1"/>')
        out.append(f'<line x1="{xz - 4:.1f}" y1="{y:.1f}" x2="{xz + 4:.1f}" y2="{y:.1f}" stroke="currentColor"/>')
        out.append(f'<text x="{xz - 8:.1f}" y="{y + 4:.1f}" text-anchor="end">{escape(lab)}</text>')
    return out


def plot(curves, xlim, ylim, xticks=None, yticks=None, xlabel="", ylabel="", *, impulses=None,
         guides=None, labels=None, dashed=None, colors=None, width=560, height=250, max_width=560,
         grid=False, title=None) -> str:
    """Piecewise-linear curves (each a list of (x, y) points, drawn as one polyline).

    impulses: list of (x, height, label) drawn as arrows from 0 (label written next to the head, may be "").
    guides:   list of ((x1, y1), (x2, y2)) thin dashed helper lines (e.g. a level or an edge).
    labels:   list of (x, y, text) free text in data coordinates.
    dashed:   list of booleans, one per curve (dashed copies, e.g. aliases drawn but not summed).
    colors:   list of colours per curve (default ACCENT).
    """
    f = _Frame(xlim, ylim, width, height, left=48, right=34, top=30, bottom=30)
    out = _open(width, height, max_width, title)
    out += _axes(f, xticks, yticks, xlabel, ylabel, grid)
    for g in guides or []:
        (a, b), (c, d) = g
        out.append(f'<line x1="{f.X(a):.1f}" y1="{f.Y(b):.1f}" x2="{f.X(c):.1f}" y2="{f.Y(d):.1f}" '
                   'stroke="currentColor" stroke-opacity="0.55" stroke-dasharray="4 4" stroke-width="1"/>')
    for i, pts in enumerate(curves or []):
        col = (colors or [ACCENT] * len(curves))[i]
        dash = ' stroke-dasharray="7 5"' if dashed and dashed[i] else ""
        d = " ".join(f"{f.X(x):.1f},{f.Y(y):.1f}" for x, y in pts)
        out.append(f'<polyline points="{d}" fill="none" stroke="{col}" stroke-width="2.4" '
                   f'stroke-linejoin="round"{dash}/>')
    for x, hgt, lab in impulses or []:
        y0, y1 = f.Y(0), f.Y(hgt)
        out.append(f'<line x1="{f.X(x):.1f}" y1="{y0:.1f}" x2="{f.X(x):.1f}" y2="{y1:.1f}" stroke="{ACCENT}" '
                   'stroke-width="2.4" marker-end="url(#arrowhead-accent)"/>')
        if lab:
            out.append(f'<text x="{f.X(x) + 6:.1f}" y="{y1 - 2:.1f}">{escape(lab)}</text>')
    for x, y, text in labels or []:
        out.append(f'<text x="{f.X(x):.1f}" y="{f.Y(y):.1f}" text-anchor="middle">{escape(text)}</text>')
    out.append("</svg>")
    return "".join(out)


def stem(n, x, xlim=None, ylim=None, xticks=None, yticks=None, xlabel="n", ylabel="", *, values=False,
         width=560, height=230, max_width=560, grid=False, title=None) -> str:
    """Stem plot of the samples x[n] at integer positions n (lists of equal length).

    xticks default: every n (or every k-th when there are many); values=True writes each nonzero
    sample's value above (below) its stem. With negative samples the index labels are written in a
    row below the lowest stem (the figure grows by 22 px). Callers choose xticks so that labels of
    neighbouring indices do not collide (two-digit labels need about 18 px per index)."""
    n = list(n)
    x = [float(v) for v in x]
    lo, hi = min(n), max(n)
    xlim = xlim or (lo - 1, hi + 1)
    top = max([0.0] + x)
    bot = min([0.0] + x)
    pad = 0.18 * max(1e-9, top - bot)
    ylim = ylim or (bot - (pad if bot < 0 else 0), top + pad)
    if xticks is None:
        step = 1 if hi - lo <= 16 else (2 if hi - lo <= 32 else 4)
        xticks = [(k, fmt_num(k)) for k in range(lo, hi + 1) if (k - lo) % step == 0]
    neg = min(x) < 0
    f = _Frame(xlim, ylim, width, height + (22 if neg else 0), left=48, right=34, top=30, bottom=30 + (22 if neg else 0))
    out = _open(width, height + (22 if neg else 0), max_width, title)
    # with negative samples the index labels go below the lowest stem, not across the stems
    out += _axes(f, xticks, yticks, xlabel, ylabel, grid, xlabel_row=(f.Y(ylim[0]) + 18) if neg else None)
    for k, v in zip(n, x):
        if v == 0:
            out.append(f'<circle cx="{f.X(k):.1f}" cy="{f.Y(0):.1f}" r="3.2" fill="{ACCENT}"/>')
            continue
        out.append(f'<line x1="{f.X(k):.1f}" y1="{f.Y(0):.1f}" x2="{f.X(k):.1f}" y2="{f.Y(v):.1f}" '
                   f'stroke="{ACCENT}" stroke-width="2"/>')
        out.append(f'<circle cx="{f.X(k):.1f}" cy="{f.Y(v):.1f}" r="4" fill="{ACCENT}"/>')
        if values:
            dy = -9 if v > 0 else 18
            out.append(f'<text x="{f.X(k):.1f}" y="{f.Y(v) + dy:.1f}" text-anchor="middle" font-size="12">'
                       f'{escape(fmt_num(v))}</text>')
    out.append("</svg>")
    return "".join(out)


def pole_zero(poles, zeros, *, lim=1.6, ticks=None, iticks=None, width=300, max_width=320, title=None,
              show_unit_circle=True) -> str:
    """z-plane plot: unit circle, poles as ×, zeros as ○ (complex numbers or (re, im) pairs).

    A repeated location is drawn once with its multiplicity written next to it ("2").
    ticks:  list of (value, label) on the real axis (default ±1, ±1/2);
    iticks: list of (value, label) on the imaginary axis (default ±j, ±j/2)."""
    def as_c(p):
        return complex(*p) if isinstance(p, (tuple, list)) else complex(p)

    def group(points):
        seen: list[tuple[complex, int]] = []
        for p in map(as_c, points):
            for i, (q, m) in enumerate(seen):
                if abs(q - p) < 1e-9:
                    seen[i] = (q, m + 1)
                    break
            else:
                seen.append((p, 1))
        return seen

    height = width
    f = _Frame((-lim, lim), (-lim, lim), width, height, left=26, right=26, top=26, bottom=26)
    out = _open(width, height, max_width, title)
    ticks = ticks if ticks is not None else [(-1, "−1"), (-0.5, "−1/2"), (0.5, "1/2"), (1, "1")]
    iticks = iticks if iticks is not None else [(-1, "−j"), (-0.5, "−j/2"), (0.5, "j/2"), (1, "j")]
    xz, yz = f.X(0), f.Y(0)
    out.append(f'<line x1="{f.X(-lim):.1f}" y1="{yz:.1f}" x2="{f.X(lim) + 8:.1f}" y2="{yz:.1f}" '
               'stroke="currentColor" stroke-width="1.1" marker-end="url(#arrowhead)"/>')
    out.append(f'<line x1="{xz:.1f}" y1="{f.Y(-lim):.1f}" x2="{xz:.1f}" y2="{f.Y(lim) - 8:.1f}" '
               'stroke="currentColor" stroke-width="1.1" marker-end="url(#arrowhead)"/>')
    out.append(f'<text x="{f.X(lim) - 6:.1f}" y="{yz - 8:.1f}" text-anchor="end" font-style="italic">Re z</text>')
    out.append(f'<text x="{xz + 8:.1f}" y="{f.Y(lim) + 4:.1f}" font-style="italic">Im z</text>')
    if show_unit_circle:
        rpx = f.X(1) - f.X(0)
        out.append(f'<circle cx="{xz:.1f}" cy="{yz:.1f}" r="{rpx:.1f}" fill="none" stroke="currentColor" '
                   'stroke-opacity="0.6" stroke-dasharray="5 4" stroke-width="1.1"/>')
    for v, lab in ticks:
        out.append(f'<line x1="{f.X(v):.1f}" y1="{yz - 3:.1f}" x2="{f.X(v):.1f}" y2="{yz + 3:.1f}" stroke="currentColor"/>')
        out.append(f'<text x="{f.X(v):.1f}" y="{yz + 21:.1f}" text-anchor="middle" font-size="12">{escape(lab)}</text>')
    for v, lab in iticks:
        out.append(f'<line x1="{xz - 3:.1f}" y1="{f.Y(v):.1f}" x2="{xz + 3:.1f}" y2="{f.Y(v):.1f}" stroke="currentColor"/>')
        out.append(f'<text x="{xz - 6:.1f}" y="{f.Y(v) + 4:.1f}" text-anchor="end" font-size="12">{escape(lab)}</text>')
    s = 7
    for p, m in group(poles):
        cx, cy = f.X(p.real), f.Y(p.imag)
        out.append(f'<path d="M{cx - s:.1f},{cy - s:.1f} L{cx + s:.1f},{cy + s:.1f} M{cx - s:.1f},{cy + s:.1f} '
                   f'L{cx + s:.1f},{cy - s:.1f}" stroke="{ACCENT2}" stroke-width="2.4"/>')
        if m > 1:
            out.append(f'<text x="{cx + 9:.1f}" y="{cy - 7:.1f}" font-size="12" fill="{ACCENT2}">{m}</text>')
    for p, m in group(zeros):
        cx, cy = f.X(p.real), f.Y(p.imag)
        out.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{s:.1f}" fill="none" stroke="{ACCENT}" stroke-width="2.4"/>')
        if m > 1:
            out.append(f'<text x="{cx + 9:.1f}" y="{cy - 7:.1f}" font-size="12" fill="{ACCENT}">{m}</text>')
    out.append("</svg>")
    return "".join(out)
