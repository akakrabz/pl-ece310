"""Block diagrams of difference equations (Lecture 5, Fig. 4) as inline SVG (writer W2).

    from ece310 import w2_diagrams
    svg = w2_diagrams.direct_form([3, -2, 1], [Fraction(1, 2), Fraction(1, 2)], form="I")

``ff`` are the feed-forward gains (on x[n-j] in direct form I, on w[n-j] in direct form II) and
``fb`` the feedback gains on y[n-i] (direct form I) or w[n-i] (direct form II), i = 1, 2, ...,
**with the Lecture 5 sign convention**: every adder input enters with a plus sign, so the drawn
system is
    direct form I :  y[n] = sum_i fb[i-1] y[n-i] + sum_j ff[j] x[n-j]
    direct form II:  w[n] = x[n] + sum_i fb[i-1] w[n-i],   y[n] = sum_j ff[j] w[n-j]
(in the Lecture 9 form y[n] + sum a_k y[n-k] = ... the feedback gains are -a_k). A zero gain
draws no branch (its delay stays in the chain if a later tap needs it). Triangles are gains (value
written above), circles with + are adders, boxes are unit delays z^-1. Strokes use currentColor
(the page's text colour) with one accent colour for the gains; text is at least 12 px; no math
typesetting inside the SVG (Unicode minus and superscript). Pure string building."""

from __future__ import annotations

from fractions import Fraction
from html import escape

ACCENT = "#d1501f"
DY = 78           # vertical distance between taps
Y0 = 46           # y of the top row


def gain_label(g) -> str:
    """Short decimal label with a Unicode minus: 2, −0.5, 0.25, 1.5 (gains must be terminating decimals)."""
    q = Fraction(g)
    if q.denominator == 1:
        s = str(q.numerator)
    else:
        s = f"{float(q):g}"
    return s.replace("-", "−")


class _Svg:
    def __init__(self, width, height):
        self.w, self.h = width, height
        self.out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" '
                    f'style="max-width:{width}px;display:block;margin:0.6em auto;color:inherit" '
                    'font-family="sans-serif" font-size="14" role="img" aria-label="block diagram">',
                    '<defs><marker id="w2arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
                    'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>']

    def line(self, x1, y1, x2, y2, arrow=False):
        self.out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="currentColor" '
                        'stroke-width="1.5"' + (' marker-end="url(#w2arr)"' if arrow else "") + "/>")

    def dot(self, x, y):
        self.out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="currentColor"/>')

    def text(self, x, y, s, anchor="middle", size=14, italic=False, color="currentColor", weight=None):
        self.out.append(f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-size="{size}" fill="{color}"'
                        + (' font-style="italic"' if italic else "") + (f' font-weight="{weight}"' if weight else "")
                        + f">{escape(s)}</text>")

    def adder(self, x, y):
        r = 13
        self.out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="none" stroke="currentColor" stroke-width="1.6"/>')
        self.line(x - 7, y, x + 7, y)
        self.line(x, y - 7, x, y + 7)

    def delay(self, x, y):
        """z^-1 box centred at (x, y)."""
        self.out.append(f'<rect x="{x - 21:.1f}" y="{y - 14:.1f}" width="42" height="28" rx="3" fill="none" '
                        'stroke="currentColor" stroke-width="1.4"/>')
        self.out.append(f'<text x="{x - 3:.1f}" y="{y + 5:.1f}" text-anchor="middle" font-size="15" font-style="italic" '
                        'fill="currentColor">z<tspan dy="-6" font-size="12" font-style="normal">−1</tspan></text>')

    def gain(self, xc, y, g, right=True):
        """Triangle centred at xc pointing right (or left) with the gain written above it."""
        a, b = (xc - 15, xc + 15) if right else (xc + 15, xc - 15)
        self.out.append(f'<polygon points="{a:.1f},{y - 13:.1f} {a:.1f},{y + 13:.1f} {b:.1f},{y:.1f}" fill="none" '
                        f'stroke="{ACCENT}" stroke-width="1.8"/>')
        self.text(xc, y - 19, gain_label(g), size=14, color=ACCENT, weight="600")

    def svg(self):
        return "".join(self.out + ["</svg>"])


def _adders(levels, top_extra):
    """Which levels of an adder column get an adder: a lower level only if another branch joins from below
    it, the top only if it has at least two inputs (branch at 0, the chain from below, top_extra others)."""
    lv = sorted(set(levels))
    lower = [lv_ for lv_ in lv if lv_ > 0]
    has = {lv_: lv_ != max(lower) for lv_ in lower}
    has[0] = (1 if 0 in lv else 0) + (1 if lower else 0) + top_extra >= 2
    return has


def _column(s: _Svg, x, has):
    """Draw the adders of a column and the vertical wires joining them (bottom to top)."""
    for lv_, on in has.items():
        if on:
            s.adder(x, Y0 + lv_ * DY)
    pts = sorted(has, reverse=True)                        # bottom first
    for lo, hi in zip(pts[:-1], pts[1:]):
        y1 = Y0 + lo * DY - (13 if has[lo] else 0)
        y2 = Y0 + hi * DY + (13 if has[hi] else 0)
        s.line(x, y1, x, y2, arrow=has[hi])


def _into(s: _Svg, x_from, y, x_col, has_adder, from_left=True):
    """Branch wire from x_from into the column at x_col: an arrow into the adder, or a plain corner."""
    if has_adder:
        s.line(x_from, y, x_col - 13 if from_left else x_col + 13, y, arrow=True)
    else:
        s.line(x_from, y, x_col, y)


def _sig(name, k):
    return f"{name}[n]" if k == 0 else f"{name}[n−{k}]"


def direct_form(ff, fb, form: str = "I") -> str:
    ff = [Fraction(v) for v in ff]
    fb = [Fraction(v) for v in fb]
    nff, nfb = len(ff) - 1, len(fb)                      # delays needed on each side
    ffl = [j for j in range(nff + 1) if ff[j] != 0]
    fbl = [i for i in range(1, nfb + 1) if fb[i - 1] != 0]
    if form == "I":
        rows = max(nff, nfb) + 1
        width, height = 600, int(Y0 + (rows - 1) * DY + 40)
        s = _Svg(width, height)
        xin, gff, xa1, xa2, gfb, xout = 74, 160, 248, 330, 420, 506
        s.text(10, Y0 + 5, "x[n]", anchor="start", italic=True)
        s.line(46, Y0, xin, Y0)
        if nff > 0 or ff[0] == 0:
            s.dot(xin, Y0)
        left = _adders(ffl, 0)
        for j in range(nff + 1):                           # input delay line and feed-forward taps
            y = Y0 + j * DY
            if j > 0:
                s.line(xin, y - DY, xin, y - DY / 2 - 14)
                s.delay(xin, y - DY / 2)
                s.line(xin, y - DY / 2 + 14, xin, y)
                s.dot(xin, y)
                s.text(xin - 9, y + 18, _sig("x", j), anchor="end", size=12, italic=True)
            if ff[j] != 0:
                s.line(xin, y, gff - 15, y)
                s.gain(gff, y, ff[j], right=True)
                _into(s, gff + 15, y, xa1, left.get(j, False))
        _column(s, xa1, left)
        x_left_out = xa1 + 13 if left[0] else xa1
        if nfb == 0:
            s.line(x_left_out, Y0, xout + 34, Y0, arrow=True)
            s.text(xout + 40, Y0 + 5, "y[n]", anchor="start", italic=True)
            return s.svg()
        right = _adders(fbl, 1)
        _into(s, x_left_out, Y0, xa2, right[0])
        _column(s, xa2, right)
        s.line(xa2 + 13 if right[0] else xa2, Y0, xout, Y0)
        s.dot(xout, Y0)
        s.line(xout, Y0, xout + 34, Y0, arrow=True)
        s.text(xout + 40, Y0 + 5, "y[n]", anchor="start", italic=True)
        for i in range(1, nfb + 1):                        # output delay line and feedback taps
            y = Y0 + i * DY
            s.line(xout, y - DY, xout, y - DY / 2 - 14)
            s.delay(xout, y - DY / 2)
            s.line(xout, y - DY / 2 + 14, xout, y)
            s.dot(xout, y)
            s.text(xout + 9, y + 18, _sig("y", i), anchor="start", size=12, italic=True)
            if fb[i - 1] != 0:
                s.line(xout, y, gfb + 15, y)
                s.gain(gfb, y, fb[i - 1], right=False)
                _into(s, gfb - 15, y, xa2, right.get(i, False), from_left=False)
        return s.svg()

    # direct form II: one delay line in the middle
    nd = max(nff, nfb)
    rows = nd + 1
    width, height = 600, int(Y0 + (rows - 1) * DY + 40)
    s = _Svg(width, height)
    xa1, gfb, xc, gff, xa2 = 120, 205, 300, 395, 480
    left = _adders(fbl, 1)
    right = _adders(ffl, 0)
    s.text(10, Y0 + 5, "x[n]", anchor="start", italic=True)
    _into(s, 46, Y0, xa1, left[0])
    _column(s, xa1, left)
    s.line(xa1 + 13 if left[0] else xa1, Y0, xc, Y0)
    s.dot(xc, Y0)
    s.text(xc + 8, Y0 - 9, "w[n]", anchor="start", size=12, italic=True)
    for k in range(nd + 1):
        y = Y0 + k * DY
        if k > 0:
            s.line(xc, y - DY, xc, y - DY / 2 - 14)
            s.delay(xc, y - DY / 2)
            s.line(xc, y - DY / 2 + 14, xc, y)
            s.dot(xc, y)
            s.text(xc + 8, y + 20, _sig("w", k), anchor="start", size=12, italic=True)
        if 1 <= k <= nfb and fb[k - 1] != 0:
            s.line(xc, y, gfb + 15, y)
            s.gain(gfb, y, fb[k - 1], right=False)
            _into(s, gfb - 15, y, xa1, left.get(k, False), from_left=False)
        if k <= nff and ff[k] != 0:
            s.line(xc, y, gff - 15, y)
            s.gain(gff, y, ff[k], right=True)
            _into(s, gff + 15, y, xa2, right.get(k, False))
    _column(s, xa2, right)
    s.line(xa2 + 13 if right[0] else xa2, Y0, xa2 + 52, Y0, arrow=True)
    s.text(xa2 + 58, Y0 + 5, "y[n]", anchor="start", italic=True)
    return s.svg()
