"""Ideal sampling of piecewise-linear analog spectra (Lecture 17), exact arithmetic.

Shared by questions/sampling/sampled-spectrum (numbers) and questions/sampling/sketch-sampled-spectrum
(pl-sketch). An even, real, bandlimited spectrum X_a(Omega) is one of

    tri   A (1 - |Omega|/Omega_B)                      on |Omega| <= Omega_B
    trap  A, flat on |Omega| <= Omega_1 = r Omega_B, then linear down to 0 at Omega_B
    rect  A                                            on |Omega| <= Omega_B   (jumps at the edges)
    vee   A |Omega|/Omega_B                            on |Omega| <= Omega_B   (jumps at the edges)

with Omega_B = 2 pi B (B in Hz). Sampling with T = 1/f_s gives (Lecture 17, notes eq. 16)

    X_d(omega) = (1/T) sum_k X_a((omega - 2 pi k)/T).

Frequencies are handled as multiples of pi, u = omega/pi: the copy centred at 0 has band edge
b = Omega_B T/pi = 2BT and height H = A/T, the copy k is centred at u = 2k. Everything is exact
(fractions.Fraction). At a jump of X_d (rect/vee edges) the value is undefined; `value` returns
None there and the questions never ask for such a point. Where all jumps cancel (touching boxes,
Lecture 17 Example 1(b)) the value is the common one-sided limit, i.e. the "flat line" of the lecture.
"""

from __future__ import annotations

import math
from fractions import Fraction as F

from . import fmt, svg

SHAPE_WORDS = {"tri": "triangle", "trap": "trapezoid", "rect": "rectangle", "vee": "V shape"}


class Spec:
    """X_a(Omega) of the given shape, sampled at f_s (Hz)."""

    def __init__(self, shape: str, A, B, fs, r=None):
        self.shape = shape
        self.A, self.B, self.fs = F(A), F(B), F(fs)
        self.r = F(r) if r is not None else None
        if shape == "trap" and not (self.r is not None and 0 < self.r < 1):
            raise ValueError("trap needs 0 < r < 1")
        self.b = 2 * self.B / self.fs           # band edge of the central copy, in multiples of pi
        self.H = self.A * self.fs               # height of every copy (A/T)
        self.c = self.r * self.b if self.r is not None else None

    # --------------------------------------------------------------------- shape
    def unit_vertices(self):
        """Vertices of X_a/A against s = Omega/Omega_B (jumps as repeated s)."""
        if self.shape == "tri":
            return [(F(-1), F(0)), (F(0), F(1)), (F(1), F(0))]
        if self.shape == "trap":
            return [(F(-1), F(0)), (-self.r, F(1)), (self.r, F(1)), (F(1), F(0))]
        if self.shape == "rect":
            return [(F(-1), F(0)), (F(-1), F(1)), (F(1), F(1)), (F(1), F(0))]
        if self.shape == "vee":
            return [(F(-1), F(0)), (F(-1), F(1)), (F(0), F(0)), (F(1), F(1)), (F(1), F(0))]
        raise ValueError(self.shape)

    @property
    def continuous(self) -> bool:
        return self.shape in ("tri", "trap")

    def central(self):
        """Vertices of the copy centred at u = 0, in (u, X_d) units."""
        return [(self.b * s, self.H * y) for s, y in self.unit_vertices()]

    # --------------------------------------------------------------------- evaluation
    @staticmethod
    def _limit(verts, x, side):
        for (x0, y0), (x1, y1) in zip(verts, verts[1:]):
            if x0 < x1 and ((side > 0 and x0 <= x < x1) or (side < 0 and x0 < x <= x1)):
                return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        return F(0)

    def _ks(self, u):
        K = math.ceil((abs(F(u)) + self.b) / 2) + 1
        return range(-K, K + 1)

    def copy_limits(self, u, k):
        v = self.central()
        x = F(u) - 2 * k
        return self._limit(v, x, -1), self._limit(v, x, 1)

    def limits(self, u):
        """(left, right) limits of X_d at omega = u*pi."""
        u = F(u)
        left = sum((self.copy_limits(u, k)[0] for k in self._ks(u)), F(0))
        right = sum((self.copy_limits(u, k)[1] for k in self._ks(u)), F(0))
        return left, right

    def value(self, u):
        """X_d(u*pi), or None at a jump."""
        left, right = self.limits(u)
        return left if left == right else None

    def copy_value(self, u, k):
        left, right = self.copy_limits(u, k)
        return left if left == right else None

    def active_copies(self, u):
        """Copies k that are nonzero near u (either side)."""
        return [k for k in self._ks(u) if any(self.copy_limits(u, k))]

    def central_only(self, u):
        """The classic mistake: only the copy centred at 0 (None at its own jumps)."""
        return self.copy_value(u, 0)

    def max_of_copies(self, u):
        """The classic mistake: copies drawn on top of each other, not added."""
        vals = [self.copy_value(u, k) for k in self._ks(u)]
        return None if any(v is None for v in vals) else max(vals)

    def overlap(self) -> bool:
        return self.b > 1

    def regime(self) -> str:
        two_b = 2 * self.B
        return "above" if self.fs > two_b else ("at" if self.fs == two_b else "below")

    # --------------------------------------------------------------------- polyline (continuous shapes)
    def breakpoints(self, lo=F(-1), hi=F(1)):
        pts = {F(lo), F(hi)}
        for k in range(-4, 5):
            for x, _ in self.central():
                if lo < x + 2 * k < hi:
                    pts.add(x + 2 * k)
        return sorted(pts)

    def polyline(self, lo=F(-1), hi=F(1)):
        """Vertices (u, X_d) of X_d on [lo, hi], collinear interior points removed (tri/trap only)."""
        if not self.continuous:
            raise ValueError("polyline() needs a continuous spectrum")
        pts = [(x, self.value(x)) for x in self.breakpoints(lo, hi)]
        out = [pts[0]]
        for i in range(1, len(pts) - 1):
            (x0, y0), (x1, y1), (x2, y2) = out[-1], pts[i], pts[i + 1]
            if (y1 - y0) * (x2 - x1) != (y2 - y1) * (x1 - x0):
                out.append(pts[i])
        out.append(pts[-1])
        return out

    def sketch_function(self, kmax=3):
        """Python expression in x = omega/pi (abs only) for pl-sketch: the sum of the copies, each written
        as sum_j (slope change at vertex j)/2 * |x - 2k - x_j| (exact for compactly supported
        continuous piecewise-linear copies)."""
        if not self.continuous:
            raise ValueError("sketch_function() needs a continuous spectrum")
        v = self.central()
        slopes = [F(0)] + [(y1 - y0) / (x1 - x0) for (x0, y0), (x1, y1) in zip(v, v[1:])] + [F(0)]
        terms = []
        for k in range(-kmax, kmax + 1):
            for j, (xj, _) in enumerate(v):
                dm = slopes[j + 1] - slopes[j]
                if dm == 0:
                    continue
                shift = 2 * k + xj
                arg = "x" if shift == 0 else f"x - ({shift.numerator}/{shift.denominator})"
                terms.append(f"({dm.numerator}/{dm.denominator})/2*abs({arg})")
        return " + ".join(terms)

    # --------------------------------------------------------------------- LaTeX / figure
    def omega_coef(self, s=F(1)):
        """Omega = s*Omega_B as a multiple of pi (Fraction)."""
        return 2 * self.B * F(s)

    def tex_omega(self, s=F(1)):
        return tex_pi(self.omega_coef(s))

    def tex_T(self):
        return r"\tfrac{1}{" + fmt.plain_num(self.fs) + "}" if self.fs.denominator == 1 else fmt.tex_num(1 / self.fs)

    def tex_xa(self):
        """X_a(Omega) as a cases block."""
        A, WB = self.A, self.tex_omega()
        Atex = "" if A == 1 else fmt.tex_num(A)
        if self.shape == "tri":
            body = rf"{Atex}\left(1 - \dfrac{{\lvert\Omega\rvert}}{{{WB}}}\right)"
            return rf"\begin{{cases}} {body}, & \lvert\Omega\rvert \le {WB} \\ 0, & \lvert\Omega\rvert \gt {WB} \end{{cases}}"
        if self.shape == "rect":
            return rf"\begin{{cases}} {fmt.tex_num(A)}, & \lvert\Omega\rvert \le {WB} \\ 0, & \lvert\Omega\rvert \gt {WB} \end{{cases}}"
        if self.shape == "vee":
            body = rf"\dfrac{{{Atex}\lvert\Omega\rvert}}{{{WB}}}"
            return rf"\begin{{cases}} {body}, & \lvert\Omega\rvert \le {WB} \\ 0, & \lvert\Omega\rvert \gt {WB} \end{{cases}}"
        W1 = self.tex_omega(self.r)
        width = tex_pi(self.omega_coef(1 - self.r))
        body = (rf"\dfrac{{{WB} - \lvert\Omega\rvert}}{{{width}}}" if A == 1
                else rf"{Atex}\cdot\dfrac{{{WB} - \lvert\Omega\rvert}}{{{width}}}")
        return (rf"\begin{{cases}} {fmt.tex_num(A)}, & \lvert\Omega\rvert \le {W1} \\ {body}, & {W1} \le \lvert\Omega\rvert \le {WB}"
                rf" \\ 0, & \lvert\Omega\rvert \gt {WB} \end{{cases}}")

    def words(self):
        """Plain-HTML description with inline math, e.g. 'a triangle of height 2 on |Omega| <= 100 pi rad/s'."""
        WB = self.tex_omega()
        A = fmt.tex_num(self.A)
        if self.shape == "tri":
            return rf"a triangle of height ${A}$ at $\Omega = 0$, reaching $0$ at $\lvert\Omega\rvert = {WB}$ rad/s"
        if self.shape == "rect":
            return rf"a rectangle of height ${A}$ on $\lvert\Omega\rvert \le {WB}$ rad/s"
        if self.shape == "vee":
            return (rf"a V shape, $0$ at $\Omega = 0$ rising linearly to ${A}$ at $\lvert\Omega\rvert = {WB}$ rad/s, "
                    r"and $0$ beyond")
        return (rf"a trapezoid of height ${A}$, flat on $\lvert\Omega\rvert \le {self.tex_omega(self.r)}$ and falling "
                rf"linearly to $0$ at $\lvert\Omega\rvert = {WB}$ rad/s")

    def figure(self):
        """Inline SVG of X_a(Omega) with the band edge(s) and the height labelled."""
        pts = [(-1.45, 0.0)] + [(float(s), float(y * self.A)) for s, y in self.unit_vertices()] + [(1.45, 0.0)]
        ticks = [(-1.0, svg_pi(-self.omega_coef())), (1.0, svg_pi(self.omega_coef()))]
        if self.shape == "trap":
            ticks += [(-float(self.r), svg_pi(-self.omega_coef(self.r))), (float(self.r), svg_pi(self.omega_coef(self.r)))]
        top = float(self.A)
        guides = [((-1.45, top), (1.45, top))] if self.shape in ("trap", "rect", "vee") else []
        return svg.plot([pts], xlim=(-1.5, 1.5), ylim=(0, 1.3 * top), xticks=sorted(ticks),
                        yticks=[(top, svg.fmt_num(self.A))], xlabel="Ω", ylabel="Xₐ(Ω)",
                        guides=guides, width=560, height=230, max_width=520,
                        title="analog spectrum X_a(Omega)")


# ------------------------------------------------------------------------- solution LaTeX
def copy_name(k):
    if k == 0:
        return r"\text{copy at }0"
    return r"\text{copy at }" + tex_pi(2 * k)


def _jumps_at(sp, u):
    ks = [k for k in sp._ks(u) if abs(u - 2 * k) <= sp.b]
    return [k for k in ks if sp.copy_limits(u, k)[0] != sp.copy_limits(u, k)[1]]


def _side_sum_tex(sp, u, side):
    """X_d(u pi^-) or X_d(u pi^+) as the sum of the copies present on that side of u."""
    ks = [k for k in sp._ks(u)
          if (side < 0 and -sp.b < u - 2 * k <= sp.b) or (side > 0 and -sp.b <= u - 2 * k < sp.b)]
    vals = [sp.copy_limits(u, k)[0 if side < 0 else 1] for k in ks]
    sup = "-" if side < 0 else "+"
    head = rf"X_d\!\left({tex_pi(u)}^{{{sup}}}\right) = "
    if not ks:
        return head + "0"
    parts = [rf"\underbrace{{{tex_dec(v)}}}_{{{copy_name(k)}}}" for k, v in zip(ks, vals)]
    return head + " + ".join(parts) + (f" = {tex_dec(sum(vals))}" if len(parts) > 1 else "")


def value_tex(sp, u, name):
    """Display math for X_d at u*pi: the sum of the copies whose (closed) band contains u. A copy that
    is zero there because u is its band edge is labelled 'ends here'/'starts here'. Where copies jump
    at u, both one-sided sums are written out (they agree at every query point); the explanation is
    in value_note_html, as HTML text after the formula."""
    u = F(u)
    val = sp.value(u)
    head = rf"{name} = "
    ks = [k for k in sp._ks(u) if abs(u - 2 * k) <= sp.b]
    if not ks:
        return head + rf"0 \quad(\text{{no copy reaches }}\omega = {tex_pi(u)})"
    if _jumps_at(sp, u):
        return (r"\begin{gathered}" + _side_sum_tex(sp, u, -1) + r" \\ " + _side_sum_tex(sp, u, 1)
                + r"\end{gathered}")
    parts = []
    for k in ks:
        lab = copy_name(k)
        if u - 2 * k == sp.b:
            lab += r"\text{ (ends here)}"
        elif u - 2 * k == -sp.b:
            lab += r"\text{ (starts here)}"
        parts.append(rf"\underbrace{{{tex_dec(sp.copy_limits(u, k)[0])}}}_{{{lab}}}")
    if len(parts) == 1:
        return head + parts[0]
    return head + " + ".join(parts) + f" = {tex_dec(val)}"


def value_note_html(sp, u, name):
    """HTML sentence after value_tex where copies jump at u ('' elsewhere). 'Flat' only when X_d has
    zero slope on both sides."""
    u = F(u)
    if not _jumps_at(sp, u):
        return ""
    val = sp.value(u)
    d = F(1, 10 ** 9)
    flat = sp.value(u - d) == val and sp.value(u + d) == val
    how = "$X_d$ is flat here" if flat else "$X_d$ is continuous here, but not flat"
    return (rf"The copies jump at $\omega = {tex_pi(u)}$, but the jumps cancel: both one-sided limits equal "
            rf"${tex_dec(val)}$, so ${name} = {tex_dec(val)}$; {how}.")


def central_tex(sp):
    """The copy centred at 0 as a function of omega."""
    H, bt = tex_dec(sp.H), tex_pi(sp.b)
    if sp.shape == "tri":
        return rf"{H}\left(1 - \dfrac{{\lvert\omega\rvert}}{{{bt}}}\right),\quad \lvert\omega\rvert \le {bt}"
    if sp.shape == "rect":
        return rf"{H},\quad \lvert\omega\rvert \le {bt}"
    if sp.shape == "vee":
        return rf"{H}\,\dfrac{{\lvert\omega\rvert}}{{{bt}}},\quad \lvert\omega\rvert \le {bt}"
    ct, wt = tex_pi(sp.c), tex_pi(sp.b - sp.c)
    return (rf"{H}\ \text{{for }}\lvert\omega\rvert \le {ct},\quad "
            rf"{H}\,\dfrac{{{bt} - \lvert\omega\rvert}}{{{wt}}}\ \text{{for }}{ct} \le \lvert\omega\rvert \le {bt}")



# ------------------------------------------------------------------------- formatting helpers
def tex_pi(q):
    """q*pi in LaTeX: 100\\pi, \\tfrac{2\\pi}{3}, -\\pi, 0."""
    q = F(q)
    if q == 0:
        return "0"
    s = "-" if q < 0 else ""
    p, d = abs(q.numerator), q.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return s + (num if d == 1 else r"\tfrac{" + num + "}{" + str(d) + "}")


def svg_pi(q):
    """Unicode tick label for q*pi: 100π, −2π/3."""
    q = F(q)
    if q == 0:
        return "0"
    s = "−" if q < 0 else ""
    p, d = abs(q.numerator), q.denominator
    num = "π" if p == 1 else f"{p}π"
    return s + (num if d == 1 else f"{num}/{d}")


def tex_dec(q):
    """Exact number as a short decimal when it terminates (37.5, 0.25), else \\tfrac."""
    q = F(q)
    d = q.denominator
    while d % 2 == 0:
        d //= 2
    while d % 5 == 0:
        d //= 5
    if d != 1 or q.denominator == 1:
        return fmt.tex_num(q)
    s = f"{float(q):.10f}".rstrip("0").rstrip(".")
    return s


def plain_dec(q):
    """What a student would type: 37.5, 2/3, 150."""
    q = F(q)
    t = tex_dec(q)
    return t if "frac" not in t else fmt.plain_num(q)


def nice_ceil(x):
    """Smallest value >= x from {1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8} x 10^k (as Fraction)."""
    x = F(x)
    e = F(1, 1000)
    while True:
        for m in (F(1), F(6, 5), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6), F(8)):
            if m * e >= x:
                return m * e
        e *= 10
