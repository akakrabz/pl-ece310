"""Helpers for the Lecture 1-2 signal questions (sampling and sequences, geometric sums, complex
exponentials): angles as exact multiples of pi, exact cos/sin at multiples of pi/4 and pi/6, and
single-term surds q*sqrt(d) for magnitudes such as 2*sqrt(2).

Angles are Fractions t meaning t*pi. Everything is exact; float() only for correct_answers."""

from __future__ import annotations

import math
from fractions import Fraction as Fr

from .fmt import Q, tex_num


# ----------------------------------------------------------------------------- angles
def principal(t) -> Fr:
    """Reduce an angle t (multiple of pi) to the principal range (-1, 1]."""
    t = Q(t)
    while t <= -1:
        t += 2
    while t > 1:
        t -= 2
    return t


def tex_pi(t) -> str:
    """t*pi in LaTeX: 0, \\pi, -\\pi, \\tfrac{3\\pi}{4}, -\\tfrac{\\pi}{6}, 2\\pi."""
    t = Q(t)
    if t == 0:
        return "0"
    s = "-" if t < 0 else ""
    p, q = abs(t.numerator), t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return s + (num if q == 1 else r"\tfrac{" + num + "}{" + str(q) + "}")


def tex_pi_n(t, var: str = "n") -> str:
    """(t*pi) var in LaTeX, written the way the notes write frequencies: \\tfrac{\\pi}{3}n, \\pi n, -\\tfrac{2\\pi}{3}n."""
    t = Q(t)
    if t == 0:
        return "0"
    return tex_pi(t) + var if t.denominator != 1 else tex_pi(t) + " " + var


def tex_phase_term(t) -> str:
    """' + \\tfrac{\\pi}{3}', ' - \\tfrac{\\pi}{4}' or '' (t = 0), to append inside cos(...)."""
    t = Q(t)
    if t == 0:
        return ""
    return (" - " if t < 0 else " + ") + tex_pi(abs(t))


def tex_ej(t) -> str:
    """e^{j t pi} with the sign in front of j: e^{j\\tfrac{\\pi}{3}}, e^{-j\\tfrac{5\\pi}{6}}; '' for t = 0."""
    t = Q(t)
    if t == 0:
        return ""
    return "e^{" + ("-" if t < 0 else "") + "j" + tex_pi(abs(t)) + "}"


def plain_frac(x) -> str:
    """3, -3, 3/4: what a student types."""
    q = Q(x)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


# ----------------------------------------------------------------------------- surds q*sqrt(d)
def _squarefree(d: int) -> tuple[int, int]:
    """d = s^2 * r with r squarefree -> (s, r)."""
    s, r, k = 1, d, 2
    while k * k <= r:
        while r % (k * k) == 0:
            r //= k * k
            s *= k
        k += 1
    return s, r


class Rad:
    """q * sqrt(d): q a Fraction, d a squarefree positive integer (d = 1 for rationals)."""

    def __init__(self, q, d: int = 1):
        q = Q(q)
        s, r = _squarefree(int(d))
        self.q, self.d = q * s, r
        if self.q == 0:
            self.d = 1

    def __mul__(self, o):
        o = o if isinstance(o, Rad) else Rad(o)
        return Rad(self.q * o.q, self.d * o.d)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, Rad) else Rad(o)
        # 1/sqrt(d) = sqrt(d)/d
        return Rad(self.q / o.q / o.d, self.d * o.d)

    def __neg__(self):
        return Rad(-self.q, self.d)

    def __pow__(self, k: int):
        out = Rad(1)
        base = self if k >= 0 else Rad(1) / self
        for _ in range(abs(k)):
            out = out * base
        return out

    def __eq__(self, o):
        o = o if isinstance(o, Rad) else Rad(o)
        return self.q == o.q and self.d == o.d

    def __hash__(self):
        return hash((self.q, self.d))

    def is_zero(self) -> bool:
        return self.q == 0

    def is_rational(self) -> bool:
        return self.d == 1

    def sq(self) -> Fr:
        """The square, a rational number."""
        return self.q * self.q * self.d

    def __float__(self):
        return float(self.q) * math.sqrt(self.d)

    def abs(self):
        return Rad(abs(self.q), self.d)

    def tex(self, paren: bool = False) -> str:
        """\\tfrac{3\\sqrt{2}}{2}, -2\\sqrt{3}, \\tfrac{1}{2}, 5."""
        if self.d == 1:
            return tex_num(self.q, paren=paren)
        p, q = abs(self.q.numerator), self.q.denominator
        top = rf"\sqrt{{{self.d}}}" if p == 1 else rf"{p}\sqrt{{{self.d}}}"
        s = top if q == 1 else r"\tfrac{" + top + "}{" + str(q) + "}"
        s = ("-" if self.q < 0 else "") + s
        if paren and self.q < 0:
            s = r"\left(" + s + r"\right)"
        return s

    def coef_tex(self) -> str:
        """As a coefficient in front of a function: '' for 1, '-' for -1, else tex()."""
        if self == Rad(1):
            return ""
        if self == Rad(-1):
            return "-"
        return self.tex()

    def plain(self) -> str:
        """A string a student could type into pl-number-input: exact fraction if rational, else 6 s.f."""
        if self.d == 1:
            return plain_frac(self.q)
        return f"{float(self):.6g}"


# exact cos(t*pi) for t a multiple of 1/4 or 1/6
_COS = {Fr(0): Rad(1), Fr(1, 6): Rad(Fr(1, 2), 3), Fr(1, 4): Rad(Fr(1, 2), 2), Fr(1, 3): Rad(Fr(1, 2)),
        Fr(1, 2): Rad(0), Fr(2, 3): Rad(Fr(-1, 2)), Fr(3, 4): Rad(Fr(-1, 2), 2), Fr(5, 6): Rad(Fr(-1, 2), 3),
        Fr(1): Rad(-1)}


def cos_pi(t) -> Rad:
    """cos(t*pi) exactly; t must be a multiple of 1/4 or 1/6."""
    t = abs(principal(t))
    if t not in _COS:
        raise ValueError(f"cos({t} pi) is not tabulated")
    return _COS[t]


def sin_pi(t) -> Rad:
    """sin(t*pi) = cos((t - 1/2) pi) exactly."""
    return cos_pi(Q(t) - Fr(1, 2))
