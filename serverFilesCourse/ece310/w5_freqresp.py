"""Exact helpers for the frequency-response questions (Lectures 15-16, questions/freqresp/*).

Angles and digital frequencies are ``Fraction`` multiples of pi (``t`` stands for ``t*pi``);
principal values lie in (-1, 1]. Magnitudes are exact SymPy expressions (pi/6, 2*sqrt(3)/3,
1/2 + sqrt(2)/2, ...), so the LaTeX in the worked solutions is exact while the graded answers
are floats. A complex value is kept in polar form ``Polar(r, t)`` = r e^{j pi t}.

    principal(t)        reduce t (multiple of pi) to (-1, 1]
    tex_pi(t)           t*pi in LaTeX: \\tfrac{3\\pi}{4}, -\\tfrac{\\pi}{6}, \\pi, 0
    tex_ej(t, var)      e^{j...}: e^{-j\\frac{5\\pi}{6}}, e^{j\\frac{\\pi}{3}n} (var appended), '' / 1 for t = 0
    plain(t)            what a student types: 3/4, -1/6, 1, 0
    cos_pi(t), sin_pi(t)  exact SymPy cos(t pi), sin(t pi)
    Polar               r e^{j pi t}: mul, conj, value(), tex(), from_real(x)
    tex_real(x), tex_r(r) LaTeX of an exact real number / magnitude (sums: positive terms first)
    is_sum(x)           whether x prints as a sum (then Polar.tex() wraps it in parentheses)
    dec(x)              short decimal for the solution text ("0.5236")
"""

from __future__ import annotations

import cmath
import math
from fractions import Fraction as F

import sympy as sp


def principal(t) -> F:
    """Reduce an angle (Fraction multiple of pi) to (-1, 1]."""
    t = F(t)
    while t <= -1:
        t += 2
    while t > 1:
        t -= 2
    return t


def tex_pi(t) -> str:
    """t*pi in LaTeX (t a Fraction): 0, \\pi, -\\pi, \\tfrac{\\pi}{6}, -\\tfrac{5\\pi}{6}, 2\\pi."""
    t = F(t)
    if t == 0:
        return "0"
    s = "-" if t < 0 else ""
    p, q = abs(t.numerator), t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return s + (num if q == 1 else r"\tfrac{" + num + "}{" + str(q) + "}")


def tex_pi_frac(t) -> str:
    """t*pi for use inside an exponent: \\frac{5\\pi}{6} (never \\tfrac), with a leading minus."""
    t = F(t)
    if t == 0:
        return "0"
    s = "-" if t < 0 else ""
    p, q = abs(t.numerator), t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return s + (num if q == 1 else r"\frac{" + num + "}{" + str(q) + "}")


def tex_ej(t, var: str = "") -> str:
    """e^{j t pi var} in LaTeX; var is appended inside the exponent ('n', 'n + ...').
    t = 0 gives '' (no factor) when var is '' and '1' never appears."""
    t = F(t)
    if t == 0:
        return ""
    s = "-" if t < 0 else ""
    p, q = abs(t.numerator), t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    body = num if q == 1 else r"\frac{" + num + "}{" + str(q) + "}"
    if var and body.endswith(r"\pi"):
        var = " " + var
    return f"e^{{{s}j{body}{var}}}"


def plain(t) -> str:
    """Plain text a student would type: 3/4, -1/6, 1, 0."""
    t = F(t)
    return str(t.numerator) if t.denominator == 1 else f"{t.numerator}/{t.denominator}"


def cos_pi(t):
    return sp.cos(sp.pi * sp.Rational(F(t).numerator, F(t).denominator))


def sin_pi(t):
    return sp.sin(sp.pi * sp.Rational(F(t).numerator, F(t).denominator))


def rat(x) -> sp.Rational:
    x = F(x)
    return sp.Rational(x.numerator, x.denominator)


def _tex_term(t) -> str:
    s = sp.latex(t)
    return s.replace(r"\frac", r"\tfrac") if len(s) < 40 else s


def _terms(x):
    x = sp.nsimplify(x) if not isinstance(x, sp.Basic) else x
    return list(sp.Add.make_args(sp.expand(x)))


def is_sum(x) -> bool:
    """True if x prints as a sum of two or more terms (needs parentheses next to a factor)."""
    return len(_terms(x)) > 1


def tex_real(x) -> str:
    """LaTeX of an exact real SymPy number. A sum is printed with its positive terms first and
    explicit signs (\\sqrt{3} - 1, \\tfrac{1}{2} - \\tfrac{\\sqrt{3}}{2}, -1 - \\sqrt{2}), never
    as '-' glued to a sum."""
    terms = _terms(x)
    if len(terms) == 1:
        return _tex_term(terms[0])
    terms.sort(key=lambda t: 0 if float(t) > 0 else 1)          # stable: positive terms first
    out = _tex_term(terms[0])
    for t in terms[1:]:
        out += (" - " + _tex_term(-t)) if float(t) < 0 else (" + " + _tex_term(t))
    return out


def tex_r(r) -> str:
    """LaTeX of an exact nonnegative SymPy magnitude (a sum: positive terms first)."""
    return tex_real(r)


def dec(x, digits: int = 4) -> str:
    """Short decimal: 0.5236, 1.155, 2."""
    x = float(x)
    if abs(x - round(x)) < 1e-12:
        return str(int(round(x)))
    return f"{x:.{digits}g}"


def is_rational(r) -> bool:
    return bool(sp.nsimplify(r).is_rational)


class Polar:
    """r e^{j pi t}: r an exact nonnegative SymPy number, t a Fraction (principal value kept)."""

    def __init__(self, r, t=0):
        self.r = sp.nsimplify(r) if not isinstance(r, sp.Basic) else sp.simplify(r)
        self.t = principal(t) if self.r != 0 else F(0)

    @staticmethod
    def from_real(x) -> "Polar":
        """A real SymPy number as r e^{j pi t} with t in {0, 1}."""
        x = sp.simplify(x)
        xv = float(x)
        if abs(xv) < 1e-15:
            return Polar(0, 0)
        return Polar(sp.Abs(x), 0 if xv > 0 else 1)

    def __mul__(self, o: "Polar") -> "Polar":
        return Polar(self.r * o.r, self.t + o.t)

    def conj(self) -> "Polar":
        return Polar(self.r, -self.t)

    def scale(self, c) -> "Polar":
        """Multiply by a real SymPy/rational number."""
        return self * Polar.from_real(sp.nsimplify(c))

    @property
    def zero(self) -> bool:
        return self.r == 0

    def mag(self) -> float:
        return float(self.r)

    def ang(self) -> float:
        return float(self.t)

    def value(self) -> complex:
        return float(self.r) * cmath.exp(1j * math.pi * float(self.t))

    def tex(self) -> str:
        """r e^{j t pi}, '0' for zero, 'r' alone for t = 0."""
        if self.zero:
            return "0"
        rt = tex_r(self.r)
        e = tex_ej(self.t)
        if not e:
            return rt
        if rt == "1":
            return e
        if is_sum(self.r):
            rt = r"\left(" + rt + r"\right)"
        return rt + r"\," + e

    def __repr__(self):  # pragma: no cover - debugging aid
        return f"Polar({self.r}, {self.t})"
