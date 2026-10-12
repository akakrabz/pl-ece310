"""Exact DFT arithmetic and LaTeX helpers for the dft/ questions (writer W8).

The N-point DFT is X[k] = sum_{n=0}^{N-1} x[n] e^{-j 2 pi k n / N} (course notation: X[k],
{x[n]}_{n=0}^{N-1}, <n>_N for n modulo N, as on the Midterm 2 keys FA2024 #8 and SP2023 #9).

    Surd(a, b, d)            a + b*sqrt(d), a and b Fractions, d in {1, 2, 3} (exact cos/sin values)
    cos_sin_pi(t)            (cos(t*pi), sin(t*pi)) as Surds, t a Fraction with denominator 1, 2, 3, 4 or 6
    dft_exact(x, k)          (Re X[k], Im X[k]) as Surds for an integer x of length N in {1,2,3,4,6,8,12}
    unit_exact(N, m)         (cos, -sin) of e^{-j 2 pi m / N}: the twiddle factor W_N^m as (re, im) Surds
    tex_cplx(re, im)         LaTeX of re + j im (Surds, Fractions or ints): 3, -2j, 1 - 2j, \\tfrac{\\sqrt{2}}{2} + ...
    plain_cplx(re, im)       what a student types into pl-matrix-input with allow-complex: 3, -2j, 1-2j
    tex_list(items)          \\{a,\\ b,\\ c\\} from LaTeX items (a sequence listed from n = 0, as on the exams)
    tex_mod(expr, N)         \\langle expr\\rangle_{N}
    tex_pi_frac(t)           t*pi in LaTeX (\\tfrac{3\\pi}{4}, -\\tfrac{\\pi}{2}, \\pi, 0)

Everything is exact (fractions.Fraction); floats only through .value()."""

from __future__ import annotations

import math
from fractions import Fraction

from .fmt import Q, tex_num


class Surd:
    """a + b*sqrt(d) with Fraction a, b and d in {1, 2, 3} (d = 1 means a plain rational)."""

    def __init__(self, a=0, b=0, d=1):
        a, b = Fraction(a), Fraction(b)
        if d == 1:
            a, b = a + b, Fraction(0)
        self.a, self.b, self.d = a, b, int(d)

    def _field(self, o: "Surd") -> int:
        if self.b == 0:
            return o.d if o.b != 0 else max(self.d, o.d)
        if o.b == 0 or o.d == self.d:
            return self.d
        raise ValueError("mixing sqrt(2) and sqrt(3)")

    def __add__(self, o):
        o = o if isinstance(o, Surd) else Surd(o)
        return Surd(self.a + o.a, self.b + o.b, self._field(o))

    def __sub__(self, o):
        o = o if isinstance(o, Surd) else Surd(o)
        return self + o.scale(-1)

    def __neg__(self):
        return self.scale(-1)

    def scale(self, c) -> "Surd":
        c = Fraction(c)
        return Surd(self.a * c, self.b * c, self.d)

    def __eq__(self, o):
        o = o if isinstance(o, Surd) else Surd(o)
        return self.a == o.a and (self.b == o.b) and (self.b == 0 or self.d == o.d)

    def __hash__(self):
        return hash((self.a, self.b, self.d if self.b else 1))

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_rational(self) -> bool:
        return self.b == 0

    def n_terms(self) -> int:
        return (self.a != 0) + (self.b != 0)

    def value(self) -> float:
        return float(self.a) + float(self.b) * math.sqrt(self.d)

    def _surd_tex(self, b: Fraction) -> str:
        p, q = abs(b.numerator), b.denominator
        top = rf"\sqrt{{{self.d}}}" if p == 1 else rf"{p}\sqrt{{{self.d}}}"
        return top if q == 1 else r"\tfrac{" + top + "}{" + str(q) + "}"

    def tex(self) -> str:
        parts = []
        if self.a != 0:
            parts.append(tex_num(self.a))
        if self.b != 0:
            t = self._surd_tex(self.b)
            if parts:
                parts.append(("- " if self.b < 0 else "+ ") + t)
            else:
                parts.append(("-" if self.b < 0 else "") + t)
        return " ".join(parts) if parts else "0"


# cos(t pi), sin(t pi) for t in [0, 1/2] (the rest by symmetry)
_H = Fraction(1, 2)
_FIRST = {
    Fraction(0): (Surd(1), Surd(0)),
    Fraction(1, 6): (Surd(0, _H, 3), Surd(_H)),
    Fraction(1, 4): (Surd(0, _H, 2), Surd(0, _H, 2)),
    Fraction(1, 3): (Surd(_H), Surd(0, _H, 3)),
    Fraction(1, 2): (Surd(0), Surd(1)),
}


def cos_sin_pi(t) -> tuple[Surd, Surd]:
    """(cos(t*pi), sin(t*pi)) exactly; t a Fraction whose denominator divides 4 or 6."""
    t = Fraction(t) % 2                      # [0, 2)
    sgn_s = 1
    if t > 1:                                # sin(t pi) = -sin((2 - t) pi), cos even
        t, sgn_s = 2 - t, -1
    sgn_c = 1
    if t > _H:                               # cos((1 - u) pi) = -cos(u pi), sin same
        t, sgn_c = 1 - t, -1
    if t not in _FIRST:
        raise ValueError(f"no exact cos/sin for {t}*pi")
    c, s = _FIRST[t]
    return c.scale(sgn_c), s.scale(sgn_s)


def unit_exact(N: int, m: int) -> tuple[Surd, Surd]:
    """e^{-j 2 pi m / N} = (re, im) exactly."""
    c, s = cos_sin_pi(Fraction(2 * m, N))
    return c, -s


def dft_exact(x, k: int) -> tuple[Surd, Surd]:
    """(Re X[k], Im X[k]) of the N-point DFT, N = len(x), exactly."""
    N = len(x)
    re, im = Surd(0), Surd(0)
    for n, v in enumerate(x):
        if v == 0:
            continue
        c, s = unit_exact(N, k * n)
        re = re + c.scale(v)
        im = im + s.scale(v)
    return re, im


def _as_surd(v) -> Surd:
    return v if isinstance(v, Surd) else Surd(Q(v))


def tex_cplx(re, im) -> str:
    """LaTeX of re + j*im: '3', '-2j', '1 - 2j', '\\tfrac{1}{2} - \\tfrac{\\sqrt{3}}{2}j', '1 + j\\left(...\\right)'."""
    re, im = _as_surd(re), _as_surd(im)
    rt = "" if re.is_zero() else re.tex()
    if im.is_zero():
        return rt or "0"
    if im.n_terms() == 1:
        neg = im.value() < 0
        mag = im.scale(-1) if neg else im
        mt = mag.tex()
        jt = "j" if mt == "1" else (mt + r"\,j" if "sqrt" in mt else mt + "j")
        if not rt:
            return ("-" if neg else "") + jt
        return f"{rt} {'-' if neg else '+'} {jt}"
    body = r"j\left(" + im.tex() + r"\right)"
    return f"{rt} + {body}" if rt else body


def plain_cplx(re, im) -> str:
    """Plain text of an exact Gaussian rational (no surds): '3', '-2j', '1-2j', '0.5+1.5j'."""
    re, im = Q(re), Q(im)

    def num(q):
        return str(q.numerator) if q.denominator == 1 else repr(float(q))

    if im == 0:
        return num(re)
    jt = ("" if abs(im) == 1 else num(abs(im))) + "j"
    if re == 0:
        return ("-" if im < 0 else "") + jt
    return f"{num(re)}{'-' if im < 0 else '+'}{jt}"


def tex_list(items) -> str:
    return r"\{" + r",\ ".join(items) + r"\}"


def tex_mod(expr: str, N: int) -> str:
    return r"\langle " + expr + r"\rangle_{" + str(N) + "}"


def tex_pi_frac(t) -> str:
    """t*pi in LaTeX: 0, \\pi, -\\pi, \\tfrac{\\pi}{4}, -\\tfrac{3\\pi}{4}."""
    t = Fraction(t)
    if t == 0:
        return "0"
    s = "-" if t < 0 else ""
    p, q = abs(t.numerator), t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return s + (num if q == 1 else r"\tfrac{" + num + "}{" + str(q) + "}")
