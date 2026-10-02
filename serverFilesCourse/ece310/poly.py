"""Polynomials in w = z^{-1} with exact Fraction (or complex) coefficients.

A polynomial is a list ``c`` with ``c[k]`` the coefficient of ``z^{-k}``.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from .fmt import Q


def trim(a: Sequence) -> list:
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a if a else [Fraction(0)]


def padd(a: Sequence, b: Sequence) -> list:
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def pscale(a: Sequence, c) -> list:
    return trim([c * x for x in a])


def pmul(a: Sequence, b: Sequence) -> list:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def from_roots(roots: Sequence) -> list:
    """prod_k (1 - r_k z^{-1})."""
    out: list = [Fraction(1)]
    for r in roots:
        out = pmul(out, [Fraction(1), -r])
    return out


def peval(a: Sequence, z) -> complex:
    """Evaluate sum_k a[k] z^{-k} at a (complex) point z != 0."""
    w = 1 / complex(z)
    return sum(complex(c) * w**k for k, c in enumerate(a))


def degree(a: Sequence) -> int:
    return len(trim(a)) - 1


def pdivmod(num: Sequence, den: Sequence) -> tuple[list, list]:
    """Long division in w = z^{-1}: num = q*den + r with deg r < deg den (Lecture 10)."""
    num, den = trim(num), trim(den)
    if len(num) < len(den):
        return [Fraction(0)], list(num)
    r = list(num)
    q = [Fraction(0)] * (len(num) - len(den) + 1)
    for i in range(len(num) - len(den), -1, -1):
        coef = r[i + len(den) - 1] / den[-1]
        q[i] = coef
        for j, d in enumerate(den):
            r[i + j] -= coef * d
    r = trim(r[: len(den) - 1]) if len(den) > 1 else [Fraction(0)]
    return trim(q), r


def as_fractions(a: Sequence) -> list[Fraction]:
    return [Q(x) for x in a]
