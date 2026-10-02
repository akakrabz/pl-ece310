"""Difference equations (Lecture 5) and their transfer functions (Lecture 9), exact arithmetic.

Coefficient convention = Lecture 9 = ``scipy.signal.lfilter(b, a, x)``:

    a[0] y[n] + a[1] y[n-1] + ... = b[0] x[n] + b[1] x[n-1] + ...      (a[0] = 1 normally)

    H(z) = (b[0] + b[1] z^{-1} + ...) / (a[0] + a[1] z^{-1} + ...)

Import as ``from ece310 import lccde`` (this module is not re-exported by ``ece310/__init__``).
Every function accepts ints / Fractions and returns Fractions or strings (LaTeX never contains a
raw ``<``, ``>`` or ``&``; SymPy strings use ``**`` and ``(p/q)`` rationals, as in ``fmt``).
"""

from __future__ import annotations

from fractions import Fraction
from typing import Callable, Sequence

from .fmt import Q, sym_num, tex_frac, tex_sum
from .poly import pdivmod, trim


# ----------------------------------------------------------------------------- terms
def sig(name: str, k: int) -> str:
    """'y[n]', 'y[n-2]', 'x[n+1]' (sample of signal ``name`` delayed by k)."""
    if k == 0:
        return f"{name}[n]"
    return f"{name}[n-{k}]" if k > 0 else f"{name}[n+{-k}]"


def tex_lccde(b: Sequence, a: Sequence, form: str = "std", scale=1) -> str:
    """LaTeX of the difference equation with coefficients (b, a), a[0] != 0.

    form="std": Lecture 9 form   a0 y[n] + a1 y[n-1] + ... = b0 x[n] + ...   (all coefficients
                multiplied by ``scale``, e.g. scale=4 turns 1/4 coefficients into integers);
    form="rec": recursion        y[n] = -(a1/a0) y[n-1] - ... + (b0/a0) x[n] + ...
    """
    a = [Q(v) for v in a]
    b = [Q(v) for v in b]
    if form == "rec":
        a0 = a[0]
        rhs = [(-c / a0, sig("y", k)) for k, c in enumerate(a) if k > 0]
        rhs += [(c / a0, sig("x", k)) for k, c in enumerate(b)]
        return "y[n] = " + tex_sum(rhs)
    s = Q(scale)
    lhs = tex_sum((s * c, sig("y", k)) for k, c in enumerate(a))
    rhs = tex_sum((s * c, sig("x", k)) for k, c in enumerate(b))
    return lhs + " = " + rhs


# ----------------------------------------------------------------------------- polynomials
def zpow_tex_pos(m: int, var: str = "z") -> str:
    """z^{m} for m >= 0: '' (m = 0), 'z', 'z^{2}'."""
    if m == 0:
        return ""
    return var if m == 1 else f"{var}^{{{m}}}"


def tex_poly_zpos(coeffs: Sequence, top: int, var: str = "z") -> str:
    """sum_k coeffs[k] z^{top-k} (positive powers), e.g. z^{2} - \\tfrac{1}{2}z."""
    return tex_sum((Q(c), zpow_tex_pos(top - k, var)) for k, c in enumerate(coeffs))


def tex_poly_zinv(coeffs: Sequence, var: str = "z") -> str:
    """sum_k coeffs[k] z^{-k}."""
    return tex_sum((Q(c), "" if k == 0 else f"{var}^{{-{k}}}") for k, c in enumerate(coeffs))


def sym_poly(coeffs: Sequence, var: str = "z") -> str:
    """'(c0 + c1*z**(-1) + ...)' as a SymPy string."""
    terms = []
    for k, c in enumerate(coeffs):
        c = Q(c)
        if c == 0:
            continue
        terms.append(sym_num(c) if k == 0 else f"{sym_num(c)}*{var}**(-{k})")
    return "(" + (" + ".join(terms) if terms else "0") + ")"


def tex_H(b: Sequence, a: Sequence, positive: bool = False, var: str = "z") -> str:
    """H(z) = B/A in LaTeX, in powers of z^{-1} (default) or of z (positive=True: numerator and
    denominator multiplied by z^N, N = max degree). A plain polynomial when A == [1]."""
    b, a = trim([Q(v) for v in b]), trim([Q(v) for v in a])
    if positive:
        top = max(len(b), len(a)) - 1
        num, den = tex_poly_zpos(b, top, var), tex_poly_zpos(a, top, var)
    else:
        num, den = tex_poly_zinv(b, var), tex_poly_zinv(a, var)
    if len(a) == 1 and a[0] == 1:
        return num
    return tex_frac(num, den)


def sym_H(b: Sequence, a: Sequence, var: str = "z") -> str:
    """SymPy string of B(z^{-1}) / A(z^{-1})."""
    return f"{sym_poly(b, var)}/{sym_poly(a, var)}"


def pgcd(p: Sequence, q: Sequence) -> list:
    """Greatest common divisor of two polynomials in z^{-1} (exact), normalised so that its
    lowest-order nonzero coefficient is 1. [1] means "coprime" (no common factor)."""
    p, q = trim([Q(v) for v in p]), trim([Q(v) for v in q])
    while not (len(q) == 1 and q[0] == 0):
        _, r = pdivmod(p, q)
        p, q = q, trim(r)
    lead = next((c for c in p if c != 0), Fraction(1))
    return trim([c / lead for c in p])


def coprime(p: Sequence, q: Sequence) -> bool:
    """True when B and A share no factor (1 - c z^{-1}) — i.e. no pole-zero cancellation away
    from z = 0 (common powers of z^{-1} are ignored)."""
    g = pgcd(p, q)
    while len(g) > 1 and g[0] == 0:
        g = g[1:]
    return len(trim(g)) == 1


# ----------------------------------------------------------------------------- recursion
def run(b: Sequence, a: Sequence, x: Callable[[int], object] | Sequence, N: int) -> list[Fraction]:
    """y[0..N-1] of the causal recursion from initial rest (x[n] = 0 and y[n] = 0 for n < 0).
    x is a callable n -> value or a list (x[n] = list[n], zero beyond its end)."""
    a = [Q(v) for v in a]
    b = [Q(v) for v in b]
    if callable(x):
        xs = [Q(x(n)) for n in range(N)]
    else:
        xs = [Q(x[n]) if n < len(x) else Fraction(0) for n in range(N)]
    y: list[Fraction] = []
    for n in range(N):
        acc = sum((b[k] * xs[n - k] for k in range(len(b)) if n - k >= 0), Fraction(0))
        acc -= sum((a[k] * y[n - k] for k in range(1, len(a)) if n - k >= 0), Fraction(0))
        y.append(acc / a[0])
    return y


def impulse(b: Sequence, a: Sequence, N: int) -> list[Fraction]:
    """h[0..N-1] of the causal LCCDE (impulse response by recursion from rest)."""
    return run(b, a, [1], N)
