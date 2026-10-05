"""Helpers for the inverse z-transform questions (ztransform/inverse-pfe, all-possible-rocs,
improper-long-division, complex-poles): nice pole sets, exact cover-up steps in LaTeX, and
pl-symbolic-input answers whose variable n is declared an INTEGER.

Why an integer n: PrairieLearn grades pl-symbolic-input with SymPy's ``equals``. For a plain
(complex) symbol n, forms that agree for every integer n but not for every complex n, such as
(-1)^n 2^n vs (-2)^n or (-2/3)^(-n) vs (-3/2)^n, are NOT recognised as equal. When the correct
answer is the SymPy JSON of an expression in Symbol("n", integer=True), the element gives the
student's n the same assumption (``_assumptions``) and these forms grade as equal.

All arithmetic is exact (``Fraction``); nothing here changes the shared modules.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from .fmt import Q, tex_factor, tex_num, tex_poly_zinv

F = Fraction

#: real poles used by the drills page (|p| != 1, small numerators and denominators)
NICE_POLES = [F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(1, 2), F(-1, 2), F(2, 3), F(-2, 3),
              F(3, 2), F(-3, 2), F(2), F(-2), F(3), F(-3)]


def nice(q, maxden: int = 6, maxnum: int = 12, bad_dens: Sequence[int] = (5,)) -> bool:
    """Exam-nice rational: small numerator and denominator (5ths excluded by default)."""
    q = Q(q)
    return q.denominator <= maxden and abs(q.numerator) <= maxnum and q.denominator not in bad_dens


def peval_w(coeffs: Sequence, w) -> Fraction:
    """Exact sum_k coeffs[k] w^k (a polynomial in w = z^{-1})."""
    return sum((Q(c) * Q(w) ** k for k, c in enumerate(coeffs)), F(0))


def coverup(num: Sequence, poles: Sequence, k: int) -> tuple[Fraction, Fraction, Fraction]:
    """(numerator value, remaining-denominator value, A_k) of the cover-up rule at pole k for
    X = num(z^{-1}) / prod_j (1 - p_j z^{-1}) (proper, distinct poles)."""
    w = 1 / Q(poles[k])
    nv = peval_w(num, w)
    dv = F(1)
    for j, pj in enumerate(poles):
        if j != k:
            dv *= 1 - Q(pj) * w
    return nv, dv, nv / dv


def coverup_tex(num: Sequence, poles: Sequence, k: int, label: str) -> str:
    """LaTeX of  A_k = [num / prod_{j != k}(1 - p_j z^{-1})] at z^{-1} = 1/p_k = n/d = A_k.

    The intermediate n/d is shown only when it says more than the result (not when d = 1, nor when
    n and d are coprime integers with d > 0, where n/d is just A_k written again); with a single
    pole there is no remaining factor and no fraction bar."""
    nv, dv, a = coverup(num, poles, k)
    others = "".join(tex_factor(pj) for j, pj in enumerate(poles) if j != k)
    w = 1 / Q(poles[k])
    if others:
        head = r"\left.\dfrac{" + tex_poly_zinv(num) + "}{" + others + r"}\right|_{z^{-1} = " + tex_num(w) + "}"
    else:
        head = r"\left.\left(" + tex_poly_zinv(num) + r"\right)\right|_{z^{-1} = " + tex_num(w) + "}"
    redundant = (dv == 1 or (nv.denominator == 1 and dv.denominator == 1 and dv > 0
                             and Fraction(nv.numerator, dv.numerator) == a
                             and abs(a.numerator) == abs(nv.numerator)))
    mid = "" if redundant else r" = \dfrac{" + tex_num(nv) + "}{" + tex_num(dv) + "}"
    return label + " = " + head + mid + " = " + tex_num(a)


def sym_json(expr: str):
    """pl-symbolic-input correct answer: SymPy JSON of ``expr`` with n an integer symbol."""
    from .answers import symbolic

    return symbolic(expr)
