"""LaTeX and SymPy-string formatting for exact rational numbers, sums, powers and sequences.

Conventions (shared by every ECE 310 question):
  * numbers are ``fractions.Fraction`` (ints are accepted everywhere);
  * LaTeX never contains a raw ``<``, ``>`` or ``&``: use ``\\lt``, ``\\gt``, ``\\le``, ``\\ge``
    so the strings are safe in HTML text and attributes;
  * SymPy strings are fully parenthesised and use ``**`` (what ``pl-symbolic-input`` parses),
    rationals as ``(3/4)``, never floats.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

Number = int | Fraction


def Q(x) -> Fraction:
    """Coerce int / str ("3/4") / Fraction to Fraction. Floats are refused (use exact values)."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, bool):
        raise TypeError("bool is not a number here")
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, str):
        return Fraction(x)
    raise TypeError(f"expected int/Fraction/str, got {type(x).__name__}: {x!r}")


# ----------------------------------------------------------------------------- numbers
def tex_num(x: Number, paren: bool = False, dfrac: bool = False) -> str:
    """3, -3, \\tfrac{3}{4}, -\\tfrac{3}{4}; paren=True wraps negatives and fractions in \\left(\\right)."""
    q = Q(x)
    f = r"\dfrac" if dfrac else r"\tfrac"
    if q.denominator == 1:
        s = str(q.numerator)
    else:
        s = ("-" if q < 0 else "") + f"{f}{{{abs(q.numerator)}}}{{{q.denominator}}}"
    if paren and (q < 0 or q.denominator != 1):
        s = r"\left(" + s + r"\right)"
    return s


def sym_num(x: Number) -> str:
    """SymPy string: 3, (-3), (3/4), (-3/4)."""
    q = Q(x)
    if q.denominator == 1:
        return str(q.numerator) if q >= 0 else f"({q.numerator})"
    return f"({q.numerator}/{q.denominator})"


def plain_num(x: Number) -> str:
    """Plain-text form a student would type: 3, -3, 3/4, -3/4."""
    q = Q(x)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


# ----------------------------------------------------------------------------- sums
def tex_sum(terms: Iterable[tuple[Number, str]]) -> str:
    """LaTeX of sum(coef * body). body == '' is a constant term. Coefficients +-1 are elided.
    Returns '0' when every coefficient is zero."""
    out = ""
    for c, body in terms:
        c = Q(c)
        if c == 0:
            continue
        a = abs(c)
        if body:
            mag = "" if a == 1 else tex_num(a)
            if mag and body[0].isdigit():
                sep = r" \cdot "  # 2 \cdot 3^{n}, never "2\,3^{n}" (reads as 23^n)
            elif mag and body.startswith((r"\tfrac", r"\dfrac", r"\left", "(")):
                sep = r"\,"
            else:
                sep = ""
            t = mag + sep + body
        else:
            t = tex_num(a)
        if not out:
            out = ("-" if c < 0 else "") + t
        else:
            out += (" - " if c < 0 else " + ") + t
    return out or "0"


def sym_sum(terms: Iterable[tuple[Number, str]]) -> str:
    """SymPy string of sum(coef * body); body is a SymPy string ('' for a constant). '0' if empty."""
    parts = []
    for c, body in terms:
        c = Q(c)
        if c == 0:
            continue
        parts.append(f"{sym_num(c)}*({body})" if body else sym_num(c))
    return " + ".join(parts) if parts else "0"


# ----------------------------------------------------------------------------- powers
def tex_pow(base: Number | complex, exponent: str = "n") -> str:
    """base^{exponent} in LaTeX: 2^{n}, (-3)^{n}, \\left(\\tfrac{1}{2}\\right)^{n}; base 1 -> '1'."""
    q = Q(base)
    if q == 1:
        return "1"
    if q.denominator == 1 and q > 0:
        return f"{q.numerator}^{{{exponent}}}"
    return tex_num(q, paren=True) + f"^{{{exponent}}}"


def sym_pow(base: Number, exponent: str = "n") -> str:
    """SymPy string: (1/2)**(n), (-3)**(n), 2**(n); base 1 -> '1'."""
    q = Q(base)
    if q == 1:
        return "1"
    return f"{sym_num(q) if (q < 0 or q.denominator != 1) else str(q.numerator)}**({exponent})"


def zpow_tex(k: int, var: str = "z") -> str:
    """z^{-k}: '' for k = 0, 'z^{-1}', 'z^{-2}', 'z' for k = -1, 'z^{3}' for k = -3."""
    if k == 0:
        return ""
    if k == -1:
        return var
    return f"{var}^{{{-k}}}"


def zpow_sym(k: int, var: str = "z") -> str:
    """SymPy string for z^{-k} ('1' for k = 0)."""
    if k == 0:
        return "1"
    return f"{var}**({-k})"


# ----------------------------------------------------------------------------- polynomials in z^-1
def tex_poly_zinv(coeffs: Sequence[Number], var: str = "z", shift: int = 0) -> str:
    """sum_k coeffs[k] z^{-(k+shift)} in LaTeX, e.g. 1 - \\tfrac{1}{2}z^{-1}."""
    return tex_sum((c, zpow_tex(k + shift, var)) for k, c in enumerate(coeffs))


def sym_poly_zinv(coeffs: Sequence[Number], var: str = "z", shift: int = 0) -> str:
    """sum_k coeffs[k] z^{-(k+shift)} as a parenthesised SymPy string."""
    terms = []
    for k, c in enumerate(coeffs):
        c = Q(c)
        if c == 0:
            continue
        e = k + shift
        terms.append(sym_num(c) if e == 0 else f"{sym_num(c)}*{zpow_sym(e, var)}")
    return "(" + (" + ".join(terms) if terms else "0") + ")"


def tex_frac(num_tex: str, den_tex: str) -> str:
    return r"\dfrac{" + num_tex + "}{" + den_tex + "}"


def tex_factor(p: Number, var: str = "z") -> str:
    """(1 - p z^{-1}) as LaTeX, e.g. \\left(1 - \\tfrac{1}{2}z^{-1}\\right), \\left(1 + 2z^{-1}\\right)."""
    return r"\left(" + tex_sum([(1, ""), (-Q(p), f"{var}^{{-1}}")]) + r"\right)"


def sym_factor(p: Number, var: str = "z") -> str:
    return f"(1 - {sym_num(p)}*{var}**(-1))"


# ----------------------------------------------------------------------------- sequences
def tex_seq(values: Sequence[Number], start: int, pad_to_zero: bool = True) -> str:
    """\\{1,\\ \\underset{\\uparrow}{2},\\ 3\\} with the arrow under n = 0.

    values[i] is x[start + i]. If n = 0 lies outside the listed samples and pad_to_zero is True,
    zeros are added so that the arrow can be shown (as the exams do)."""
    vals = [Q(v) for v in values]
    if pad_to_zero:
        if start > 0:
            vals = [Fraction(0)] * start + vals
            start = 0
        end = start + len(vals) - 1
        if end < 0:
            vals = vals + [Fraction(0)] * (-end)
    items = []
    for i, v in enumerate(vals):
        t = tex_num(v)
        if start + i == 0:
            t = r"\underset{\uparrow}{" + t + "}"
        items.append(t)
    return r"\{" + r",\ ".join(items) + r"\}"


def plain_seq(values: Sequence[Number]) -> str:
    """'1, -2, 3/4' — what a student types into a list answer."""
    return ", ".join(plain_num(v) for v in values)
