"""Correct-answer helpers for PrairieLearn elements.

`symbolic` is the recommended way to store the correct answer of a `pl-symbolic-input` whose answer
is a function of the time index n: the SymPy expression is serialized with PrairieLearn's own JSON
format and with n declared an *integer*. The element carries that assumption over to the student's
expression, so SymPy can prove equalities that only hold on the integers, e.g.

    (-3)**(-n) == (-1/3)**n        (-1)**n * 2**n == (-2)**n        2**(-n) == (1/2)**n

With a plain string correct answer n is a complex symbol and such correct submissions score 0.
"""

from __future__ import annotations

from typing import Iterable


def symbolic(expr: str, integer: Iterable[str] = ("n",), extra_locals: dict | None = None):
    """PrairieLearn JSON for a pl-symbolic-input correct answer; `expr` is a SymPy string (e.g. from
    ece310.fmt). Symbols listed in `integer` are declared integer=True; other names parse as usual
    (pi, sqrt, cos, sin, exp, ... are SymPy functions)."""
    import prairielearn as pl
    import sympy

    loc = {name: sympy.Symbol(name, integer=True) for name in integer}
    loc.update(extra_locals or {})
    return pl.to_json(sympy.sympify(expr, locals=loc))
