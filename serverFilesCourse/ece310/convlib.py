"""Small helpers shared by the convolution questions (questions/convolution/*, Writer B).

Import as ``from ece310 import convlib`` (it is deliberately NOT imported by ``ece310/__init__``,
so the existing modules are unchanged). Everything is exact (``fractions.Fraction``) and every
LaTeX string follows the ``fmt`` conventions (no raw ``<``, ``>``, ``&``).
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from .fmt import Q, sym_num, sym_pow, tex_num, tex_pow


# ----------------------------------------------------------------------------- indices
def tex_arg(k: int, var: str = "n") -> str:
    """Argument of a shifted signal: 'n', 'n-3', 'n+1' (k is the delay)."""
    if k == 0:
        return var
    return f"{var}-{k}" if k > 0 else f"{var}+{-k}"


def tex_u(k: int, var: str = "n") -> str:
    """u[n - k] in LaTeX."""
    return f"u[{tex_arg(k, var)}]"


def tex_delta(k: int, var: str = "n") -> str:
    r"""\delta[n - k] in LaTeX."""
    return rf"\delta[{tex_arg(k, var)}]"


def tex_coef(c, body: str) -> str:
    """c * body with +-1 elided; a digit-led body gets a \\cdot, a parenthesis-led one a thin space."""
    c = Q(c)
    if c == 1:
        return body
    if c == -1:
        return "-" + body
    mag = tex_num(c)
    if body[:1].isdigit():
        return mag + r"\cdot " + body
    if body.startswith((r"\left", "(", r"\tfrac", r"\dfrac")):
        return mag + r"\," + body
    return mag + body


def tex_signed_terms(terms: Sequence[tuple]) -> str:
    """sum of (coef, body) pairs with proper signs ('0' if all coefficients vanish)."""
    out = ""
    for c, body in terms:
        c = Q(c)
        if c == 0:
            continue
        t = tex_coef(abs(c), body) if body else tex_num(abs(c))
        if not out:
            out = ("-" if c < 0 else "") + t
        else:
            out += (" - " if c < 0 else " + ") + t
    return out or "0"


# ----------------------------------------------------------------------------- exponentials
def tex_exp(base, expo: str = "n") -> str:
    """base^{expo}, with parentheses for negative or fractional bases."""
    return tex_pow(base, expo)


def tex_exp_step(A, a, k: int) -> str:
    """A a^n u[n - k] in LaTeX (A = 1 elided)."""
    return tex_coef(A, tex_pow(a) + r"\," + tex_u(k))


def sym_recip_pow(base, expo: str = "n") -> str:
    """The same power written with the reciprocal base: (1/2)^n -> 2**(-(n)), (-3)^n -> (-1/3)**(-(n))."""
    r = 1 / Q(base)
    s = sym_num(r) if (r < 0 or r.denominator != 1) else str(r.numerator)
    return f"{s}**(-({expo}))"


def sym_shift_pow(base, shift: int) -> str:
    """base**(n - shift) as a SymPy string."""
    return sym_pow(base, f"n - ({shift})") if shift else sym_pow(base)


def plain_row(values: Sequence, sep: str = ", ") -> str:
    """'[1, -2, 3]' (or MATLAB style with sep=' ') for integer/terminating-decimal values."""
    def one(v):
        q = Q(v)
        if q.denominator == 1:
            return str(q.numerator)
        return repr(float(q))
    return "[" + sep.join(one(v) for v in values) + "]"


def shift_add_table_html(rows: Sequence[tuple[str, Sequence, int]], total_label: str,
                         total: Sequence, total_start: int) -> str:
    """Bootstrap table: one row per (label_tex, values, start) aligned on n, then the column sums.

    values[i] sits at n = start + i. Columns run over the support of ``total`` (and of the rows)."""
    lo = min([total_start] + [s for _, _, s in rows])
    hi = max([total_start + len(total) - 1] + [s + len(v) - 1 for _, v, s in rows])
    ns = list(range(lo, hi + 1))
    head = "".join(f"<th>$n={n}$</th>" for n in ns)
    body = []
    for label, vals, start in rows:
        cells = []
        for n in ns:
            j = n - start
            cells.append(f"<td>${tex_num(vals[j])}$</td>" if 0 <= j < len(vals) else "<td></td>")
        body.append(f"<tr><th>${label}$</th>{''.join(cells)}</tr>")
    tot = []
    for n in ns:
        j = n - total_start
        v = total[j] if 0 <= j < len(total) else Fraction(0)
        tot.append(f"<td><b>${tex_num(v)}$</b></td>")
    return ('<table class="table table-sm table-bordered text-center" style="width:auto">'
            f"<tr><th></th>{head}</tr>{''.join(body)}<tr><th>${total_label}$</th>{''.join(tot)}</tr></table>")


def sym_answer_json(sym: str):
    """PrairieLearn JSON for a pl-symbolic-input answer in n, with n declared an INTEGER.

    With an integer n, SymPy accepts every correct way of writing a power of a negative base,
    e.g. (-3)**(-n) for (-1/3)**n; with a plain symbol n those compare as different functions."""
    from .answers import symbolic

    return symbolic(sym)
