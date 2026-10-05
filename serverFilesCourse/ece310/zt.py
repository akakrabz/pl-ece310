"""z-transform machinery: regions of convergence, partial fractions, inverse transforms, and the
standard signal terms (with their transforms) used by the generators.

Everything is exact (``Fraction``) for real rational poles. Complex poles work numerically
(``complex``) in ``pfe`` and ``InverseZ.value``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Sequence

from .fmt import (
    Q,
    sym_factor,
    sym_num,
    sym_pow,
    sym_sum,
    tex_factor,
    tex_frac,
    tex_num,
    tex_pow,
    tex_sum,
    zpow_sym,
    zpow_tex,
)
from .poly import from_roots, pdivmod, trim


def mag(p) -> Fraction | float:
    """|p| — exact for real rationals, float otherwise."""
    if isinstance(p, (int, Fraction)):
        return abs(Q(p))
    return abs(complex(p))


# ============================================================================ ROC
@dataclass(frozen=True)
class ROC:
    """inner < |z| < outer. outer=None means infinity. has0 (only if inner == 0) says whether z = 0
    itself belongs to the ROC; hasinf (only if outer is None) whether z = infinity does."""

    inner: Fraction = Fraction(0)
    outer: Fraction | None = None
    has0: bool = False
    hasinf: bool = True
    empty: bool = False

    # ---------------------------------------------------------------- constructors
    @staticmethod
    def right(r, starts_at: int = 0) -> "ROC":
        """Right-sided term with largest pole magnitude r whose first sample is at n = starts_at."""
        r = Q(r) if isinstance(r, (int, Fraction)) else r
        return ROC(inner=r, outer=None, has0=False, hasinf=starts_at >= 0)

    @staticmethod
    def left(r, ends_at: int = -1) -> "ROC":
        """Left-sided term with smallest pole magnitude r whose last sample is at n = ends_at."""
        r = Q(r) if isinstance(r, (int, Fraction)) else r
        return ROC(inner=Fraction(0), outer=r, has0=ends_at <= 0, hasinf=False)

    @staticmethod
    def finite(n_first: int, n_last: int) -> "ROC":
        return ROC(inner=Fraction(0), outer=None, has0=n_last <= 0, hasinf=n_first >= 0)

    @staticmethod
    def everything() -> "ROC":
        return ROC(inner=Fraction(0), outer=None, has0=True, hasinf=True)

    # ---------------------------------------------------------------- algebra
    def intersect(self, other: "ROC") -> "ROC":
        if self.empty or other.empty:
            return ROC(empty=True)
        inner = max(self.inner, other.inner)
        if self.outer is None:
            outer = other.outer
        elif other.outer is None:
            outer = self.outer
        else:
            outer = min(self.outer, other.outer)
        has0 = inner == 0 and (self.has0 if self.inner == 0 else False) and (other.has0 if other.inner == 0 else False)
        hasinf = outer is None and (self.hasinf if self.outer is None else False) and (other.hasinf if other.outer is None else False)
        if outer is not None and inner >= outer:
            return ROC(empty=True)
        return ROC(inner=inner, outer=outer, has0=has0, hasinf=hasinf)

    # ---------------------------------------------------------------- predicates
    def contains_radius(self, r) -> bool:
        if self.empty:
            return False
        if r == 0:
            return self.inner == 0 and self.has0
        return self.inner < r and (self.outer is None or r < self.outer)

    def contains_unit_circle(self) -> bool:
        return self.contains_radius(1)

    def is_causal_system(self) -> bool:
        """Causal LTI system <=> ROC is the exterior of a circle INCLUDING infinity."""
        return not self.empty and self.outer is None and self.hasinf

    def kind(self) -> str:
        """'right', 'left', 'two-sided', 'finite' or 'empty' (shape of the ROC)."""
        if self.empty:
            return "empty"
        if self.inner == 0 and self.outer is None:
            return "finite"
        if self.outer is None:
            return "right"
        if self.inner == 0:
            return "left"
        return "two-sided"

    def test_point(self) -> float:
        """A radius strictly inside the ROC (for numerical checks)."""
        if self.empty:
            raise ValueError("empty ROC")
        lo = float(self.inner)
        if self.outer is None:
            return lo * 1.5 + 0.7 if lo > 0 else 1.3
        hi = float(self.outer)
        if lo == 0:
            return hi * 0.6
        return math.sqrt(lo * hi)

    # ---------------------------------------------------------------- output
    def tex(self) -> str:
        """LaTeX (uses \\lt / \\gt, never raw < >)."""
        z = r"\lvert z\rvert"
        if self.empty:
            return r"\text{empty (no z-transform)}"
        if self.inner == 0 and self.outer is None:
            if self.has0 and self.hasinf:
                return r"\text{all } z"
            if self.hasinf:
                return z + r" \gt 0"
            if self.has0:
                return z + r" \lt \infty"
            return r"0 \lt " + z + r" \lt \infty"
        if self.outer is None:
            r = _tex_radius(self.inner)
            return (z + r" \gt " + r) if self.hasinf else (r + r" \lt " + z + r" \lt \infty")
        R = _tex_radius(self.outer)
        if self.inner == 0:
            return (z + r" \lt " + R) if self.has0 else (r"0 \lt " + z + r" \lt " + R)
        return _tex_radius(self.inner) + r" \lt " + z + r" \lt " + R

    def key(self) -> str:
        """Canonical string, handy for de-duplicating multiple-choice options."""
        return self.tex()

    def to_json(self) -> dict:
        return {
            "inner": str(self.inner),
            "outer": None if self.outer is None else str(self.outer),
            "has0": self.has0,
            "hasinf": self.hasinf,
            "empty": self.empty,
        }

    @staticmethod
    def from_json(d: dict) -> "ROC":
        return ROC(
            inner=Fraction(d["inner"]),
            outer=None if d["outer"] is None else Fraction(d["outer"]),
            has0=d["has0"],
            hasinf=d["hasinf"],
            empty=d["empty"],
        )


def _tex_radius(r) -> str:
    if isinstance(r, (int, Fraction)):
        return tex_num(r)
    return f"{float(r):.4g}"


def rocs_for_poles(poles: Sequence, finite_right: bool = True) -> list[ROC]:
    """All ROCs of a rational X(z) with these poles (proper in z^{-1}, no poles at 0/infinity
    besides the trivial ones): the annuli between consecutive distinct pole magnitudes."""
    mags = sorted({mag(p) for p in poles})
    out = [ROC(inner=Fraction(0), outer=mags[0], has0=True, hasinf=False)]
    for a, b in zip(mags, mags[1:]):
        out.append(ROC(inner=a, outer=b, has0=False, hasinf=False))
    out.append(ROC(inner=mags[-1], outer=None, has0=False, hasinf=finite_right))
    return out


# ============================================================================ partial fractions
def pfe(num: Sequence, poles: Sequence) -> tuple[list, list]:
    """X(z) = num(z^{-1}) / prod_k (1 - p_k z^{-1}) with DISTINCT nonzero poles.

    Returns (C, A): X = sum_k C[k] z^{-k} + sum_k A[k] / (1 - p_k z^{-1}).
    C is [0] when X is proper (deg num < number of poles). Exact for Fractions."""
    if len(set(poles)) != len(poles):
        raise ValueError("repeated poles: use a different expansion")
    if any(p == 0 for p in poles):
        raise ValueError("pole at z = 0 is a shift, not a PFE pole")
    den = from_roots(poles)
    num = trim(num)
    if len(num) >= len(den):
        C, R = pdivmod(num, den)
    else:
        C, R = [Fraction(0)], list(num)
    A = []
    for k, p in enumerate(poles):
        w = 1 / p
        numv = sum(c * w**i for i, c in enumerate(R))
        denv = 1
        for j, pj in enumerate(poles):
            if j != k:
                denv *= 1 - pj * w
        A.append(numv / denv)
    return C, A


# ============================================================================ inverse transform
@dataclass
class InverseZ:
    """x[n] = sum_k C[k] delta[n-k] + sum_{right} A p^n u[n] - sum_{left} A p^n u[-n-1]."""

    C: list
    poles: list
    A: list
    sides: list  # 'R' or 'L' per pole

    def value(self, n: int):
        v = self.C[n] if 0 <= n < len(self.C) else 0
        for p, a, s in zip(self.poles, self.A, self.sides):
            if s == "R" and n >= 0:
                v += a * p**n
            elif s == "L" and n <= -1:
                v -= a * p**n
        return v

    def right_terms(self):
        return [(a, p) for p, a, s in zip(self.poles, self.A, self.sides) if s == "R"]

    def left_terms(self):
        return [(a, p) for p, a, s in zip(self.poles, self.A, self.sides) if s == "L"]

    def has_direct_terms(self) -> bool:
        return any(c != 0 for c in self.C)

    def sym_nonneg(self) -> str:
        """SymPy string of x[n] for n >= len(C) (all n >= 0 when there are no direct terms)."""
        return sym_sum((a, sym_pow(p, "n")) for a, p in self.right_terms())

    def sym_neg(self) -> str:
        """SymPy string of x[n] for n <= -1."""
        return sym_sum((-a, sym_pow(p, "n")) for a, p in self.left_terms())

    def tex(self) -> str:
        terms = []
        for k, c in enumerate(self.C):
            terms.append((c, r"\delta[n]" if k == 0 else rf"\delta[n-{k}]"))
        for a, p in self.right_terms():
            terms.append((a, exp_body_tex(p, "n")))
        for a, p in self.left_terms():
            terms.append((-a, exp_body_tex(p, "-n-1")))
        return tex_sum(terms)


def inverse_pfe(C: Sequence, poles: Sequence, A: Sequence, roc: ROC) -> InverseZ:
    """Assign each pole to the right- or left-sided pair according to the ROC."""
    sides = []
    for p in poles:
        m = mag(p)
        if m <= roc.inner:
            sides.append("R")
        elif roc.outer is not None and m >= roc.outer:
            sides.append("L")
        else:
            raise ValueError(f"pole {p} lies inside the ROC {roc.tex()}")
    return InverseZ(list(C), list(poles), list(A), sides)


# ============================================================================ signal terms
def u_arg_right(k: int) -> str:
    """n-k, n, n+2."""
    if k == 0:
        return "n"
    return f"n-{k}" if k > 0 else f"n+{-k}"


def u_arg_left(m: int) -> str:
    """argument of u[-n+m]: -n, -n+2, -n-1."""
    if m == 0:
        return "-n"
    return f"-n+{m}" if m > 0 else f"-n-{-m}"


def exp_body_tex(a, u_arg: str) -> str:
    """a^n u[arg] (just u[arg] when a = 1)."""
    if Q(a) == 1:
        return f"u[{u_arg}]"
    return tex_pow(a, "n") + r"\,u[" + u_arg + "]"


@dataclass
class ZTerm:
    """One additive term of a signal together with its z-transform."""

    x_tex: str            # time-domain LaTeX of the term (with its sign/coefficient)
    X_tex: str            # its z-transform, LaTeX
    X_sym: str            # its z-transform, SymPy string in z
    roc: ROC
    poles: list = field(default_factory=list)
    kind: str = ""
    spec: dict = field(default_factory=dict)  # JSON-able description for checkers


def _signed_frac(c, num_body: str, den_tex: str) -> str:
    """c * num_body / den with the sign in front: '-\\dfrac{2z^{-1}}{...}', never '\\dfrac{-2z^{-1}}{...}'."""
    c = Q(c)
    frac = tex_frac(tex_sum([(abs(c), num_body)]), den_tex)
    return "-" + frac if c < 0 else frac


def term_right_exp(A, a, k: int = 0) -> ZTerm:
    """x[n] = A a^n u[n-k]  <->  A a^k z^{-k} / (1 - a z^{-1}),  |z| > |a| (and z != inf if k < 0)."""
    A, a = Q(A), Q(a)
    c = A * a**k
    x_tex = tex_sum([(A, exp_body_tex(a, u_arg_right(k)))])
    X_tex = _signed_frac(c, zpow_tex(k), tex_sum([(1, ""), (-a, "z^{-1}")]))
    X_sym = f"{sym_num(c)}*{zpow_sym(k)}/{sym_factor(a)}"
    return ZTerm(x_tex, X_tex, X_sym, ROC.right(abs(a), starts_at=k), [a], "right_exp",
                 {"type": "right_exp", "A": str(A), "a": str(a), "k": k})


def term_left_exp(B, b, m: int = -1) -> ZTerm:
    """x[n] = B b^n u[-n+m] (nonzero for n <= m)  <->  -B b^{m+1} z^{-(m+1)} / (1 - b z^{-1}),
    |z| < |b| (and z != 0 if m >= 1)."""
    B, b = Q(B), Q(b)
    c = -B * b ** (m + 1)
    x_tex = tex_sum([(B, exp_body_tex(b, u_arg_left(m)))])
    X_tex = _signed_frac(c, zpow_tex(m + 1), tex_sum([(1, ""), (-b, "z^{-1}")]))
    X_sym = f"{sym_num(c)}*{zpow_sym(m + 1)}/{sym_factor(b)}"
    return ZTerm(x_tex, X_tex, X_sym, ROC.left(abs(b), ends_at=m), [b], "left_exp",
                 {"type": "left_exp", "B": str(B), "b": str(b), "m": m})


def term_n_exp(A, a) -> ZTerm:
    """x[n] = A n a^n u[n]  <->  A a z^{-1} / (1 - a z^{-1})^2,  |z| > |a|."""
    A, a = Q(A), Q(a)
    rest = exp_body_tex(a, "n") if a != 1 else "u[n]"
    body = ("n \\cdot " if rest[0].isdigit() else "n\\,") + rest  # n \cdot 3^{n}, never "n\,3^{n}"
    x_tex = tex_sum([(A, body)])
    X_tex = _signed_frac(A * a, "z^{-1}", tex_factor(a) + "^{2}")
    X_sym = f"{sym_num(A * a)}*z**(-1)/{sym_factor(a)}**2"
    return ZTerm(x_tex, X_tex, X_sym, ROC.right(abs(a), 0), [a, a], "n_exp",
                 {"type": "n_exp", "A": str(A), "a": str(a)})


def term_finite(values: Sequence, start: int) -> ZTerm:
    """x[n] = sum values[i] delta[n - (start+i)]  <->  sum values[i] z^{-(start+i)}."""
    vals = [Q(v) for v in values]
    while vals and vals[0] == 0:
        vals.pop(0)
        start += 1
    while vals and vals[-1] == 0:
        vals.pop()
    if not vals:
        raise ValueError("all-zero finite sequence")
    terms_x, terms_X, sym = [], [], []
    for i, v in enumerate(vals):
        n = start + i
        if v == 0:
            continue
        terms_x.append((v, r"\delta[n]" if n == 0 else (rf"\delta[n-{n}]" if n > 0 else rf"\delta[n+{-n}]")))
        terms_X.append((v, zpow_tex(n)))
        sym.append(f"{sym_num(v)}*{zpow_sym(n)}")
    return ZTerm(tex_sum(terms_x), tex_sum(terms_X), "(" + " + ".join(sym) + ")",
                 ROC.finite(start, start + len(vals) - 1), [], "finite",
                 {"type": "finite", "values": [str(v) for v in vals], "start": start})


_COS = {"pi/2": Fraction(0), "pi/3": Fraction(1, 2), "2*pi/3": Fraction(-1, 2), "pi": Fraction(-1)}
_W_TEX = {"pi/2": r"\tfrac{\pi}{2}", "pi/3": r"\tfrac{\pi}{3}", "2*pi/3": r"\tfrac{2\pi}{3}", "pi": r"\pi"}


def term_cos(A, r, w: str) -> ZTerm:
    """x[n] = A r^n cos(w n) u[n], w in {'pi/2','pi/3','2*pi/3','pi'}
    <->  A (1 - r cos w z^{-1}) / (1 - 2 r cos w z^{-1} + r^2 z^{-2}),  |z| > r."""
    A, r = Q(A), Q(r)
    c = _COS[w]
    rn = "" if r == 1 else tex_pow(r, "n") + r"\,"
    x_tex = tex_sum([(A, rn + r"\cos\!\left(" + _W_TEX[w] + r"n\right)u[n]")])
    num = [A, -A * r * c]
    den = [Fraction(1), -2 * r * c, r * r]
    X_tex = tex_frac(tex_sum([(num[0], ""), (num[1], "z^{-1}")]),
                     tex_sum([(den[0], ""), (den[1], "z^{-1}"), (den[2], "z^{-2}")]))
    X_sym = (f"({sym_num(num[0])} + {sym_num(num[1])}*z**(-1))/"
             f"(1 + {sym_num(den[1])}*z**(-1) + {sym_num(den[2])}*z**(-2))")
    import cmath
    th = {"pi/2": math.pi / 2, "pi/3": math.pi / 3, "2*pi/3": 2 * math.pi / 3, "pi": math.pi}[w]
    poles = [float(r) * cmath.exp(1j * th), float(r) * cmath.exp(-1j * th)]
    return ZTerm(x_tex, X_tex, X_sym, ROC.right(r, 0), poles, "cos",
                 {"type": "cos", "A": str(A), "r": str(r), "w": w})


def sum_terms(terms: Sequence[ZTerm]) -> tuple[str, str, str, ROC]:
    """(x_tex, X_tex, X_sym, roc) of a sum of terms; X is left as a sum of the term transforms."""
    x_tex = terms[0].x_tex
    for t in terms[1:]:
        x_tex += (" " + t.x_tex) if t.x_tex.startswith("-") else (" + " + t.x_tex)
    X_tex = terms[0].X_tex
    for t in terms[1:]:
        X_tex += (" " + t.X_tex) if t.X_tex.startswith("-") else (" + " + t.X_tex)
    X_sym = " + ".join(f"({t.X_sym})" for t in terms)
    roc = terms[0].roc
    for t in terms[1:]:
        roc = roc.intersect(t.roc)
    return x_tex, X_tex, X_sym, roc


# ============================================================================ numerical helpers
def x_value(spec: dict, n: int):
    """Exact sample of a term described by ZTerm.spec (used by checkers and solutions)."""
    t = spec["type"]
    if t == "right_exp":
        return Q(spec["A"]) * Q(spec["a"]) ** n if n >= spec["k"] else Fraction(0)
    if t == "left_exp":
        return Q(spec["B"]) * Q(spec["b"]) ** n if n <= spec["m"] else Fraction(0)
    if t == "n_exp":
        return Q(spec["A"]) * n * Q(spec["a"]) ** n if n >= 0 else Fraction(0)
    if t == "finite":
        i = n - spec["start"]
        return Q(spec["values"][i]) if 0 <= i < len(spec["values"]) else Fraction(0)
    if t == "cos":
        th = {"pi/2": math.pi / 2, "pi/3": math.pi / 3, "2*pi/3": 2 * math.pi / 3, "pi": math.pi}[spec["w"]]
        return float(Q(spec["A"])) * float(Q(spec["r"])) ** n * round(math.cos(th * n), 12) if n >= 0 else 0.0
    raise ValueError(t)
