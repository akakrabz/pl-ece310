"""Helpers for the independent answer checks in tools/checks/ (NOT used by PrairieLearn).

Import from a checker as ``from checklib import *``; tools/ is on sys.path when the test runner
loads checkers, and serverFilesCourse/ too (so ``from ece310 import ...`` works)."""

from __future__ import annotations

import cmath
import math
import pathlib
import sys
from fractions import Fraction
from typing import Callable

ROOT = pathlib.Path(__file__).resolve().parents[1]
for p in (str(ROOT / "serverFilesCourse"), str(ROOT / "tools")):
    if p not in sys.path:
        sys.path.insert(0, p)

import plsim  # noqa: E402,F401  (makes sympy importable from the harness cache)
import numpy as np  # noqa: E402
import sympy  # noqa: E402
from scipy import signal  # noqa: E402,F401

Z, N = sympy.symbols("z n")


def S(expr: str, **extra):
    """Parse a correct-answer string the way a SymPy user would (variables z, n, plus extras)."""
    loc = {"z": Z, "n": N, "j": sympy.I, "pi": sympy.pi}
    loc.update({k: sympy.Symbol(k) for k in extra})
    return sympy.sympify(expr.replace("^", "**"), locals=loc)


def eval_z(expr: str, zval: complex) -> complex:
    return complex(S(expr).subs(Z, zval).evalf(30))


def eval_n(expr: str, nval: int) -> complex:
    return complex(S(expr).subs(N, nval).evalf(30))


def _log_abs(x) -> float:
    if isinstance(x, Fraction):
        return math.log(abs(x.numerator)) - math.log(x.denominator)
    if isinstance(x, int):
        return math.log(abs(x))
    return math.log(abs(complex(x)))


def zsum(x: Callable[[int], object], zval: complex, nmin: int = -800, nmax: int = 800) -> complex:
    """sum_{n=nmin}^{nmax} x[n] z^{-n}, computed in log space so that huge x[n] times tiny z^{-n}
    (or the reverse) does not overflow. Terms below e^-60 relative to 1 are dropped."""
    lz, az = math.log(abs(zval)), cmath.phase(zval)
    total = 0j
    for m in range(nmin, nmax + 1):
        v = x(m)
        if v == 0:
            continue
        la = _log_abs(v) - m * lz
        if la < -60:
            continue
        if la > 700:
            raise OverflowError(f"term n={m} too large: the point is probably outside the ROC")
        ph = (0.0 if (isinstance(v, (int, Fraction)) and v > 0) else (math.pi if isinstance(v, (int, Fraction)) else cmath.phase(complex(v)))) - m * az
        total += math.exp(la) * cmath.exp(1j * ph)
    return total


def impulse_response(b, a, N: int) -> np.ndarray:
    """First N samples of h for the causal LCCDE y[n] + sum a_k y[n-k] = sum b_k x[n-k] (a[0] = 1)."""
    x = np.zeros(N)
    x[0] = 1.0
    return signal.lfilter([float(v) for v in b], [float(v) for v in a], x)


def close(a, b, tol: float = 1e-9) -> bool:
    return abs(complex(a) - complex(b)) <= tol * max(1.0, abs(complex(b)))
