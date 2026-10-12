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


# --------------------------------------------------------------------------- pl-sketch submissions
def sketch_raw(params: dict, name: str, drawings: dict[str, list[list[tuple[float, float]]]]) -> dict[str, str]:
    """Raw submission for a pl-sketch element, to use as an ``overrides`` entry in submissions().

    ``drawings`` maps a tool id to a list of drawings in GRAPH coordinates: for polyline / spline /
    free-draw / line tools each drawing is a list of (x, y) points, for point tools a list with one
    (x, y), for horizontal-line [(None, y)], for vertical-line [(x, None)]. Tools not mentioned
    keep their initial state (read-only helpers, given curves). The canvas transform and the
    gradeable format are taken from the element itself (``params[name]``, written by pl-sketch's
    prepare, and pl-sketch's own ``_solution_to_gradeable``). Polyline and spline drawings are then
    re-fitted the way the browser does it (sketchresponse's ``fitCurve(points, 0)``: one cubic
    Bezier per segment with centre tangents at interior vertices, rounded to 0.01 px), so a
    polyline is graded as a SMOOTH curve through its vertices, not as straight segments — exactly
    what a student's browser sends (PrairieLearn's own Test button sends straight segments).
    Example::

        sketch_raw(params, "xd", {"pl": [[(-3.14, 0), (-1, 0), (0, 2), (1, 0), (3.14, 0)]]})
        # -> {"xd-sketchresponse-submission": "<base64>"}
    """
    import base64
    import copy as _copy
    import importlib.util
    import json as _json

    el = plsim.PL_APP / "elements" / "pl-sketch"
    if str(el) not in sys.path:
        sys.path.insert(0, str(el))
    spec = importlib.util.spec_from_file_location("_pl_sketch_element", el / "pl-sketch.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]

    p = params[name]
    r = p["ranges"]

    def sx(x: float) -> float:
        return (x - r["x_start"]) * r["width"] / (r["x_end"] - r["x_start"])

    def sy(y: float) -> float:
        return (y - r["y_end"]) * r["height"] / (r["y_start"] - r["y_end"])

    state = _copy.deepcopy(p["initial_state"])
    for tool_id, items in drawings.items():
        kind = p["tool_data"][tool_id]["name"]
        if kind == "point":
            state[tool_id] = [{"x": sx(it[0][0]), "y": sy(it[0][1])} for it in items]
        elif kind == "horizontal-line":
            state[tool_id] = [{"y": sy(it[0][1])} for it in items]
        elif kind == "vertical-line":
            state[tool_id] = [{"x": sx(it[0][0])} for it in items]
        else:
            state[tool_id] = [[{"x": sx(x), "y": sy(y)} for (x, y) in it] for it in items]
    gradeable = mod._solution_to_gradeable(state, p["tool_data"], p["config"]["width"], p["config"]["height"])
    for tool_id in drawings:
        if p["tool_data"][tool_id]["name"] in ("polyline", "spline"):
            gradeable[tool_id] = [{"spline": _browser_fit([(pt["x"], pt["y"]) for pt in curve])}
                                  for curve in state[tool_id] if len(curve) > 1]
    payload = base64.b64encode(_json.dumps({"data": state, "gradeable": gradeable}).encode("utf-8")).decode("utf-8")
    return {f"{name}-sketchresponse-submission": payload}


def _browser_fit(pts: list[tuple[float, float]]) -> list[list[float]]:
    """sketchresponse's fitCurve(points, 0) as the polyline/spline plugins call it (sketch_tool/lib/
    plugins/polyline.js splineData): each segment p_i -> p_{i+1} becomes a cubic Bezier with control
    points p_i + t_i*L/3 and p_{i+1} + s_{i+1}*L/3 (L = |p_{i+1} - p_i|), where the tangents are the
    chord directions at the two ends and unit(p_{i+1} - p_{i-1}) at interior vertices; coordinates
    rounded to 0.01 px. Checked against the JavaScript on 200 random polylines (0 mismatches)."""
    def unit(a, b):
        dx, dy = b[0] - a[0], b[1] - a[1]
        n = math.hypot(dx, dy) or 1.0
        return (dx / n, dy / n)

    n = len(pts)
    out = [list(pts[0])]
    for i in range(n - 1):
        a, b = pts[i], pts[i + 1]
        L = math.hypot(b[0] - a[0], b[1] - a[1]) / 3
        t1 = unit(a, b) if i == 0 else unit(pts[i - 1], pts[i + 1])
        t2 = unit(b, a) if i + 1 == n - 1 else unit(pts[i + 2], pts[i])
        out += [[a[0] + t1[0] * L, a[1] + t1[1] * L], [b[0] + t2[0] * L, b[1] + t2[1] * L], list(b)]
    return [[round(100 * x) / 100, round(100 * y) / 100] for x, y in out]
