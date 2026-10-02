"""Independent check of questions/signals/index-transformations: brute-force y[n] = x[a n + b]."""

import numpy as np


def _brute(x, nx, a, b, lo=-30, hi=30):
    xs = dict(zip(range(nx, nx + len(x)), x))
    ns = np.arange(lo, hi + 1)
    y = np.array([xs.get(int(a * n + b), 0) for n in ns], dtype=float)
    nz = np.flatnonzero(y)
    return y[nz[0]:nz[-1] + 1], int(ns[nz[0]])


def _row(vals, sep=", "):
    return "[" + sep.join(str(int(v)) for v in vals) + "]"


def check(params, correct):
    probs = []
    x, nx, a, b = params["x"], params["nx"], params["a"], params["b"]
    y, ny = _brute(x, nx, a, b)
    got = np.array(correct["y"]["_value"][0])
    if correct["ny"] != ny:
        probs.append(f"first index {correct['ny']} != brute force {ny}")
    if got.shape != y.shape or not np.allclose(got, y):
        probs.append(f"y {got} != brute force {y}")
    if x[0] == 0 or x[-1] == 0 or not 4 <= len(x) <= 6:
        probs.append("x must have 4-6 samples with nonzero ends")
    if a not in (1, -1, 2, -2) or b == 0 or abs(b) > 3:
        probs.append(f"unexpected a={a}, b={b}")
    if abs(a) == 2:
        support = range(nx, nx + len(x))
        used = {a * n + b for n in range(-30, 31)}
        dropped = [m for m in support if m not in used and x[m - nx] != 0]
        survivors = np.count_nonzero(y)
        if not dropped:
            probs.append("decimation drops no nonzero sample")
        if survivors < 2:
            probs.append("fewer than two nonzero samples survive decimation")
    return probs


def submissions(params, correct):
    x, nx, a, b = params["x"], params["nx"], params["a"], params["b"]
    y, ny = _brute(x, nx, a, b)
    yl = [int(v) for v in y]
    cases = [
        ({"ny": str(ny + 1)}, {"ny": 0, "y": 1}),                          # right values, start off by one
        ({"ny": str(ny - 1)}, {"ny": 0, "y": 1}),
        ({"y": _row(yl, " ")}, {"y": 1}),                                   # MATLAB-style [1 2 3]
        ({"y": "[" + ", ".join(f"{v}.0" for v in yl) + "]"}, {"y": 1}),    # decimals
        ({"y": _row(yl + [0])}, {"y": 0}),                                 # trailing zero
        ({"ny": "0.5"}, {"ny": "invalid"}),
    ]
    if yl != yl[::-1]:
        cases.append(({"y": _row(yl[::-1])}, {"y": 0}))                    # reversed vector
    if a != 1:                                                              # x[a(n+b)]: wrong order
        yw, nyw = _brute(x, nx, a, a * b)
        exp = {}
        if nyw != ny:
            exp["ny"] = 0
        if len(yw) != len(y) or not np.allclose(yw, y):
            exp["y"] = 0
        if exp:
            cases.append(({"ny": str(nyw), "y": _row(yw)}, exp))
    yb, nyb = _brute(x, nx, a, -b)                                          # sign of the shift flipped
    exp = {}
    if nyb != ny:
        exp["ny"] = 0
    if len(yb) != len(y) or not np.allclose(yb, y):
        exp["y"] = 0
    if exp:
        cases.append(({"ny": str(nyb), "y": _row(yb)}, exp))
    return cases
