"""Independent check of questions/dtft/parseval-inner-product: X_d and Y_d are evaluated with numpy
from the definition on a uniform grid over one period and the integrals are computed with the
periodic rectangle rule (exact for trigonometric polynomials of degree below the grid size)."""

import numpy as np

NGRID = 2048


def _dtft(v, s, w):
    n = np.arange(s, s + len(v))
    return np.exp(-1j * np.outer(w, n)) @ np.asarray(v, dtype=float)


def _num(v):
    return str(int(v)) if float(v).is_integer() else repr(float(v))


def _vals(params):
    w = -np.pi + 2 * np.pi * np.arange(NGRID) / NGRID
    X = _dtft(params["x"], params["sx"], w)
    Y = _dtft(params["y"], params["sy"], w)
    n0 = params["n0"]
    dw = 2 * np.pi / NGRID
    return {
        "A": (X * np.conj(Y)).sum() * dw / (2 * np.pi),
        "B": (X * np.exp(1j * w * n0)).sum() * dw / np.pi,
        "C": (X * Y).sum() * dw / (2 * np.pi),
        "D": (np.abs(X) ** 2).sum() * dw / np.pi,
    }


def check(params, correct):
    probs = []
    for name, val in _vals(params).items():
        if abs(val.imag) > 1e-9:
            probs.append(f"{name}: numerical value {val} is not real")
        if abs(val.real - correct[name]) > 1e-8 * max(1.0, abs(val.real)):
            probs.append(f"{name}: answer {correct[name]} != numerical {val.real}")
    if abs(correct["A"] - correct["C"]) < 0.5:
        probs.append("(a) and (c) coincide: the missing conjugate is not tested")
    if abs(correct["A"]) < 0.5 or abs(correct["C"]) < 0.5:
        probs.append("(a) or (c) is zero")
    x, y, sx, sy, n0 = params["x"], params["y"], params["sx"], params["sy"], params["n0"]
    for v in (x, y):
        if not (3 <= len(v) <= 4) or v[0] == 0 or v[-1] == 0 or max(abs(t) for t in v) > 3:
            probs.append(f"sequence {v} violates the design")
    dx = {sx + i: v for i, v in enumerate(x)}
    if n0 == 0 or dx.get(n0, 0) == 0 or dx.get(-n0, 0) == dx[n0] or dx.get(0, 0) == dx[n0]:
        probs.append(f"n0 = {n0} does not separate x[n0] from x[-n0] / x[0]")
    return probs


def submissions(params, correct):
    x, y, sx, sy, n0 = params["x"], params["y"], params["sx"], params["sy"], params["n0"]
    A, B, C, D = (int(round(correct[k])) for k in "ABCD")
    dx = {sx + i: v for i, v in enumerate(x)}
    dy = {sy + i: v for i, v in enumerate(y)}
    cases = [
        ({"A": f"{2 * A}/2", "C": f"{C}.0"}, {"A": 1, "C": 1}),       # equivalent forms
        ({"B": f"{3 * B}/3", "D": f"{D}.000"}, {"B": 1, "D": 1}),
        ({"C": _num(A)}, {"C": 0}),                                  # (c) with the conjugate anyway
        ({"A": _num(C)}, {"A": 0}),                                  # (a) with y reversed
        ({"B": _num(B // 2)}, {"B": 0}),                             # forgot 2 pi: x[n0]
        ({"D": _num(D // 2)}, {"D": 0}),                             # forgot 2 pi: energy
        ({"A": _num(2 * A)}, {"A": 0}),                              # 2 pi too many
        ({"B": "pi"}, {"B": "invalid"}),
    ]
    if 2 * dx.get(-n0, 0) != B:
        cases.append(({"B": _num(2 * dx.get(-n0, 0))}, {"B": 0}))   # e^{j w n0} read as x[-n0]
    if 2 * dx.get(0, 0) != B:
        cases.append(({"B": _num(2 * dx.get(0, 0))}, {"B": 0}))     # the n = 0 shortcut
    pos = sum(a * b for a, b in zip(x, y))                           # aligned by list position
    if pos != A:
        cases.append(({"A": _num(pos)}, {"A": 0}))
    if 2 * sum(x) ** 2 != D:
        cases.append(({"D": _num(2 * sum(x) ** 2)}, {"D": 0}))       # (sum x)^2
    if -C != C and -C != A:
        cases.append(({"C": _num(-C)}, {"C": 0}))
    return cases
