"""Independent check of questions/dtft/dtft-quantities: the DTFT is evaluated numerically on a
fine frequency grid and integrated over one period (the periodic rectangle rule is exact for
trigonometric polynomials of degree below the number of grid points)."""

import numpy as np

NGRID = 4096


def _dtft(x, start, w):
    n = np.arange(start, start + len(x))
    return np.exp(-1j * np.outer(w, n)) @ np.asarray(x, dtype=float)


def _num(v):
    """Plain number string for a submission (ints without a trailing .0)."""
    return str(int(v)) if float(v).is_integer() else repr(float(v))


def check(params, correct):
    probs = []
    x, start = params["x"], params["start"]
    w = -np.pi + 2 * np.pi * np.arange(NGRID) / NGRID
    X = _dtft(x, start, w)
    X0 = _dtft(x, start, np.array([0.0]))[0]
    Xpi = _dtft(x, start, np.array([np.pi]))[0]
    intX = (2 * np.pi / NGRID) * X.sum()
    intX2 = (2 * np.pi / NGRID) * (np.abs(X) ** 2).sum()
    for name, val in (("X0", X0), ("Xpi", Xpi), ("intX", intX / np.pi), ("intX2", intX2 / np.pi)):
        if abs(val.imag) > 1e-9:
            probs.append(f"{name}: numerical value {val} is not real")
        if abs(val.real - correct[name]) > 1e-8 * max(1.0, abs(val.real)):
            probs.append(f"{name}: correct answer {correct[name]} != numerical {val.real}")
    if abs(correct["X0"] - correct["Xpi"]) < 1e-12:
        probs.append("X_d(0) == X_d(pi): forgetting (-1)^n would not be detected")
    if abs(correct["intX"]) < 1e-12:
        probs.append("x[0] == 0: the integral of X_d is trivially 0")
    if not (3 <= len(x) <= 6) or x[0] == 0 or x[-1] == 0 or not (start <= 0 <= start + len(x) - 1):
        probs.append(f"sequence {x} starting at {start} violates the design")
    return probs


def submissions(params, correct):
    x = params["x"]
    X0, Xpi, iX, iX2 = (int(round(correct[k])) for k in ("X0", "Xpi", "intX", "intX2"))
    x0 = iX // 2
    cases = [
        # equivalent forms (score 1)
        ({"X0": f"{2 * X0}/2"}, {"X0": 1}),                       # fraction format
        ({"Xpi": f"{Xpi}.0"}, {"Xpi": 1}),                        # decimal format
        ({"intX": f"{4 * x0}/2"}, {"intX": 1}),
        ({"intX2": f"{3 * iX2}/3"}, {"intX2": 1}),
        # classic mistakes (score 0)
        ({"intX2": _num(iX2 // 2)}, {"intX2": 0}),               # Parseval without the 2*pi: sum |x|^2
        ({"intX": _num(x0)}, {"intX": 0}),                        # inverse DTFT without the 2*pi: x[0]
        ({"Xpi": _num(X0)}, {"Xpi": 0}),                          # sum x instead of sum (-1)^n x
        ({"X0": "abc"}, {"X0": "invalid"}),
    ]
    abs_sum = sum(abs(v) for v in x)
    if abs_sum != abs(X0) and abs_sum != X0:
        cases.append(({"X0": _num(abs_sum)}, {"X0": 0}))            # summed magnitudes
    wrong_parity = sum((-1) ** i * v for i, v in enumerate(x))      # (-1)^i with i = list position
    if wrong_parity != Xpi:
        cases.append(({"Xpi": _num(wrong_parity)}, {"Xpi": 0}))
    sq_of_sum = 2 * X0 ** 2                                          # 2*pi*(sum x)^2
    if sq_of_sum != iX2:
        cases.append(({"intX2": _num(sq_of_sum)}, {"intX2": 0}))
    if 2 * X0 != iX:
        cases.append(({"intX": _num(2 * X0)}, {"intX": 0}))          # 2*pi*X_d(0) instead of 2*pi*x[0]
    return cases
