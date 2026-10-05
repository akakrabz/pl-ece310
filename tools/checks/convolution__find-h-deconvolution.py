"""Independent check of questions/convolution/find-h-deconvolution: convolve the answer back."""

import numpy as np


def _h(correct):
    return np.array(correct["h"]["_value"][0])


def check(params, correct):
    probs = []
    p = params
    h, nh = _h(correct), correct["nh"]
    if h[0] == 0 or h[-1] == 0:
        probs.append(f"h has a zero end sample: {h}")
    y, ny = np.array(p["y"], dtype=float), p["ny"]
    if p["template"] == "fir":
        x, nx = np.array(p["x"], dtype=float), p["nx"]
        yy = np.convolve(x, h)
        if nx + nh != ny:
            probs.append(f"start: nx + nh = {nx + nh} != ny = {ny}")
        if yy.shape != y.shape or not np.allclose(yy, y):
            probs.append(f"x * h = {yy} != y = {y}")
    else:
        a, M = float(p["a"]), 20
        x = a ** np.arange(M)                       # x[n] = a^n u[n], n = 0..M-1 (exact up to n = M-1)
        yy = np.convolve(x, h)[:M]                  # starts at n = 0 + nh
        want = np.zeros(M)
        off = ny - nh
        if off < 0:
            probs.append("y starts before h")
        else:
            want[off:off + len(y)] = y
        if not np.allclose(yy, want):
            probs.append(f"(a^n u[n]) * h = {yy[:8]} != y = {want[:8]}")
    return probs


def _row(v, sep=", "):
    return "[" + sep.join(str(int(round(t))) for t in v) + "]"


def submissions(params, correct):
    h, nh = [int(round(t)) for t in _h(correct)], correct["nh"]
    cases = [
        ({"nh": str(nh + 1)}, {"nh": 0, "h": 1}),                       # correct samples, start one late
        ({"nh": str(nh - 1)}, {"nh": 0}),                               # ... or one early
        ({"h": _row(h, " ")}, {"h": 1}),                                # MATLAB style
        ({"h": "[" + _row(h) + "]"}, {"h": 1}),                         # Python style [[...]]
        ({"h": _row(h + [0])}, {"h": 0}),                               # trailing zero appended
        ({"nh": "0.5"}, {"nh": "invalid"}),
    ]
    if h != h[::-1]:
        cases.append(({"h": _row(h[::-1])}, {"h": 0}))                 # reversed
    p = params
    if p["template"] == "fir" and p["nx"] != 0:
        cases.append(({"nh": str(p["ny"] + p["nx"])}, {"nh": 0}))       # added the start indices
    if p["template"] == "geo":
        cases.append(({"h": _row(p["y"])}, {"h": 0}))                   # took y itself as h
        alt = [p["y"][0] - 0] + [p["y"][j] - p["y"][j - 1] for j in range(1, len(p["y"]))] + [-p["y"][-1]]
        if alt != h:
            cases.append(({"h": _row(alt)}, {"h": 0}))                  # used the step rule y[n]-y[n-1] for a != 1
    return cases
