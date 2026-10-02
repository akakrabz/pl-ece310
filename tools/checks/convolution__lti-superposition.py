"""Independent check of questions/convolution/lti-superposition.

Superposition with numpy index shifts on a common grid (not the server's list code), plus a
second route: recover h = deconvolve(y1, x1) and compute x2 * h with np.convolve."""

import numpy as np
from scipy import signal

LO, LEN = -12, 40                       # common grid n = LO .. LO + LEN - 1


def _grid(vals, start):
    g = np.zeros(LEN)
    g[start - LO: start - LO + len(vals)] = vals
    return g


def _shift(g, d):
    """g[n - d] on the same grid (d > 0 delays)."""
    return np.roll(g, d) if d >= 0 else np.roll(g, d)


def _trim(g):
    nz = np.nonzero(np.abs(g) > 1e-9)[0]
    return g[nz[0]: nz[-1] + 1], LO + nz[0]


def check(params, correct):
    probs = []
    p = params
    y1g = _grid(p["y1"], p["ny1"])
    y2g = p["alpha"] * _shift(y1g, p["k"]) + p["beta"] * _shift(y1g, p["m"])
    want, start = _trim(y2g)
    got = np.array(correct["y2"]["_value"][0])
    if correct["ny2"] != start:
        probs.append(f"start {correct['ny2']} != {start}")
    if got.shape != want.shape or not np.allclose(got, want):
        probs.append(f"y2 {got} != superposition {want}")
    # route 2: h from deconvolution, then x2 * h
    x1, y1 = np.array(p["x1"], float), np.array(p["y1"], float)
    hq, rem = signal.deconvolve(y1, x1)
    if np.max(np.abs(rem)) > 1e-9:
        probs.append("y1 is not x1 * (finite h)")
    nh = p["ny1"] - p["n1"]
    x2g = p["alpha"] * _shift(_grid(p["x1"], p["n1"]), p["k"]) + p["beta"] * _shift(_grid(p["x1"], p["n1"]), p["m"])
    x2, nx2 = _trim(x2g)
    if nx2 != p["nx2"] or x2.shape != (len(p["x2"]),) or not np.allclose(x2, p["x2"]):
        probs.append(f"x2 shown {p['x2']}@{p['nx2']} != {x2}@{nx2}")
    y2b = np.convolve(x2, hq)
    if nx2 + nh != start or y2b.shape != want.shape or not np.allclose(y2b, want):
        probs.append(f"x2 * h = {y2b} (start {nx2 + nh}) != {want} (start {start})")
    if p["template"] == "seq" and abs(p["k"] - p["m"]) < len(p["x1"]):
        probs.append("seq template with overlapping copies")
    return probs


def _row(v, sep=", "):
    return "[" + sep.join(str(int(round(t))) for t in v) + "]"


def submissions(params, correct):
    p = params
    y2 = [int(round(t)) for t in correct["y2"]["_value"][0]]
    ny2 = correct["ny2"]
    y1g = _grid(p["y1"], p["ny1"])
    cases = [
        ({"y2": _row(y2, " ")}, {"y2": 1}),                         # MATLAB style
        ({"y2": "[" + _row(y2) + "]"}, {"y2": 1}),                  # Python style
        ({"y2": _row([float(v) for v in y2])}, {"y2": 1}),
        ({"ny2": str(ny2 + 1)}, {"ny2": 0, "y2": 1}),               # start off by one
        ({"y2": _row(y2 + [0])}, {"y2": 0}),                        # trailing zero
        ({"ny2": "x"}, {"ny2": "invalid"}),
    ]
    mistakes = {
        "no shift": p["alpha"] * y1g + p["beta"] * y1g,
        "alpha dropped": (1 if p["alpha"] != 1 else p["alpha"]) * _shift(y1g, p["k"]) + (p["beta"] if p["alpha"] != 1 else 1) * _shift(y1g, p["m"]),
        "wrong direction": p["alpha"] * _shift(y1g, -p["k"]) + p["beta"] * _shift(y1g, -p["m"]),
    }
    for g in mistakes.values():
        if np.any(np.abs(g) > 1e-9):
            vals, _ = _trim(g)
            vals = [int(round(t)) for t in vals]
            if vals != y2:
                cases.append(({"y2": _row(vals)}, {"y2": 0}))
        else:
            cases.append(({"y2": "[0]"}, {"y2": 0}))
    conv_y1 = [int(round(t)) for t in np.convolve(p["x2"], p["y1"])]   # convolved x2 with y1 instead of h
    cases.append(({"y2": _row(conv_y1)}, {"y2": 0}))
    if p["ny1"] != ny2:
        cases.append(({"ny2": str(p["ny1"])}, {"ny2": 0}))           # kept the start of y1
    return cases
