"""Independent check of questions/convolution/finite-convolution."""

import numpy as np


def check(params, correct):
    probs = []
    x, h = params["x"], params["h"]
    y = np.convolve(x, h)
    ny = params["nx"] + params["nh"]
    got = np.array(correct["y"]["_value"][0])
    if correct["ny"] != ny:
        probs.append(f"start index {correct['ny']} != {ny}")
    if got.shape != y.shape or not np.allclose(got, y):
        probs.append(f"y {got} != np.convolve {y}")
    if y[0] == 0 or y[-1] == 0:
        probs.append("y has a zero end sample (ambiguous length)")
    return probs


def submissions(params, correct):
    y = [int(v) for v in correct["y"]["_value"][0]]
    ny = correct["ny"]
    off = list(y)
    off[len(off) // 2] += 1
    flipped_h = np.convolve(params["x"], params["h"][::-1]).astype(int).tolist()   # forgot to flip (correlation)
    cases = [
        ({"ny": str(ny + 1)}, {"ny": 0, "y": 1}),                               # start index off by one
        ({"y": "[" + " ".join(map(str, y)) + "]"}, {"y": 1}),                     # MATLAB-style row vector
        ({"y": "[" + ", ".join(map(str, off)) + "]"}, {"y": 0}),                  # one wrong sample
        ({"y": "[" + ", ".join(map(str, y + [0])) + "]"}, {"y": 0}),              # extra trailing zero
        ({"ny": "1.5"}, {"ny": "invalid"}),
    ]
    if flipped_h != y:
        cases.append(({"y": "[" + ", ".join(map(str, flipped_h)) + "]"}, {"y": 0}))
    return cases
