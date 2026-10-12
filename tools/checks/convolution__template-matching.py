"""Independent check of questions/convolution/template-matching.

h = t[-n] by reversing the array (start -(a + 2)); y = np.convolve(x, h) with start index
0 + (-(a + 2)); the peak by argmax (unique, on the copy); rho from numpy norms of the template and
of the window x[n0 + a .. n0 + a + 2]; Cauchy-Schwarz bound checked."""

import math

import numpy as np


def _y(params):
    t, a, x = params["t"], params["a"], params["x"]
    h = np.array(t[::-1], float)
    hs = -(a + 2)
    y = np.convolve(np.array(x, float), h)
    return y, hs, h          # y[i] is y at n = i + hs (x starts at 0)


def _window(params, n):
    x, a = params["x"], params["a"]
    return np.array([x[k] if 0 <= k < len(x) else 0 for k in range(n + a, n + a + 3)], float)


def check(params, correct):
    probs = []
    y, ys, h = _y(params)
    if correct["hfirst"] != ys:
        probs.append(f"h starts at {ys}, answer {correct['hfirst']}")
    got = np.array(correct["hrow"]["_value"][0], float)
    if got.shape != h.shape or not np.allclose(got, h):
        probs.append(f"h row {got} != {h}")
    i = int(np.argmax(y))
    srt = np.sort(y)
    if srt[-1] - srt[-2] < 2:
        probs.append(f"peak not unique by a margin of 2: {y}")
    if correct["npeak"] != i + ys or correct["ypeak"] != int(y[i]):
        probs.append(f"peak ({correct['npeak']}, {correct['ypeak']}) != numpy ({i + ys}, {y[i]})")
    n1 = params["n1"]
    v1 = y[n1 - ys] if 0 <= n1 - ys < len(y) else 0.0
    if correct["yn1"] != int(v1) or n1 == i + ys:
        probs.append(f"y[{n1}] = {v1}, answer {correct['yn1']} (peak at {i + ys})")
    t = np.array(params["t"], float)
    w = _window(params, i + ys)
    if abs(float(w @ t) - y[i]) > 1e-9:
        probs.append("window inner product differs from the convolution at the peak")
    rho = y[i] / (np.linalg.norm(t) * np.linalg.norm(w))
    if not abs(rho - correct["rho"]) <= 1e-9:
        probs.append(f"rho {correct['rho']} != numpy {rho}")
    if not (0 < rho <= 1 + 1e-12):
        probs.append(f"rho = {rho} outside (0, 1]")
    if list(params["t"]) == list(params["t"])[::-1]:
        probs.append("symmetric template: the flip probe is meaningless")
    if len(_rho_mistakes(params, correct)) < 2:
        probs.append("fewer than two distinct wrong rho values")
    return probs


def _row(v, sep=", "):
    return "[" + sep.join(str(int(round(float(t)))) for t in v) + "]"


def _far(a, b):
    return abs(float(a) - float(b)) > 0.02 * max(1.0, abs(float(b)))


def _rho_mistakes(params, correct):
    y0 = correct["ypeak"]
    t = np.array(params["t"], float)
    et = float(t @ t)
    w = _window(params, correct["npeak"])
    ew = float(w @ w)
    out = []
    for v in (y0, y0 / et, y0 / (et * ew), y0 / et ** 0.5):     # unnormalized, by E_t only, no sqrt, one factor
        if _far(v, correct["rho"]) and all(_far(v, o) for o in out):
            out.append(v)
    return out


def submissions(params, correct):
    y, ys, h = _y(params)
    t, a = list(params["t"]), params["a"]
    n0, y0 = correct["npeak"], correct["ypeak"]
    cases = [
        ({"hrow": _row(h, " ")}, {"hrow": 1}),                                  # MATLAB style
        ({"hrow": _row(t)}, {"hrow": 0}),                                       # forgot the flip
        ({"hfirst": str(a)}, {"hfirst": 0}) if a != ys else ({"hfirst": str(ys + 1)}, {"hfirst": 0}),
        ({"npeak": str(n0 + 2)}, {"npeak": 0}),                                 # end of the copy / causal delay
        ({"npeak": str(n0 - 2)}, {"npeak": 0}),
        ({"npeak": str(n0 + 3)}, {"npeak": 0}),                                 # off by the template length
        ({"ypeak": f" {y0} "}, {"ypeak": 1}),
        ({"rho": f"{correct['rho']:.4g}"}, {"rho": 1}),                         # 4 significant digits
        ({"rho": f"{correct['rho']:.6f}"}, {"rho": 1}),
        ({"npeak": "2.5"}, {"npeak": "invalid"}),
    ]
    if -a != ys:
        cases.append(({"hfirst": str(-a)}, {"hfirst": 0}))                    # reversed t but forgot its length
    corr_t = np.convolve(np.array(params["x"], float), np.array(t, float))     # x * t: no flip
    if int(np.max(corr_t)) != y0:
        cases.append(({"ypeak": str(int(np.max(corr_t)))}, {"ypeak": 0}))
    if params["c"] * sum(v * v for v in t) != y0:
        cases.append(({"ypeak": str(params["c"] * sum(v * v for v in t))}, {"ypeak": 0}))   # ignored the noise
    for m in _rho_mistakes(params, correct):
        cases.append(({"rho": repr(float(m))}, {"rho": 0}))
    n1 = params["n1"]
    v1 = int(y[n1 - ys]) if 0 <= n1 - ys < len(y) else 0
    nb = [int(y[n1 - ys + d]) for d in (-1, 1) if 0 <= n1 - ys + d < len(y)]
    for v in nb:
        if v != v1:
            cases.append(({"yn1": str(v)}, {"yn1": 0}))
            break
    return cases
