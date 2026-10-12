"""Independent check of questions/dft/circular-convolution.

x, h are read from the displayed LaTeX. The N-point result is computed twice, with numpy's FFT
(ifft(fft(x, N) * fft(h, N))) and with the explicit modular sum sum_m x[m] h[(n - m) mod N]; the
smallest N and the radix-2 padding are found by brute force (try N upward until the FFT product
reproduces np.convolve, for x, h and for random sequences of the lengths P, Q of part (d))."""

import re

import numpy as np


def _parse(tex):
    m = re.fullmatch(r"\\\{(\w)\[n\]\\\}_\{n=0\}\^\{(\d+)\} = \\\{(.*)\\\}", tex.strip())
    assert m, tex
    vals = [int(v) for v in m.group(3).split(r",\ ")]
    assert len(vals) == int(m.group(2)) + 1
    return vals


def _circ_fft(x, h, N):
    return np.real(np.fft.ifft(np.fft.fft(x, N) * np.fft.fft(h, N)))


def _circ_sum(x, h, N):
    xp = np.zeros(N)
    hp = np.zeros(N)
    xp[:len(x)] = x
    hp[:len(h)] = h
    return np.array([sum(xp[m] * hp[(n - m) % N] for m in range(N)) for n in range(N)])


def _equals_linear(x, h, N):
    y = np.convolve(x, h)
    if N < max(len(x), len(h)):
        return False
    yl = np.zeros(max(N, len(y)))
    yl[:len(y)] = y
    return bool(np.allclose(_circ_fft(x, h, N), yl[:N]) and not np.any(np.abs(yl[N:]) > 1e-9))


def check(params, correct):
    probs = []
    x, h = _parse(params["x_tex"]), _parse(params["h_tex"])
    if x != params["x"] or h != params["h"]:
        probs.append("displayed sequences differ from the stored ones")
    N = params["Nc"]
    got = np.array(correct["yN"]["_value"][0])
    c1, c2 = _circ_fft(x, h, N), _circ_sum(x, h, N)
    if got.shape != (N,) or not (np.allclose(got, c1) and np.allclose(got, c2)):
        probs.append(f"y_N {got} != fft {np.round(c1, 6)} / modular sum {c2}")
    y = np.convolve(x, h)
    if y[0] == 0 or y[-1] == 0:
        probs.append("linear convolution has a zero end sample")
    if not (max(len(x), len(h)) <= N < len(y)):
        probs.append(f"N = {N} outside [max(L1, L2), L1 + L2 - 1)")
    nmin = next(n for n in range(1, 64) if _equals_linear(x, h, n))
    if correct["Nmin"] != nmin:
        probs.append(f"smallest N is {nmin}, answer {correct['Nmin']}")
    differ = [n for n in range(N) if abs(c2[n] - y[n]) > 1e-9]
    marked = sorted(int(o["html"].strip().strip("$").split("=")[1]) for o in correct["diff"])
    if differ != marked:
        probs.append(f"indices that differ {differ}, marked {marked}")
    if not differ:
        probs.append("no index differs (checkbox needs a correct option)")
    if [int(o["html"].strip().strip("$").split("=")[1]) for o in params["diff"]] != list(range(N)):
        probs.append("checkbox options are not n = 0..N-1 in order")
    P, Q = params["P"], params["Q"]
    rng = np.random.default_rng(P * 100 + Q)
    u, v = rng.integers(1, 5, P), rng.integers(1, 5, Q)
    n2 = next(2 ** k for k in range(12) if _equals_linear(u, v, 2 ** k))
    if (correct["zu"], correct["zv"]) != (n2 - P, n2 - Q):
        probs.append(f"radix-2 zeros should be {(n2 - P, n2 - Q)}, answer {(correct['zu'], correct['zv'])}")
    return probs


def submissions(params, correct):
    yN = [int(round(v)) for v in correct["yN"]["_value"][0]]
    N = params["Nc"]
    y = [int(v) for v in np.convolve(params["x"], params["h"])]
    L = len(y)
    cases = [
        ({"yN": "[" + " ".join(map(str, yN)) + "]"}, {"yN": 1}),                          # MATLAB style
        ({"yN": "[" + ", ".join(f"{v}.0" for v in yN) + "]"}, {"yN": 1}),                  # decimals
        ({"yN": "[" + ", ".join(map(str, y[:N])) + "]"}, {"yN": 0}),                       # truncated, not wrapped
    ]
    wrong_dir = list(y[:N])                                                              # overflow added at the end
    for i, v in enumerate(y[N:]):
        wrong_dir[N - (L - N) + i] += v
    if wrong_dir != yN:
        cases.append(({"yN": "[" + ", ".join(map(str, wrong_dir)) + "]"}, {"yN": 0}))
    xs, hs = params["x"], params["h"]
    hp = list(hs) + [0] * (N - len(hs))
    corr = [sum(xs[m] * hp[(m - n) % N] for m in range(len(xs))) for n in range(N)]   # shifted h the wrong way
    if corr != yN:
        cases.append(({"yN": "[" + ", ".join(map(str, corr)) + "]"}, {"yN": 0}))
    cases.append(({"yN": "[" + ", ".join(map(str, y)) + "]"}, {"yN": 0}))               # the linear result itself
    nmin = correct["Nmin"]
    cases.append(({"Nmin": f"0{nmin}"}, {"Nmin": 1}))
    cases.append(({"Nmin": f" {nmin} "}, {"Nmin": 1}))
    cases.append(({"Nmin": str(max(params["L1"], params["L2"]))}, {"Nmin": 0}))         # longer sequence only
    cases.append(({"Nmin": str(params["L1"] + params["L2"])}, {"Nmin": 0}))             # forgot the -1
    cases.append(({"Nmin": "7.5"}, {"Nmin": "invalid"}))
    # (c) checkbox, net-correct
    opts = params["diff"]
    ck = [o["key"] for o in correct["diff"]]
    c = len(ck)
    wrong = [o["key"] for o in opts if o["key"] not in ck]
    cases.append(({"diff": list(reversed(ck))}, {"diff": 1}))
    last = [o["key"] for o in opts][-c:]                                                # wrapped onto the end
    hit = sum(1 for k in last if k in ck)
    cases.append(({"diff": last}, {"diff": max(0.0, (hit - (c - hit)) / c)}))
    cases.append(({"diff": ck + wrong[:1]}, {"diff": (c - 1) / c}))
    first = [o["key"] for o in opts][:params["W"]]                                     # "the first L - N indices"
    hit = sum(1 for k in first if k in ck)
    cases.append(({"diff": first}, {"diff": max(0.0, (hit - (len(first) - hit)) / c)}))
    cases.append(({"diff": []}, {"diff": "invalid"}))
    # (d)
    P, Q, M, N2 = params["P"], params["Q"], params["M"], params["N2"]
    zu, zv = correct["zu"], correct["zv"]
    cases.append(({"zu": f"0{zu}", "zv": f"+{zv}"}, {"zu": 1, "zv": 1}))
    cases.append(({"zu": f" {zu} ", "zv": f"0{zv}"}, {"zu": 1, "zv": 1}))
    cases.append(({"zu": str(N2), "zv": str(N2)}, {"zu": 0, "zv": 0}))                  # the padded length, not the zeros
    cases.append(({"zu": str(M - P), "zv": str(M - Q)}, {"zu": 0, "zv": 0}))             # plain DFT length P + Q - 1
    big = 1 << (max(P, Q) - 1).bit_length()                                             # power of 2 for the longer one only
    if big != N2:
        cases.append(({"zu": str(big - P), "zv": str(big - Q)}, {"zu": 0, "zv": 0}))
    return cases
