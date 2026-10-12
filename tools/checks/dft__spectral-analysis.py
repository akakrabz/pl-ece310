"""Independent check of questions/dft/spectral-analysis.

(a) the spacing and the bin frequency from numpy.fft.fftfreq; the alias by brute force: cosines at
    every multiple of the bin spacing in (f_s/2, 3 f_s] are sampled and the first one whose N-point DFT
    peaks at k_lo is taken.
(b) the plotted magnitudes must equal |numpy.fft.fft| of x[n] rebuilt from the ANSWERS (A1, A2, Omega1,
    Omega2, with random phases), the frequencies must be below the Nyquist frequency and Omega1 < Omega2,
    and every peak height must be written in the figure."""

import re

import numpy as np

RNG = np.random.default_rng(9)


def differs(wrong, right):
    return abs(float(wrong) - float(right)) > 2 * (1e-6 + 1e-3 * abs(float(right)))


def check(params, correct):
    probs = []
    N, fs = params["Na"], params["fsa"]
    freqs = np.fft.fftfreq(N, d=1.0 / fs)
    if abs(correct["df"] - fs / N) > 1e-9:
        probs.append(f"df {correct['df']} != fs/N = {fs / N}")
    if abs(correct["f_hi"] - freqs[params["k_hi"]]) > 1e-6:
        probs.append(f"f_hi {correct['f_hi']} != fftfreq {freqs[params['k_hi']]}")
    if not params["k_hi"] > N / 2 > params["k_lo"] > 0:
        probs.append("k_lo / k_hi on the wrong sides of N/2")
    n = np.arange(N)
    alias = None
    for f in np.arange(fs / N, 3 * fs + 1, fs / N):
        if f <= fs / 2 + 1e-9:
            continue
        X = np.abs(np.fft.fft(np.cos(2 * np.pi * f * n / fs)))
        if abs(X[params["k_lo"]] - N / 2) < 1e-6 and int(np.argmax(X[: N // 2 + 1])) == params["k_lo"]:
            alias = f
            break
    if alias is None or abs(alias - correct["f_alias"]) > 1e-6:
        probs.append(f"f_alias {correct['f_alias']} but brute force gives {alias}")
    # (b)
    Nb, fsb = params["Nb"], params["fsb"]
    T = 1.0 / fsb
    if abs(Nb - fsb * float(params["Db_tex"])) > 1e-9:
        probs.append(f"N = {Nb} but fs * D = {fsb * float(params['Db_tex'])}")
    t = np.arange(Nb) * T
    A1, A2, W1, W2 = (correct[k] for k in ("A1", "A2", "W1", "W2"))
    for _ in range(2):
        ph = RNG.uniform(0, 2 * np.pi, 2)
        x = A1 * np.cos(np.pi * W1 * t + ph[0]) + A2 * np.cos(np.pi * W2 * t + ph[1])
        if not np.allclose(np.abs(np.fft.fft(x)), params["mag"], atol=1e-6):
            probs.append("|fft| of the signal rebuilt from the answers does not match the plot")
            break
    if not (0 < W1 < W2 < fsb):                                  # Omega/pi < f_s  <=>  Omega < pi f_s
        probs.append(f"frequencies {W1}, {W2} not ordered / above Nyquist (f_s = {fsb})")
    if A1 == A2 or A1 != int(A1) or A2 != int(A2):
        probs.append("amplitudes not distinct integers")
    ticks = sorted((float(xp), lab) for xp, lab in
                   re.findall(r'<text x="([\d.]+)" y="[\d.]+" text-anchor="middle">(\d+)</text>', params["fig_html"]))
    for (x0, a), (x1, b) in zip(ticks, ticks[1:]):       # x-axis index labels must stay apart
        if Nb >= 24 and int(b) - int(a) < 2:
            probs.append(f"tick labels {a} and {b} are on adjacent bins")
        if x1 - x0 < 3.8 * (len(a) + len(b)) + 6:        # ~7.5 px per digit at 13 px, half widths + margin
            probs.append(f"tick labels {a} and {b} are only {x1 - x0:.1f} px apart")
    for k in (params["k1"], params["k2"], Nb - params["k1"], Nb - params["k2"]):
        if str(k) not in [lab for _, lab in ticks]:
            probs.append(f"peak bin {k} is not labelled")
    for h in set(params["mag"]) - {0}:
        if f">{int(h)}</text>" not in params["fig_html"]:
            probs.append(f"peak height {h} is not written in the figure")
    return probs


def _num(v):
    v = float(v)
    return str(int(round(v))) if abs(v - round(v)) < 1e-12 else f"{v:.6g}"


def submissions(params, correct):
    N, fs = params["Na"], params["fsa"]
    df, fhi, fal = correct["df"], correct["f_hi"], correct["f_alias"]
    Nb, fsb = params["Nb"], params["fsb"]
    A1, A2, W1, W2 = (correct[k] for k in ("A1", "A2", "W1", "W2"))
    k1, k2 = params["k1"], params["k2"]
    cases = [
        ({"df": f"{fs}/{N}", "f_hi": f"{_num(fhi)}.0"}, {"df": 1, "f_hi": 1}),
        ({"df": f"{df:.4g}", "f_alias": f"{2 * fal:g}/2"}, {"df": 1, "f_alias": 1}),
        ({"df": _num(N / fs)}, {"df": 0}),                                   # N / f_s
        ({"df": _num(2 * np.pi * df)}, {"df": 0}),                           # rad/s instead of Hz
        ({"f_hi": _num(params["k_hi"] * df)}, {"f_hi": 0}),                  # forgot the upper half is negative
        ({"f_hi": _num(-fhi)}, {"f_hi": 0}),                                 # sign
        ({"f_alias": _num(fs + params["k_lo"] * df)}, {"f_alias": 0}),       # the next alias, not the smallest
        ({"f_alias": _num(fs / 2 + params["k_lo"] * df)}, {"f_alias": 0}),
        ({"A1": f"{2 * A1:g}/2", "A2": f"{A2:.1f}"}, {"A1": 1, "A2": 1}),
        ({"W1": f"{W1:g}", "W2": f"{2 * W2:g}/2"}, {"W1": 1, "W2": 1}),
        ({"A1": f"{2 * params['h1']}/{Nb}", "A2": f"{2 * params['h2']}/{Nb}"}, {"A1": 1, "A2": 1}),  # 2|X|/N as typed
        ({"A1": _num(params["h1"] / Nb), "A2": _num(params["h2"] / Nb)}, {"A1": 0, "A2": 0}),       # |X|/N
        ({"A1": _num(params["h1"]), "A2": _num(params["h2"])}, {"A1": 0, "A2": 0}),                 # raw height
        ({"A1": _num(A2), "A2": _num(A1)}, {"A1": 0, "A2": 0}),                                     # swapped
        ({"W1": _num(2 * k1 / Nb), "W2": _num(2 * k2 / Nb)}, {"W1": 0, "W2": 0}),                   # rad/sample (k/N, no f_s)
        ({"W1": _num(k1 * fsb / Nb), "W2": _num(k2 * fsb / Nb)}, {"W1": 0, "W2": 0}),               # hertz typed in the pi box
        ({"W1": _num(np.pi * W1), "W2": _num(np.pi * W2)}, {"W1": 0, "W2": 0}),                     # rad/s typed in the pi box
        ({"W1": _num(2 * (Nb - k1) * fsb / Nb)}, {"W1": 0}),                                        # mirror bin
        ({"A1": "3 apples"}, {"A1": "invalid"}),
    ]
    return cases
