"""Independent check of questions/dft/dft-plot-matching.

x[n] is rebuilt from the ANSWERS (A1, A2, and Ω_i T = π·W_i·T, or ω_i = π·W_i for the lecture type),
the four DFTs are recomputed with numpy (numpy.hamming(32) for the window), and every plotted value
list must equal the DFT its dropdown names. The plots must be unmistakable: DFT I alone has 16 stems,
DFT II's even bins repeat DFT I, DFT III's peaks are >= 1.6 times DFT II's, DFT IV has nothing beyond
two bins of its peaks while DFT II has visible side lobes there, and DFT I leaks. The SVGs must carry
16 / 32 stems, y tick labels at both DFT III peak heights, and the bins must be readable (labels at
even k). Number ranges and wording per variant type are checked against the source exams."""

import re

import numpy as np


def _rebuild(params, correct):
    lec = params["is_lecture"]
    T = 1.0 if lec else 1.0 / params["inv_T"]
    n = np.arange(32)
    x = correct["A1"] * np.cos(np.pi * correct["W1"] * T * n) + correct["A2"] * np.cos(np.pi * correct["W2"] * T * n)
    return {"I": np.abs(np.fft.fft(x[:16])), "II": np.abs(np.fft.fft(x[:16], 32)),
            "III": np.abs(np.fft.fft(x)), "IV": np.abs(np.fft.fft(np.hamming(32) * x))}


def check(params, correct):
    probs = []
    D = _rebuild(params, correct)
    named = {}
    for L in "abcd":
        name = correct["plot_" + L]["html"].strip().replace("DFT ", "")
        named[name] = np.array(params["vals_" + L])
        ref = D[name]
        if len(ref) != len(params["vals_" + L]) or np.max(np.abs(ref - named[name])) > 6e-3:
            probs.append(f"plot ({L}) does not match DFT {name} of the signal rebuilt from the answers")
        svg = params[f"fig_{L}_html"]
        stems = svg.count("<circle")
        if stems != len(ref):
            probs.append(f"plot ({L}) draws {stems} stems, expected {len(ref)}")
    if sorted(named) != ["I", "II", "III", "IV"]:
        probs.append(f"dropdown answers are not a permutation: {sorted(named)}")
        return probs
    I, II, III, IV = (named[k] for k in ("I", "II", "III", "IV"))
    # unmistakable by design
    if np.max(np.abs(II[::2] - I)) > 2e-3:
        probs.append("DFT II even bins differ from DFT I")
    if not III.max() >= 1.6 * II.max():
        probs.append(f"DFT III peak {III.max()} not clearly above DFT II's {II.max()}")
    pk = [k for k in range(32) if III[k] > 1e-6]
    if len(pk) != 4:
        probs.append(f"DFT III has {len(pk)} nonzero bins, expected 4 (on-grid tones)")
    far = [k for k in range(32) if min(min(abs(k - q), 32 - abs(k - q)) for q in pk) > 2]
    if far and IV[far].max() > 0.03 * IV.max():
        probs.append(f"DFT IV side lobes up to {IV[far].max():.3f} (peak {IV.max():.3f})")
    if not far or II[far].max() < 0.1 * II.max():
        probs.append("DFT II shows no visible side lobes away from the peaks")
    near = [q + d for q in pk for d in (-1, 1)]
    if IV[[k % 32 for k in near]].min() < 0.25 * IV[pk].min():
        probs.append("DFT IV main lobes are not visibly three stems wide")
    if np.sum(I > 0.05 * I.max()) <= 4:
        probs.append("DFT I does not leak (both tones on the 16-point grid)")
    # tones, grid and ranges
    k_of = {}
    T = 1.0 if params["is_lecture"] else 1.0 / params["inv_T"]
    for i in (1, 2):
        k = correct[f"W{i}"] * T * 16          # ω = 2πk/32 = πk/16
        if abs(k - round(k)) > 1e-9:
            probs.append(f"tone {i} not on the 32-point grid (k = {k})")
        k_of[i] = int(round(k))
        h = III[k_of[i]]
        if abs(h - 16 * correct[f"A{i}"]) > 1e-3:
            probs.append(f"DFT III height at k = {k_of[i]} is {h}, not 16 A{i}")
        if not 3 <= k_of[i] <= 13:
            probs.append(f"tone bin {k_of[i]} outside [3, 13]")
        lab = f'text-anchor="end">{int(round(16 * correct[f"A{i}"]))}</text>'
        svg_iii = params[f"fig_{params['plot_of']['III']}_html"]
        if abs(16 * correct[f"A{i}"] - round(16 * correct[f"A{i}"])) > 1e-9 or lab not in svg_iii:
            probs.append(f"DFT III peak height {16 * correct[f'A{i}']} has no y tick label")
    gap = abs(k_of[1] - k_of[2])
    if gap < (3 if params["is_fa25"] else 4):
        probs.append(f"bins only {gap} apart")
    if correct["A1"] == correct["A2"]:
        probs.append("equal amplitudes: the numbering would be ambiguous")
    if params["is_fa24"]:
        if not correct["A1"] > correct["A2"] or correct["A1"] != int(correct["A1"]) or correct["A2"] != int(correct["A2"]):
            probs.append("FA2024 type needs integer A1 > A2")
        if abs(params["fs"] * float(params["D_tex"]) - 32) > 1e-9:
            probs.append("f_s times the duration is not 32 samples")
    else:
        if not correct["W1"] < correct["W2"]:
            probs.append("tones not numbered by frequency")
    if params["is_fa25"] and params["inv_T"] != 100:
        probs.append("FA2025 type must keep T = 1/100 s")
    if not params["is_lecture"]:
        for i in (1, 2):
            if not 0 < correct[f"W{i}"] < params["inv_T"]:     # Ω < π/T
                probs.append(f"Ω{i} above the Nyquist frequency")
    # x labels: every even bin labelled on 32-bin plots, every bin on the 16-bin plot
    for L in "abcd":
        svg = params[f"fig_{L}_html"]
        labs = re.findall(r'text-anchor="middle">(\d+)</text>', svg)
        M = len(params["vals_" + L])
        want = [str(k) for k in range(0, M, 1 if M == 16 else 2)]
        if [s for s in labs if s in want] != want:
            probs.append(f"plot ({L}) x labels {labs[:20]}")
    # answer panel text: each line names the right DFT
    for it in params["sol_items"]:
        if correct["plot_" + it["letter"]]["html"].strip() != "DFT " + it["dft"]:
            probs.append(f"solution line for plot ({it['letter']}) names DFT {it['dft']}")
    return probs


def _key(params, name, text):
    return next(o["key"] for o in params[name] if o["html"].strip() == text)


def submissions(params, correct):
    cases = []
    by = {correct["plot_" + L]["html"].strip(): L for L in "abcd"}
    keys = {L: correct["plot_" + L]["key"] for L in "abcd"}
    cases.append(({"plot_" + L: keys[L] for L in "abcd"}, {"plot_" + L: 1 for L in "abcd"}))
    # zero padding mistaken for more data: II <-> III swapped
    a2, a3, a4, a1 = by["DFT II"], by["DFT III"], by["DFT IV"], by["DFT I"]
    cases.append(({"plot_" + a2: _key(params, "plot_" + a2, "DFT III"), "plot_" + a3: _key(params, "plot_" + a3, "DFT II")},
                  {"plot_" + a2: 0, "plot_" + a3: 0}))
    # Hamming mistaken for zero padding: II <-> IV
    cases.append(({"plot_" + a4: _key(params, "plot_" + a4, "DFT II"), "plot_" + a2: _key(params, "plot_" + a2, "DFT IV")},
                  {"plot_" + a4: 0, "plot_" + a2: 0}))
    cases.append(({"plot_" + a1: _key(params, "plot_" + a1, "DFT II")}, {"plot_" + a1: 0}))
    A1, A2, W1, W2 = (correct[k] for k in ("A1", "A2", "W1", "W2"))
    from fractions import Fraction
    fr = {k: Fraction(v).limit_denominator(64) for k, v in (("A1", A1), ("A2", A2), ("W1", W1), ("W2", W2))}
    cases += [
        ({"A1": f"{fr['A1'].numerator}/{fr['A1'].denominator}", "W1": f"{fr['W1'].numerator}/{fr['W1'].denominator}"},
         {"A1": 1, "W1": 1}),
        ({"A2": f"{A2:.4g}", "W2": f"{W2:.6g}"}, {"A2": 1, "W2": 1}),
        ({"A1": f"{A1 * 1.03:.6g}", "W2": f"{W2 * 0.97:.6g}"}, {"A1": 1, "W2": 1}),      # careful reading, within 5 %
        ({"A1": f"{A1 / 2:.6g}", "A2": f"{A2 / 2:.6g}"}, {"A1": 0, "A2": 0}),           # |X|/32 (forgot the 2)
        ({"A1": f"{2 * A1:.6g}", "A2": f"{2 * A2:.6g}"}, {"A1": 0, "A2": 0}),           # 2|X|/16 on DFT III
        ({"A1": f"{A1 * 16:.6g}"}, {"A1": 0}),                                           # raw height
        ({"A1": f"{A1 * 0.53:.6g}"}, {"A1": 0}),                                         # Hamming peak read with 2/32
        ({"A1": f"{A2:.6g}", "A2": f"{A1:.6g}"}, {"A1": 0, "A2": 0}),                    # swapped numbering
        ({"W1": f"{2 * W1:.6g}", "W2": f"{2 * W2:.6g}"}, {"W1": 0, "W2": 0}),            # 2πk/16 on a 32-point plot
        ({"W1": f"{W1 / 2:.6g}"}, {"W1": 0}),                                            # hertz typed in the π box
        ({"W1": f"{W1 * np.pi:.6g}"}, {"W1": 0}),                                        # rad/s typed with π included
        ({"A2": "two"}, {"A2": "invalid"}),
    ]
    if not params["is_lecture"]:
        T = 1.0 / params["inv_T"]
        cases.append(({"W1": f"{W1 * T:.6g}", "W2": f"{W2 * T:.6g}"}, {"W1": 0, "W2": 0}))   # ω instead of Ω
    k1 = round(W1 * (1.0 if params["is_lecture"] else 1.0 / params["inv_T"]) * 16)
    step = W1 / k1
    cases.append(({"W1": f"{(k1 + 1) * step:.6g}"}, {"W1": 0}))                         # one bin off
    cases.append(({"W1": f"{(32 - k1) * step:.6g}"}, {"W1": 0}))                        # mirror bin
    return cases
