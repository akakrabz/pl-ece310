"""Independent check of questions/dft/window-functions.

Window facts are measured, not assumed: the DTFT of the rectangular window and of numpy.hamming(N)
(0.54 - 0.46 cos(2πn/(N-1))) is evaluated on a dense grid (zero-padded FFT); the first null gives the
main-lobe width and the largest value beyond it the tallest side lobe.
(a) 10^(-L/20) and 20 log10(1/m) with numpy;
(b) the two tones are parsed back from the displayed x[n]; a window is "visible" when A_weak/A_strong
    exceeds that window's measured tallest side lobe; the ratio must be a factor >= 2 away from both
    thresholds and the tones farther apart than the Hamming main lobe;
(c) the Hamming width from the measured first null; N by brute force (smallest N whose main lobe,
    4π/N or 8π/N, is narrower than Δω);
(d) SP2025 #1(f) statements against the measured widths and side lobes; FA2019 #1(c) by summing the
    finite cosine's DTFT; the "which window" scenario by applying both criteria to both windows."""

import math
import re

import numpy as np

NFFT = 1 << 16


def lobe_stats(w):
    """(main-lobe width / pi, tallest side lobe relative to the peak) of window w."""
    W = np.abs(np.fft.fft(w, NFFT))[: NFFT // 2]
    W = W / W[0]
    i = 1
    while W[i + 1] < W[i]:
        i += 1
    return 2 * (2 * i / NFFT), float(W[i:].max())


def windows(N):
    return lobe_stats(np.ones(N)), lobe_stats(np.hamming(N))


def parse_tones(tex):
    """'0.08\\cos(0.3\\pi n) + 4\\cos(0.65\\pi n)' -> [(amp, w/pi), ...]."""
    out = []
    for a, w in re.findall(r"(\\tfrac\{1\}\{10\}|[\d.]*)\\cos\((?:([\d.]+)\\pi|\\pi) n\)", tex):
        amp = 0.1 if a == r"\tfrac{1}{10}" else (float(a) if a else 1.0)
        out.append((amp, float(w) if w else 1.0))
    return out


def _key(params, name, text):
    return next(o["key"] for o in params[name] if o["html"].strip() == text)


def check(params, correct):
    probs = []
    # (a)
    if abs(correct["ratio"] - 10 ** (-params["aL"] / 20)) > 1e-12:
        probs.append(f"ratio {correct['ratio']} != 10^(-{params['aL']}/20)")
    if abs(correct["db"] - 20 * np.log10(1 / params["am"])) > 1e-9:
        probs.append(f"db {correct['db']} != 20 log10(1/{params['am']})")
    if params["a_lect"] and (params["a_win"], params["aL"]) not in (("rectangular", 13), ("Hamming", 42)):
        probs.append(f"lecture level for {params['a_win']} is {params['aL']} dB")
    if not 10 <= abs(correct["db"]) < 100:
        probs.append(f"dB answer {correct['db']} outside [10, 100): 3 significant digits would need 2 decimals")
    # (b)
    spec = params["b_spec"]
    tones = parse_tones(params["b_sig"])
    want = sorted([(spec["A_s"], spec["w_s"]), (spec["A_w"], spec["w_w"])], key=lambda t: t[1])
    if len(tones) != 2 or any(abs(a - b) > 1e-12 or abs(u - v) > 1e-12 for (a, u), (b, v) in zip(tones, want)):
        probs.append(f"displayed x[n] {params['b_sig']} parses to {tones}, spec {want}")
    N = spec["N"]
    (ml_r, sl_r), (ml_h, sl_h) = windows(N)
    ratio = spec["A_w"] / spec["A_s"]
    if abs(ratio - params["b_ratio"]) > 1e-12:
        probs.append("b_ratio does not match the amplitudes")
    vis = {"rectangular window": ratio > sl_r, "Hamming window": ratio > sl_h}
    vis["neither window"] = not (vis["rectangular window"] or vis["Hamming window"])
    keys = {o["key"] for o in correct["vis"]}
    for o in params["vis"]:
        txt = o["html"].strip()
        if vis[txt] != (o["key"] in keys):
            probs.append(f"(b) option {txt}: measured {vis[txt]}, marked {o['key'] in keys}")
    for name, sl in (("rect", sl_r), ("ham", sl_h)):
        if 1 / 2 < ratio / sl < 2:
            probs.append(f"(b) ratio {ratio} within a factor 2 of the {name} side lobe {sl:.4f}")
    sep = abs(spec["w_w"] - spec["w_s"])
    if params["b_lect"]:
        if sep < 0.95 * ml_h:                       # L23 eq. 16: Hamming lobes just touch (Δω = 8π/N)
            probs.append("(b) lecture tones: Hamming main lobes overlap")
    elif sep < 2 * ml_h:
        probs.append("(b) tones closer than two Hamming main lobes, but the question says the lobes do not merge")
    if not 0.2 < sl_r < 0.23 or not 0.006 < sl_h < 0.009:
        probs.append(f"measured side lobes {sl_r:.4f}, {sl_h:.5f} differ from -13 / -42 dB")
    # (c)
    Nc = params["cN"]
    ml = windows(Nc)[1][0]
    if abs(correct["W_ham"] - 8 / Nc) > 1e-12:
        probs.append(f"W_ham {correct['W_ham']} != 8/{Nc}")
    if abs(ml - 8 / Nc) > 0.2 * 8 / Nc:
        probs.append(f"measured Hamming main lobe {ml:.4f} pi far from 8/N")
    cs = params["c_spec"]
    if cs["kind"] == "digital":
        dw = (cs["w2"] - cs["w1"]) * math.pi
        if abs(float(params["c_w1"].replace(r"\pi", "")) - cs["w1"]) > 1e-12:
            probs.append("c_w1 text mismatch")
        if abs(float(params["c_w2"].replace(r"\pi", "")) - cs["w2"]) > 1e-12:
            probs.append("c_w2 text mismatch")
    else:
        dw = 2 * math.pi * (cs["f2"] - cs["f1"]) / cs["fs"]
        if (params["c_f1"], params["c_f2"], params["c_fs"]) != (cs["f1"], cs["f2"], cs["fs"]):
            probs.append("notes text mismatch")
    for name, width in (("N_rect", 4 * math.pi), ("N_ham", 8 * math.pi)):
        n = 1
        while not width / n < dw:
            n += 1
        if correct[name] != n:
            probs.append(f"{name} {correct[name]} but brute force gives {n}")
        b = width / dw
        if abs(b - round(b)) < 0.02:
            probs.append(f"{name}: bound {b:.4f} nearly an integer (strict/non-strict ambiguity)")
    # (d)
    ds = params["d_spec"]
    ans = correct["d"]["html"].strip()
    if ds["kind"] == "sp25":
        (ml_r, sl_r), (ml_h, sl_h) = windows(64)
        truth = ((ds["side"] == "lower") == (sl_h < sl_r)) and ((ds["main"] == "narrower") == (ml_h < ml_r))
        if ans != ("True" if truth else "False"):
            probs.append(f"(d) SP2025 statement {ds}: measured {truth}, answer {ans}")
        if f"{ds['side']} sidelobes and {ds['main']} main lobes" not in params["d_statement"]:
            probs.append("(d) statement text does not match its spec")
    elif ds["kind"] == "fa19":
        n = np.arange(ds["N"])
        w0 = ds["w0"] * math.pi
        om = np.linspace(-math.pi, math.pi, 2001)
        X = np.exp(-1j * np.outer(om, n)) @ np.cos(w0 * n)

        def R(o):
            return np.exp(-1j * np.outer(o, n)).sum(axis=1)
        shifted = 0.5 * R(om - w0) + 0.5 * R(om + w0)
        if np.max(np.abs(X - shifted)) > 1e-9:
            probs.append("(d) finite cosine DTFT differs from the shifted-window form")
        away = np.abs(om - w0) > 4 * math.pi / ds["N"]
        away &= np.abs(om + w0) > 4 * math.pi / ds["N"]
        if not np.max(np.abs(X[away])) > 0.05:
            probs.append("(d) no visible side lobes away from +-w0")
        truth = not ds["impulses"]
        if ans != ("True" if truth else "False"):
            probs.append(f"(d) FA2019 statement impulses={ds['impulses']}: answer {ans}")
        if (r"\delta" in params["d_formula"]) != ds["impulses"]:
            probs.append("(d) formula text does not match its spec")
    else:
        N = ds["N"]
        (ml_r, sl_r), (ml_h, sl_h) = windows(N)
        dw = abs(ds["w2"] - ds["w1"])
        ratio = min(ds["A1"], ds["A2"]) / max(ds["A1"], ds["A2"])
        rect_ok = dw > ml_r and ratio > sl_r
        ham_ok = dw > ml_h and ratio > sl_h
        if rect_ok == ham_ok:
            probs.append(f"(d) scenario not decisive: rect {rect_ok}, Hamming {ham_ok}")
        win = "Rectangular" if rect_ok else "Hamming"
        reason = "main lobe" if rect_ok else "side lobes"
        if not (ans.startswith(win) and reason in ans):
            probs.append(f"(d) scenario answer {ans!r}, expected {win} / {reason}")
        u = dw * N
        if min(abs(u - 4), abs(u - 8)) < 0.9:
            probs.append(f"(d) Δω·N/π = {u} too close to 4 or 8")
        if min(abs(math.log(ratio / sl_r)), abs(math.log(ratio / sl_h))) < math.log(2):
            probs.append(f"(d) ratio {ratio} within a factor 2 of a side-lobe level")
        tones = parse_tones(params["d_sig"])
        if [(a, w) for a, w in tones] != [(ds["A1"], ds["w1"]), (ds["A2"], ds["w2"])]:
            probs.append(f"(d) displayed {params['d_sig']} parses to {tones}")
    # display slips in generated text
    for k in ("b_sig", "d_sig", "d_formula"):
        s = params.get(k, "")
        if re.search(r"(^|[^\d.])1\\cos|\\cos\(1\\pi|\+ -|- -", s):
            probs.append(f"display slip in {k}: {s}")
    return probs


def _num(v):
    return f"{v:.6g}"


def submissions(params, correct):
    r, db = correct["ratio"], correct["db"]
    eps_r = 0.51 * 10 ** (math.floor(math.log10(abs(r))) - 2)
    cases = [
        ({"ratio": f"{r:.3g}", "db": f"{db:.3g}"}, {"ratio": 1, "db": 1}),
        ({"ratio": f"{r:.5g}", "db": f"{db:.6g}"}, {"ratio": 1, "db": 1}),
        ({"ratio": f"{10 ** (-params['aL'] / 10):.3g}"}, {"ratio": 0}),          # power ratio (10 log10)
        ({"ratio": f"{params['aL'] / 20:.3g}"}, {"ratio": 0}),                    # L/20 without the power
        ({"db": f"{-db:.3g}"}, {"db": 0}),                                        # sign
        ({"db": f"{db / 2:.3g}"}, {"db": 0}),                                     # 10 log10
        ({"db": f"{-20 * math.log(params['am']):.3g}"}, {"db": 0}),              # natural log
        ({"ratio": "abc"}, {"ratio": "invalid"}),
    ]
    two = float(f"{r:.2g}")
    if abs(two - r) > 1.2 * eps_r:
        cases.append(({"ratio": f"{two:g}"}, {"ratio": 0}))                       # 2 significant digits
    # (b) checkbox
    keys = sorted(o["key"] for o in correct["vis"])
    rect = _key(params, "vis", "rectangular window")
    ham = _key(params, "vis", "Hamming window")
    nei = _key(params, "vis", "neither window")
    allk = sorted([rect, ham, nei])
    cases.append(({"vis": keys}, {"vis": 1}))
    if rect in keys:                                    # both windows correct
        cases.append(({"vis": [ham]}, {"vis": 0.5}))
        cases.append(({"vis": allk}, {"vis": 0.5}))
        cases.append(({"vis": [nei]}, {"vis": 0}))
    else:                                               # Hamming only
        cases.append(({"vis": [rect, ham]}, {"vis": 0}))
        cases.append(({"vis": [rect]}, {"vis": 0}))
        cases.append(({"vis": [nei]}, {"vis": 0}))
    # (c)
    Nc = params["cN"]
    cases += [
        ({"W_ham": f"8/{Nc}"}, {"W_ham": 1}),
        ({"W_ham": f"{8 / Nc:.6g}", "N_rect": f" {correct['N_rect']} "}, {"W_ham": 1, "N_rect": 1}),
        ({"W_ham": f"4/{Nc}"}, {"W_ham": 0}),                                     # rectangular width / Hamming half-width
        ({"W_ham": f"16/{Nc}"}, {"W_ham": 0}),
        ({"W_ham": "8"}, {"W_ham": 0}),
        ({"N_rect": str(correct["N_rect"] - 1)}, {"N_rect": 0}),                  # rounded down
        ({"N_ham": str(correct["N_rect"])}, {"N_ham": 0}),                        # rectangular N for Hamming
        ({"N_rect": str(correct["N_rect"] // 2 + 1)}, {"N_rect": 0}),             # about the half-lobe bound
        ({"N_ham": str(correct["N_ham"] - 1)}, {"N_ham": 0}),
        ({"N_rect": "12.5"}, {"N_rect": "invalid"}),
    ]
    if 2 * correct["N_rect"] != correct["N_ham"]:
        cases.append(({"N_ham": str(2 * correct["N_rect"])}, {"N_ham": 0}))      # doubled the rounded N
    # (d)
    dk = correct["d"]["key"]
    cases.append(({"d": dk}, {"d": 1}))
    for o in params["d"]:
        if o["key"] != dk:
            cases.append(({"d": o["key"]}, {"d": 0}))
    return cases
