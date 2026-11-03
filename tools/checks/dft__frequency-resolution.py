"""Independent check of questions/dft/frequency-resolution.

l23   the tones are parsed back from the displayed LaTeX; numpy converts them to digital frequencies and
      takes the smallest separation over ALL pairs; the bounds 4pi/dw, 2pi/dw, 8pi/dw and the time
      (4pi/dw) T are recomputed; brute force: the smallest integer N whose main lobes 4pi/N do not
      overlap for any pair is floor(bound) + 1; the marked MC option must say how N_full scales with
      the new f_s (recomputed) and that the time is unchanged.
fa23  |X[k]| = |numpy.fft.fft(cos(lambda0 n))| must equal the plotted values (and the three written ones);
      the tallest bin, the estimate and the bracket are checked against the true lambda0, and the
      bracket must contain every frequency on a fine grid whose DFT shows the same picture (same
      tallest bin, same taller neighbour); the friend's plot must be numpy's 2N-point DFT of the
      zero-padded record or of 2N samples, as the marked option says, and the two must look different
      (even bins equal the old ones vs. a peak about twice as tall)."""

import html
import re

import numpy as np


def _tones(params):
    return [int(f) for f in re.findall(r"(\d+)\\ \\text\{Hz\}", params["tones_tex"])]


def _fs(params):
    t = params["fs_tex"]
    m = re.fullmatch(r"([\d.]+)\\ \\text\{(kHz|Hz)\}", t)
    return float(m.group(1)) * (1000 if m.group(2) == "kHz" else 1)


def check(params, correct):
    probs = []
    if params["l23"]:
        fr = np.array(_tones(params), float)
        fs = _fs(params)
        if abs(fs - params["fs"]) > 1e-9:
            probs.append(f"f_s shown {fs} vs {params['fs']}")
        if len(fr) not in (2, 3) or not np.all(np.diff(fr) > 0) or fr[-1] >= fs / 2:
            probs.append(f"tones {fr}")
        w = 2 * np.pi * fr / fs
        dw = min(abs(a - b) for i, a in enumerate(w) for b in w[i + 1:])
        for name, val in (("Nf", 4 * np.pi / dw), ("Nh", 2 * np.pi / dw), ("Nham", 8 * np.pi / dw), ("Tf", 4 * np.pi / dw / fs)):
            if abs(correct[name] - val) > 1e-9 * val:
                probs.append(f"{name} {correct[name]} vs {val}")
        Nmin = next(N for N in range(2, 10 ** 6)
                    if all(abs(a - b) > 4 * np.pi / N * (1 + 1e-12) for i, a in enumerate(w) for b in w[i + 1:]))
        if Nmin != int(np.floor(correct["Nf"] + 1e-9)) + 1:
            probs.append(f"brute-force smallest N {Nmin} vs bound {correct['Nf']}")
        gaps = np.diff(fr)
        if len(gaps) == 2 and gaps[0] == gaps[1]:
            probs.append("two equal gaps: the closest pair is ambiguous")
        marked = html.unescape(correct["fsmc"]["html"])
        word = "doubles" if params["fs_change"] == "doubled" else "halves"
        if not (marked.startswith(r"$N_{\text{full}}$ " + word) and "is unchanged" in marked):
            probs.append(f"marked fs option wrong: {marked}")
        if params["Nf2_tex"] != str(int(round(4 * np.pi / dw * (2 if word == "doubles" else 0.5)))):
            probs.append(f"Nf at the new rate {params['Nf2_tex']}")
        return probs
    if params["rec"]:
        fs, L, Np = params["fs"], params["L"], params["Np"]
        D = L / fs
        if abs(float(eval(params["D_tex"].replace("\\tfrac{", "(").replace("}{", ")/(").replace("}", ")"))) - D) > 1e-9:
            probs.append(f"duration shown {params['D_tex']} vs L / f_s = {D}")
        # main-lobe half-widths measured on a 64x zero-padded spectrum of the windows themselves
        M = 64 * L
        rect = np.abs(np.fft.fft(np.ones(L), M))
        ham = np.abs(np.fft.fft(np.hamming(L), M))
        z_rect = int(np.argmax(np.diff(rect[: M // 4]) > 0))       # first local minimum
        z_ham = int(np.argmax(np.diff(ham[: M // 4]) > 0))
        hw_rect, hw_ham = 2 * np.pi * z_rect / M, 2 * np.pi * z_ham / M
        if abs(hw_rect - 2 * np.pi / L) > 2 * np.pi / M * 1.01:
            probs.append("rectangular half-width")
        if abs(hw_ham - 4 * np.pi / L) > 0.05 * 4 * np.pi / L:
            probs.append(f"Hamming half-width {hw_ham * L / np.pi:.3f} pi/L, not about 4 pi / L")
        dw_full = 2 * (2 * np.pi / L)
        for name, val in (("dff", dw_full * fs / (2 * np.pi)), ("dfh", dw_full / 2 * fs / (2 * np.pi)),
                          ("dwf", dw_full / np.pi), ("dfham", 2 * (4 * np.pi / L) * fs / (2 * np.pi))):
            if abs(correct[name] - val) > 1e-9 * val:
                probs.append(f"{name} {correct[name]} vs {val}")
        x = np.ones(L)
        if not np.allclose(np.abs(np.fft.fft(x, 8 * Np))[:: 8 * Np // Np], np.abs(np.fft.fft(x, Np)), atol=1e-9):
            probs.append("padding check")
        marked = html.unescape(correct["pad"]["html"])
        if not marked.startswith("The same as without padding"):
            probs.append(f"marked padding option: {marked}")
        if Np <= L or Np & (Np - 1):
            probs.append(f"N = {Np} not a power of two above L = {L}")
        return probs
    N, ks = params["N"], params["ks"]
    d = params["d"][0] / params["d"][1]
    lam = 2 * np.pi * (ks + d) / N
    n = np.arange(N)
    X = np.abs(np.fft.fft(np.cos(lam * n)))
    if not np.allclose(X, params["vals1"], atol=1e-3):
        probs.append("plotted |X[k]| differ from numpy")
    for k in (ks - 1, ks, ks + 1):
        if f">{X[k]:.1f}</text>" not in params["fig1_html"]:
            probs.append(f"value of bin {k} ({X[k]:.1f}) not written in the figure")
    half = X[: N // 2 + 1]
    if int(np.argmax(half)) != ks:
        probs.append(f"tallest bin {int(np.argmax(half))} vs {ks}")
    if abs(correct["lam"] - 2 * ks / N) > 1e-12:
        probs.append("estimate != 2 ks / N")
    a, b = correct["a"], correct["b"]
    if not (a * np.pi <= lam <= b * np.pi) or abs((b - a) - 1 / N) > 1e-12:
        probs.append(f"bracket [{a}, {b}] pi misses lambda0 = {lam / np.pi:.4f} pi or is not half a bin")
    side = X[ks + 1] > X[ks - 1]
    grid = np.linspace(0.02, 0.98, 4801) * np.pi
    same = []
    for g in grid:
        Y = np.abs(np.fft.fft(np.cos(g * n)))[: N // 2 + 1]
        if int(np.argmax(Y)) == ks and (Y[ks + 1] > Y[ks - 1]) == side:
            same.append(g)
    if not same or min(same) < (a - 0.02 / N) * np.pi or max(same) > (b + 0.02 / N) * np.pi:
        probs.append(f"frequencies with the same picture span [{min(same) / np.pi:.4f}, {max(same) / np.pi:.4f}] pi, "
                     f"bracket [{a:.4f}, {b:.4f}]")
    M = 2 * N
    pad = np.abs(np.fft.fft(np.cos(lam * n), M))
    long = np.abs(np.fft.fft(np.cos(lam * np.arange(M))))
    marked = html.unescape(correct["proc"]["html"])
    want = pad if "zero-padded" in marked else long if f"used ${M}$ samples" in marked else None
    if want is None or not np.allclose(want, params["vals2"], atol=1e-3):
        probs.append(f"friend's plot does not match the marked processing: {marked[:50]}")
    if not np.allclose(pad[::2], X, atol=1e-9):
        probs.append("padded even bins differ from the N-point DFT")
    if long.max() < 1.25 * pad.max():
        probs.append("padded and longer-record plots too alike")
    k2 = params["k2"]
    if f">{params['vals2'][k2]:.1f}</text>" not in params["fig2_html"] or k2 != int(np.argmax(np.array(params["vals2"])[: N + 1])):
        probs.append("friend's tallest stem not labelled")
    return probs


def _num(v):
    v = float(v)
    return str(int(round(v))) if abs(v - round(v)) < 1e-12 else f"{v:.6g}"


def submissions(params, correct):
    if params["l23"]:
        fs, df = params["fs"], params["df"]
        Nf, Nh, Tf, Nham = (correct[k] for k in ("Nf", "Nh", "Tf", "Nham"))
        cases = [
            ({"Nf": f"{2 * fs}/{df}", "Nh": f"{fs}/{df}"}, {"Nf": 1, "Nh": 1}),
            ({"Nf": _num(Nf + 1), "Nh": _num(Nh + 1), "Nham": _num(Nham + 1)}, {"Nf": 1, "Nh": 1, "Nham": 1}),   # smallest integer N
            ({"Tf": f"2/{df}"}, {"Tf": 1}),
            ({"Tf": f"{Tf:.4g}"}, {"Tf": 1}),
            ({"Nf": _num(Nh), "Nh": _num(Nf)}, {"Nf": 0, "Nh": 0}),          # criteria swapped
            ({"Nf": _num(Nf + 2)}, {"Nf": 1 if 2 <= 1e-3 * Nf else 0}),
            ({"Nf": _num(2 * np.pi * Nf)}, {"Nf": 0}),                      # Hz used as rad/sample
            ({"Tf": _num(Tf / 2)}, {"Tf": 0}),                             # half-lobe time
            ({"Tf": _num(Nf * fs)}, {"Tf": 0}),                            # N f_s
            ({"Nham": _num(Nf)}, {"Nham": 0}),                             # window ignored
            ({"Nham": _num(4 * Nf)}, {"Nham": 0}),
            ({"Nf": "nan"}, {"Nf": "invalid"}),
        ]
        fr = _tones(params)
        span = fr[-1] - fr[0]
        if span != df:
            cases.append(({"Nf": _num(2 * fs / span)}, {"Nf": 0}))           # extreme pair instead of the closest
        for o in params["fsmc"]:
            cases.append(({"fsmc": o["key"]}, {"fsmc": 1 if o["key"] == correct["fsmc"]["key"] else 0}))
        return cases
    if params["rec"]:
        fs, L, Np = params["fs"], params["L"], params["Np"]
        dff, dfh, dwf, dfham = (correct[k] for k in ("dff", "dfh", "dwf", "dfham"))
        cases = [
            ({"dff": f"{2 * fs}/{L}", "dfh": f"{fs}/{L}"}, {"dff": 1, "dfh": 1}),
            ({"dwf": f"4/{L}", "dfham": f"{4 * fs}/{L}"}, {"dwf": 1, "dfham": 1}),
            ({"dwf": f"{dwf:.4g}"}, {"dwf": 1}),
            ({"dff": _num(dfh), "dfh": _num(dff)}, {"dff": 0, "dfh": 0}),     # criteria swapped
            ({"dff": _num(2 * fs / Np)}, {"dff": 0}),                       # padded length used
            ({"dff": _num(2 * np.pi * dff)}, {"dff": 0}),                   # rad/s
            ({"dwf": _num(2 / L)}, {"dwf": 0}),                             # half-width
            ({"dwf": _num(4 * np.pi / L)}, {"dwf": 0}),                     # radians in the pi box
            ({"dfham": _num(dff)}, {"dfham": 0}),
            ({"dfham": _num(dff / 2)}, {"dfham": 0}),                       # "Hamming is narrower"
        ]
        for o in params["pad"]:
            cases.append(({"pad": o["key"]}, {"pad": 1 if o["key"] == correct["pad"]["key"] else 0}))
        return cases
    N, ks = params["N"], params["ks"]
    lam, a, b = correct["lam"], correct["a"], correct["b"]
    cases = [
        ({"lam": f"{2 * ks}/{N}"}, {"lam": 1}),
        ({"lam": f"{lam:.4g}", "a": f"{a:.4g}", "b": f"{b:.4g}"}, {"lam": 1, "a": 1, "b": 1}),
        ({"lam": _num(ks / N)}, {"lam": 0}),                               # pi k / N grid
        ({"lam": _num(2 * ks / (N - 1))}, {"lam": 0}),                     # 2 pi k / (N - 1)
        ({"lam": _num(2 * (N - ks) / N)}, {"lam": 0}),                     # the mirror bin
        ({"lam": _num(np.pi * lam)}, {"lam": 0}),                          # radians in the pi box
        ({"a": _num(a - 1 / N), "b": _num(b + 1 / N)}, {"a": 0, "b": 0}),   # a whole bin either side
        ({"a": _num(2 * ks / N - 1 / N), "b": _num(2 * ks / N + 1 / N)}, {"a": 1 if a != lam else 0, "b": 1 if b != lam else 0}),
    ]
    for o in params["proc"]:
        cases.append(({"proc": o["key"]}, {"proc": 1 if o["key"] == correct["proc"]["key"] else 0}))
    return cases
