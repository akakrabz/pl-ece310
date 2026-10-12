"""Independent check of questions/sampling/true-false.

Fixed statements: answer table typed from the official keys (study-site Midterm 2 True/False bank,
Q17-Q43) and from the 2019 course notes (Singer & Munson) and Lecture 17, with the source of each
entry in a comment. Parameterized statements: recomputed numerically from the stored parameters --
FFTs of densely sampled test signals (Nyquist rates, aliases), sums of spectral copies, Poisson sums
of actual samples, SymPy's DiracDelta, sinc interpolation, simulated ideal digital filtering, numpy
DFTs on random sequences, and an instrumented radix-2 FFT -- never through server.py's code path."""

import re

import numpy as np
import sympy as sp

OFFICIAL = {   # sid -> answer
    "bandlimited_def": True,          # SP2023 MT2 #1(e), bank Q22 (key: True, definition)
    "every_signal_nyquist": False,    # FA2023 MT2 #1(d), bank Q23 (key: False, only bandlimited signals)
    "larger_T": False,                # FA2021 MT2 #1(e), bank Q24 (key: False, "the DTFT would be expanded"); reworded
    "recon_causal": False,            # FA2021 MT2 #1(c), bank Q29 (key: False, sinc is two-sided)
    "da_through_samples": True,       # 2019 notes eq. (3.3): ideal reconstruction formula interpolates the samples
    "zoh_flat_lpf": False,            # 2019 notes (ZOH chapter): F_a must be e^{jWT/2}/sinc(WT/2), not flat
    "adda_lti_bandlimited": True,     # 2019 notes section 5.3: bandlimited input, ideal converters -> H_d(WT)
    "adda_always_ti": False,          # 2019 notes section 5.3: "linear ... in general not time-invariant"
    "dft_same": False,                # FA2021 MT2 #1(b), bank Q33 (key: False, DFT invertible)
    "zp_mainlobe": False,             # 2019 notes section 2.6.2: zero padding samples the same DTFT more densely
    "more_samples_narrower": True,    # 2019 notes section 2.6.1: main lobe width 4 pi / N
    "fft_is_dft": True,               # FA2019 MT2 #1(j), bank Q43 (key: True)
    "fft_power2": True,               # 2019 notes section 2.6.2: pad to a length with a fast FFT (power of 2)
}

RNG = np.random.default_rng(17)


def f(pq):
    return pq[0] / pq[1]


def highest_freq(sig, fs=4096):
    """Highest frequency (Hz) present in sig(t), from a 1-second FFT at fs (integer-Hz tones land on bins)."""
    t = np.arange(fs) / fs
    X = np.abs(np.fft.rfft(sig(t)))
    k = np.nonzero(X > 1e-6 * X.max())[0]
    return int(k.max())


def lpf_out_freq(f0, fs, wc):
    """Simulate A/D at fs (1 s of samples), an ideal digital low-pass (cutoff wc rad/sample) in the DFT domain and the
    ideal D/A: frequency (Hz) of the output tone, or None if the output is zero."""
    n = np.arange(fs)
    X = np.fft.fft(np.cos(2 * np.pi * f0 * n / fs))
    k = np.arange(fs)
    wk = 2 * np.pi * np.minimum(k, fs - k) / fs
    X[wk > wc] = 0
    y = np.fft.ifft(X).real
    if np.abs(y).max() < 1e-6:
        return None
    Y = np.abs(np.fft.rfft(y))
    return int(np.argmax(Y))           # bin index = Hz, in [0, fs/2]: the bandlimited tone the ideal D/A outputs


def fft_structure(N):
    """Radix-2 decimation-in-time FFT with counters: (number of stages, butterflies per stage)."""
    x = RNG.standard_normal(N) + 1j * RNG.standard_normal(N)
    bits = N.bit_length() - 1
    rev = [int(format(i, f"0{bits}b")[::-1], 2) for i in range(N)]
    a = x[rev].astype(complex)
    stages, per_stage = 0, []
    size = 2
    while size <= N:
        half, count = size // 2, 0
        tw = np.exp(-2j * np.pi * np.arange(half) / size)
        for start in range(0, N, size):
            for j in range(half):
                u, v = a[start + j], a[start + j + half] * tw[j]
                a[start + j], a[start + j + half] = u + v, u - v
                count += 1
        stages += 1
        per_stage.append(count)
        size *= 2
    assert np.allclose(a, np.fft.fft(x))
    assert len(set(per_stage)) == 1
    return stages, per_stage[0]


_FFT_CACHE = {}


def recompute(sid, prm, probs):
    if sid in OFFICIAL:
        return OFFICIAL[sid]
    if sid == "nyq_scale":                    # test signal with tones at 36 and 24 Hz: Nyquist rate F = 72 Hz
        c, claim = f(prm["c"]), f(prm["claim"])
        F = 2 * highest_freq(lambda t: np.cos(2 * np.pi * 36 * t) + 0.5 * np.cos(2 * np.pi * 24 * t))
        Fy = 2 * highest_freq(lambda t: np.cos(2 * np.pi * 36 * c * t) + 0.5 * np.cos(2 * np.pi * 24 * c * t))
        return bool(abs(Fy - claim * F) < 1e-9)
    if sid == "nyq_product":
        f1, f2, form = prm["f1"], prm["f2"], prm["form"]
        sig = {"prod": lambda t: np.cos(2 * np.pi * f1 * t) * np.sin(2 * np.pi * f2 * t),
               "square": lambda t: np.cos(2 * np.pi * f1 * t) ** 2,
               "sum": lambda t: np.cos(2 * np.pi * f1 * t) + 3 * np.sin(2 * np.pi * f2 * t)}[form]
        return 2 * highest_freq(sig) == prm["R"]
    if sid == "cos_to_sequence":
        fn = np.cos if prm["kind"] == "cos" else np.sin
        n = np.arange(41)
        x = fn(2 * np.pi * prm["f0"] * n * f(prm["T"]))
        return bool(np.allclose(x, fn(np.pi * f(prm["p"]) * n), atol=1e-9))
    if sid == "alias_reconstruct":
        fs, f0, fc = prm["fs"], prm["f0"], prm["fc"]
        n = np.arange(fs)
        X = np.abs(np.fft.rfft(np.cos(2 * np.pi * f0 * n / fs)))
        fout = int(np.argmax(X))              # the unique tone below fs/2 with these samples
        return fout == fc and fc < fs / 2
    if sid == "guard_band":                   # add the spectral copies on a grid of (0, pi)
        B, fs = prm["B"], prm["fs"]
        w = np.linspace(0, np.pi, 4001)[1:-1]
        Xa = lambda W: np.clip(1 - np.abs(W) / (2 * np.pi * B), 0, None)   # noqa: E731
        Xd = fs * sum(Xa((w - 2 * np.pi * k) * fs) for k in range(-4, 5))
        return bool(np.any(Xd < 1e-12))
    if sid == "height_scaling":               # X_d(0) = sum_n x_a(n/fs), x_a(t) = A B sinc^2(B t)
        A, B, fs = prm["A"], prm["B"], prm["fs"]
        n = np.arange(-200000, 200001)
        Xd0 = np.sum(A * B * np.sinc(B * n / fs) ** 2)
        return bool(abs(Xd0 - f(prm["V"])) < 1e-3 * max(1.0, abs(Xd0)))
    if sid == "sifting":
        t = sp.symbols("t", real=True)
        b = sp.oo if prm["b"] is None else prm["b"]
        val = sp.integrate((t + prm["c"]) * sp.DiracDelta(t - prm["t0"]), (t, prm["a"], b))
        return bool(sp.simplify(val - prm["V"]) == 0)
    if sid == "da_impulse":                   # ideal reconstruction (T = 1) vs the claimed waveform at off-sample times
        A, k = prm["A"], prm["k"]
        tt = k + np.array([0.25, 0.5, -0.6, 1.3])
        xa = A * np.sinc(tt - k)
        claim = {"sinc": A * np.sinc(tt - k), "impulse": 0 * tt, "rect": np.where((tt >= k) & (tt < k + 1), A, 0.0)}[prm["form"]]
        return bool(np.allclose(xa, claim))
    if sid == "da_value":
        t0 = f(prm["t0"])
        xa = prm["a"] * np.sinc(t0 - prm["n1"]) + prm["b"] * np.sinc(t0 - prm["n2"])
        if t0 % 1 and abs(xa) < 0.05:
            probs.append("da_value: half-sample value too close to 0")
        return bool(abs(xa - prm["V"]) < 1e-9)
    if sid == "zoh_value":
        y, Tms, t0 = prm["y"], prm["T"], f(prm["t0"])
        n = int(np.floor(t0 / Tms))
        val = y[n] if 0 <= n < len(y) else 0
        return bool(abs(val - f(prm["V"])) < 1e-12)
    if sid == "adda_lowpass":
        fout = lpf_out_freq(prm["f0"], prm["fs"], np.pi * f(prm["wc"]))
        return fout == prm["fc"]
    if sid == "adda_cutoff":                  # tones 2 % below and above the claimed cutoff
        fs, wc, fc = prm["fs"], np.pi * f(prm["wc"]), prm["fc"]
        if fc >= fs / 2:
            return False
        lo, hi = int(np.floor(0.98 * fc)), int(np.ceil(1.02 * fc))
        return lpf_out_freq(lo, fs, wc) == lo and lpf_out_freq(hi, fs, wc) is None
    if sid == "dft_X0":
        X = np.fft.fft(np.array(prm["x"], float))
        return bool(abs(X[prm["k"]] - prm["V"]) < 1e-9)
    if sid == "dft_real_k":
        N, k1, k2 = prm["N"], prm["k1"], prm["k2"]
        ok = True
        for _ in range(6):
            X = np.fft.fft(RNG.standard_normal(N))
            ok &= bool(abs(X[k1].imag) < 1e-9) if prm["kind"] == "real" else bool(abs(X[k1] - np.conj(X[k2])) < 1e-9)
        return ok
    if sid == "idft_periodic":
        Xk = np.array(prm["X"], float)
        N = len(Xk)
        n = np.arange(-8, 12)
        x = np.array([np.mean(Xk * np.exp(2j * np.pi * np.arange(N) * m / N)) for m in n])
        xn = dict(zip(n.tolist(), x))
        form = prm["form"]
        if form == "zero":
            return bool(all(abs(xn[m]) < 1e-12 for m in n.tolist() if m < 0 or m > N - 1))
        if form == "periodic":
            return bool(all(abs(xn[m + N] - xn[m]) < 1e-12 for m in range(-8, 12 - N)))
        return bool(abs(xn[N] - f(prm["V"])) < 1e-12)
    if sid == "dft_dtft_sample":
        N, k0, th = prm["N"], prm["k0"], np.pi * f(prm["th"])
        ok = True
        for _ in range(3):
            x = RNG.standard_normal(N) + 1j * RNG.standard_normal(N)
            ok &= bool(abs(np.fft.fft(x)[k0] - np.sum(x * np.exp(-1j * th * np.arange(N)))) < 1e-9)
        return ok
    if sid == "zp_index":
        L, M = prm["L"], prm["M"]
        ok = True
        for _ in range(3):
            x = RNG.standard_normal(L) + 1j * RNG.standard_normal(L)
            ok &= bool(abs(np.fft.fft(x)[prm["k0"]] - np.fft.fft(x, M)[prm["m0"]]) < 1e-9)
        return ok
    if sid == "cos_bins":
        N, c = prm["N"], f(prm["c"])
        x = np.cos(2 * np.pi * c * np.arange(N) / N)
        count = int(np.sum(np.abs(np.fft.fft(x)) > 1e-8 * N))
        return count == prm["m"]
    if sid == "hamming":                      # main-lobe width and peak sidelobe of both windows, N = 64
        N = 64
        w = np.linspace(1e-6, np.pi, 20001)

        def lobes(win):
            W = np.abs(np.exp(-1j * np.outer(w, np.arange(N))) @ win)
            i = 1
            while W[i] <= W[i - 1]:
                i += 1                         # first local minimum = edge of the main lobe
            return w[i - 1], W[i:].max() / W[0]
        null_r, side_r = lobes(np.ones(N))
        null_h, side_h = lobes(np.hamming(N))
        lower, wider = side_h < side_r, null_h > null_r
        return prm["side"] == ("lower" if lower else "higher") and prm["main"] == ("wider" if wider else "narrower")
    if sid == "rect_mainlobe":
        N, c = prm["N"], np.pi * f(prm["c"])

        def amp(v):                            # W_d(w) e^{jw(N-1)/2}: real, changes sign at each zero
            v = np.atleast_1d(v)
            return (np.exp(-1j * np.outer(v, np.arange(N))).sum(axis=1) * np.exp(1j * v * (N - 1) / 2)).real
        A = amp(np.linspace(1e-4, c * 0.999, 4000))
        return bool(abs(amp(c)[0]) < 1e-9 and np.all(A > 0))
    if sid == "circ_shift":
        x, m = np.array(prm["x"], float), prm["m"]
        N = len(x)
        y = np.fft.ifft(np.exp(-2j * np.pi * np.arange(N) * m / N) * np.fft.fft(x)).real
        return bool(np.allclose(y, prm["claim"]))
    if sid == "circ_conv":
        y = np.fft.ifft(np.fft.fft(prm["x"]) * np.fft.fft(prm["h"])).real
        return bool(np.allclose(y, prm["claim"]))
    if sid == "circ_equals_linear":
        Lx, Lh, N = prm["Lx"], prm["Lh"], prm["N"]
        ok = True
        for _ in range(3):
            x, h = RNG.standard_normal(Lx), RNG.standard_normal(Lh)
            lin = np.convolve(x, h)
            circ = np.fft.ifft(np.fft.fft(x, N) * np.fft.fft(h, N)).real
            ok &= bool(len(lin) <= N and np.allclose(circ[:len(lin)], lin) and np.allclose(circ[len(lin):], 0))
        return ok
    if sid == "fft_count":
        N = 2 ** prm["m"]
        if N not in _FFT_CACHE:
            _FFT_CACHE[N] = fft_structure(N)
        return _FFT_CACHE[N] == (prm["s"], prm["b"])
    raise KeyError(sid)


# ----------------------------------------------------------------------------- display slips
# Typesetting slips that KaTeX accepts but students would misread; checked in every generated variant.
SLIPS = [
    (re.compile(r"e\^\{\s*[+-]?\s*j\s*[-+\u2212]"), "sign right after j in an exponent (e^{j-...})"),
    (re.compile(r"(?:\+|(?<![<>])-|\u2212)\s*(?:-|\u2212)(?![-\u2212>])"), "double sign (+ -, - -)"),
    (re.compile(r"(?<![\d.\w}\]])1\s*(?:\\[,;:!]\s*)?(?:e\^\{|\\cos|\\sin|\\delta|\\pi|\\omega)"), "coefficient 1 written out"),
    (re.compile(r"\}1\s*[a-z]\b"), "factor 1 glued to a symbol"),
    (re.compile(r"(?<![\d.\w])[01]\s*T(?![\w])"), "0T or 1T"),
    (re.compile(r"(?<![\d.\w])-\s*0(?![\d.])"), "minus zero"),
    (re.compile(r"\\left\(\s*-?\d+\s*\\right\)\s*\\(?:cos|sin)"), "parenthesised integer before cos/sin"),
    (re.compile(r"e\^\{-?j0\}"), "e^{j0}"),
]
TEXT_SLIPS = [(re.compile(r"(?<![\d.])1\$?\s+samples\b"), "1 samples")]
MATH_RX = re.compile(r"\$\$(.+?)\$\$|\$(.+?)\$", re.S)


def display_slips(sid, html):
    out = []
    for a, b in MATH_RX.findall(html):
        e = a or b
        for rx, what in SLIPS:
            m = rx.search(e)
            if m:
                out.append(f"{sid}: display slip ({what}): ...{e[max(0, m.start() - 20):m.end() + 20]}...")
    for rx, what in TEXT_SLIPS:
        if rx.search(html):
            out.append(f"{sid}: display slip ({what})")
    return out


def check(params, correct):
    probs = []
    items = params["statements"]
    if len(items) != 6:
        probs.append(f"{len(items)} statements instead of 6")
    if len({it["sid"] for it in items}) != len(items):
        probs.append("duplicate statement")
    for it in items:
        probs += display_slips(it["sid"], it["stmt_html"]) + display_slips(it["sid"], it["why"])
        want = recompute(it["sid"], it["param"], probs)
        graded = correct[it["name"]]["html"].strip()
        if graded != ("True" if want else "False"):
            probs.append(f"{it['name']} ({it['sid']}, {it['param']}): graded {graded}, recomputed {want}")
        if it["ans"] != graded:
            probs.append(f"{it['name']}: answer panel says {it['ans']} but grading uses {graded}")
    return probs


def _key(params, name, text):
    for opt in params[name]:
        if opt["html"].strip() == text:
            return opt["key"]
    raise KeyError(name)


def submissions(params, correct):
    names = [it["name"] for it in params["statements"]]
    right = {nm: correct[nm]["html"].strip() for nm in names}
    flip = {"True": "False", "False": "True"}
    cases = []
    for nm in names[:3]:                    # one statement flipped: that one scores 0, the other five 1 (total 5/6)
        cases.append(({nm: _key(params, nm, flip[right[nm]])}, {m: (0 if m == nm else 1) for m in names}))
    cases.append(({nm: _key(params, nm, flip[right[nm]]) for nm in names}, {nm: 0 for nm in names}))    # all flipped
    cases.append(({nm: _key(params, nm, right[nm]) for nm in names}, {nm: 1 for nm in names}))         # all correct
    cases.append(({names[0]: correct[names[0]]["key"]}, {names[0]: 1}))
    cases.append(({nm: _key(params, nm, "False") for nm in names},                                    # "all False" guess
                  {nm: (1 if right[nm] == "False" else 0) for nm in names}))
    return cases
