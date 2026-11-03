"""Independent check of questions/dft/true-false.

Every statement is decided again from its stored parameters (and, where the claimed number is only in the
text, from the number parsed out of the statement) with numpy, never through server.py's code path:
DFTs and inverse DFTs of the given or random test sequences, direct DTFT sums, brute-force searches
(which analog tone puts peaks in the stated bins), numerically measured main-lobe widths, peak heights and
side-lobe levels of the rectangular and Hamming windows, simulated windowing for linearity and time
invariance, and partial sums for absolute summability. The conceptual SP2021 #6 statements and two FA2021
#2(c) generalisations come from a table typed from the official key / Lecture 14, each with its source."""

import re

import numpy as np

RNG = np.random.default_rng(310)

TABLE = {   # which_transform ids and hd_substitute fixed forms -> answer
    "dft_cont": False,          # SP2021 MT2 #6: the DFT has N values (discrete k); continuous frequency is CTFT/DTFT
    "dtft_periodic": True,      # SP2021 #6(c) key: "unlike DTFT, [the CTFT] is not always periodic"; checked numerically below
    "dft_computer": True,       # SP2021 #6(b) key: DFT
    "ctft_periodic": False,     # SP2021 #6(c) key: CTFT "not always periodic"
    "dtft_from_z": True,        # SP2021 #6(a) key: DTFT (absolutely summable -> ROC contains the unit circle)
    "dtft_only_N": False,       # Lecture 13: the DTFT is defined for every real w
    "dft_periodic_k": True,     # Lecture 20: X[k + N] = X[k]; checked numerically below
    "every": False,             # FA2021 MT2 #2(c) key: H_d != H(e^{jw}) when the ROC misses the unit circle
    "stable": True,             # Lecture 14: BIBO stable <=> ROC contains |z| = 1 (causal or not)
}


def f(pq):
    return pq[0] / pq[1]


def dtft(x, n0, w):
    n = n0 + np.arange(len(x))
    return np.exp(-1j * np.outer(np.atleast_1d(w), n)) @ np.asarray(x, dtype=complex)


def parse_num(t):
    """A number as typeset by the question: 7, -3, \\tfrac{5}{2}, -\\tfrac{3}{4}, 0.05."""
    t = t.strip()
    m = re.fullmatch(r"(-?)\\tfrac\{(\d+)\}\{(\d+)\}", t)
    if m:
        return (-1 if m.group(1) else 1) * int(m.group(2)) / int(m.group(3))
    return float(t)


def win_dtft(w, grid):
    return np.abs(dtft(w, 0, grid))


def first_null(w, M=20000):
    """First local minimum of |W_d| on (0, pi]."""
    g = np.linspace(0, np.pi, M + 1)[1:]
    a = win_dtft(w, g)
    i = next(i for i in range(1, len(a) - 1) if a[i] <= a[i - 1] and a[i] <= a[i + 1])
    return g[i]


def peak_sidelobe(w):
    g = np.linspace(0, np.pi, 40001)
    a = win_dtft(w, g)
    z = first_null(w)
    return 20 * np.log10(a[g > z * 1.0001].max() / a[0])


def hamming(N):
    n = np.arange(N)
    return 0.54 - 0.46 * np.cos(2 * np.pi * n / (N - 1))


def summable(h):
    s1, s2 = np.abs(h[: len(h) // 2]).sum(), np.abs(h).sum()
    return bool(np.isfinite(s2) and s2 - s1 < 1e-6 * max(1.0, s1))


# ----------------------------------------------------------------------------- recomputation
def recompute(it, probs):
    sid, prm, text = it["sid"], it["param"], it["stmt_html"]
    if sid == "x0_zero":
        X = np.fft.fft(np.array(prm["x"], float))
        return bool(abs(X[prm["k"]]) < 1e-9)
    if sid == "real_bins":
        N, k1, k2 = prm["N"], prm["k1"], prm["k2"]
        ok = True
        for _ in range(4):
            X = np.fft.fft(RNG.standard_normal(N))
            ok &= bool(abs(X[k1] - (np.conj(X[k1]) if prm["kind"] == "real" else np.conj(X[k2]))) < 1e-9)
        return ok
    if sid == "fa25_dtft":
        a, b, c = prm["a"], prm["b"], prm["c"]
        L = 2 * b + 1
        x, n0 = np.array([c / 2] + [0] * (b - 1) + [a] + [0] * (b - 1) + [c / 2]), -b
        w = RNG.uniform(-np.pi, np.pi, 7)
        if not np.allclose(dtft(x, n0, w), a + c * np.cos(b * w)):
            probs.append("fa25_dtft: x does not have the stated DTFT")
        k = np.arange(L)
        X = np.array([np.sum(x * np.exp(-2j * np.pi * kk * (n0 + np.arange(len(x))) / L)) for kk in k])
        if not np.allclose(X, np.fft.fft(np.roll(np.concatenate([x[b:], x[:b]]), 0))):
            probs.append("fa25_dtft: DFT over the support differs from the DFT of the circular arrangement")
        claim = {"ok": a + c * np.cos(2 * np.pi * b * k / L), "drop": a + c * np.cos(2 * np.pi * k / L),
                 "halfgrid": a + c * np.cos(np.pi * b * k / L), "nosub": a + c * np.cos(b * k)}[prm["form"]]
        return bool(np.allclose(X, claim))
    if sid == "parseval_dft":
        X = np.fft.fft(np.array(prm["x"], float))
        return bool(abs(np.sum(np.abs(X) ** 2) - f(prm["V"])) < 1e-9)
    if sid == "idft_periodic":
        X = np.array(prm["X"], complex)
        N = len(X)
        n = np.arange(-2 * N, 3 * N)
        x = np.array([np.mean(X * np.exp(2j * np.pi * np.arange(N) * k / N)) for k in n])
        at = dict(zip(n.tolist(), x))
        form = prm["form"]
        if prm["kind"] == "list":
            if form == "zero":
                return bool(all(abs(at[k]) < 1e-9 for k in at if k < 0 or k > N - 1))
            if form == "periodic":
                return bool(all(abs(at[k + N] - at[k]) < 1e-9 for k in at if k + N in at))
            v = parse_num(re.search(r"\$x\[\d+\] = (.+?)\$", text).group(1))
            return bool(abs(at[N] - v) < 1e-9)
        if prm["kind"] == "fa25":
            c = X[0].real / N
            if form == "window":
                ok = all(abs(at[k] - c) < 1e-9 for k in range(N)) and all(abs(at[k]) < 1e-9 for k in (-2, -1, N, N + 1))
                return bool(ok)
            return bool(all(abs(at[k] - c) < 1e-9 for k in at))
        B = X[0].real
        if form == "delta_one":
            return bool(abs(at[0] - B) < 1e-9 and all(abs(at[k]) < 1e-9 for k in range(1, N)))
        if form == "spread":
            return bool(all(abs(at[k] - B / N) < 1e-9 for k in range(N)))
        return bool(all(abs(at[k] - (B if k == 0 else 0)) < 1e-9 for k in at))
    if sid == "dft_same":
        N = 6
        F = np.exp(-2j * np.pi * np.outer(np.arange(N), np.arange(N)) / N)
        invertible = np.linalg.matrix_rank(F) == N
        if prm["form"] == "same":
            return not invertible
        if prm["form"] == "unique":
            return bool(invertible)
        x = RNG.standard_normal(N)
        A, B = np.fft.fft(x), np.fft.fft(x, 2 * N)
        return bool(A.shape == B.shape and np.allclose(A, B))
    if sid == "dup_rule":
        N = 5
        x = RNG.standard_normal(N) + 1j * RNG.standard_normal(N)
        X = np.fft.fft(x)
        if prm["form"].startswith("dup"):
            y = np.fft.ifft(np.concatenate([X, X]))
            inter = np.zeros(2 * N, complex)
            inter[0::2] = x
            claim = inter if prm["form"] == "dup_ok" else np.tile(x, 2)
        else:
            Y = np.zeros(2 * N, complex)
            Y[0::2] = X
            y = np.fft.ifft(Y)
            claim = np.tile(x, 2) / 2 if prm["form"] == "inter_ok" else np.tile(x, 2)
        return bool(np.allclose(y, claim))
    if sid == "bin_freq":
        N, k0, th = prm["N"], prm["k0"], np.pi * f(prm["th"])
        ok = True
        for _ in range(3):
            x = RNG.standard_normal(N)
            v = dtft(x, 0, th)[0]
            v = np.conj(v) if prm["conj"] else v
            ok &= bool(abs(np.fft.fft(x)[k0] - v) < 1e-9)
        return ok
    if sid == "zp_scale":
        c, N = prm["c"], 5
        x = RNG.standard_normal(N) + 1j * RNG.standard_normal(N)
        X, Y = np.fft.fft(x), np.fft.fft(x, c * N)
        w = RNG.uniform(-np.pi, np.pi, 6)
        same_dtft = np.allclose(dtft(x, 0, w), dtft(np.concatenate([x, np.zeros((c - 1) * N)]), 0, w))
        form = prm["form"]
        if form == "ok":
            return bool(np.allclose(X, Y[c * np.arange(N)]) and same_dtft)
        if form == "inv":
            return bool(abs(X[c] - Y[1]) < 1e-9 and same_dtft)
        if form == "scaled":
            return bool(np.allclose(X, Y[c * np.arange(N)]) and np.allclose(dtft(x, 0, w), dtft(x, 0, c * w)))
        return bool(np.allclose(X, Y[:N]) and same_dtft)
    if sid == "zp_index":
        L, M, k0, m0 = prm["L"], prm["M"], prm["k0"], prm["m0"]
        ok = True
        for _ in range(3):
            x = RNG.standard_normal(L) if prm["conj"] else RNG.standard_normal(L) + 1j * RNG.standard_normal(L)
            X, Y = np.fft.fft(x), np.fft.fft(x, M)
            ok &= bool(abs(X[k0] - (np.conj(Y[m0]) if prm["conj"] else Y[m0])) < 1e-9)
        return ok
    if sid == "which_transform":
        sid2 = prm["id"]
        if sid2 == "dtft_periodic":
            x = RNG.standard_normal(7)
            w = RNG.uniform(-np.pi, np.pi, 5)
            return bool(np.allclose(dtft(x, -2, w), dtft(x, -2, w + 2 * np.pi)))
        if sid2 == "dft_periodic_k":
            x = RNG.standard_normal(6)
            n = np.arange(6)
            Xk = lambda k: np.sum(x * np.exp(-2j * np.pi * k * n / 6))       # noqa: E731
            return bool(all(abs(Xk(k) - Xk(k + 6)) < 1e-9 for k in range(-3, 9)))
        if sid2 == "dtft_from_z":
            n = np.arange(400)
            h = 0.5 ** n
            w = RNG.uniform(-np.pi, np.pi, 5)
            if not np.allclose(dtft(h, 0, w), 1 / (1 - 0.5 * np.exp(-1j * w))):
                probs.append("dtft_from_z: example fails")
        return TABLE[sid2]
    if sid == "hd_substitute":
        form = prm["form"]
        if form in ("every", "stable"):
            return TABLE[form]
        n = np.arange(600)
        if form == "fa21":                                    # y[n] = y[n-2] + x[n] - x[n-1], impulse response by recursion
            h = np.zeros(len(n))
            for k in n:
                h[k] = (h[k - 2] if k >= 2 else 0) + (1 if k == 0 else 0) - (1 if k == 1 else 0)
            if not np.allclose(h[:6], [1, -1, 1, -1, 1, -1]):
                probs.append("hd_substitute fa21: recursion does not give (-1)^n u[n]")
            return summable(h)
        a = f(prm["a"])
        with np.errstate(all="ignore"):
            h = a ** n.astype(float)
        if not summable(h):
            return False
        w = RNG.uniform(-np.pi, np.pi, 5)
        return bool(np.allclose(dtft(h, 0, w), 1 / (1 - a * np.exp(-1j * w))))
    if sid == "half_shift":
        N, form = prm["N"], prm["form"]
        x = RNG.standard_normal(N) + 1j * RNG.standard_normal(N)
        X = np.fft.fft(x)
        k = np.arange(N)
        Y = (-1.0) ** k * X if form.startswith("alt") else X[(k - N // 2) % N]
        y = np.roll(x, N // 2) if form in ("alt_shift", "kshift_shift") else (-1.0) ** k * x
        return bool(np.allclose(np.fft.ifft(Y), y))
    if sid == "phase_ramp":
        x = np.array(prm["x"], float)
        k = np.arange(8)
        s = 1 if prm["adv"] else -1
        y = np.fft.ifft(np.exp(s * 2j * np.pi * prm["m"] * k / 8) * np.fft.fft(x))
        return bool(np.allclose(y, prm["claim"]))
    if sid == "cos_bins":
        N, c = prm["N"], f(prm["c"])
        n = np.arange(N)
        X = np.abs(np.fft.fft(np.cos(2 * np.pi * c * n / N)))
        count = int(np.sum(X > 1e-6 * X.max()))
        if prm["kind"] == "analog":
            m = re.search(r"\\cos\((\d+)\\pi t\)", text)
            fs = int(re.search(r"f_s = (\d+)\$ Hz", text).group(1))
            xa = np.cos(int(m.group(1)) * np.pi * n / fs)
            X2 = np.abs(np.fft.fft(xa))
            count = int(np.sum(X2 > 1e-6 * X2.max()))
        return count == prm["m"]
    if sid == "peak_height":
        N, m, A = prm["N"], prm["m"], f(prm["A"])
        n = np.arange(N)
        x = A * (np.cos if prm["trig"] == "cos" else np.sin)(2 * np.pi * m * n / N)
        return bool(abs(np.abs(np.fft.fft(x))[m] - f(prm["V"])) < 1e-9)
    if sid == "peak_freq":
        N, fs, k0 = prm["N"], prm["fs"], prm["k0"]
        n = np.arange(N)
        hits = []
        for fr in np.arange(1, 8 * N // 2) * fs / (8 * N):        # candidate tones on a grid 8 times finer than the bins
            X = np.abs(np.fft.fft(np.cos(2 * np.pi * fr * n / fs)))
            nz = set(np.nonzero(X > 1e-6 * X.max())[0].tolist())
            if nz == {k0, N - k0}:
                hits.append(fr)
        if len(hits) != 1:
            probs.append(f"peak_freq: {len(hits)} tones fit the stated bins")
            return None
        f0 = hits[0]
        m = re.search(r"Then \$f_0 = (.+?)\$ Hz", text)
        if m:
            return bool(abs(parse_num(m.group(1)) - f0) < 1e-6)
        m = re.search(r"Then \$\\Omega_0 = (\d+)\\pi\$ rad/s", text)
        return bool(abs(int(m.group(1)) * np.pi - 2 * np.pi * f0) < 1e-6)
    if sid in ("zp_resolution", "full_lobe"):
        N, d = prm["N"], np.pi * f(prm["d"])
        if sid == "full_lobe":
            z = first_null(np.ones(N))
            return bool(d > (2 * z if prm["kind"] == "full" else z))
        c, form = prm["c"], prm["form"]
        rect = np.ones(N)
        if not d < 2 * first_null(rect):
            probs.append("zp_resolution: main lobes do not overlap to begin with")
        if form == "pad":
            padded = np.concatenate([rect, np.zeros((c - 1) * N)])
            return bool(d > 2 * first_null(padded))
        if form == "more":
            return bool(d > 2 * first_null(np.ones(c * N)))
        padded = np.concatenate([rect, np.zeros((c - 1) * N)])
        g = np.linspace(-np.pi, np.pi, 777)
        return bool(np.allclose(win_dtft(rect, g), win_dtft(padded, g)) and abs(first_null(padded) - 2 * np.pi / N) < 1e-3)
    if sid == "record_time":
        df = prm["df"]
        if prm["form"] == "count":
            dw = 2 * np.pi * df / prm["fs"]
            z = first_null(np.ones(64))                          # 2 pi / 64: main-lobe half width scales as 1/N
            thr = 2 * (z * 64) / dw                              # N with 2 * (2 pi / N) = dw
            return bool(abs(thr - prm["V"]) < 1e-2 * thr)
        if prm["form"] == "double":
            t = [4 * np.pi / (2 * np.pi * df / fs) / fs for fs in (prm["fs"], 2 * prm["fs"])]
            return bool(abs(t[1] - t[0] / 2) < 1e-9)
        fs = 8000.0
        t = 4 * np.pi / (2 * np.pi * df / fs) / fs
        return bool(abs(t - f(prm["t"])) < 1e-9)
    if sid == "fa25_plots":
        N = prm["N"]
        rect, rect2, ham2 = np.ones(N), np.ones(2 * N), hamming(2 * N)
        padded = np.concatenate([rect, np.zeros(N)])
        g = np.linspace(-np.pi, np.pi, 999)
        form = prm["form"]
        if form == "zp_same":
            return bool(np.allclose(win_dtft(rect, g), win_dtft(padded, g)))
        if form == "more_taller":
            return bool(abs(win_dtft(rect2, [0.0])[0] / win_dtft(rect, [0.0])[0] - 2) < 1e-9
                        and abs(first_null(rect2) / first_null(rect) - 0.5) < 1e-3)
        if form == "zp_taller":
            return bool(abs(win_dtft(padded, [0.0])[0] / win_dtft(rect, [0.0])[0] - 2) < 0.1)
        return bool(first_null(ham2) < first_null(rect2))
    if sid == "hamming":
        N = 64
        lower = peak_sidelobe(hamming(N)) < peak_sidelobe(np.ones(N))
        wider = first_null(hamming(N)) > first_null(np.ones(N))
        return (prm["side"] == "lower") == lower and (prm["main"] == "wider") == wider
    if sid == "window_lti":
        N = 16
        w = hamming(N)

        def sys(x):                                              # x on n = 0..N-1+pad; window fixed at 0..N-1
            out = np.zeros_like(x)
            out[:N] = w * x[:N]
            return out
        x1, x2 = RNG.standard_normal(2 * N), RNG.standard_normal(2 * N)
        lin = np.allclose(sys(2 * x1 - 3 * x2), 2 * sys(x1) - 3 * sys(x2))
        ti = np.allclose(sys(np.roll(x1, 1) * (np.arange(2 * N) > 0)), np.roll(sys(x1), 1) * (np.arange(2 * N) > 0))
        return {"linear": lin, "lti": lin and ti, "ti": ti, "notti": lin and not ti}[prm["form"]]
    if sid == "window_product":
        x, w = RNG.standard_normal(8), RNG.standard_normal(8)
        om = RNG.uniform(-np.pi, np.pi, 5)
        Y = dtft(x * w, 0, om)
        M = 64
        th = 2 * np.pi * np.arange(M) / M
        conv = np.array([np.mean(dtft(x, 0, th) * dtft(w, 0, o - th)) for o in om])
        if prm["form"] == "conv":
            return bool(np.allclose(Y, conv))
        return bool(np.allclose(Y, dtft(x, 0, om) * dtft(w, 0, om)))
    if sid == "mainlobe_width":
        N, form = prm["N"], prm["form"]
        m = re.search(r"width (?:about )?\$(.+?)\$", text)
        tw = m.group(1)
        mm = re.fullmatch(r"\\tfrac\{(\d*)\\pi\}\{(\d+)\}", tw)
        claimed_w = np.pi * (int(mm.group(1) or 1) / int(mm.group(2))) if mm else np.pi * float(tw.replace(r"\pi", "") or 1)
        if form.startswith("rect"):
            hm = re.search(r"height \$(.+?)\$", text)
            h_claim = parse_num(hm.group(1))
            rect = np.ones(N)
            return bool(abs(2 * first_null(rect) - claimed_w) < 1e-3 and abs(win_dtft(rect, [0.0])[0] - h_claim) < 1e-9)
        width = 2 * first_null(hamming(N))
        return bool(abs(width - claimed_w) < 0.2 * claimed_w)
    raise KeyError(sid)


# ----------------------------------------------------------------------------- display slips
SLIPS = [
    (re.compile(r"e\^\{\s*[+-]?\s*j\s*[-+\u2212]"), "sign right after j in an exponent (e^{j-...})"),
    (re.compile(r"(?:\+|(?<![<>])-|\u2212)\s*(?:-|\u2212)(?![-\u2212>])"), "double sign (+ -, - -)"),
    (re.compile(r"(?<![\d.\w}\]])1\s*(?:\\[,;:!]\s*)?(?:e\^\{|\\cos|\\sin|\\delta|\\pi|\\omega)"), "coefficient 1 written out"),
    (re.compile(r"\}1\s*[a-z]\b"), "factor 1 glued to a symbol"),
    (re.compile(r"(?<![\d.\w])-\s*0(?![\d.])"), "minus zero"),
    (re.compile(r"e\^\{-?j0\}"), "e^{j0}"),
    (re.compile(r"k\s*\\cdot\s*1(?![\d.])"), "factor k times 1"),
]
MATH_RX = re.compile(r"\$\$(.+?)\$\$|\$(.+?)\$", re.S)


def display_slips(sid, html):
    out = []
    for a, b in MATH_RX.findall(html):
        e = a or b
        for rx, what in SLIPS:
            m = rx.search(e)
            if m:
                out.append(f"{sid}: display slip ({what}): ...{e[max(0, m.start() - 20):m.end() + 20]}...")
    return out


def check(params, correct):
    probs = []
    items = params["statements"]
    if len(items) != 5:
        probs.append(f"{len(items)} statements instead of 5")
    if len({it["sid"] for it in items}) != len(items):
        probs.append("duplicate statement")
    answers = []
    for it in items:
        probs += display_slips(it["sid"], it["stmt_html"]) + display_slips(it["sid"], it["why_html"])
        want = recompute(it, probs)
        graded = correct[it["name"]]["html"].strip()
        answers.append(graded)
        if want is None or graded != ("True" if want else "False"):
            probs.append(f"{it['name']} ({it['sid']}, {it['param']}): graded {graded}, recomputed {want}")
        if it["ans"] != graded:
            probs.append(f"{it['name']}: answer panel says {it['ans']} but grading uses {graded}")
    if len(set(answers)) < 2:
        probs.append("all five answers equal")
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
    for nm in names[:3]:                    # one statement flipped: that one scores 0, the other four 1
        cases.append(({nm: _key(params, nm, flip[right[nm]])}, {m: (0 if m == nm else 1) for m in names}))
    cases.append(({nm: _key(params, nm, flip[right[nm]]) for nm in names}, {nm: 0 for nm in names}))    # all flipped
    cases.append(({nm: _key(params, nm, right[nm]) for nm in names}, {nm: 1 for nm in names}))         # all correct
    cases.append(({names[0]: correct[names[0]]["key"]}, {names[0]: 1}))
    cases.append(({nm: _key(params, nm, "True") for nm in names},                                     # "all True" guess
                  {nm: (1 if right[nm] == "True" else 0) for nm in names}))
    return cases
