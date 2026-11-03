"""Independent check of questions/dft/leakage-check.

(a) every displayed statement is parsed back from its TEXT (amplitudes, cos/sin, frequencies p pi/q,
    N; or the analog cosine, f_s and duration), the signal is rebuilt and numpy.fft.fft counts the
    nonzero bins (|X[k]| > 1e-9 N); the statement is true iff the count is the claimed one (two, four).
    The marked set must match, with 1-3 true statements.
(b) numpy.fft.fft of the record: the nonzero bins must be exactly {k1, k2} and |X[k1]| = Xk.
(c) the off-bin record (l22: the same cosine, N2 samples; fa21: the L samples zero-padded with
    numpy.fft.fft(x, N)): kappa = N w0 / (2 pi) (l22: also N2 / period found by brute force), the
    tallest bin in 0..N/2 must be kmax, at least 5 % above the runner-up, and the record must leak
    (more than two nonzero bins)."""

import html
import re

import numpy as np

TOL = 1e-9


def nonzero(X, N):
    return [k for k in range(len(X)) if abs(X[k]) > TOL * N]


def parse_stmt(text):
    t = html.unescape(text)
    claim = 4 if "exactly four" in t else 2
    m = re.search(r"\\cos\((\d*)\\pi t\)", t)
    if m:
        W0 = int(m.group(1)) if m.group(1) else 1
        fs = int(re.search(r"f_s = (\d+)\$ Hz", t).group(1))
        D = float(re.search(r"for \$([\d.]+)\$ s", t).group(1))
        N = int(round(fs * D))
        if int(re.search(r"\}_\{n=0\}\^\{ (\d+) \}", t).group(1)) != N - 1:
            raise ValueError("analog statement: N inconsistent")
        n = np.arange(N)
        return np.cos(np.pi * W0 * n / fs), N, claim
    m = re.search(r"\$(\d+)\$-point DFT", t)
    if m:
        N = int(m.group(1))
        if f"0 \\le n \\le {N - 1}$" not in t:
            raise ValueError("discrete statement: range inconsistent")
    else:
        N = int(re.search(r"\}_\{n=0\}\^\{ (\d+) \}", t).group(1)) + 1
    n = np.arange(N)
    x = np.zeros(N)
    terms = re.findall(r"(\d*)\\(cos|sin)\\!\\left\(\\tfrac\{(\d*)\\pi\}\{(\d+)\}n\\right\)", t)
    for a, fn, p, q in terms:
        a = int(a) if a else 1
        p = int(p) if p else 1
        x += a * (np.cos if fn == "cos" else np.sin)(np.pi * p / int(q) * n)
    if len(terms) != claim // 2:
        raise ValueError(f"parsed {len(terms)} terms for claim {claim}")
    return x, N, claim


def stmt_truth(params):
    out = {}
    for o in params["stmt"]:
        x, N, claim = parse_stmt(o["html"])
        cnt = len(nonzero(np.fft.fft(x), N))
        out[o["key"]] = (cnt == claim, cnt, claim)
    return out


def records(params):
    """(x_on, N_on), (x_off, N_off, Npts), w0."""
    A = params["A"]
    w0 = np.pi * params["t"][0] / params["t"][1]
    if params["l22"]:
        N1, N2 = params["N1"], params["N2"]
        return (A * np.cos(w0 * np.arange(N1)), N1), (A * np.cos(w0 * np.arange(N2)), N2, N2), w0
    m = re.search(r"\\cos\((\d*)\\pi t\)", params["xc_tex"])
    W0 = int(m.group(1)) if m.group(1) else 1
    fs, L = params["fs"], params["L"]
    if abs(W0 * np.pi / fs - w0) > 1e-12:
        raise ValueError("w0 inconsistent with x_c(t)")
    x = A * np.cos(np.pi * W0 * np.arange(L) / fs)
    return (x, L), (x, L, params["Np"]), w0


def check(params, correct):
    probs = []
    truth = stmt_truth(params)
    marked = {o["key"] for o in correct["stmt"]}
    for key, (ok, cnt, claim) in truth.items():
        if ok != (key in marked):
            probs.append(f"(a) statement {key}: {cnt} nonzero bins, claim {claim}, marked {key in marked}")
    nt = sum(ok for ok, _, _ in truth.values())
    if not 1 <= nt <= 3 or len(truth) != 4:
        probs.append(f"(a) {nt} true of {len(truth)}")
    (x1, N1), (x2, N2, Np), w0 = records(params)
    X1 = np.fft.fft(x1)
    nz = nonzero(X1, N1)
    if nz != [correct["k1"], correct["k2"]]:
        probs.append(f"(b) nonzero bins {nz} vs {correct['k1']}, {correct['k2']}")
    if abs(abs(X1[correct["k1"]]) - correct["Xk"]) > 1e-6 * N1:
        probs.append(f"(b) |X[k1]| {abs(X1[correct['k1']]):.6f} vs {correct['Xk']}")
    X2 = np.abs(np.fft.fft(x2, Np))
    if len(nonzero(X2, N2)) <= 2:
        probs.append("(c) the record does not leak")
    kap = w0 * Np / (2 * np.pi)
    if abs(kap - correct["kap"]) > 1e-9 * Np:
        probs.append(f"(c) kappa {correct['kap']} vs {kap}")
    if params["l22"]:
        per = next(P for P in range(1, 200) if np.allclose(np.cos(w0 * (np.arange(50) + P)), np.cos(w0 * np.arange(50))))
        if abs(N2 / per - correct["kap"]) > 1e-9:
            probs.append(f"(c) periods {N2}/{per} != {correct['kap']}")
    half = X2[: Np // 2 + 1]
    order = np.argsort(-half)
    if order[0] != correct["kmax"]:
        probs.append(f"(c) argmax {order[0]} vs kmax {correct['kmax']}")
    if half[order[0]] < 1.05 * half[order[1]]:
        probs.append(f"(c) tallest bin only {half[order[0]] / half[order[1]]:.3f} x the runner-up")
    frac = kap - np.floor(kap)
    if abs(frac - 0.5) < 0.12 or min(frac, 1 - frac) < 0.08 or 2e-3 * kap > min(frac, 1 - frac):
        probs.append(f"(c) fractional part {frac:.3f}")
    return probs


def _num(v):
    v = float(v)
    return str(int(round(v))) if abs(v - round(v)) < 1e-12 else f"{v:.6g}"


def submissions(params, correct):
    truth = stmt_truth(params)
    T = [k for k, v in truth.items() if v[0]]
    F = [k for k, v in truth.items() if not v[0]]
    cases = [({"stmt": T}, {"stmt": 1}), ({"stmt": F}, {"stmt": 0})]
    if F:
        cases.append(({"stmt": T + F[:1]}, {"stmt": (len(T) - 1) / len(T)}))
    for k in F:
        cases.append(({"stmt": [k]}, {"stmt": 0}))
    if len(T) > 1:
        cases.append(({"stmt": T[:1]}, {"stmt": 1 / len(T)}))
    A = params["A"]
    k1, k2, Xk, kap, kmax = (correct[n] for n in ("k1", "k2", "Xk", "kap", "kmax"))
    Non = params["N1"] if params["l22"] else params["L"]
    Npts = params["N2"] if params["l22"] else params["Np"]
    cases += [
        ({"k1": f"+{k1}", "k2": f" {k2} "}, {"k1": 1, "k2": 1}),
        ({"Xk": f"{A * Non}/2"}, {"Xk": 1}),
        ({"Xk": f"{Xk:.6g}"}, {"Xk": 1}),
        ({"k2": str(-k1)}, {"k2": 0}),                                 # -k0 instead of N - k0
        ({"k1": str(2 * k1)}, {"k1": 0}),                              # w0 N / pi
        ({"Xk": _num(A * Non)}, {"Xk": 0}),                            # AN
        ({"Xk": _num(A / 2)}, {"Xk": 0}),
        ({"kap": f"{kap:.4g}"}, {"kap": 1}),
        ({"kap": _num(2 * kap)}, {"kap": 0}),                          # w0 N / pi
        ({"kap": _num(round(kap))}, {"kap": 0}),                       # rounded
        ({"kmax": str(kmax + (1 if kap > kmax else -1))}, {"kmax": 0}),  # the other neighbour
        ({"kmax": str(Npts - kmax)}, {"kmax": 0}),                     # the mirror (outside 0..N/2)
        ({"k1": f"{k1}.0"}, {"k1": "invalid"}),
    ]
    if A != 1:
        cases.append(({"Xk": _num(Non / 2)}, {"Xk": 0}))               # forgot the amplitude
    if params["fa21"]:
        cases.append(({"kap": _num(k1)}, {"kap": 0}))                  # the unpadded bin
    else:
        cases.append(({"kap": _num(Npts / kap)}, {"kap": 0}))          # period / ... inverted
    return cases
