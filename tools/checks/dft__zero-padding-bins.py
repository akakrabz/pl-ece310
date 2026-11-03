"""Independent check of questions/dft/zero-padding-bins.

A random real test signal of length N is used throughout; its DTFT is summed directly at the displayed
frequencies and compared with numpy.fft.fft (zero-padded where asked):
(a) every displayed relation X[k] = X_d(theta) / X_d^*(theta) is parsed from the option HTML and decided
    numerically; it must match the marked set;
(b) for every candidate L the padded DFT must contain both DTFT values (brute-force bin search);
(c) the smallest such L >= 0 by brute force;
(d) the bins of the padded DFT that hold X_d(w_b) and X_d^*(w_b)."""

import re
from fractions import Fraction as Fr

import numpy as np

RNG = np.random.default_rng(329)


def _dtft(x, w):
    return complex(np.sum(np.asarray(x) * np.exp(-1j * w * np.arange(len(x)))))


def _parse_pi(t):
    """LaTeX multiple of pi -> Fraction (0, \\pi, -\\pi, 2\\pi, \\tfrac{3\\pi}{8}, -\\tfrac{\\pi}{2})."""
    t = t.strip()
    sgn = -1 if t.startswith("-") else 1
    t = t.lstrip("-")
    if t == "0":
        return Fr(0)
    m = re.fullmatch(r"\\tfrac\{(\d*)\\pi\}\{(\d+)\}", t)
    if m:
        return sgn * Fr(int(m.group(1) or 1), int(m.group(2)))
    m = re.fullmatch(r"(\d*)\\pi", t)
    if m:
        return sgn * Fr(int(m.group(1) or 1))
    raise ValueError(t)


def _relations(params, x):
    """key -> (holds, largest gap over three real test signals, k, conj, theta)."""
    xs = [x] + [RNG.standard_normal(len(x)) for _ in range(2)]
    out = {}
    for o in params["rel"]:
        m = re.fullmatch(r"\$X\[(\d+)\] = X_d(\^\*)?\\left\((.+)\\right\)\$", o["html"].strip())
        k, conj, th = int(m.group(1)), bool(m.group(2)), _parse_pi(m.group(3))
        gap = 0.0
        for xx in xs:
            v = _dtft(xx, np.pi * float(th))
            v = np.conj(v) if conj else v
            gap = max(gap, abs(np.fft.fft(xx)[k] - v))
        out[o["key"]] = (gap < 1e-9, gap, k, conj, th)
    return out


def _bins(x, M, w):
    Y = np.fft.fft(x, M)
    v = _dtft(x, w)
    return [k for k in range(M) if abs(Y[k] - v) < 1e-9]


def _freqs(params):
    if params["hz"]:
        return [2 * np.pi * params["F1"] / params["fs"], 2 * np.pi * params["F2"] / params["fs"]]
    return [np.pi * float(_parse_pi(params["w1_tex"])), np.pi * float(_parse_pi(params["w2_tex"]))]


def _works(x, L, ws):
    return all(_bins(x, len(x) + L, w) for w in ws)


def _x(params):
    return RNG.standard_normal(params["N"])


def check(params, correct):
    probs = []
    N = params["N"]
    x = _x(params)
    if N not in (16, 20, 24, 32):
        probs.append(f"N = {N}")
    # (a)
    truth = _relations(params, x)
    marked = {o["key"] for o in correct["rel"]}
    for key, (ok, gap, k, conj, th) in truth.items():
        if ok != (key in marked):
            probs.append(f"(a) X[{k}] = X_d{'*' if conj else ''}({th} pi) is {ok} numerically but marked {key in marked}")
        if not ok and gap < 0.1:
            probs.append(f"(a) false relation at k = {k} differs by only {gap:.3g}")
    nt = sum(v[0] for v in truth.values())
    if len(truth) != 5 or not 2 <= nt <= 4:
        probs.append(f"(a) {len(truth)} relations, {nt} true")
    # (b)
    ws = _freqs(params)
    if params["hz"]:
        if params["F1"] >= params["fs"] / 2 or params["F2"] >= params["fs"] / 2:
            probs.append("(b) a frequency is not below f_s/2")
    marked = {o["key"] for o in correct["Lset"]}
    for o in params["Lset"]:
        L = int(re.fullmatch(r"\$L = (\d+)\$", o["html"].strip()).group(1))
        if _works(x, L, ws) != (o["key"] in marked):
            probs.append(f"(b) L = {L}: numerically {_works(x, L, ws)}, marked {o['key'] in marked}")
    if not 1 <= len(marked) <= 3:
        probs.append(f"(b) {len(marked)} correct options")
    # (c)
    best = next(L for L in range(0, 400) if _works(x, L, ws))
    if best != correct["Lmin"]:
        probs.append(f"(c) brute force L = {best}, answer {correct['Lmin']}")
    if best == 0:
        probs.append("(c) no padding needed")
    if len(_mistakes_c(params, correct)) < 2:
        probs.append("(c) fewer than two classic mistakes")
    # (d)
    M = N + best
    wb = ws[params["b"] - 1]
    kb = _bins(x, M, wb)
    v = np.conj(_dtft(x, wb))
    Y = np.fft.fft(x, M)
    kc = [k for k in range(M) if abs(Y[k] - v) < 1e-9]
    if kb != [correct["k"]] or kc != [correct["kc"]]:
        probs.append(f"(d) numpy bins {kb}, {kc}; answers {correct['k']}, {correct['kc']}")
    return probs


def _mistakes_c(params, correct):
    N, q, d, q1, q2 = params["N"], params["q"], params["d"], params["q1"], params["q2"]
    cands = [(-N) % d, q, (-N) % max(q1, q2), (-N) % min(q1, q2), (-N) % (q1 * q2)]
    out = []
    for c in cands:
        if c != correct["Lmin"] and c not in out:
            out.append(c)
    return out


def submissions(params, correct):
    Lmin, k, kc = correct["Lmin"], correct["k"], correct["kc"]
    M = params["N"] + Lmin
    cases = [
        ({"Lmin": f"+{Lmin}", "k": f" {k} ", "kc": f"0{kc}"}, {"Lmin": 1, "k": 1, "kc": 1}),
        ({"Lmin": f" {Lmin}", "k": f"+{k}", "kc": f"+{kc}"}, {"Lmin": 1, "k": 1, "kc": 1}),
        ({"Lmin": f"{Lmin}.0"}, {"Lmin": "invalid"}),
    ]
    for w in _mistakes_c(params, correct):
        cases.append(({"Lmin": str(w)}, {"Lmin": 0}))
    for w in {2 * k % M, kc, k + 1}:                       # w/pi grid; swapped; off by one
        if w != k:
            cases.append(({"k": str(w)}, {"k": 0}))
    for w in {k, M - 1 - k, -k}:                            # no conjugate; plain flip; negative index
        if w != kc:
            cases.append(({"kc": str(w)}, {"kc": 0}))
    x = _x(params)
    truth = _relations(params, x)
    T = [kk for kk, v in truth.items() if v[0]]
    F = [kk for kk, v in truth.items() if not v[0]]
    cases += [({"rel": T}, {"rel": 1}),
              ({"rel": T + F[:1]}, {"rel": (len(T) - 1) / len(T)}),
              ({"rel": F}, {"rel": 0})]
    for kk in F:
        cases.append(({"rel": [kk]}, {"rel": 0}))
    ws = _freqs(params)
    good, bad = [], []
    for o in params["Lset"]:
        L = int(re.fullmatch(r"\$L = (\d+)\$", o["html"].strip()).group(1))
        (good if _works(x, L, ws) else bad).append(o["key"])
    cases += [({"Lset": good}, {"Lset": 1}),
              ({"Lset": good + bad[:1]}, {"Lset": (len(good) - 1) / len(good)}),
              ({"Lset": bad}, {"Lset": 0}),
              ({"Lset": good + bad}, {"Lset": max(0, (len(good) - len(bad)) / len(good))})]
    return cases
