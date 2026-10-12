"""Independent check of questions/dft/dft-and-dtft.

(a), (b): a random complex test signal; its DTFT is summed directly at the displayed frequencies and
compared with numpy.fft.fft of the zero-padded signal: (a) the matching bin must be the answer,
(b) the smallest N (brute force, N = 4..200) whose DFT contains both DTFT values must be the answer.
(c) every displayed relation X[a] = Y[b] / X[a] = Y^*[b] is decided with numpy.fft.fft(x) and
numpy.fft.fft(x, M) and must match the marked set. (d) the marked option is the true statement."""

import html
import re

import numpy as np

TRUE_PREFIXES = ("The DTFT is unchanged", "It gives more closely spaced samples", r"$X_d(0) = \sum_n x[n]$ and the energy",
                 "The appended zeros only add terms", "The width of a sinusoid's main lobe")
RNG = np.random.default_rng(310)
XT = RNG.normal(size=4) + 1j * RNG.normal(size=4)          # generic test signal (all values distinct)


def _dtft(x, w):
    return complex(np.sum(x * np.exp(-1j * w * np.arange(len(x)))))


def _bins_with(x, N, w):
    X = np.fft.fft(x, N)
    v = _dtft(x, w)
    return [k for k in range(N) if abs(X[k] - v) < 1e-9 * max(1, abs(v))]


def _frac(pq):
    return pq[0] / pq[1]


def differs(wrong, right):
    return abs(float(wrong) - float(right)) > 2 * (1e-6 + 1e-3 * abs(float(right)))


def check(params, correct):
    probs = []
    # (a)
    N = params["Na"]
    w = np.pi * _frac(params["wa"])
    hits = _bins_with(XT, N, w)
    if hits != [correct["ka"]]:
        probs.append(f"(a) bins holding X_d({w:.4f}) for N={N}: {hits}, answer {correct['ka']}")
    # (b)
    w1, w2 = (np.pi * _frac(t) for t in params["wb"])
    best = next((n for n in range(4, 201) if _bins_with(XT, n, w1) and _bins_with(XT, n, w2)), None)
    if best != correct["Nb"]:
        probs.append(f"(b) brute force gives N = {best}, answer {correct['Nb']}")
    if len(_mistakes_b(params, correct)) < 2:
        probs.append("(b) fewer than two distinct classic mistakes")
    # (c)
    truth = _relations(params)
    marked = {o["key"] for o in correct["rel"]}
    for key, (ok, gap, a, b, conj) in truth.items():
        if ok != (key in marked):
            probs.append(f"(c) X[{a}] = Y{'*' if conj else ''}[{b}] is {ok} by numpy but marked {key in marked}")
        if not ok and gap < 0.5:
            probs.append(f"(c) false relation X[{a}] vs Y[{b}] differs by only {gap:.3g}")
    n_true = sum(ok for ok, *_ in truth.values())
    if len(truth) != 6 or not 1 <= n_true <= 4:
        probs.append(f"(c) {len(truth)} relations, {n_true} true")
    L_, M_ = len(params["xc"]), params["M"]

    def shortcut(a, b, conj):                       # "plain iff b = a M/L, conjugate iff b != a M/L"
        on_grid = (a * M_) % L_ == 0 and (a * M_ // L_) % M_ == b
        return (not on_grid) if conj else on_grid

    if all(shortcut(a, b, conj) == ok for ok, gap, a, b, conj in truth.values()):
        probs.append("(c) the index-only shortcut (plain iff b = aM/L, conjugate iff b != aM/L) answers every relation")
    L, M = len(params["xc"]), params["M"]
    if params["Mz"] != M - L or M <= L or (2 * M) % L:
        probs.append(f"(c) L = {L}, M = {M}, Mz = {params['Mz']} inconsistent")
    if M % L and not any((a * M) % L for ok, gap, a, b, conj in truth.values()):
        probs.append("(c) M is not a multiple of L but no relation uses a bin without a partner")
    # (d)
    marked = html.unescape(correct["zp"]["html"]).strip()
    if not marked.startswith(TRUE_PREFIXES):
        probs.append(f"(d) marked option is not a true statement: {marked[:60]}")
    if sum(html.unescape(o["html"]).strip().startswith(TRUE_PREFIXES) for o in params["zp"]) != 1:
        probs.append("(d) not exactly one true statement among the options")
    return probs


def _relations(params):
    """key -> (holds, |difference|, a, b, conj) for each displayed relation."""
    x = np.array(params["xc"], float)
    X, Y = np.fft.fft(x), np.fft.fft(x, params["M"])
    out = {}
    for o in params["rel"]:
        m = re.fullmatch(r"\$X\[(\d+)\] = Y(\^\*)?\[(\d+)\]\$", o["html"].strip())
        a, conj, b = int(m.group(1)), bool(m.group(2)), int(m.group(3))
        v = np.conj(Y[b]) if conj else Y[b]
        gap = abs(X[a] - v)
        out[o["key"]] = (gap < 1e-9, gap, a, b, conj)
    return out


def _mistakes_b(params, correct):
    from math import gcd
    q1, q2 = params["qb"]
    d1, d2 = params["wb"][0][1], params["wb"][1][1]           # denominators of w/pi
    cands = [d1 * d2 // gcd(d1, d2), q1 * q2, max(q1, q2)]
    out = []
    for c in cands:
        if c != correct["Nb"] and c not in out:
            out.append(c)
    return out


def _num(v):
    v = float(v)
    return str(int(round(v))) if abs(v - round(v)) < 1e-12 else f"{v:.6g}"


def submissions(params, correct):
    ka, Na, Nb = correct["ka"], params["Na"], correct["Nb"]
    cases = [
        ({"ka": f"+{ka}"}, {"ka": 1}),
        ({"ka": f" {ka} "}, {"ka": 1}),
        ({"Nb": f"+{Nb}"}, {"Nb": 1}),
        ({"Nb": f"0{Nb}"}, {"Nb": 1}),
        ({"ka": f"{ka}.0"}, {"ka": "invalid"}),
    ]
    for wrong in {2 * ka, Na - ka, ka - Na}:                   # pi k/N grid; mirror bin; negative index
        if wrong != ka:
            cases.append(({"ka": str(wrong)}, {"ka": 0}))
    for wrong in _mistakes_b(params, correct):
        cases.append(({"Nb": str(wrong)}, {"Nb": 0}))
    truth = _relations(params)
    T = [k for k, v in truth.items() if v[0]]
    F = [k for k, v in truth.items() if not v[0]]
    cases += [
        ({"rel": T}, {"rel": 1}),
        ({"rel": T + F[:1]}, {"rel": (len(T) - 1) / len(T)}),        # net-correct: one wrong tick costs one right tick
        ({"rel": T[:1]}, {"rel": 1 / len(T)}),
        ({"rel": F}, {"rel": 0}),
        ({"rel": T[:1] + F[:2]}, {"rel": 0}),
    ]
    for k, (ok, gap, a, b, conj) in truth.items():                  # each trap alone scores 0
        if not ok:
            cases.append(({"rel": [k]}, {"rel": 0}))
    for o in params["zp"]:
        cases.append(({"zp": o["key"]}, {"zp": 1 if o["key"] == correct["zp"]["key"] else 0}))
    return cases
