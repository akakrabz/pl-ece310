"""Independent check of questions/lccde/recursion: scipy.signal.lfilter on an impulse for h[0..3];
FIR/IIR from an exact (Fraction) recursion run far past the input terms — h finite <=> FIR, whatever
the equation looks like (recursive FIR variants have feedback that cancels)."""

from fractions import Fraction

import numpy as np
from scipy import signal


def _ba(params):
    b = [float(Fraction(v)) for v in params["b"]]
    a = [float(Fraction(v)) for v in params["a"]]
    return b, a


def _imp(b, a, N):
    x = np.zeros(N)
    x[0] = 1.0
    return signal.lfilter(b, a, x)


def _imp_exact(params, N):
    b = [Fraction(v) for v in params["b"]]
    a = [Fraction(v) for v in params["a"]]
    h = []
    for n in range(N):
        acc = b[n] if n < len(b) else Fraction(0)
        acc -= sum(a[k] * h[n - k] for k in range(1, len(a)) if n - k >= 0)
        h.append(acc / a[0])
    return h


def _plain(v: float) -> str:
    f = Fraction(v).limit_denominator(1000)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def check(params, correct):
    probs = []
    b, a = _ba(params)
    h = _imp(b, a, 4)
    for n in range(4):
        if not np.isclose(correct[f"h{n}"], h[n], rtol=1e-12, atol=1e-12):
            probs.append(f"h[{n}] = {correct[f'h{n}']} but lfilter gives {h[n]}")
    he = _imp_exact(params, 60)
    fir_true = all(v == 0 for v in he[len(b):])
    want = "FIR" if fir_true else "IIR"
    if want not in correct["firiir"]["html"]:
        probs.append(f"marked answer {correct['firiir']['html']!r}, expected {want}")
    kind = params["kind"]
    if (kind == "iir") == fir_true:
        probs.append(f"kind {kind} but the exact impulse response is {'finite' if fir_true else 'infinite'}")
    rb = np.roots(b) if len(b) > 1 else np.array([])
    ra = np.roots(a) if len(a) > 1 else np.array([])
    shared = any(np.min(np.abs(ra - r)) < 1e-9 for r in rb) if len(ra) and len(rb) else False
    if kind == "iir" and shared:
        probs.append("IIR variant whose B and A share a root (cancellation)")
    if kind == "fir_rec" and not shared:
        probs.append("recursive FIR variant without a pole-zero cancellation")
    if kind == "fir_plain" and len(a) != 1:
        probs.append("plain FIR variant with feedback terms")
    for v in h:
        f = Fraction(float(v)).limit_denominator(64)
        if abs(float(f) - v) > 1e-12 or f.denominator > 32 or abs(f.numerator) > 99:
            probs.append(f"h value not exam-nice: {v}")
    return probs


def _key(params, text):
    for o in params["firiir"]:
        if o["html"].strip() == text:
            return o["key"]
    raise KeyError(text)


def submissions(params, correct):
    b, a = _ba(params)
    h = _imp(b, a, 6)
    names = [f"h{n}" for n in range(4)]
    cases = []
    # equivalent forms: exact fractions, full decimals, and 4-significant-digit decimals all score 1
    cases.append(({k: _plain(h[i]) for i, k in enumerate(names)}, {k: 1 for k in names}))
    cases.append(({k: repr(float(h[i])) for i, k in enumerate(names)}, {k: 1 for k in names}))
    cases.append(({k: f"{float(h[i]):.4g}" for i, k in enumerate(names)}, {k: 1 for k in names}))
    # classic mistake 1: sign error from moving the feedback terms (a_k -> -a_k)
    if len(a) > 1:
        a_bad = [a[0]] + [-v for v in a[1:]]
        hb = _imp(b, a_bad, 6)
        ov = {k: _plain(hb[i]) for i, k in enumerate(names)}
        cases.append((ov, {k: (1 if np.isclose(hb[i], h[i], rtol=1e-3, atol=1e-6) else 0) for i, k in enumerate(names)}))
    # classic mistake 2: the step response instead of the impulse response
    s = np.cumsum(h)
    ov = {k: _plain(s[i]) for i, k in enumerate(names)}
    cases.append((ov, {k: (1 if np.isclose(s[i], h[i], rtol=1e-3, atol=1e-6) else 0) for i, k in enumerate(names)}))
    # classic mistake 3: the wrong FIR/IIR call (for recursive FIR variants: "it has feedback, so IIR")
    fir = params["kind"] != "iir"
    cases.append(({"firiir": _key(params, "IIR" if fir else "FIR")}, {"firiir": 0}))
    cases.append(({"firiir": _key(params, "FIR" if fir else "IIR")}, {"firiir": 1}))
    cases.append(({"h2": "h[2]"}, {"h2": "invalid"}))
    return cases
