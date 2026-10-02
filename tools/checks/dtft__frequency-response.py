"""Independent check of questions/dtft/frequency-response with scipy.signal.freqz: gains at
0, pi/2, pi, stability (poles inside the unit circle), and the filter type read off the shape of
|H_d| on a dense grid (monotone -> LP/HP, single interior peak/dip -> BP/BS, constant -> allpass)."""

from fractions import Fraction

import numpy as np
from scipy import signal

TYPES = ["lowpass", "highpass", "bandpass", "bandstop", "allpass"]


def _ba(params):
    return ([float(Fraction(v)) for v in params["b"]], [float(Fraction(v)) for v in params["a"]])


def _gains(params, w):
    b, a = _ba(params)
    _, h = signal.freqz(b, a, worN=np.asarray(w, float))
    return np.abs(h)


def _shape(mag):
    tol = 1e-9 * max(1.0, mag.max())
    d = np.diff(mag)
    if mag.max() - mag.min() < tol:
        return "allpass"
    if np.all(d <= tol):
        return "lowpass"
    if np.all(d >= -tol):
        return "highpass"
    i = int(np.argmax(mag))
    if 0 < i < len(mag) - 1 and np.all(d[:i] >= -tol) and np.all(d[i:] <= tol):
        return "bandpass"
    i = int(np.argmin(mag))
    if 0 < i < len(mag) - 1 and np.all(d[:i] <= tol) and np.all(d[i:] >= -tol):
        return "bandstop"
    return "other"


def _num(v):
    f = Fraction(v).limit_denominator(1000)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def check(params, correct):
    probs = []
    b, a = _ba(params)
    if len(a) > 1 and np.any(np.abs(np.roots(a)) >= 1 - 1e-12):
        probs.append(f"unstable: poles {np.roots(a)}")
    g0, g2, gpi = _gains(params, [0, np.pi / 2, np.pi])
    for name, val in (("H0", g0), ("Hpi", gpi), ("B", params["A"] * g2)):
        if abs(val - correct[name]) > 1e-9 * max(1, abs(val)):
            probs.append(f"{name}: correct {correct[name]} != freqz {val}")
    shape = _shape(_gains(params, np.linspace(0, np.pi, 4001)))
    if shape != params["kind"]:
        probs.append(f"type {params['kind']!r} but the magnitude response looks {shape!r}")
    corr = [c["text"] for c in params["type_choices"] if c["correct"]]
    if corr != [shape]:
        probs.append(f"MC marks {corr} correct, freqz shape is {shape}")
    for v in (correct["H0"], correct["Hpi"], correct["B"]):
        f = Fraction(v).limit_denominator(1000)
        if abs(float(f) - v) > 1e-12 or f.denominator > 30:
            probs.append(f"answer {v} is not a nice fraction")
    if abs(correct["B"] - round(correct["B"])) > 1e-12:
        probs.append(f"B = {correct['B']} is not an integer")
    return probs


def _differs(val, ref):
    """A wrong value must be far outside rtol = 1e-3 / atol = 1e-6 to be a meaningful probe."""
    return abs(val - ref) > 1e-2 * max(abs(ref), 1e-3)


def _key(params, text):
    return next(o["key"] for o in params["type"] if o["html"].strip() == text)


def submissions(params, correct):
    b, a = _ba(params)
    H0, Hpi, B, A = correct["H0"], correct["Hpi"], correct["B"], params["A"]
    kind = params["kind"]
    wrong_kind = {"lowpass": "highpass", "highpass": "lowpass", "bandpass": "bandstop",
                  "bandstop": "bandpass", "allpass": "lowpass"}[kind]
    cases = [
        ({"H0": _num(H0)}, {"H0": 1}),                                   # fraction format
        ({"Hpi": f"{Hpi:.12f}"}, {"Hpi": 1}),                            # long decimal
        ({"B": _num(B)}, {"B": 1}),
        ({"type": _key(params, kind)}, {"type": 1}),
        ({"type": _key(params, wrong_kind)}, {"type": 0}),               # mixed up the type
        ({"H0": "1/0"}, {"H0": "invalid"}),
    ]
    mist = []
    s0 = sum(b) / sum(a)                                                 # H(1) before taking |.|
    if s0 < 0:
        mist.append(("H0", s0))                                          # forgot the absolute value
    if len(a) > 1:                                                       # LCCDE sign slip: a_k -> -a_k
        a_flip = [a[0]] + [-v for v in a[1:]]
        if np.all(np.abs(np.roots(a_flip)) < 1):
            _, h = signal.freqz(b, a_flip, worN=[0.0, np.pi])
            mist += [("H0", abs(h[0])), ("Hpi", abs(h[1]))]
        mist.append(("H0", abs(sum(b))))                                 # forgot to divide by A(1)
    mist.append(("Hpi", abs(sum(b) / sum(a))))                           # forgot (-1)^k
    mist.append(("B", float(A)))                                         # forgot to scale by the gain
    mist.append(("B", A * H0))                                           # used the gain at w = 0
    for name in ("H0", "Hpi", "B"):                                      # 4 significant digits suffice
        v = correct[name]
        if abs(v - round(v)) > 1e-9:
            cases.append(({name: f"{v:.4g}"}, {name: 1}))
    seen = set()
    for name, val in mist:
        if _differs(val, correct[name]) and (name, round(val, 9)) not in seen:
            seen.add((name, round(val, 9)))
            cases.append(({name: _num(val)}, {name: 0}))
    return cases
