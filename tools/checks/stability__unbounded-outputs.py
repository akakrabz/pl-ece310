"""Independent check of questions/stability/unbounded-outputs.

Builds b, a of H(z) from this checker's own reading of each displayed denominator factor (not from
the server's pole list), synthesises every input sample by sample, runs scipy.signal.lfilter for
N = 8000 samples and applies a growth test: max|y| over the whole record versus the first quarter
(resonance grows linearly, ratio ~ 4; a bounded output saturates early, ratio ~ 1)."""

import cmath
import html
import math
import re
from fractions import Fraction

import numpy as np
from scipy.signal import lfilter

N = 8000
FACTORS = {  # coefficients of z^0, z^-1, z^-2 of the displayed unit-circle factor
    "1": [1, -1], "-1": [1, 1], "j": [1, -1j], "±j": [1, 0, 1],
    "pi4": [1, -cmath.exp(1j * math.pi / 4)], "±pi4": [1, -math.sqrt(2), 1],
    "pi3": [1, -cmath.exp(1j * math.pi / 3)], "±pi3": [1, -1, 1],
    "2pi3": [1, -cmath.exp(2j * math.pi / 3)], "±2pi3": [1, 1, 1],
    "rad23": [1, -cmath.exp(2j / 3)], "±1": [1, 0, -1],
}
DISPLAY = {  # what the factor must look like in H_tex
    "1": r"\left(1 - z^{-1}\right)", "-1": r"\left(1 + z^{-1}\right)", "±j": r"\left(1 + z^{-2}\right)",
    "±pi3": r"\left(1 - z^{-1} + z^{-2}\right)", "±2pi3": r"\left(1 + z^{-1} + z^{-2}\right)", "±1": r"\left(1 - z^{-2}\right)",
    "rad23": r"e^{j2/3}z^{-1}", "2pi3": r"e^{j2\pi/3}z^{-1}", "pi3": r"e^{j\pi/3}z^{-1}", "pi4": r"e^{j\pi/4}z^{-1}",
    "±pi4": r"\sqrt{2}\,z^{-1} + z^{-2}", "j": r"1 - j\,z^{-1}",
}


def _theta(ang):
    kind, num, den = ang
    return math.pi * num / den if kind == "pi" else num / den


def _x(spec):
    n = np.arange(N)
    k = spec["kind"]
    if k == "exp":
        return np.exp(1j * _theta(spec["ang"]) * n) * (n >= spec.get("k", 0))
    if k == "combo":
        return np.exp(1j * _theta(spec["ang"]) * n) - float(Fraction(spec["r"])) ** n
    if k == "cos":
        return np.cos(_theta(spec["ang"]) * n)
    if k == "sin":
        return np.sin(_theta(spec["ang"]) * n)
    if k == "geo":
        return float(Fraction(spec["a"])) ** n
    if k == "nexp":
        return n * float(Fraction(spec["a"])) ** n
    if k == "fir":
        x = np.zeros(N)
        x[: len(spec["c"])] = spec["c"]
        return x
    raise ValueError(k)


def _ba(params):
    a = np.array(FACTORS[params["group"]], dtype=complex)
    if params["inside"] is not None:
        a = np.convolve(a, [1, -float(Fraction(params["inside"]))])
    b = np.array([1.0]) if params["zero"] is None else np.array([1, -float(Fraction(params["zero"]))])
    return b, a


def _growth(params, spec):
    b, a = _ba(params)
    y = np.abs(lfilter(b, a, _x(spec)))
    return y.max() / max(y[: N // 4].max(), 1e-300)


def _norm(s):
    return re.sub(r"\s+", "", html.unescape(s or ""))


def check(params, correct):
    probs = []
    if DISPLAY[params["group"]] not in params["H_tex"]:
        probs.append(f"H_tex does not show the factor of group {params['group']}")
    want = set()
    for inp in params["inputs"]:
        r = _growth(params, inp["spec"])
        unb = r > 2.5
        if 1.3 < r <= 2.5:
            probs.append(f"ambiguous growth ratio {r:.2f} for {inp['tex']}")
        if unb != inp["unbounded"]:
            probs.append(f"{inp['tex']}: lfilter growth ratio {r:.2f} but server says unbounded={inp['unbounded']}")
        if unb:
            want.add(_norm(f"$x[n] = {inp['tex']}$"))
    got = {_norm(o["html"]) for o in correct["unb"]}
    if got != want:
        probs.append(f"correct option set {sorted(got)} != lfilter set {sorted(want)}")
    n_opt = len(params["unb"])
    if not (5 <= n_opt <= 6) or not (1 <= len(want) <= min(3, n_opt - 2)):
        probs.append(f"{n_opt} options with {len(want)} unbounded")
    return probs


def submissions(params, correct):
    """net-correct: score = max(0, (#correct ticked - #wrong ticked) / #correct); a blank submission is invalid."""
    opts = params["unb"]
    ck = [o["key"] for o in correct["unb"]]
    c = len(ck)
    wrong = [o["key"] for o in opts if o["key"] not in ck]
    cases = []
    # equivalent forms: another order of the same keys; a single key as a plain string
    cases.append(({"unb": list(reversed(ck))}, {"unb": 1}))
    cases.append((({"unb": ck[0]} if c == 1 else {"unb": ck[1:] + ck[:1]}), {"unb": 1}))
    # classic mistakes, with the exact net-correct score
    if c >= 2:
        cases.append(({"unb": ck[:-1]}, {"unb": (c - 1) / c}))                  # missed one resonant input
    else:
        cases.append(({"unb": wrong[:1]}, {"unb": 0}))                          # ticked a bounded input instead
    trap = None
    texts = {_norm(o["html"]): o["key"] for o in opts}
    for inp in params["inputs"]:
        spec = inp["spec"]
        if not inp["unbounded"] and spec["kind"] in ("exp", "cos", "sin") and (
                spec["ang"][0] == "rad" or params["group"] == "rad23"):
            trap = texts[_norm(f"$x[n] = {inp['tex']}$")]
    extra = trap if trap is not None else wrong[0]
    cases.append(({"unb": ck + [extra]}, {"unb": (c - 1) / c}))                  # fell for 2/3 rad vs 2*pi/3 (or one bounded input)
    cases.append(({"unb": ck + wrong}, {"unb": max(0.0, (c - len(wrong)) / c)}))  # ticked everything
    cases.append(({"unb": wrong}, {"unb": 0}))                                   # every decision wrong
    cases.append(({"unb": []}, {"unb": "invalid"}))                              # blank is not a valid answer
    cases.append(({"unb": ["z"]}, {"unb": "invalid"}))
    return cases
