"""Independent check of questions/stability/bibo-from-impulse-response.

Evaluates h[n] sample by sample from the term list in params["spec"] (floats, no closed forms) and
forms partial sums of |h[n]| over |n| <= 150 and |n| <= 300: a stable h has converged (the two
agree), an unstable one keeps growing (harmonic: +ln 2; bounded periodic: linear; exponential)."""

import math
from fractions import Fraction

import numpy as np


def _h(spec, n):
    v = 0.0
    for t in spec:
        k = t["t"]
        if k == "rexp" and n >= t["k"]:
            v += float(Fraction(t["A"])) * float(Fraction(t["a"])) ** n
        elif k == "lexp" and n <= t["m"]:
            v += float(Fraction(t["B"])) * float(Fraction(t["b"])) ** n
        elif k == "abs":
            v += float(Fraction(t["A"])) * float(Fraction(t["a"])) ** abs(n)
        elif k == "nexp" and n >= 0:
            v += float(Fraction(t["A"])) * n * float(Fraction(t["a"])) ** n
        elif k == "harm" and n >= t["s"]:
            v += float(Fraction(t["A"])) / (n + t["c"])
        elif k == "trig" and n >= 0:
            w = math.pi * t["w"][0] / t["w"][1]
            f = math.cos if t["fn"] == "cos" else math.sin
            v += float(Fraction(t["A"])) * float(Fraction(t["r"])) ** n * round(f(w * n), 12)
        elif k == "fin":
            i = n - t["start"]
            if 0 <= i < len(t["values"]):
                v += float(Fraction(t["values"][i]))
        elif k == "delta" and n == t["k"]:
            v += float(Fraction(t["A"]))
    return v


def _partial(spec, N):
    tot = 0.0
    for n in range(-N, N + 1):
        try:
            tot += abs(_h(spec, n))
        except OverflowError:
            return math.inf
    return tot


def _signed(spec, N=300):
    return sum(_h(spec, n) for n in range(-N, N + 1))


def check(params, correct):
    probs = []
    spec = params["spec"]
    s1, s2 = _partial(spec, 150), _partial(spec, 300)
    stable = math.isfinite(s2) and abs(s2 - s1) < 1e-7 * max(1.0, s2)
    if stable != params["is_stable"]:
        probs.append(f"partial sums {s1}, {s2} say stable={stable}, server says {params['is_stable']}")
    want = "Yes" if stable else "No"
    if correct["stable"]["html"].strip() != want:
        probs.append(f"MC answer {correct['stable']['html']!r}, expected {want}")
    if stable:
        if not isinstance(correct["sumabs"], float) or not np.isclose(correct["sumabs"], s2, rtol=1e-9, atol=1e-12):
            probs.append(f"sum |h| = {correct['sumabs']} but partial sum gives {s2}")
        if abs(float(Fraction(params["signed_str"])) - _signed(spec)) > 1e-9 * max(1, abs(_signed(spec))):
            probs.append("signed-sum probe value is wrong")
    elif correct["sumabs"] != "":
        probs.append("unstable variant must have a blank correct answer for (b)")
    return probs


def _key(params, text):
    return next(o["key"] for o in params["stable"] if o["html"].strip() == text)


def submissions(params, correct):
    stable = params["is_stable"]
    yes, no = _key(params, "Yes"), _key(params, "No")
    cases = []
    if stable:
        S = Fraction(params["S_str"])
        frac = f"{S.numerator}/{S.denominator}" if S.denominator != 1 else str(S.numerator)
        cases += [
            ({"stable": yes, "sumabs": frac}, {"stable": 1, "sumabs": 1}),                       # exact fraction
            ({"sumabs": f"{float(S):.10f}"}, {"sumabs": 1}),                                     # decimal
            ({"sumabs": f"{2 * S.numerator}/{2 * S.denominator}"}, {"sumabs": 1}),               # unreduced fraction
            ({"stable": no, "sumabs": ""}, {"stable": 0, "sumabs": 0}),                          # "unstable", box left empty
        ]
        if S.denominator != 1:                                                                # 4 significant digits suffice
            cases.append(({"sumabs": f"{float(S):.4g}"}, {"sumabs": 1}))
        signed = Fraction(params["signed_str"])
        if signed != S:                                                                       # forgot the absolute values
            cases.append(({"sumabs": str(float(signed))}, {"sumabs": 0}))
        for w in params["wrong_strs"]:                                                        # wrong start index / n = 0 twice ...
            if Fraction(w) != S:
                cases.append(({"sumabs": str(Fraction(w))}, {"sumabs": 0}))
    else:
        cases += [
            ({"stable": no, "sumabs": ""}, {"stable": 1, "sumabs": 1}),                         # empty box
            ({"sumabs": "   "}, {"sumabs": 1}),                                                  # whitespace counts as empty
            ({"stable": yes, "sumabs": ""}, {"stable": 0, "sumabs": 0}),                         # "stable" + empty box: no credit for (b)
            ({"stable": yes, "sumabs": "2"}, {"stable": 0, "sumabs": 0}),                        # "bounded, so stable" + a number
            ({"stable": no, "sumabs": "0"}, {"stable": 1, "sumabs": 0}),                         # contradiction: a number for an unstable h
        ]
    cases.append(({"sumabs": "abc"}, {"sumabs": "invalid"}))
    return cases
