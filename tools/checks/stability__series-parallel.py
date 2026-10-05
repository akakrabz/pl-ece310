"""Independent check of questions/stability/series-parallel: impulse responses of the subsystems
(exact, C delta[n] + A p^n u[n]) combined by convolution (series) or addition (parallel); the overall
H(z) answer is compared with the truncated z-transform of that h, and the stability answers with
the decay of h."""

from fractions import Fraction

import numpy as np
from checklib import eval_z

L = 70


def _h(spec):
    K, q, p = (Fraction(v) for v in spec)
    b0, b1 = K, -K * q
    C = -b1 / p
    A = b0 - C
    return [C * (n == 0) + A * p ** n for n in range(L)]


def _overall(params):
    s = params["spec"]
    h1, h2 = _h(s["s1"]), _h(s["s2"])
    if s["conn"] == "series":
        h = [sum(h1[k] * h2[n - k] for k in range(n + 1)) for n in range(L)]
    else:
        h = [x + y for x, y in zip(h1, h2)]
    return h1, h2, h


def _stable(h):
    """Absolutely summable <=> the tail has died out (all poles here have |p| <= 3/4 or >= 1)."""
    return max(abs(float(v)) for v in h[L - 15:]) < 1e-5


def _choice_text(params, correct, name):
    return correct[name]["html"].strip()


def check(params, correct):
    probs = []
    h1, h2, h = _overall(params)
    z0 = 3.6 * np.exp(0.7j)
    want = sum(complex(float(v)) * z0 ** (-n) for n, v in enumerate(h))
    got = eval_z(correct["H"], z0)
    if abs(got - want) > 1e-8 * max(1, abs(want)):
        probs.append(f"H(z0) = {got} but sum h[n] z0^-n = {want}")
    for name, seq in (("st1", h1), ("st2", h2), ("st", h)):
        truth = "Yes" if _stable(seq) else "No"
        if _choice_text(params, correct, name) != truth:
            probs.append(f"{name}: marked {correct[name]['html']!r}, decay test says {truth}")
    return probs


def _key(params, name, text):
    for o in params[name]:
        if o["html"].strip() == text:
            return o["key"]
    raise KeyError(text)


def _sym(spec, form="inv"):
    K, q, p = (Fraction(v) for v in spec)
    if form == "pos":       # multiply by z: K (z - q)/(z - p)
        return f"({K})*(z - ({q}))/(z - ({p}))"
    if form == "badpole":   # sign error in the pole
        return f"({K})*(1 - ({q})*z^(-1))/(1 + ({p})*z^(-1))"
    return f"({K})*(1 - ({q})*z^(-1))/(1 - ({p})*z^(-1))"


def submissions(params, correct):
    s = params["spec"]
    H1, H2 = _sym(s["s1"]), _sym(s["s2"])
    op, wrong_op = ("*", "+") if s["conn"] == "series" else ("+", "*")
    cases = [
        ({"H": f"({H1}){op}({H2})"}, {"H": 1}),                                           # unsimplified combination
        ({"H": f"({_sym(s['s1'], 'pos')}){op}({_sym(s['s2'], 'pos')})"}, {"H": 1}),          # positive powers of z
        ({"H": f"({H1}){wrong_op}({H2})"}, {"H": 0}),                                      # series <-> parallel mix-up
        ({"H": f"({_sym(s['s1'], 'badpole')}){op}({H2})"}, {"H": 0}),                      # sign error in a pole
        ({"H": "1/(1-0.5*z^(-1))"}, {"H": "invalid"}),
    ]
    st = correct["st"]["html"].strip()
    cases.append(({"st": _key(params, "st", "No" if st == "Yes" else "Yes")}, {"st": 0}))
    # "unstable subsystem => unstable cascade" is the classic wrong call when a pole cancels
    h1, h2, _ = _overall(params)
    both = "Yes" if (_stable(h1) and _stable(h2)) else "No"
    cases.append(({"st": _key(params, "st", both)}, {"st": 1 if both == st else 0}))
    return cases
