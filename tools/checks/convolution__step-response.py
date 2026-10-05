"""Independent check of questions/convolution/step-response.

fir: s = h * u via np.convolve with a long step (running sum), read off the window and the tail.
exp: h = s[n] - s[n-1] via np.diff of the sampled step response, against h[0] and the closed form."""

from fractions import Fraction as F

import numpy as np
from checklib import S, N


def _s_exp(p, n):
    a, c, d = F(p["a"]), F(p["c"]), F(p["d"])
    return float(c + d * a ** n) if n >= 0 else 0.0


def check(params, correct):
    probs = []
    p = params
    if p["template"] == "fir":
        h, nh, L = np.array(p["h"], dtype=float), p["nh"], p["L"]
        lo = nh - 3
        u = np.ones(40)                                    # u[n] for n = lo .. lo+39 (starts before h)
        hh = np.zeros(40)
        hh[nh - lo: nh - lo + L] = h
        s = np.convolve(hh, u)[:40]                        # s[lo + i], exact for i < 40
        sval = lambda n: s[n - lo]                         # noqa: E731
        win = np.array(correct["win"]["_value"][0])
        want = np.array([sval(n) for n in range(p["w0"], p["w_last"] + 1)])
        if win.shape != want.shape or not np.allclose(win, want):
            probs.append(f"window {win} != {want}")
        tail = [sval(n) for n in range(p["n_last"], p["n_last"] + 10)]
        if not np.allclose(tail, correct["steady"]):
            probs.append(f"steady {correct['steady']} != tail {tail[:3]}")
        if np.isclose(sval(p["n_last"] - 1), correct["steady"]):
            probs.append("steady state reached before the stated index (ambiguous window)")
        if p["w0"] != nh - 1 or sval(p["w0"]) != 0:
            probs.append("window does not start one sample before h")
    else:
        ns = np.arange(-3, 30)
        s = np.array([_s_exp(p, int(n)) for n in ns])
        h = np.diff(s, prepend=0.0)                         # h[n] = s[n] - s[n-1]  (s[-4] = 0)
        hval = lambda n: h[n + 3]                          # noqa: E731
        if any(abs(hval(n)) > 1e-12 for n in range(-3, 0)):
            probs.append("h nonzero for n < 0")
        if abs(hval(0) - correct["h0"]) > 1e-9:
            probs.append(f"h[0] = {correct['h0']} vs {hval(0)}")
        e = S(correct["hn"]["_value"])
        for n in range(1, 20):
            if abs(float(e.subs(N, n)) - hval(n)) > 1e-9 * max(1.0, abs(hval(n))):
                probs.append(f"closed form at n={n}: {float(e.subs(N, n))} vs {hval(n)}")
                break
        if "." in correct["hn"]["_value"]:
            probs.append("closed form contains a float")
    return probs


def _pb(q):
    q = F(q)
    s = str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return s if (q > 0 and q.denominator == 1) else f"({s})"


def _pn(q):
    q = F(q)
    return f"({q.numerator})" if q.denominator == 1 else f"({q.numerator}/{q.denominator})"


def _row(v, sep=", "):
    return "[" + sep.join(str(int(round(t))) for t in v) + "]"


def submissions(params, correct):
    p = params
    cases = []
    if p["template"] == "fir":
        win = [int(round(t)) for t in correct["win"]["_value"][0]]
        steady = int(round(correct["steady"]))
        h = p["h"]
        cases += [({"win": _row(win, " ")}, {"win": 1}),                     # MATLAB style
                  ({"win": "[" + _row(win) + "]"}, {"win": 1}),              # Python style
                  ({"steady": f"{steady}.0"}, {"steady": 1})]
        hwin = [0] + h[:-1]                                                 # entered h instead of s
        if hwin != win:
            cases.append(({"win": _row(hwin)}, {"win": 0}))
        early = [int(t) for t in np.cumsum([0] + h)[1:len(win) + 1]]       # running sum one sample early
        if early != win:
            cases.append(({"win": _row(early)}, {"win": 0}))
        if _far(h[-1], steady):
            cases.append(({"steady": str(h[-1])}, {"steady": 0}))           # last sample of h, not the sum
        cases.append(({"win": _row(win + [steady])}, {"win": 0}))           # one entry too many
    else:
        a, c, d = F(p["a"]), F(p["c"]), F(p["d"])
        K = d * (a - 1) / a
        form = correct["hn"]["_value"]
        cases += [({"hn": f"{_pn(d)}*({_pb(a)}-1)*{_pb(a)}^(n-1)"}, {"hn": 1}),   # unsimplified difference
                  ({"hn": f"{_pn(K)}*{_pb(1 / a)}^(-n)"}, {"hn": 1}),              # reciprocal base
                  ({"hn": f"({_pn(c)} + {_pn(d)}*{_pb(a)}^n) - ({_pn(c)} + {_pn(d)}*{_pb(a)}^(n-1))"}, {"hn": 1}),
                  ({"h0": p["h0_frac"]}, {"h0": 1}),
                  ({"hn": f"{_pn(d * (a - 1))}*{_pb(a)}^n"}, {"hn": 0}),          # forward difference s[n+1]-s[n]
                  ({"hn": f"{_pn(d)}*{_pb(a)}^n"}, {"hn": 0}),                     # kept only the exponential of s
                  ({"hn": f"-({form})"}, {"hn": 0})]                               # sign error
        if _far(c, c + d):                                                         # dropped the exponential at n = 0
            cases.append(({"h0": f"{c.numerator}/{c.denominator}"}, {"h0": 0}))
        h0 = c + d
        if h0.denominator != 1:                                                    # 4 significant digits (rtol 1e-3)
            cases.append(({"h0": f"{float(h0):.4g}"}, {"h0": 1}))
        if _far(K, h0):                                                            # used the n >= 1 formula at n = 0
            cases.append(({"h0": f"{K.numerator}/{K.denominator}"}, {"h0": 0}))
    return cases


def _far(wrong, right):
    """A mistake probe is only meaningful if it is clearly outside the rtol = 1e-3 tolerance."""
    return abs(float(wrong) - float(right)) > 1e-2 * max(abs(float(right)), 1e-4)
