"""Independent check of questions/dtft/dtft-at-frequency: X_d(w0) by a direct numpy sum of
x[n] e^{-j w0 n} (true indices), |X_d| = np.abs, principal phase = np.angle (compared on the
circle, so pi and -pi are the same answer); classic mistakes are probed against the grading:
e^{+jwn}, indices counted from the first listed sample, |X|^2, |Re| + |Im|, arctan without the
quadrant correction, an angle off by 2 pi, degrees."""

import math
import re

import numpy as np


def _X(x, start, w):
    n = np.arange(start, start + len(x))
    return complex(np.sum(np.asarray(x, float) * np.exp(-1j * w * n)))


def _ws(params):
    return [a / b * math.pi for a, b in params["ts"]]


def _wrap(a):
    return (a + math.pi) % (2 * math.pi) - math.pi


def check(params, correct):
    probs = []
    x, start = params["x"], params["start"]
    ws = _ws(params)
    if len(ws) != (3 if params["has3"] else 2):
        probs.append("number of frequencies does not match has3")
    allowed = {0.0, 1.0, 0.5, -0.5, 1 / 3, 2 / 3}
    if any(min(abs(w / math.pi - a) for a in allowed) > 1e-12 for w in ws) or len(set(ws)) != len(ws):
        probs.append(f"frequencies {params['ts']} not from the list / repeated")
    Xs = [_X(x, start, w) for w in ws]
    for i, X in enumerate(Xs, start=1):
        for name, v in ((f"re{i}", X.real), (f"im{i}", X.imag)):
            if abs(v - correct[name]) > 1e-9 * max(1.0, abs(v)):
                probs.append(f"{name}: answer {correct[name]} != numerical {v}")
    if all(abs(X.imag) < 1e-9 for X in Xs):
        probs.append("no asked frequency has Im X_d != 0 (sign error e^{+jwn} undetected)")
    Xm = Xs[params["k"]]
    if abs(abs(Xm) - correct["mag"]) > 1e-9 * max(1.0, abs(Xm)):
        probs.append(f"|X| = {abs(Xm)}, answer {correct['mag']}")
    if abs(Xm) < 1 - 1e-9:
        probs.append(f"|X_d| = {abs(Xm)} < 1 at the phase frequency")
    if abs(Xm.imag) < 1e-9 and Xm.real > 0:
        probs.append("phase frequency has X_d > 0 (phase 0 is trivial)")
    ph = correct["ph"]
    if abs(_wrap(np.angle(Xm) - ph)) > 1e-9:
        probs.append(f"phase {np.angle(Xm)} != answer {ph}")
    if not (-math.pi < ph <= math.pi + 1e-12):
        probs.append(f"phase answer {ph} not principal")
    shown = re.findall(r"\\approx (-?[0-9.]+)", params["phase_line_tex"])
    if not shown or abs(float(shown[-1]) - ph) > 1e-4:
        probs.append(f"solution phase {shown} != {ph}")
    if params["mag_exact_is_dec"] is False:
        if abs(float(params["mag_dec"]) - correct["mag"]) > 1e-4:
            probs.append(f"solution magnitude {params['mag_dec']} != {correct['mag']}")
    d = {start + i: v for i, v in enumerate(x)}
    if all(d.get(n, 0) == d.get(-n, 0) for n in d):
        probs.append("sequence is even about n = 0")
    if not (3 <= len(x) <= 5) or x[0] == 0 or x[-1] == 0:
        probs.append(f"sequence {x} violates the design")
    for s in (params["phase_line_tex"], params["Xm_tex"]):
        if "+ -" in s or "- -" in s or "j\\,1 " in s:
            probs.append(f"display slip in {s!r}")
    return probs


def _band(got, want, rtol=1e-3, atol=1e-6):
    """1 if clearly accepted, 0 if clearly rejected, None if too close to the tolerance edge."""
    tol = atol + rtol * abs(want)
    dd = abs(got - want)
    if dd <= 0.4 * tol:
        return 1
    if dd > 2 * tol:
        return 0
    return None


def _ph_expect(got, ans):
    tol = 1e-3 * abs(ans) + 1e-3
    dd = abs(_wrap(got - ans))
    if dd <= 0.4 * tol and abs(got) <= math.pi + 0.4 * tol:
        return 1
    if dd > 2 * tol or abs(got) > math.pi + 2 * tol:
        return 0
    return None


def submissions(params, correct):
    x, start = params["x"], params["start"]
    ws = _ws(params)
    cases = [({"re1": "abc"}, {"re1": "invalid"}), ({"ph": "pi/2"}, {"ph": "invalid"})]
    for i, w in enumerate(ws, start=1):
        X = _X(x, start, w)
        cases.append(({f"re{i}": f"{X.real:.6g}", f"im{i}": f"{X.imag:.4g}"}, {f"re{i}": 1, f"im{i}": 1}))
        probes = [X.conjugate(), _X(x, 0, w)]                  # e^{+jwn}; indices from the first sample
        for pv in probes:
            exp = {f"re{i}": _band(pv.real, X.real), f"im{i}": _band(pv.imag, X.imag)}
            if None not in exp.values() and 0 in exp.values():
                cases.append(({f"re{i}": f"{pv.real:.6g}", f"im{i}": f"{pv.imag:.6g}"}, exp))
    # magnitude
    Xm = _X(x, start, ws[params["k"]])
    mag = abs(Xm)
    cases.append(({"mag": f"{mag:.4g}"}, {"mag": 1}))
    for wrong in (mag * mag, abs(Xm.real) + abs(Xm.imag)):
        if _band(wrong, mag) == 0:
            cases.append(({"mag": f"{wrong:.6g}"}, {"mag": 0}))
    # phase
    ph = correct["ph"]
    cases.append(({"ph": f"{ph:.4g}"}, {"ph": 1}))
    cases.append(({"ph": f"{ph:.6f}"}, {"ph": 1}))
    if abs(abs(ph) - math.pi) < 1e-9:                          # the other spelling of +-pi
        cases.append(({"ph": "-3.1416"}, {"ph": 1}))
        cases.append(({"ph": "3.142"}, {"ph": 1}))
    other = ph + 2 * math.pi if ph < 0 or abs(ph - math.pi) < 1e-9 else ph - 2 * math.pi   # not principal
    cases.append(({"ph": f"{other:.6g}"}, {"ph": 0}))
    wrongs = [math.degrees(ph), ph / math.pi]
    if Xm.real < -1e-9:
        wrongs.append(math.atan(Xm.imag / Xm.real))             # no quadrant correction
    for wv in wrongs:
        e = _ph_expect(wv, ph)
        if e == 0:
            cases.append(({"ph": f"{wv:.6g}"}, {"ph": 0}))
    return cases
