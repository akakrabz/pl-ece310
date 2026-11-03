"""Independent check of questions/sampling/max-bandwidth-filter.

The chain is rebuilt in numpy in physical units (rad/s and rad, not the multiples of pi the server
uses): X_a is the triangle of height 1 on |Omega| < 2 pi B, X_d(w) = (1/T) sum_{|k|<=8} X_a((w - 2 pi k)/T),
Y_d = H_d X_d with the ideal LPF/HPF (cutoff included in the passband), Y_a(Omega) = T Y_d(Omega T) for
|Omega| < pi/T. Checked: every value, the edge of Y_a found by scanning, the Yes/No answer by comparing
Y_a with the ideal analog filter applied to X_a on a grid, and the largest B (or T) of part (c) by moving the bandlimit
(or the period) by bisection, testing on a frequency grid (passband edges included) whether a
neighbouring copy is nonzero inside the passband of H_d (aliasing in Y_d) or anywhere in [-pi, pi]
(aliasing in X_d)."""

import math
from fractions import Fraction

import numpy as np

K = np.arange(-8, 9)


def _setup(params):
    fs = float(Fraction(params["fs"]))
    B = float(Fraction(params["B"]))
    c = float(Fraction(params["c"])) * math.pi
    return fs, B, c, params["lp"]


def xa(W, B):
    return np.maximum(0.0, 1 - np.abs(W) / (2 * math.pi * B))


def xd(w, B, fs, skip_central=False):
    w = np.atleast_1d(np.asarray(w, float))
    ks = K[K != 0] if skip_central else K
    return fs * xa((w[:, None] - 2 * math.pi * ks[None, :]) * fs, B).sum(axis=1)


def hd(w, c, lp):
    wp = np.abs((np.asarray(w, float) + math.pi) % (2 * math.pi) - math.pi)
    return (wp <= c + 1e-12).astype(float) if lp else (wp >= c - 1e-12).astype(float)


def ya(W, B, fs, c, lp):
    W = np.asarray(W, float)
    T = 1 / fs
    return np.where(np.abs(W) < math.pi * fs, T * hd(W * T, c, lp) * xd(W * T, B, fs), 0.0)


def _wgrid(c):
    return np.concatenate([np.linspace(-math.pi, math.pi, 2001), [c, -c, math.pi, -math.pi]])


def aliased_yd(B, fs, c, lp):
    """Is a neighbouring copy nonzero somewhere in the passband of H_d (edges included in the grid)?"""
    w = _wgrid(c)
    return bool(np.max(hd(w, c, lp) * xd(w, B, fs, skip_central=True)) > 1e-9 * fs)


def aliased_xd(B, fs, c=0.0):
    w = _wgrid(c)
    return bool(np.max(xd(w, B, fs, skip_central=True)) > 1e-9 * fs)


def boundary(pred, lo, hi):
    """pred(lo) False, pred(hi) True, monotone: the switching point by bisection."""
    assert not pred(lo) and pred(hi), (lo, hi)
    for _ in range(50):
        mid = (lo + hi) / 2
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def _close(a, b, scale=1.0):
    return abs(a - b) <= 1e-7 * max(1.0, abs(b), scale)


def check(params, correct):
    probs = []
    fs, B, c, lp = _setup(params)
    u1 = float(Fraction(params["u1"])) * math.pi
    if not _close(correct["Xd0"], xd(0.0, B, fs)[0], fs):
        probs.append(f"X_d(0) {correct['Xd0']} vs {xd(0.0, B, fs)[0]}")
    w = np.linspace(0, 2 * math.pi, 200001)
    cen = fs * xa(w * fs, B)
    edge = w[cen > 1e-12].max() / math.pi
    if abs(edge - correct["wB"]) > 2e-4:
        probs.append(f"band edge {correct['wB']} vs scan {edge}")
    if not _close(correct["Xd1"], xd(u1, B, fs)[0], fs):
        probs.append(f"X_d(u1) {correct['Xd1']} vs {xd(u1, B, fs)[0]}")
    for d in (-1e-7, 1e-7):
        if abs(xd(u1 + d, B, fs)[0] - xd(u1, B, fs)[0]) > 1e-3 * fs:
            probs.append("X_d jumps at u1")
    yd_c = hd(c, c, lp) * xd(c, B, fs)[0]
    if not _close(correct["Ydc"], yd_c, fs):
        probs.append(f"Y_d(w_c) {correct['Ydc']} vs {yd_c}")
    W = np.linspace(0, math.pi * fs * 0.999999, 40001)
    Y = ya(W, B, fs, c, lp)
    nz = W[Y > 1e-12]
    e = (nz.max() if lp else nz.min()) / math.pi
    if abs(e - correct["Omc"]) > 1e-4 * fs + 1e-12:
        probs.append(f"Y_a edge {correct['Omc']} pi vs scan {e} pi")
    yac = ya(np.array([correct["Omc"] * math.pi]), B, fs, c, lp)[0]
    if not _close(correct["Yac"], yac):
        probs.append(f"Y_a(Omega_c) {correct['Yac']} vs {yac}")
    if not _close(correct["Yac"] * fs, correct["Ydc"], fs):
        probs.append("Y_a(Omega_c) != T Y_d(w_c)")
    # LTI on this input? compare with the ideal analog filter of cutoff w_c / T = c fs rad/s (no band limit)
    Wall = np.linspace(-6 * math.pi * B, 6 * math.pi * B, 60001)
    pas = (np.abs(Wall) <= c * fs) if lp else (np.abs(Wall) >= c * fs)
    ideal = pas * xa(Wall, B)
    far = (np.abs(np.abs(Wall) - c * fs) > 1e-6 * fs) & (np.abs(np.abs(Wall) - math.pi * fs) > 1e-6 * fs)
    lti = bool(np.max(np.abs(ya(Wall, B, fs, c, lp) - ideal)[far]) < 1e-9)
    marked = correct["lti"]["html"].strip()
    if marked != ("Yes" if lti else "No"):
        probs.append(f"LTI answer {marked}, simulation says {lti}")
    # (c)
    if params["kind"] == "fa19" and params["Omax_tex"] != f"{int(2 * B)}\\pi":
        probs.append(f"displayed Omega_max {params['Omax_tex']} != {2 * B} pi")
    if params["kind"] in ("fa25", "hpf"):
        by = boundary(lambda x: aliased_yd(x, fs, c, lp), 1e-3 * fs, 1.2 * fs)
        bx = boundary(lambda x: aliased_xd(x, fs), 1e-3 * fs, 1.2 * fs)
        if abs(by - correct["c_y"]) > 1e-6 * fs:
            probs.append(f"(c) largest B for Y_d: bisection {by}, answer {correct['c_y']}")
        if abs(bx - correct["c_x"]) > 1e-6 * fs:
            probs.append(f"(c) largest B for X_d: bisection {bx}, answer {correct['c_x']}")
        for f in np.linspace(0.02, 0.98, 25):                    # monotone: clean below, aliased above
            if aliased_yd(f * correct["c_y"], fs, c, lp) or not aliased_yd(correct["c_y"] * (2 - f), fs, c, lp):
                probs.append("(c) aliasing in Y_d is not monotone around the answer")
                break
        if float(params["B_tex"]) != B:
            probs.append(f"displayed B {params['B_tex']} != {B}")
    else:
        TN = 1 / (2 * B)
        ty = boundary(lambda T: aliased_yd(B, 1 / T, c, lp), 0.01 * TN, 2.5 * TN)
        tx = boundary(lambda T: aliased_xd(B, 1 / T), 0.01 * TN, 2.5 * TN)
        if abs(ty - correct["c_y"]) > 1e-6 * TN:
            probs.append(f"(c) largest T: bisection {ty}, answer {correct['c_y']}")
        if abs(tx - correct["c_x"]) > 1e-6 * TN:
            probs.append(f"(c) Nyquist period: bisection {tx}, answer {correct['c_x']}")
    b = 2 * B / fs
    if not (float(Fraction(params["c"])) < b < 2) or abs(b - (2 - float(Fraction(params["c"])))) < 1e-12:
        probs.append(f"b = {b} outside the allowed range")
    return probs


def _s(v):
    q = Fraction(v).limit_denominator(100000)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def _differs(u, v, atol=1e-6):
    return abs(float(u) - float(v)) > 2 * (atol + 1e-3 * abs(float(v)))


def _twice(v):
    q = Fraction(v).limit_denominator(100000)
    return f"{2 * q.numerator}/{2 * q.denominator}"


def _key(params, text):
    return next(o["key"] for o in params["lti"] if o["html"].strip() == text)


def submissions(params, correct):
    fs, B, c, lp = _setup(params)
    names = ["Xd0", "wB", "Xd1", "Ydc", "Omc", "Yac", "c_y", "c_x"]
    yes = correct["lti"]["html"].strip() == "Yes"
    cases = [
        ({n: _s(correct[n]) for n in names}, {n: 1 for n in names}),
        ({n: f"{correct[n]:.5g}" for n in names}, {n: 1 for n in names}),
        ({"Xd0": f"{2 * int(correct['Xd0'])}/2", "Ydc": _twice(correct["Ydc"])}, {"Xd0": 1, "Ydc": 1}),
        ({"lti": _key(params, "Yes")}, {"lti": 1 if yes else 0}),
        ({"lti": _key(params, "No")}, {"lti": 0 if yes else 1}),
        ({"Xd0": "1"}, {"Xd0": 0}),                                               # forgot 1/T
        ({"wB": _s(Fraction(correct["wB"]) / 2)}, {"wB": 0}),                     # B T instead of 2 pi B T
        ({"wB": _s(2 * B * fs)}, {"wB": 0}),                                      # Omega_B / T
        ({"Omc": _s(Fraction(params["c"]))}, {"Omc": 0}),                         # forgot the 1/T
        ({"Omc": _s(Fraction(params["c"]) * Fraction(params["fs"]) / 2)}, {"Omc": 0}),   # hertz
        ({"Yac": "1/0"}, {"Yac": "invalid"}),
    ]
    if _differs(correct["Ydc"], correct["Yac"]):
        cases.append(({"Yac": _s(correct["Ydc"])}, {"Yac": 0}))                  # forgot the D/A gain T
        cases.append(({"Ydc": _s(correct["Yac"])}, {"Ydc": 0}))                  # forgot 1/T
    if _differs(1.0, correct["Yac"]):
        cases.append(({"Yac": "1"}, {"Yac": 0}))
    for name, u in (("Xd1", float(Fraction(params["u1"]))), ("Ydc", float(Fraction(params["c"])))):
        cu = fs * max(0.0, 1 - u / (2 * B / fs))                                  # the central copy alone
        if _differs(cu, correct[name]):
            cases.append(({name: repr(cu)}, {name: 0}))                           # copies not added
    if params["kind"] == "fa25":
        cases += [({"c_y": _s(fs / 2)}, {"c_y": 0}),                             # Nyquist answered for Y_d
                  ({"c_y": _s(Fraction(params["c"]) * Fraction(params["fs"]) / 2)}, {"c_y": 0}),   # analog cutoff in Hz
                  ({"c_y": _s(2 * correct["c_y"])}, {"c_y": 0}),
                  ({"c_x": _s(fs)}, {"c_x": 0})]                                  # forgot the 2
    elif params["kind"] == "hpf":
        c_ = Fraction(params["c"])
        cases += [({"c_y": _s((2 - c_) * Fraction(params["fs"]) / 2)}, {"c_y": 0}),   # the low-pass trick
                  ({"c_x": _s(fs)}, {"c_x": 0})]
    else:
        cases += [({"c_y": f"1/{int(2 * B)}"}, {"c_y": 0}),                      # Nyquist period
                  ({"c_x": f"1/{int(B)}"}, {"c_x": 0}),                          # forgot the 2
                  ({"c_y": _s(2 * Fraction(correct["c_y"]).limit_denominator(100000))}, {"c_y": 0}),
                  ({"c_y": f"{correct['c_y']:.4g}", "c_x": f"{correct['c_x']:.4g}"}, {"c_y": 1, "c_x": 1})]
    return cases
