"""Independent check of questions/freqresp/noncausal-fir.

The system is read back from the displayed LaTeX (difference equation, impulse response or closed
form) and must match params["spec"]; H_d(w) is then evaluated with numpy (DTFT sum over the taps, or
the closed form). Checked: |H_d(w0)|; the principal phases with np.angle (never +-pi, H_d != 0), with
at least one point where the real amplitude is negative and one whose angle needs a 2 pi wrap; for
the zero-phase type the real values H_d(w1), H_d(w2) (one negative); the number of +-pi jumps from
the phase itself: on a dense grid (never on a zero) a +-pi jump is a step whose wrapped size is
about pi, a 2 pi wrap wraps to about 0; the group delay from a centred difference of the phase.
Probes: decimals, fractions; the forgotten 2 of 2cos/2j sin; the sign kept in |H_d|; the missing
+-pi of a negative amplitude; the missing pi/2 of a factor j; the wrong sign of the ramp; the
unreduced angle; |H_d| instead of a negative H_d; zeros counted as jumps where A only touches 0;
2 pi wraps counted as jumps; only (0, pi) counted; -tau."""

import cmath
import math
import re
from fractions import Fraction

import numpy as np

PH_CANDS = [k / 12 for k in (2, 3, 4, 6, 8, 9, 10)]


def _f(s):
    return float(Fraction(s))


def H(params, w):
    w = np.atleast_1d(np.asarray(w, float))
    spec = params["spec"]
    if "taps" in spec:
        n = np.array([t[0] for t in spec["taps"]], float)
        v = np.array([t[1] for t in spec["taps"]], float)
        return np.exp(-1j * np.outer(w, n)) @ v
    c = spec["closed"]
    return c["g"] * np.exp(-1j * c["M"] * w) * np.sin(c["K"] * w)


def _parse_terms(tex, sym):
    """[(coef, shift)] from 'x[n+3] - 2x[n-1]' style sums (sym = 'x' or '\\delta')."""
    pat = re.compile(r"([+-]?)\s*(\d*)\s*" + re.escape(sym) + r"\[n(?:([+-])(\d+))?\]")
    out = []
    for sg, c, s2, k in pat.findall(tex):
        coef = (-1 if sg == "-" else 1) * (int(c) if c else 1)
        shift = 0 if not k else (int(k) if s2 == "+" else -int(k))
        out.append((coef, shift))
    return out


def displayed(params):
    """The system as displayed: {'taps': {n: v}} or {'closed': (g, M, K)}."""
    tex = params["sys_tex"]
    if tex.startswith("y[n] ="):
        return {"taps": {-s: c for c, s in _parse_terms(tex[6:], "x")}}
    if tex.startswith("h[n] ="):
        return {"taps": {-s: c for c, s in _parse_terms(tex[6:], r"\delta")}}
    m = re.fullmatch(r"H_d\(\\omega\) = (\d*)e\^\{-j(\d*)\\omega\}\\sin ?(\d*)\\omega", tex)
    assert m, tex
    return {"closed": (int(m.group(1) or 1), int(m.group(2) or 1), int(m.group(3) or 1))}


def centre_and_amp(params, w):
    """(m, J, A(w)): H = j^J e^{-jmw} A, A real."""
    spec = params["spec"]
    if "taps" in spec:
        ns = [t[0] for t in spec["taps"]]
        m = (min(ns) + max(ns)) / 2
    else:
        m = spec["closed"]["M"]
    v = H(params, w) * np.exp(1j * m * np.asarray(w, float))
    if np.max(np.abs(v.imag)) > np.max(np.abs(v.real)):
        return m, 1, (v / 1j)
    return m, 0, v


def count_jumps(params, n=20000):
    w = np.pi * (-1 + (2 * np.arange(n) + 1) / n)
    ang = np.angle(H(params, w))
    d = np.angle(np.exp(1j * np.diff(ang)))
    jumps = int(np.sum(np.abs(d) > np.pi / 2))
    wraps = int(np.sum((np.abs(np.diff(ang)) > np.pi * 1.5) & (np.abs(d) < np.pi / 2)))
    return jumps, wraps


def zeros_of_mag(params, n=20000):
    """Local minima of |H| that are ~0 on (-pi, pi) (crossings and touches)."""
    w = np.pi * (-1 + (2 * np.arange(n) + 1) / n)
    mag = np.abs(H(params, w))
    return [w[i] for i in range(1, n - 1) if mag[i] <= mag[i - 1] and mag[i] <= mag[i + 1] and mag[i] < 1e-2]


def _ang(z):
    a = float(np.angle(z) / np.pi)
    return 1.0 if a <= -1 + 1e-12 else a


def check(params, correct):
    probs = []
    shown = displayed(params)
    spec = params["spec"]
    if "taps" in spec:
        want = {t[0]: t[1] for t in spec["taps"]}
        if shown.get("taps") != want:
            probs.append(f"displayed system {shown} != spec {want}")
    elif shown.get("closed") != (spec["closed"]["g"], spec["closed"]["M"], spec["closed"]["K"]):
        probs.append(f"displayed closed form {shown} != spec {spec}")
    # (a)
    m0 = float(np.abs(H(params, math.pi * _f(params["f0"])))[0])
    if m0 < 1e-9 or abs(m0 - correct["mag"]) > 1e-9 * max(1, m0):
        probs.append(f"|H({params['f0']} pi)| = {m0}, answer {correct['mag']}")
    # (b)
    if params["is_ramp"]:
        negs = wraps = 0
        for i in (1, 2):
            w = math.pi * _f(params[f"f{i}"])
            z = H(params, w)[0]
            if abs(z) < 1e-9:
                probs.append(f"H_d = 0 at f{i}")
                continue
            a = _ang(z)
            if abs(abs(a) - 1) < 1e-9:
                probs.append(f"phase at f{i} is +-pi")
            if abs(a - correct[f"ph{i}"]) > 1e-9:
                probs.append(f"phase at {params[f'f{i}']} pi = {a}, answer {correct[f'ph{i}']}")
            m, J, A = centre_and_amp(params, [w])
            neg = A[0].real < 0
            negs += neg
            raw = -m * w / math.pi + J / 2 + (1 if neg else 0)
            wraps += not (-1 < raw <= 1 + 1e-12)
        grid = [s * t for t in PH_CANDS for s in (1, -1)]
        m, J, A = centre_and_amp(params, np.pi * np.array(grid))
        if np.any(A.real < -1e-9) and not negs:
            probs.append("no phase point with a negative amplitude although A < 0 occurs")
        if not wraps:
            probs.append("no phase point needs a 2 pi wrap")
        if abs(correct["ph1"] - correct["ph2"]) < 1e-9:
            probs.append("the two phase answers coincide")
    else:
        vals = []
        for i in (1, 2):
            z = H(params, math.pi * _f(params[f"f{i}"]))[0]
            if abs(z.imag) > 1e-9 or abs(z.real - correct[f"hv{i}"]) > 1e-9 or abs(z) < 1e-9:
                probs.append(f"H_d({params[f'f{i}']} pi) = {z}, answer {correct[f'hv{i}']}")
            vals.append(z.real)
        if not any(v < 0 for v in vals):
            probs.append("no negative value in (b)")
    # (c)
    nj, _ = count_jumps(params)
    if nj != correct["njumps"]:
        probs.append(f"{nj} pi jumps on the grid, answer {correct['njumps']}")
    # (d)
    for w0 in (0.31, 1.07, 2.39, -0.77, -1.9):
        if abs(H(params, w0)[0]) < 1e-2:
            continue
        eps = 1e-6
        tau = -np.angle(H(params, w0 + eps)[0] / H(params, w0 - eps)[0]) / (2 * eps)
        if abs(tau - correct["tau"]) > 1e-5:
            probs.append(f"group delay at {w0} is {tau}, answer {correct['tau']}")
    # display slips in generated LaTeX
    texts = [params.get(k, "") for k in ("sys_tex", "H_tex", "A_tex", "def_tex", "pair_tex", "h_tex", "A0_sub_tex")]
    texts += [r.get("steps", "") + r.get("A_line", "") for r in params.get("ph_rows", [])]
    texts += [r.get("sub", "") for r in params.get("val_rows", [])]
    slips = [r"(?<![\d.{])1\\omega", r"\+ -", r"- -", r"e\^\{-?j-", r"e\^\{-?j1\\omega", r"\{\}",
             r"(?<![\d.{])1(?:x\[|\\delta|\\sin|\\cos)", r"\\cos ?0\\omega", r"\\sin\\!\\left\(0\\right\)"]
    for t in texts:
        for pat in slips:
            if re.search(pat, t):
                probs.append(f"display slip /{pat}/ in {t!r}")
    return probs


def _num(x):
    q = Fraction(x).limit_denominator(1000)
    if abs(float(q) - x) < 1e-12:
        return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return f"{x:.6g}"


def _differs(a, b):
    return abs(a - b) > 1e-2 * max(abs(b), 1e-3) + 1e-6


def _principal(x):
    while x <= -1:
        x += 2
    while x > 1:
        x -= 2
    return x


def submissions(params, correct):
    cases = []
    mag = correct["mag"]
    cases.append(({"mag": f"{mag:.4g}" if abs(mag - round(mag)) > 1e-9 else str(int(round(mag)))}, {"mag": 1}))
    cases.append(({"mag": f"{mag:.6f}"}, {"mag": 1}))
    if _differs(mag / 2, mag):
        cases.append(({"mag": f"{mag / 2:.6f}"}, {"mag": 0}))                     # dropped the 2 of 2cos / 2j sin
    w0 = math.pi * _f(params["f0"])
    m, J, A = centre_and_amp(params, [w0])
    if A[0].real < 0:
        cases.append(({"mag": f"{A[0].real:.6f}"}, {"mag": 0}))                  # sign kept in the magnitude
    if params["is_ramp"]:
        for i in (1, 2):
            ph = correct[f"ph{i}"]
            w = math.pi * _f(params[f"f{i}"])
            cases.append(({f"ph{i}": _num(ph)}, {f"ph{i}": 1}))
            cases.append(({f"ph{i}": f"{ph:.4f}"}, {f"ph{i}": 1}))
            wrong = [ph + 2 if ph <= 0 else ph - 2]                               # not the principal value
            m, J, A = centre_and_amp(params, [w])
            neg = A[0].real < 0
            t = w / math.pi
            raw = -m * t + J / 2 + (1 if neg else 0)
            if neg:
                wrong.append(_principal(raw - 1))                                  # forgot the +-pi of A < 0
            if J:
                wrong.append(_principal(raw - 0.5))                                # forgot the j (pi/2)
            if m:
                wrong.append(_principal(m * t + J / 2 + (1 if neg else 0)))        # ramp with the wrong sign
            if not (-1 < raw <= 1):
                wrong.append(raw)                                                  # not reduced
            seen = []
            for v in wrong:
                if _differs(v, ph) and all(abs(v - s) > 1e-9 for s in seen):
                    seen.append(v)
                    cases.append(({f"ph{i}": _num(v)}, {f"ph{i}": 0}))
    else:
        for i in (1, 2):
            v = correct[f"hv{i}"]
            cases.append(({f"hv{i}": _num(v)}, {f"hv{i}": 1}))
            cases.append(({f"hv{i}": f"{v:.5g}"}, {f"hv{i}": 1}))
            if v < 0:
                cases.append(({f"hv{i}": _num(-v)}, {f"hv{i}": 0}))                # |H_d| instead of H_d
    nj = correct["njumps"]
    cases.append(({"njumps": str(nj)}, {"njumps": 1}))
    nz = len(zeros_of_mag(params))
    if nz != nj:
        cases.append(({"njumps": str(nz)}, {"njumps": 0}))                        # touched zeros counted too
    _, nw = count_jumps(params)
    if nw:
        cases.append(({"njumps": str(nj + nw)}, {"njumps": 0}))                   # 2 pi wraps counted too
    half = len([z for z in zeros_of_mag(params) if z > 0])
    if half != nj:
        cases.append(({"njumps": str(half)}, {"njumps": 0}))                      # only (0, pi)
    cases.append(({"njumps": "1.5"}, {"njumps": "invalid"}))
    tau = correct["tau"]
    cases.append(({"tau": _num(tau)}, {"tau": 1}))
    cases.append(({"tau": f"{tau:.4f}"}, {"tau": 1}))
    if tau != 0:
        cases.append(({"tau": _num(-tau)}, {"tau": 0}))                           # delay vs advance
    return cases
