"""Independent check of questions/sampling/two-rate-system.

Verification by samples, not by principal values: the ideal D/A output is the unique signal band-limited
to |Omega| < pi/T whose samples at nT are x[n]. So a claimed output y(t) is right iff every term has
Omega' < pi/T and y(nT) = x[n] for all n (checked for n = 0..199 with numpy).
(a) x[n] = x_a(n T1) is compared with the answered digital frequencies; y_a(n T2) with x[n].
(b) each answered term b_i f_i(Omega_i' t) must reproduce the samples of its own input term.
(c) every period in the stored list must reproduce the displayed y_a(t) (parsed from the LaTeX), the
    answer must be the largest and the count the list length; candidate periods from the aliasing
    equations Omega_i T = +-Omega_g T + 2 pi k and Omega_i T = k pi, plus 3000 random periods, must not
    add a solution."""

import math
import random
import re
from fractions import Fraction

import numpy as np

N = np.arange(200)


def _pi(s):
    """'\\tfrac{3\\pi}{4}' / '20\\pi' / '\\pi' -> Fraction multiple of pi."""
    s = s.strip()
    m = re.fullmatch(r"\\tfrac\{(\d*)\\pi\}\{(\d+)\}", s)
    if m:
        return Fraction(int(m.group(1) or 1), int(m.group(2)))
    m = re.fullmatch(r"(\d*)\\pi", s)
    assert m, s
    return Fraction(int(m.group(1) or 1))


TERM = re.compile(r"([+-]?)\s*(\d*)\\(sin|cos)\\!\\left\(([^ ]+) t(?: ([+-]) (\\tfrac\{\\pi\}\{\d+\}|\\pi))?\\right\)|([+-])\s*(\d+)$")


def parse_sum(tex):
    """'y_a(t) = 2\\cos\\!\\left(100\\pi t + \\tfrac{\\pi}{4}\\right) + 1' -> [(c, fn, W, ph)]."""
    body = tex.split("=", 1)[1].strip()
    if not body.startswith(("-", "+")):
        body = "+" + body
    out = []
    for m in TERM.finditer(body):
        if m.group(3):
            c = int(m.group(2) or 1) * (-1 if m.group(1) == "-" else 1)
            ph = Fraction(0)
            if m.group(5):
                ph = _pi(m.group(6)) * (-1 if m.group(5) == "-" else 1)
            out.append((c, m.group(3), _pi(m.group(4)), ph))
        else:
            out.append((int(m.group(8)) * (-1 if m.group(7) == "-" else 1), "const", Fraction(0), Fraction(0)))
    return out


def ev(terms, t):
    """Sum of c fn(W pi t + ph pi) (W in multiples of pi rad/s) at times t (array)."""
    y = np.zeros_like(np.asarray(t, float))
    for c, fn, W, ph in terms:
        if fn == "const":
            y = y + c
        else:
            arg = math.pi * float(W) * t + math.pi * float(ph)
            y = y + c * (np.sin(arg) if fn == "sin" else np.cos(arg))
    return y


def _comps(params):
    return [(c["A"], c["fn"], 2 * Fraction(c["f"]), Fraction(c["ph"])) for c in params["comps"]]


def produces(comps, target, fs):
    """Does A/D -> D/A at rate fs (Hz) turn the input into the target exactly?"""
    fs = float(fs)
    if any(t[1] != "const" and float(t[2]) >= fs for t in target):      # Omega' >= pi/T: not a D/A output
        return False
    n = N / fs
    return bool(np.max(np.abs(ev(comps, n) - ev(target, n))) < 1e-7)


def check(params, correct):
    probs = []
    comps = _comps(params)
    shown = parse_sum(params["xa_tex"])
    if [(c, fn, W, ph) for c, fn, W, ph in shown] != comps:
        probs.append(f"displayed x_a {params['xa_tex']} != stored components {comps}")
    target = parse_sum(params["yc_tex"])
    stored_t = [(t["c"], t["fn"], Fraction(t["W"]), Fraction(t["ph"])) for t in params["target"]]
    if target != stored_t:
        probs.append(f"displayed y_a {params['yc_tex']} != stored target {stored_t}")
    fs1, fs2, fsb = (Fraction(params["fs"][k]) for k in ("a1", "a2", "b"))
    if fs1 == fs2:
        probs.append("(a) needs T1 != T2")
    # (a)
    x = ev(comps, N / float(fs1))
    w = [Fraction(correct["w1"]).limit_denominator(10000), Fraction(correct["w2"]).limit_denominator(10000)]
    if not all(0 <= v <= 1 for v in w):
        probs.append(f"(a) digital frequencies {w} not in [0, pi]")
    xw = ev([(c, fn, wv, ph) for (c, fn, W, ph), wv in zip(comps, w)], N.astype(float))
    if np.max(np.abs(x - xw)) > 1e-7:
        probs.append("(a) x[n] does not match the answered digital frequencies")
    if any(W >= fs1 for c, fn, W, ph in comps):
        probs.append("(a) a tone aliases at T1")
    Wa = [correct["Wa1"], correct["Wa2"]]
    ya = [(c, fn, Fraction(v).limit_denominator(10000), ph) for (c, fn, W, ph), v in zip(comps, Wa)]
    if not produces(_scaled(comps, fs1, fs2), ya, fs2):
        probs.append("(a) y_a(n T2) != x[n] or y_a not band-limited to pi/T2")
    # (b)
    if not (float(comps[0][2]) < float(fsb) < float(comps[1][2])):
        probs.append("(b) T is not below the Nyquist rate of exactly the second tone")
    for i in (0, 1):
        c, fn, W, ph = comps[i]
        b, Wb = correct[f"b{i + 1}"], Fraction(correct[f"Wb{i + 1}"]).limit_denominator(10000)
        term_ph = ph if (i == 0 or ph == 0) else None
        if term_ph is None:
            probs.append("(b) second tone with a phase is not supported by the form")
            continue
        if not produces([comps[i]], [(b, fn, Wb, ph)], fsb):
            probs.append(f"(b) term {i + 1}: {b} {fn}({Wb} pi t) does not reproduce the samples of its input term")
        if Wb <= 0:
            probs.append(f"(b) term {i + 1}: Omega' = {Wb} pi not positive")
    if not produces(comps, [(correct["b1"], comps[0][1], Fraction(correct["Wb1"]).limit_denominator(10000), comps[0][3]),
                            (correct["b2"], comps[1][1], Fraction(correct["Wb2"]).limit_denominator(10000), comps[1][3])], fsb):
        probs.append("(b) the whole output does not reproduce x[n]")
    # (c)
    sols = [Fraction(s) for s in params["sols"]]
    for T in sols:
        if not produces(comps, target, 1 / T):
            probs.append(f"(c) T = {T} does not produce y_a")
    if abs(max(float(T) for T in sols) - correct["Tmax"]) > 1e-12 or correct["count"] != len(sols) or len(set(sols)) != len(sols):
        probs.append(f"(c) answer {correct['Tmax']}, {correct['count']} vs list {sols}")
    gs = {Fraction(t[2]) / 2 for t in target} | {Fraction(0)}
    fi = [Fraction(W) / 2 for c, fn, W, ph in comps]
    gmax = max(gs)
    cands = set()
    for f in fi:
        for g in gs | set(fi):
            for num in (f + g, f - g, 2 * f):
                k = 1
                while num > 0 and num / k > 2 * gmax and k < 400:
                    cands.add(num / k)
                    k += 1
    found = sorted(1 / fs for fs in cands if produces(comps, target, fs))
    if found != sorted(sols):
        probs.append(f"(c) candidate scan finds {found}, list says {sorted(sols)}")
    rng = random.Random(1)
    for _ in range(3000):
        T = rng.uniform(1e-4, 2 / float(2 * gmax))
        if produces(comps, target, 1 / T) and not any(abs(T - float(s)) < 1e-9 for s in sols):
            probs.append(f"(c) random T = {T} also works")
            break
    if params["kind"] == "kill" and len(target) != 1:
        probs.append("kill variant must keep one term")
    return probs


def _scaled(comps, fs1, fs2):
    """x[n] = x_a(n T1) rewritten as a signal sampled at T2: x_a(t T2/T1) (time scaling by fs2/fs1)."""
    k = Fraction(fs2) / Fraction(fs1)
    return [(c, fn, W * k, ph) for c, fn, W, ph in comps]


def _s(v):
    q = Fraction(v).limit_denominator(100000)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def _differs(u, v, atol=1e-6):
    return abs(float(u) - float(v)) > 2 * (atol + 1e-3 * abs(float(v)))


def submissions(params, correct):
    comps = _comps(params)
    W1, W2 = comps[0][2], comps[1][2]
    fs1, fs2, fsb = (Fraction(params["fs"][k]) for k in ("a1", "a2", "b"))
    names = ["w1", "w2", "Wa1", "Wa2", "b1", "Wb1", "b2", "Wb2", "Tmax"]
    cases = [
        ({n: _s(correct[n]) for n in names} | {"count": str(correct["count"])}, {n: 1 for n in names} | {"count": 1}),
        ({n: f"{correct[n]:.4g}" for n in names} | {"count": f" {correct['count']} "}, {n: 1 for n in names} | {"count": 1}),
        ({"Wa2": f"{2 * int(correct['Wa2'])}/2", "b2": f"{correct['b2']:.1f}"}, {"Wa2": 1, "b2": 1}),
        ({"count": "2.5"}, {"count": "invalid"}),
    ]
    # (a)
    cases.append(({"w1": repr(correct["w1"] * math.pi)}, {"w1": 0}))               # radians, not multiples of pi
    cases.append(({"Wa1": _s(W1), "Wa2": _s(W2)}, {"Wa1": 0, "Wa2": 0}))           # D/A with T1: the input back
    cases.append(({"Wa2": _s(Fraction(correct["Wa2"]) / 2)}, {"Wa2": 0}))          # hertz instead of rad/s
    if _differs(correct["w2"] / float(fs2), correct["Wa2"]):
        cases.append(({"Wa2": repr(correct["w2"] / float(fs2))}, {"Wa2": 0}))      # omega * T2
    if _differs(float(W2 / fs2), correct["w2"]):
        cases.append(({"w2": _s(W2 / fs2)}, {"w2": 0}))                              # sampled with T2
    # (b)
    cases.append(({"b2": _s(-Fraction(correct["b2"]))}, {"b2": 0}))                 # wrong sign
    if _differs(float(W2), correct["Wb2"]):
        cases.append(({"Wb2": _s(W2)}, {"Wb2": 0}))                                  # "no aliasing"
    fold = abs(2 * fsb - W2)                                                         # Omega_s - Omega_2
    if _differs(float(fold), correct["Wb2"]):
        cases.append(({"Wb2": _s(fold)}, {"Wb2": 0}))
    cases.append(({"b1": _s(-Fraction(correct["b1"]))}, {"b1": 0}))
    # (c)
    sols = sorted(Fraction(s) for s in params["sols"])
    K = len(sols)
    if K > 1:
        cases.append(({"Tmax": _s(sols[0])}, {"Tmax": 0}))                         # the smallest instead
    cases.append(({"count": str(K + 1)}, {"count": 0}))
    cases.append(({"count": str(K - 1)}, {"count": 0}))
    nyq1 = 1 / W1                                                                    # T = pi / Omega_1
    if _differs(float(nyq1), correct["Tmax"], 1e-9):
        cases.append(({"Tmax": _s(nyq1)}, {"Tmax": 0}))                              # Nyquist of the kept tone
    if params["kind"] == "kill":
        nxt = Fraction(K + 1) / W2                                                   # k = K + 1 aliases the first tone
        cases.append(({"Tmax": _s(nxt)}, {"Tmax": 0}))
    elif params["kind"] == "dc":
        cases.append(({"Tmax": _s(Fraction(2 * K + 1) / W2)}, {"Tmax": 0}))         # odd multiple of pi: (-1)^n
    else:
        f2, f3 = W2 / 2, None
        t = [Fraction(x["W"]) / 2 for x in params["target"]]
        f3 = [g for g in t if g != W1 / 2][0]
        if _differs(float(1 / (f2 + f3)), correct["Tmax"], 1e-9):
            cases.append(({"Tmax": _s(1 / (f2 + f3))}, {"Tmax": 0}))                 # folded: wrong sign
    return cases
