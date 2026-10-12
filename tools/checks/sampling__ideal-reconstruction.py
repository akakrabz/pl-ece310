"""Independent check of questions/sampling/ideal-reconstruction.

(a) The ideal D/A output is evaluated by brute force: the interpolation sum
    sum_{|n| <= M} x[n] sinc(t/T - n) (numpy's normalized sinc, M = 2e5 terms) at off-grid times,
    with x[n] read from the displayed LaTeX, and compared with C cos(2 pi f0 t + phi pi).
(b) x[n] is parsed from the displayed LaTeX and the times from the displayed t_1, t_2, t_3 and T;
    the ZOH value is x[floor(t1/T)], the ideal values are finite sinc sums.
(c) Every displayed option is sampled at f_s (400 samples) and compared with the samples of the
    displayed output, which is band-limited below f_s/2; an option is correct iff they agree."""

import math
import re
from fractions import Fraction

import numpy as np

M = 200_000


def _frac_pi(s):
    """LaTeX multiple of pi -> Fraction: '\\tfrac{3\\pi}{4}' -> 3/4, '\\pi' -> 1, '2\\pi' -> 2, '0' -> 0."""
    s = s.strip()
    neg = s.startswith("-")
    s = s.lstrip("-").strip()
    if s == "0":
        return Fraction(0)
    m = re.fullmatch(r"\\tfrac\{(\d*)\\pi\}\{(\d+)\}", s)
    if m:
        v = Fraction(int(m.group(1) or 1), int(m.group(2)))
    else:
        m = re.fullmatch(r"(\d*)\\pi", s)
        assert m, s
        v = Fraction(int(m.group(1) or 1))
    return -v if neg else v


def _parse_cos(tex):
    """'A\\cos\\left(W n + P\\right)' (or with t) -> (A, W/pi, P/pi), W, P as Fractions of pi."""
    m = re.fullmatch(r"(\d*)\\cos\\left\((.*?)([nt])(?: ([+-]) (.*))?\\right\)", tex.strip())
    assert m, tex
    A = int(m.group(1) or 1)
    W = _frac_pi(m.group(2))
    P = Fraction(0)
    if m.group(4):
        P = _frac_pi(m.group(5)) * (1 if m.group(4) == "+" else -1)
    return A, W, P


def _parse_x(params):
    """x[n] as {n: value} from the displayed LaTeX (sum of deltas or a sequence with an arrow)."""
    b = params["b"]
    tex = b["x_tex"]
    if tex.startswith(r"\{"):
        items = tex[2:-2].split(r",\ ")
        zero = next(i for i, it in enumerate(items) if r"\underset{\uparrow}" in it)
        out = {}
        for i, it in enumerate(items):
            it = re.sub(r"\\underset\{\\uparrow\}\{(.*)\}", r"\1", it)
            v = int(it)
            if v:
                out[i - zero] = v
        return out
    out = {}
    for sgn, coef, k in re.findall(r"([+-]?)\s*(\d*)\\delta\[n(?:-(\d+))?\]", tex):
        v = int(coef or 1) * (-1 if sgn == "-" else 1)
        out[int(k or 0)] = v
    return out


def _num_after_eq(tex):
    m = re.search(r"=\s*([\d.]+)\\", tex)
    return float(m.group(1))


def _ideal(xs, u):
    """Ideal D/A output at t = uT for a finite sequence."""
    return sum(v * np.sinc(u - n) for n, v in xs.items())


def check(params, correct):
    probs = []
    # ---------------------------------------------------------------- (a)
    a = params["a"]
    A, W, P = _parse_cos(a["x_tex"])
    fs = a["fs"]
    T = 1.0 / fs
    n = np.arange(-M, M + 1)
    xn = A * np.cos(math.pi * float(W) * n + math.pi * float(P))
    C, f0, phi = correct["amp"], correct["f0"], correct["phi"]
    for u in (0.31, 1.73, -2.46, 7.5):
        y = float(np.sum(xn * np.sinc(u - n)))
        want = C * math.cos(2 * math.pi * f0 * u * T + math.pi * phi)
        if abs(y - want) > 2e-3 * A:
            probs.append(f"(a) sinc sum at t = {u}T is {y:.5f}, answer gives {want:.5f}")
            break
    if not (0 < f0 < fs / 2 and C > 0 and -1 < phi <= 1):
        probs.append(f"(a) answer out of range: f0={f0}, C={C}, phi={phi}")
    yA, yW, yP = _parse_cos(a["y_tex"])
    if (yA, float(yW) / 2, float(yP)) != (C, f0, phi):
        probs.append(f"(a) displayed y_a {a['y_tex']} disagrees with the answers")
    # ---------------------------------------------------------------- (b)
    b = params["b"]
    xs = _parse_x(params)
    if xs != {int(k): v for k, v in b["xs"].items()}:
        probs.append(f"(b) displayed x {b['x_tex']} != stored {b['xs']}")
    Tb = _num_after_eq(b["Tb_tex"])
    u1, u2, u3 = (_num_after_eq(b[k]) / Tb for k in ("t1_tex", "t2_tex", "t3_tex"))
    zoh = xs.get(math.floor(u1 + 1e-12), 0)
    if abs(zoh - correct["zoh"]) > 1e-12:
        probs.append(f"(b) ZOH at {u1}T should be {zoh}, answer {correct['zoh']}")
    if abs(u1 - round(u1)) < 1e-9:
        probs.append("(b) t1 is a sample instant (ZOH value ambiguous at a jump)")
    v2, v3 = _ideal(xs, u2), _ideal(xs, u3)
    if abs(v2 - correct["ideal_s"]) > 1e-9:
        probs.append(f"(b) ideal at {u2}T is {v2}, answer {correct['ideal_s']}")
    if abs(v3 - correct["ideal_h"]) > 1e-9 * max(1, abs(v3)):
        probs.append(f"(b) ideal at {u3}T is {v3}, answer {correct['ideal_h']}")
    if abs(u2 - round(u2)) > 1e-9 or abs(u3 - math.floor(u3) - 0.5) > 1e-9:
        probs.append(f"(b) t2/T = {u2} should be an integer and t3/T = {u3} a half-integer")
    if abs(v3) < 0.05:
        probs.append(f"(b) ideal value at t3 too small: {v3}")
    # ---------------------------------------------------------------- (c)
    c = params["c"]
    fsc = c["fsc"]
    B, W0, P0 = _parse_cos(c["yc_tex"])
    f_out = float(W0) / 2
    if not (0 < f_out < fsc / 2):
        probs.append(f"(c) output frequency {f_out} Hz not inside (0, f_s/2)")
    nn = np.arange(400)
    y_s = B * np.cos(2 * math.pi * f_out * nn / fsc + math.pi * float(P0))
    keys = sorted(o["key"] for o in correct["inputs"])
    want = []
    for o in params["inputs"]:
        tex = o["html"].strip().strip("$")
        Bi, Wi, Pi = _parse_cos(tex.split("=", 1)[1])
        xs_c = Bi * np.cos(2 * math.pi * (float(Wi) / 2) * nn / fsc + math.pi * float(Pi))
        if np.allclose(xs_c, y_s, atol=1e-7):
            want.append(o["key"])
        d = (float(Wi) / 2) / fsc * 2 % 1       # omega/pi mod 1: 0 means the samples sit on 0 or pi
        if abs(d) < 1e-12 or abs(d - 1) < 1e-12:
            probs.append(f"(c) option {tex} lands on omega = 0 or pi")
    if sorted(want) != keys:
        probs.append(f"(c) correct options {keys} but sampling says {sorted(want)}")
    if not want or len(want) == len(params["inputs"]):
        probs.append("(c) needs at least one correct and one wrong option")
    if len({o["html"] for o in params["inputs"]}) != len(params["inputs"]):
        probs.append("(c) duplicate options")
    return probs


def _num(v):
    v = float(v)
    return str(int(v)) if v.is_integer() else repr(v)


def _differs(w, r):
    return abs(float(w) - float(r)) > 2 * (1e-6 + 1e-3 * abs(float(r)))


def submissions(params, correct):
    a, b = params["a"], params["b"]
    cases = []
    # (a) frequency, amplitude, phase
    f0, C, phi = correct["f0"], correct["amp"], correct["phi"]
    fs = a["fs"]
    A, W, P = _parse_cos(a["x_tex"])
    cases.append(({"f0": _num(f0) + ".0", "amp": f"{2 * int(C)}/2"}, {"f0": 1, "amp": 1}))
    cases.append(({"f0": f"{int(2 * f0)}/2", "phi": f"{phi:.4g}"}, {"f0": 1, "phi": 1}))
    cases.append(({"f0": f"{2 * math.pi * f0:.6g}"}, {"f0": 0}))                       # rad/s instead of Hz
    cases.append(({"f0": f"{f0 / fs:.6g}"}, {"f0": 0}))                                 # omega_0/(2 pi), no T
    raw = float(W) * fs / 2
    if _differs(raw, f0):
        cases.append(({"f0": _num(raw)}, {"f0": 0}))                                    # the non-principal frequency
    cases.append(({"amp": _num(C * fs)}, {"amp": 0}))                                   # kept the sampler's 1/T
    cases.append(({"amp": _num(C / 2)}, {"amp": 0}))                                    # half of each exponential
    if _differs(float(P), phi):
        cases.append(({"phi": str(Fraction(P))}, {"phi": 0}))                           # phase not flipped
    else:
        cases.append(({"phi": str(-Fraction(P))}, {"phi": 0}))                          # flipped when it should not be
    fr = Fraction(phi).limit_denominator(100)
    cases.append(({"amp": _num(C) + ".0", "phi": f"{2 * fr.numerator}/{2 * fr.denominator}"}, {"amp": 1, "phi": 1}))
    cases.append(({"phi": f"{phi * math.pi:.4g}"}, {"phi": 0}))                         # radians, not a multiple of pi
    # (b)
    xs = {int(k): v for k, v in b["xs"].items()}
    m1, m2, m3 = b["m1"], b["m2"], b["m3"]
    cases.append(({"zoh": _num(correct["zoh"]) + ".0", "ideal_s": _num(correct["ideal_s"]) + ".0"},
                  {"zoh": 1, "ideal_s": 1}))
    cases.append(({"zoh": f"{2 * int(correct['zoh'])}/2", "ideal_s": f"{2 * int(correct['ideal_s'])}/2"},
                  {"zoh": 1, "ideal_s": 1}))
    cases.append(({"zoh": _num(xs.get(m1 + 1, 0))}, {"zoh": 0}))                        # held the next sample
    avg = (xs.get(m1, 0) + xs.get(m1 + 1, 0)) / 2
    if _differs(avg, correct["zoh"]):
        cases.append(({"zoh": _num(avg)}, {"zoh": 0}))                                  # interpolated instead of held
    if xs.get(m1 - 1, 0) != xs.get(m1, 0):
        cases.append(({"zoh": _num(xs.get(m1 - 1, 0))}, {"zoh": 0}))                    # the previous sample
    if xs.get(m2 + 1, 0) != xs.get(m2, 0):
        cases.append(({"ideal_s": _num(xs.get(m2 + 1, 0))}, {"ideal_s": 0}))            # sinc centred one sample late
    if xs.get(m2 - 1, 0) != xs.get(m2, 0):
        cases.append(({"ideal_s": _num(xs.get(m2 - 1, 0))}, {"ideal_s": 0}))            # sinc centred one sample early
    unnorm = sum(v * (1.0 if n == m2 else math.sin(m2 - n) / (m2 - n)) for n, v in xs.items())
    if _differs(unnorm, correct["ideal_s"]):
        cases.append(({"ideal_s": f"{unnorm:.6g}"}, {"ideal_s": 0}))                    # sin(x)/x read in radians
    v3 = correct["ideal_h"]
    cases.append(({"ideal_h": f"{v3:.4g}"}, {"ideal_h": 1}))
    cases.append(({"ideal_h": f"{v3:.7f}"}, {"ideal_h": 1}))
    lin = (xs.get(m3, 0) + xs.get(m3 + 1, 0)) / 2
    if _differs(lin, v3):
        cases.append(({"ideal_h": _num(lin)}, {"ideal_h": 0}))                          # linear interpolation
    near = 2 / math.pi * (xs.get(m3, 0) + xs.get(m3 + 1, 0))
    if _differs(near, v3):
        cases.append(({"ideal_h": f"{near:.6g}"}, {"ideal_h": 0}))                      # only the two nearest sincs
    cases.append(({"ideal_h": "0"}, {"ideal_h": 0}))
    # (c) checkbox, net-correct
    opts = params["inputs"]
    ck = [o["key"] for o in correct["inputs"]]
    wrong = [o["key"] for o in opts if o["key"] not in ck]
    c = len(ck)
    cases.append(({"inputs": list(reversed(ck))}, {"inputs": 1}))
    if c >= 2:
        cases.append(({"inputs": ck[:-1]}, {"inputs": (c - 1) / c}))                    # missed one alias
    cases.append(({"inputs": ck + wrong[:1]}, {"inputs": (c - 1) / c}))                 # one wrong input ticked
    cases.append(({"inputs": wrong}, {"inputs": 0}))
    cases.append(({"inputs": []}, {"inputs": "invalid"}))
    return cases
