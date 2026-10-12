"""Independent check of questions/signals/sampling-and-sequences.

(a) the displayed LaTeX of x[n] is translated to Python (u, delta as functions) and evaluated
    sample by sample; (b) the displayed x_a(t) and sampling rate/period are parsed and x[k] =
    x_a(kT) is computed with numpy; (c) the classification of each bank description is typed
    here from Lecture 1's definitions (continuous/discrete time x continuous/discrete values)."""

import math
import re
from fractions import Fraction

import numpy as np


# ----------------------------------------------------------------------------- (a)
def u(m):
    return 1 if m >= 0 else 0


def d(m):
    return 1 if m == 0 else 0


def seq_to_py(tex):
    s = tex.replace(r"\left(", "(").replace(r"\right)", ")").replace(r"\,", "").replace(r"\delta[", "d[")
    s = s.replace("^{n}", "**n").replace(" ", "")
    s = s.replace("[", "(").replace("]", ")")
    s = re.sub(r"(?<=[\dn)])(?=[(udn])", "*", s)         # implicit products: 2u( , 2n , n( , )u( , )(
    return s


def eval_seq(tex, lo=-20, hi=20):
    code = seq_to_py(tex)
    vals = {}
    for n in range(lo, hi + 1):
        vals[n] = eval(code, {"__builtins__": {}}, {"u": u, "d": d, "n": n})
    return vals


def trim(vals):
    nz = sorted(n for n, v in vals.items() if v != 0)
    return nz[0], [vals[n] for n in range(nz[0], nz[-1] + 1)]


# ----------------------------------------------------------------------------- (b)
def tex_frac_pi(s):
    """' + \\tfrac{\\pi}{3}' / ' - \\tfrac{2\\pi}{3}' / ' + \\pi' / '' -> multiple of pi (Fraction)."""
    s = s.replace(" ", "")
    if not s:
        return Fraction(0)
    sign = -1 if s[0] == "-" else 1
    s = s[1:]
    m = re.fullmatch(r"\\tfrac\{(\d*)\\pi\}\{(\d+)\}", s)
    if m:
        return sign * Fraction(int(m.group(1) or 1), int(m.group(2)))
    m = re.fullmatch(r"(\d*)\\pi", s)
    if m:
        return sign * Fraction(int(m.group(1) or 1))
    raise ValueError(f"cannot parse phase {s!r}")


def parse_xa(tex):
    """'4\\cos\\!\\left(2\\pi(1500)t + \\tfrac{\\pi}{3}\\right)' -> (A, Omega0 in rad/s, theta in rad)."""
    m = re.fullmatch(r"(\d+)\\cos\\!\\left\((?:2\\pi\((\d+)\)t|(\d+)\\pi t)(.*)\\right\)", tex)
    if not m:
        raise ValueError(f"cannot parse x_a {tex!r}")
    A = int(m.group(1))
    Om = 2 * math.pi * int(m.group(2)) if m.group(2) else math.pi * int(m.group(3))
    th = float(tex_frac_pi(m.group(4))) * math.pi
    return A, Om, th


def parse_rate(tex):
    """'f_s = 8\\text{ kHz}' / 'T = 125\\ \\mu\\text{s}' / 'T = \\tfrac{1}{6000}\\text{ s}' -> f_s in Hz."""
    m = re.fullmatch(r"f_s = ([\d.]+)\\text\{ (k?)Hz\}", tex)
    if m:
        return float(m.group(1)) * (1000 if m.group(2) else 1)
    m = re.fullmatch(r"T = (\d+)\\ \\mu\\text\{s\}", tex)
    if m:
        return 1e6 / int(m.group(1))
    m = re.fullmatch(r"T = (\d+)\\text\{ ms\}", tex)
    if m:
        return 1e3 / int(m.group(1))
    m = re.fullmatch(r"T = \\tfrac\{1\}\{(\d+)\}\\text\{ s\}", tex)
    if m:
        return float(m.group(1))
    raise ValueError(f"cannot parse rate {tex!r}")


# ----------------------------------------------------------------------------- (c)
EXPECTED = {   # description start -> (continuous time?, continuous values?)
    "The voltage $v(t)$ produced by a microphone": (True, True),
    "The air pressure at your eardrum": (True, True),
    "The height of the tide": (True, True),
    "The temperature of a cooling cup": (True, True),
    "The voltage on a digital logic line": (True, False),
    "The number of people inside a store": (True, False),
    "The floor number of an elevator": (True, False),
    "The number of cars waiting at a red light": (True, False),
    "The samples $x[n] = x_a(nT)$ of a microphone": (False, True),
    "The samples of a seismometer voltage": (False, True),
    "The exact air temperature at noon": (False, True),
    "The exact average speed": (False, True),
    "Hourly temperature readings": (False, False),
    "A song on a CD": (False, False),
    "The daily closing price": (False, False),
    "The number of cars passing a sensor": (False, False),
}


def expected_option(desc):
    hits = [v for k, v in EXPECTED.items() if desc.startswith(k)]
    assert len(hits) == 1, desc
    ct, cv = hits[0]
    t = "continuous-time" if ct else "discrete-time"
    v = "continuous-valued" if cv else "discrete-valued"
    return t, v


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(b))


def check(params, correct):
    probs = []
    # (a)
    vals = eval_seq(params["seq_tex"])
    if vals[-20] != 0 or vals[20] != 0:
        probs.append("x[n] is not finite on [-20, 20]")
    n0, ys = trim(vals)
    got = [float(v) for v in correct["xs"]["_value"][0]]
    if correct["nfirst"] != n0:
        probs.append(f"first index {correct['nfirst']} != {n0} (from {seq_to_py(params['seq_tex'])})")
    if got != [float(v) for v in ys]:
        probs.append(f"samples {got} != {ys}")
    if len(ys) < 2 or sum(v != 0 for v in ys) < 2 or len(ys) > 8:
        probs.append(f"unexpected length / nonzero count: {ys}")
    if any(v != int(v) for v in ys):
        probs.append("non-integer samples")
    m = re.search(r"On those samples \$x\[n\] = ([^$]+)\$", params["seq_expl_html"])
    if params["seq_family"] == "ramp":                       # the printed linear formula a n + b
        if not m:
            probs.append("ramp explanation without the formula for x[n]")
        else:
            code = re.sub(r"(?<=\d)n", "*n", m.group(1).replace(" ", ""))
            for n in range(n0, n0 + len(ys)):
                if vals[n] != 0 and eval(code, {"__builtins__": {}}, {"n": n}) != vals[n]:
                    probs.append(f"printed formula x[n] = {m.group(1)} gives {eval(code, {}, {'n': n})} at n = {n}, "
                                 f"but x[{n}] = {vals[n]}")
                    break
    listed = params["seq_samples_tex"].split(r"\}")[0].count(r",\ ") + 1
    if listed != len(ys):
        probs.append(f"the solution lists {listed} samples but the graded row has {len(ys)}")
    # (b)
    A, Om, th = parse_xa(params["xa_tex"])
    fs = parse_rate(params["rate_tex"])
    if not math.isclose(fs, params["fs"], rel_tol=1e-12):
        probs.append(f"rate text {params['rate_tex']} != fs {params['fs']}")
    w0 = Om / fs
    if not 0 < w0 < math.pi:
        probs.append(f"w0 = {w0} not in (0, pi): aliasing would be involved")
    if not close(w0 / math.pi, correct["w0"]):
        probs.append(f"w0/pi {correct['w0']} != {w0 / math.pi}")
    k = params["k"]
    xk = A * math.cos(Om * k / fs + th)
    if not close(xk, correct["xk"], 1e-9):
        probs.append(f"x[{k}] {correct['xk']} != x_a(kT) = {xk}")
    if abs(xk) < 1e-9:
        probs.append("x[k] = 0")
    if Fraction(correct["xk"]).limit_denominator(4) != Fraction(correct["xk"]):
        probs.append(f"x[k] = {correct['xk']} is not a simple rational number")
    if len(_xk_mistakes(params, correct)) < 2:
        probs.append("fewer than two classic mistakes for x[k]")
    # (c)
    t, v = expected_option(params["sig_desc_html"])
    html = correct["sig"]["html"]
    if not (html.startswith(f"{t} and {v}")):
        probs.append(f"classification {html!r} but expected {t}, {v}")
    if len(params["sig"]) != 4:
        probs.append("expected 4 options in (c)")
    return probs


def _num(x):
    return f"{x:.6g}"


def _row(v, sep=", "):
    return "[" + sep.join(str(int(a)) for a in v) + "]"


def _xk_mistakes(params, correct):
    A, Om, th = parse_xa(params["xa_tex"])
    fs = parse_rate(params["rate_tex"])
    k, xk = params["k"], correct["xk"]
    cands = [
        -xk,                                          # sign slip
        A * math.cos(Om * k / fs - th),               # phase sign flipped
        A * math.cos(Om * k + th),                    # x_a(k): t = k seconds instead of t = kT
        A * math.cos(Om * k / fs / 2 + th),           # w0 = f0/fs*pi (forgot the factor 2)
        A * math.cos(Om * k / fs),                    # phase dropped
    ]
    out = []
    for c in cands:
        if abs(c - xk) > 2e-2 * max(1.0, abs(xk)) and all(abs(c - o) > 1e-6 for o in out):
            out.append(c)
    return out


def submissions(params, correct):
    vals = eval_seq(params["seq_tex"])
    n0, ys = trim(vals)
    cases = [
        ({"nfirst": str(n0 + 1)}, {"nfirst": 0, "xs": 1}),
        ({"nfirst": str(n0 - 1)}, {"nfirst": 0, "xs": 1}),
        ({"nfirst": "0.5"}, {"nfirst": "invalid"}),
        ({"nfirst": f" {n0} "}, {"nfirst": 1}),                                # surrounding spaces
        ({"nfirst": f"+{n0}" if n0 >= 0 else f"-0{-n0}"}, {"nfirst": 1}),      # explicit sign / leading zero
        ({"xs": _row(ys, " ")}, {"xs": 1}),                                   # MATLAB style
        ({"xs": "[" + ", ".join(f"{int(v)}.0" for v in ys) + "]"}, {"xs": 1}),  # decimals
        ({"xs": _row(list(ys) + [0])}, {"xs": 0}),                             # trailing zero
    ]
    if list(ys) != list(ys)[::-1]:
        cases.append(({"xs": _row(list(ys)[::-1])}, {"xs": 0}))
    ef, eys = params["seq_ext"]                                              # off-by-one window edge
    if eys and (ef != n0 or list(eys) != list(ys)):
        exp = {}
        if ef != n0:
            exp["nfirst"] = 0
        if list(eys) != list(ys):
            exp["xs"] = 0
        cases.append(({"nfirst": str(ef), "xs": _row(eys)}, exp))
    if params["seq_lead_zero"]:                                              # kept the zero sample at the start
        a = params["seq_a"]
        lead = [vals[n] for n in range(a, n0)] + list(ys)
        cases.append(({"nfirst": str(a), "xs": _row(lead)}, {"nfirst": 0, "xs": 0}))

    w = Fraction(*params["w0_pq"])
    cases += [
        ({"w0": f"{w.numerator}/{w.denominator}"}, {"w0": 1}),
        ({"w0": f"{2 * w.numerator}/{2 * w.denominator}"}, {"w0": 1}),
        ({"w0": f"{float(w):.4g}"}, {"w0": 1}),
        ({"w0": _num(float(w) / 2)}, {"w0": 0}),                      # f0/fs: forgot the 2
        ({"w0": _num(2 * params["f0"] * params["fs"])}, {"w0": 0}),    # 2 pi f0 / T: T used as a rate
        ({"w0": _num(float(w) * math.pi)}, {"w0": 0}),                 # radians, not a multiple of pi
        ({"w0": _num(2 * float(w))}, {"w0": 0}),                       # 2 pi applied twice (Omega read as Hz)
    ]
    xk = correct["xk"]
    q = Fraction(xk).limit_denominator(4)
    cases.append(({"xk": f"{q.numerator}/{q.denominator}"}, {"xk": 1}))
    cases.append(({"xk": f"{xk:.4g}"}, {"xk": 1}))
    for m in _xk_mistakes(params, correct):
        cases.append(({"xk": _num(m)}, {"xk": 0}))

    right = correct["sig"]["key"]
    for opt in params["sig"]:
        cases.append(({"sig": opt["key"]}, {"sig": 1 if opt["key"] == right else 0}))
    return cases
