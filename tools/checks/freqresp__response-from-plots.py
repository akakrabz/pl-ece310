"""Independent check of questions/freqresp/response-from-plots.

H_d is rebuilt from params["filt"] (levels, edges, slope) with the frequency reduced to (-pi, pi]
by atan2. The displayed input x_tex is turned into a numpy expression and must equal the input
described in params["inp"]; the true output is computed term by term with complex exponentials
(the sinusoid split into its two halves, each multiplied by H_d at its own frequency) and must
equal D + E trig(w1 n + phi) + F (exp term) + G (-1)^n built from the answers, at n = -12..12; the
printed y_tex must be the same sequence. Variant rules: no input frequency on a band edge, phases
never +-pi, E >= 0, removed terms have amplitude 0 and phase 0, the figures carry the edge, level
and phase-end labels. Probes: decimals and fractions; the gain forgotten; the other step read;
the phase sign flipped; the unreduced angle; the phase of a removed term; (-1)^n kept although
removed, or without the sign (-1)^M."""

import cmath
import math
import re
from fractions import Fraction

import numpy as np

N = np.arange(-12, 13)


def _f(s):
    return float(Fraction(s))


def Hd(params, w):
    fl = params["filt"]
    w = math.atan2(math.sin(w), math.cos(w))
    if abs(w + math.pi) < 1e-12:
        w = math.pi
    a = abs(w) / math.pi
    if a < _f(fl["w1"]) - 1e-12:
        mag = fl["L1"]
    elif a < _f(fl["w2"]) - 1e-12 or (fl["to_pi"] and a <= 1 + 1e-12):
        mag = fl["L2"]
    else:
        mag = 0
    return mag * cmath.exp(-1j * fl["M"] * w)


def components(params):
    """[(coefficient, frequency in rad)] of x[n]."""
    inp = params["inp"]
    w1, th = math.pi * _f(inp["w1"]), math.pi * _f(inp["th"])
    if inp["trig"] == "sin":
        th -= math.pi / 2
    w2 = {"jn": math.pi / 2, "mjn": -math.pi / 2}.get(inp["exp_form"], math.pi * _f(inp["w2"]))
    out = [(inp["A0"], 0.0), (inp["A1"] / 2 * cmath.exp(1j * th), w1), (inp["A1"] / 2 * cmath.exp(-1j * th), -w1),
           (inp["B"], w2)]
    if inp["has_alt"]:
        out.append((inp["C"], math.pi))
    return out, w1, w2


def tex_to_fn(tex):
    s = tex.replace(r"\left(", "(").replace(r"\right)", ")").replace(r"\,", " ")
    s = re.sub(r"\\[td]?frac\{(\d*)\\pi\}\{(\d+)\}", lambda m: f"({m.group(1) or 1}*pi/{m.group(2)})", s)
    s = re.sub(r"e\^\{(-?)j([^{}]*)\}", lambda m: f"exp({m.group(1)}1j*({m.group(2)}))", s)
    s = s.replace("(-j)^n", "((-1j)**n)").replace("(-1)^n", "((-1.0)**n)").replace("j^n", "((1j)**n)")
    s = s.replace(r"\cos", "cos").replace(r"\sin", "sin").replace(r"\pi", "pi")
    s = re.sub(r"(\d)\s*pi", r"\1*pi", s)
    s = re.sub(r"(\d|\))\s*(n|\(|cos|sin|exp)", r"\1*\2", s)
    s = re.sub(r"pi\s*(n)", r"pi*\1", s)
    s = re.sub(r"\)\s*(n|\()", r")*\1", s)
    assert re.fullmatch(r"[\d\s()+\-*/.npijcosxe]+", s), (tex, s)
    return lambda n: eval(s, {"__builtins__": {}}, {"n": np.asarray(n, float), "pi": math.pi, "cos": np.cos, "sin": np.sin,
                                                     "exp": np.exp})


def y_true(params):
    comps, _, _ = components(params)
    return sum(c * Hd(params, w) * np.exp(1j * w * N) for c, w in comps)


def y_from(params, ans):
    inp = params["inp"]
    _, w1, w2 = components(params)
    trig = np.cos if inp["trig"] == "cos" else np.sin
    y = ans["D"] + ans["E"] * trig(w1 * N + math.pi * ans["phi"]) \
        + ans["Fmag"] * np.exp(1j * math.pi * ans["Fang"]) * np.exp(1j * w2 * N)
    if inp["has_alt"]:
        y = y + ans["G"] * (-1.0) ** N
    return y


def check(params, correct):
    probs = []
    inp, fl = params["inp"], params["filt"]
    comps, w1, w2 = components(params)
    x_disp = np.asarray(tex_to_fn(params["x_tex"])(N), complex) * np.ones(len(N))
    x_true = sum(c * np.exp(1j * w * N) for c, w in comps)
    if np.max(np.abs(x_disp - x_true)) > 1e-9:
        probs.append(f"displayed x[n] {params['x_tex']!r} differs from params['inp']")
    yt = y_true(params)
    ya = y_from(params, correct)
    if np.max(np.abs(yt - ya)) > 1e-9:
        probs.append(f"answers give a different output: max err {np.max(np.abs(yt - ya))}")
    y_disp = np.asarray(tex_to_fn(params["y_tex"])(N), complex) * np.ones(len(N))
    if np.max(np.abs(y_disp - yt)) > 1e-9:
        probs.append(f"printed y[n] {params['y_tex']!r} is wrong")
    for w in (w1, w2):
        a = abs(math.atan2(math.sin(w), math.cos(w))) / math.pi
        if min(abs(a - _f(fl["w1"])), abs(a - _f(fl["w2"]))) < 1e-9:
            probs.append(f"input frequency {a} pi on a band edge")
    for name in ("phi", "Fang"):
        if not -1 < correct[name] < 1:
            probs.append(f"{name} = {correct[name]} not in (-1, 1)")
    if correct["E"] < 0:
        probs.append("E < 0")
    g1, g2 = abs(Hd(params, w1)), abs(Hd(params, w2))
    if g1 == 0 and (correct["E"] != 0 or correct["phi"] != 0):
        probs.append("removed sinusoid must have E = phi = 0")
    if g2 == 0 and (correct["Fmag"] != 0 or correct["Fang"] != 0):
        probs.append("removed exponential must have F = 0, angle 0")
    if g1 == 0 and g2 == 0:
        probs.append("both the sinusoid and the exponential are removed")
    if g1:
        want = (_f(inp["th"]) - fl["M"] * w1 / math.pi + 1) % 2 - 1
        if abs(want - correct["phi"]) > 1e-9 and abs(abs(want) - 1) > 1e-9:
            probs.append(f"phi {correct['phi']} vs {want}")
    # figures
    labels = set(re.findall(r">([^<>]+)</text>", params["mag_html"]))
    want_l = {_fmt_pi(_f(fl["w1"])), _fmt_pi(-_f(fl["w1"])), "π", "−π", str(fl["L1"]), str(fl["L2"])}
    if not fl["to_pi"]:
        want_l |= {_fmt_pi(_f(fl["w2"])), _fmt_pi(-_f(fl["w2"]))}
    if not want_l <= labels:
        probs.append(f"magnitude figure misses labels {want_l - labels}")
    plabels = set(re.findall(r">([^<>]+)</text>", params["phase_html"]))
    top = fl["M"] * _f(fl["w2"])
    if not {_fmt_pi(top), _fmt_pi(-top)} <= plabels:
        probs.append(f"phase figure misses the end labels {_fmt_pi(top)}")
    return probs


def _fmt_pi(c):
    c = Fraction(c).limit_denominator(100)
    if c == 0:
        return "0"
    a, b = abs(c.numerator), c.denominator
    num = "π" if a == 1 else f"{a}π"
    return ("−" if c < 0 else "") + (num if b == 1 else f"{num}/{b}")


def _num(x):
    q = Fraction(x).limit_denominator(1000)
    if abs(float(q) - x) < 1e-12:
        return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return f"{x:.6g}"


def _pr(x):
    return (x + 1) % 2 - 1 if (x + 1) % 2 != 0 else 1.0


def _add(cases, name, v, correct):
    if abs(v - correct[name]) > 1e-2 * max(abs(correct[name]), 1e-3) + 1e-6:
        cases.append(({name: _num(v)}, {name: 0}))


def submissions(params, correct):
    inp, fl = params["inp"], params["filt"]
    _, w1, w2 = components(params)
    L1, L2, M = fl["L1"], fl["L2"], fl["M"]
    g1, g2 = abs(Hd(params, w1)), abs(Hd(params, w2))
    cases = []
    for name in ("D", "E", "phi", "Fmag", "Fang") + (("G",) if inp["has_alt"] else ()):
        v = correct[name]
        cases.append(({name: _num(v)}, {name: 1}))
        cases.append(({name: f"{v:.5f}"}, {name: 1}))
    _add(cases, "D", inp["A0"], correct)                                    # gain forgotten
    _add(cases, "D", inp["A0"] * L2, correct)                               # the outer step read
    if g1:
        _add(cases, "E", inp["A1"], correct)
        _add(cases, "E", inp["A1"] * (L1 + L2 - g1), correct)               # the other step
        th = _f(inp["th"])
        _add(cases, "phi", _pr(th + M * w1 / math.pi), correct)             # phase sign flipped
        raw = th - M * w1 / math.pi
        if not -1 < raw <= 1:
            _add(cases, "phi", raw, correct)                                # not reduced
    else:
        _add(cases, "E", inp["A1"] * L2, correct)                           # removed term kept
        _add(cases, "phi", _pr(_f(inp["th"]) - M * w1 / math.pi), correct)
    if g2:
        _add(cases, "Fmag", inp["B"] * (L1 + L2 - g2), correct)
        _add(cases, "Fmag", inp["B"], correct)
        _add(cases, "Fang", _pr(M * w2 / math.pi), correct)
        raw = -M * w2 / math.pi
        if not -1 < raw <= 1:
            _add(cases, "Fang", raw, correct)
    else:
        _add(cases, "Fmag", inp["B"] * L2, correct)
    if inp["has_alt"]:
        _add(cases, "G", inp["C"], correct)
        _add(cases, "G", inp["C"] * L2, correct)
    cases.append(({"D": "abc"}, {"D": "invalid"}))
    return cases
