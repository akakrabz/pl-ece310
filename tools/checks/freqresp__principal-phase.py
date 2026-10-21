"""Independent check of questions/freqresp/principal-phase.

(a), (b): the displayed a + jb is parsed and np.angle (numpy's atan2) gives the principal angle.
(c): H_d(w) = e^{-jMw} A(w) is rebuilt from the displayed LaTeX (lin_tex, A_tex, the frequency
labels) and np.angle(H_d) is compared with the answer. Variant rules: A(w) != 0, no angle of exactly
+-pi in (c), (a) or (b) has Re < 0 off the real axis, (c) has an A < 0 point, an A > 0 point and a
point whose raw angle -Mw (+ pi) needs a 2 pi wrap. The final value printed in each solution row is
compared with the answer. Probes: decimals (4 and 6 digits), pi vs -pi on the negative real axis,
+-2 pi off, the sign, tan^{-1}(b/a) without the quadrant, degrees, multiples of pi, the unwrapped
angle and the angle without the pi of A < 0."""

import math
import re

import numpy as np


def _cplx(tex):
    s = tex.replace(" ", "")
    s = re.sub(r"(?<![\d])j", "1j", s)
    return complex(s)


def _pi(tex):
    t = tex.strip()
    m = re.fullmatch(r"(-?)\\tfrac\{(\d*)\\pi\}\{(\d+)\}", t)
    if m:
        return (-1 if m.group(1) else 1) * int(m.group(2) or 1) / int(m.group(3)) * math.pi
    m = re.fullmatch(r"(-?)(\d*)\\pi", t)
    if m:
        return (-1 if m.group(1) else 1) * int(m.group(2) or 1) * math.pi
    if t == "0":
        return 0.0
    raise ValueError(t)


def _A_fn(tex):
    s = tex.replace(r"\tfrac", r"\frac")
    s = re.sub(r"\\frac\{(\d+)\}\{(\d+)\}", r"(\1/\2)", s)
    s = s.replace(r"\cos(2\omega)", "cos(2*w)").replace(r"\cos\omega", "cos(w)")
    s = re.sub(r"(\d)\s*cos", r"\1*cos", s)
    s = re.sub(r"\)\s*cos", r")*cos", s)
    assert re.fullmatch(r"[\d\s()+\-*/cosw.]+", s), (tex, s)
    return lambda w: eval(s, {"__builtins__": {}}, {"cos": math.cos, "w": w})


def _M(lin_tex):
    m = re.fullmatch(r"e\^\{-j(\d*)\\omega\}", lin_tex)
    assert m, lin_tex
    return int(m.group(1) or 1)


def _vals(params):
    """[(name, true angle, H value, extra)] for all five inputs."""
    out = []
    for name in ("za", "zb"):
        z = _cplx(params[f"{name}_tex"])
        out.append((name, float(np.angle(z)), z, None))
    A = _A_fn(params["A_tex"])
    M = _M(params["lin_tex"])
    for i in (1, 2, 3):
        w = _pi(params[f"f{i}_tex"])
        a = A(w)
        H = np.exp(-1j * M * w) * a
        out.append((f"w{i}", float(np.angle(H)), H, (w, a, M)))
    return out


def check(params, correct):
    probs = []
    neg_a = wrap = pos_a = 0
    for name, ang, H, extra in _vals(params):
        if abs(ang - correct[name]) > 1e-9:
            probs.append(f"{name}: np.angle = {ang}, answer {correct[name]}")
        if extra is None:
            if [int(H.real), int(H.imag)] != params[f"{name}_ab"]:
                probs.append(f"{name}: displayed {params[name + '_tex']} != {params[name + '_ab']}")
            continue
        w, a, M = extra
        if abs(a) < 1e-3:
            probs.append(f"{name}: A({w}) = {a} ~ 0")
        if abs(abs(ang) - math.pi) < 1e-6:
            probs.append(f"{name}: angle is +-pi")
        raw = -M * w + (math.pi if a < 0 else 0)
        neg_a += a < 0
        pos_a += a > 0
        wrap += not (-math.pi < raw <= math.pi + 1e-12)
    zr = [_cplx(params["za_tex"]), _cplx(params["zb_tex"])]
    if not any(z.real < 0 and z.imag != 0 for z in zr):
        probs.append("(a)/(b): no value with Re < 0 off the real axis")
    if not (neg_a and pos_a and wrap):
        probs.append(f"(c): A<0 at {neg_a}, A>0 at {pos_a}, wraps {wrap} points")
    rows = params["zrows"] + params["wrows"]
    for (name, *_), row in zip(_vals(params), rows):
        last = float(row["steps"].rsplit(r"\approx", 1)[1].strip(" )")) if r"\approx" in row["steps"] else math.pi
        if abs(last - correct[name]) > 6e-5:
            probs.append(f"{name}: solution row ends with {last}, answer {correct[name]}")
    if params["wrows"] and any("1\\omega" in r["steps"] for r in params["wrows"]):
        probs.append("display: 1\\omega")
    return probs


def _differs(a, b):
    return abs(a - b) > 1e-2 * max(abs(b), 0.1) + 2e-3


def submissions(params, correct):
    cases = []
    for name, ang, H, extra in _vals(params):
        ans = correct[name]
        cases.append(({name: f"{ans:.4g}"}, {name: 1}))
        cases.append(({name: f"{ans:.6f}"}, {name: 1}))
        if abs(abs(ans) - math.pi) < 1e-9:                       # negative real axis: both spellings
            cases.append(({name: "-3.1416"}, {name: 1}))
            cases.append(({name: "3.1416"}, {name: 1}))
            cases.append(({name: "3.1"}, {name: 0}))
            cases.append(({name: "0"}, {name: 0}))
            continue
        wrong = [ans + 2 * math.pi, ans - 2 * math.pi, -ans, ans * 180 / math.pi, ans / math.pi]
        if extra is None:
            if H.real != 0:
                wrong.append(math.atan(H.imag / H.real))            # quadrant ignored
        else:
            w, a, M = extra
            raw = -M * w + (math.pi if a < 0 else 0)
            wrong.append(raw)                                        # not wrapped
            if a < 0:
                p = -((M * w + math.pi) % (2 * math.pi) - math.pi)     # principal(-Mw)
                wrong.append(p)                                      # forgot the pi of A < 0
        seen = []
        for v in wrong:
            if _differs(v, ans) and all(abs(v - s) > 1e-6 for s in seen):
                seen.append(v)
                cases.append(({name: f"{v:.6f}"}, {name: 0}))
    return cases
