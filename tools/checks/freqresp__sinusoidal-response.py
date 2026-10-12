"""Independent check of questions/freqresp/sinusoidal-response.

H_d is re-implemented with numpy from the system description (closed forms evaluated after reducing
the frequency into (-pi, pi]; FIR by the DTFT sum; LCCDE by scipy.signal.freqz). Checked: the
Hermitian test H_d(-w) = H_d*(w) on a grid against the Yes/No answer; D = c H_d(0); B, C+ and C-
from the eigenfunction rule with numpy (principal angles, nonzero magnitudes); for FIR and LCCDE
systems also in the time domain: the input is built sample by sample and convolved with h[n]
(or run through the recursion from n = -400), and the output must equal
D + B E[n] + C+ e^{j w2 n} + C- e^{-j w2 n} built from the answers."""

import cmath
import math
import re
from fractions import Fraction

import numpy as np
from scipy import signal


def _f(s):
    return float(Fraction(s))


def Hd(spec, w, reduce=True):
    w = np.atleast_1d(np.asarray(w, float))
    wr = np.angle(np.exp(1j * w)) if reduce else w
    k = spec["kind"]
    if k == "absw":
        return np.abs(wr) * np.exp(1j * spec["s"] * np.pi * np.sin(wr))
    if k == "jw":
        return 1j * wr * np.exp(1j * spec["s"] * np.pi * np.sin(wr))
    if k == "wcos":
        return wr * np.exp(1j * spec["s"] * np.pi * np.cos(wr))
    if k == "cosdelay":
        return (1 + np.cos(wr)) * np.exp(-1j * spec["k"] * wr)
    if k == "halfcos":
        return (_f(spec["b"]) + np.cos(wr)) * np.exp(1j * spec["s"] * np.abs(wr))
    if k == "fir":
        n = np.arange(spec["start"], spec["start"] + len(spec["h"]))
        return np.exp(-1j * np.outer(w, n)) @ np.asarray(spec["h"], float)
    g, q, a = _f(spec["g"]), spec["q"], _f(spec["a"])
    _, h = signal.freqz([g, g * q], [1, -a], worN=w)
    return h


def coeffs(params):
    """(D, B, C+, C-) recomputed with numpy."""
    spec, inp = params["spec"], params["inp"]
    r1, t1 = inp["a1"][0], _f(inp["a1"][1])
    a1 = r1 * np.exp(1j * np.pi * t1)
    D = inp["c"] * Hd(spec, 0.0)[0]
    B = a1 * Hd(spec, np.pi * _f(inp["fw"]))[0]
    A, th, w2 = inp["A"], np.pi * _f(inp["theta"]), np.pi * _f(inp["f2"])
    if inp["trig"] == "cos":
        cp, cm = A / 2 * np.exp(1j * th), A / 2 * np.exp(-1j * th)
    else:
        cp, cm = A / 2j * np.exp(1j * th), -A / 2j * np.exp(-1j * th)
    return D, B, cp * Hd(spec, w2)[0], cm * Hd(spec, -w2)[0], a1, cp, cm


def is_real(spec):
    w = np.linspace(-np.pi + 1e-3, np.pi - 1e-3, 1001)
    return bool(np.max(np.abs(Hd(spec, -w) - np.conj(Hd(spec, w)))) < 1e-9)


def ang(z):
    a = np.angle(z) / np.pi
    return 1.0 if a <= -1 + 1e-12 else float(a)


def same_angle(a, b, tol=1e-9):
    return abs(a - b) <= tol or (abs(abs(a) - 1) <= tol and abs(abs(b) - 1) <= tol)


def E_of_n(params, n):
    inp = params["inp"]
    return np.exp(1j * np.pi * _f(inp["fw"]) * n)


def x_of_n(params, n):
    inp = params["inp"]
    r1, t1 = inp["a1"][0], _f(inp["a1"][1])
    a1 = r1 * np.exp(1j * np.pi * t1)
    w2, th = np.pi * _f(inp["f2"]), np.pi * _f(inp["theta"])
    trig = np.cos if inp["trig"] == "cos" else np.sin
    return inp["c"] + a1 * E_of_n(params, n) + inp["A"] * trig(w2 * n + th)


# ---------------------------------------------------------------- evaluate printed LaTeX numbers
def _match_brace(s, i):
    depth = 0
    for k in range(i, len(s)):
        if s[k] == "{":
            depth += 1
        elif s[k] == "}":
            depth -= 1
            if depth == 0:
                return k
    raise ValueError(f"unbalanced braces in {s!r}")


def _sub_cmd(s, cmd, nargs, fn):
    while cmd + "{" in s:
        i = s.index(cmd + "{")
        j = i + len(cmd)
        args = []
        for _ in range(nargs):
            k = _match_brace(s, j)
            args.append(s[j + 1:k])
            j = k + 1
        s = s[:i] + fn(*args) + s[j:]
    return s


def tex_eval(s):
    """Numeric value of a printed constant: \\tfrac, \\sqrt, \\pi, e^{j...}, j, products, sums."""
    s = (s.replace(r"\left", "").replace(r"\right", "").replace(r"\,", " ").replace(r"\!", "")
         .replace(r"\cdot", "*").replace(r"\dfrac", r"\frac").replace(r"\tfrac", r"\frac"))
    s = _sub_cmd(s, r"\frac", 2, lambda a, b: f"(({a})/({b}))")
    s = _sub_cmd(s, r"\sqrt", 1, lambda a: f"sqrt({a})")
    s = _sub_cmd(s, "e^", 1, lambda a: f"exp({a})")
    s = s.replace(r"\pi", " pi ").replace("{", "(").replace("}", ")")
    toks = re.findall(r"\d+\.\d+|\d+|[A-Za-z]+|[()+\-*/^]", s)
    out = []
    for t in toks:
        if out:
            prev = out[-1]
            if ((re.fullmatch(r"[\d.]+", prev) or prev in ("pi", "1j", ")"))
                    and (re.fullmatch(r"[\d.]+", t) or t in ("pi", "j", "exp", "sqrt", "("))):
                out.append("*")
        out.append("1j" if t == "j" else ("**" if t == "^" else t))
    return complex(eval("".join(out), {"__builtins__": {}}, {"exp": cmath.exp, "sqrt": cmath.sqrt, "pi": math.pi}))


def _printed(params, correct, probs):
    """Every printed value in the worked solution equals the number it claims to be."""
    spec, inp = params["spec"], params["inp"]
    D, B, Cp, Cm, *_ = coeffs(params)
    want = {"B_tex": B, "Cp_tex": Cp, "Cm_tex": Cm, "D_tex": D,
            "H0_tex": Hd(spec, 0.0)[0], "H1_tex": Hd(spec, np.pi * _f(inp["fw"]))[0],
            "Hp_tex": Hd(spec, np.pi * _f(inp["f2"]))[0], "Hm_tex": Hd(spec, -np.pi * _f(inp["f2"]))[0],
            "a1_polar_tex": coeffs(params)[4], "Cp_in_tex": coeffs(params)[5], "Cm_in_tex": coeffs(params)[6]}
    for key, z in want.items():
        try:
            got = tex_eval(params[key])
        except Exception as exc:
            probs.append(f"cannot evaluate {key} = {params[key]!r}: {exc}")
            continue
        if abs(got - z) > 1e-9 * max(1, abs(z)):
            probs.append(f"printed {key} = {params[key]!r} is {got}, true value {z}")
    for name in ("B", "Cp", "Cm"):
        if abs(tex_eval(params[f"{name}_mag_tex"]) - correct[f"{name}_mag"]) > 1e-9 * max(1, correct[f"{name}_mag"]):
            probs.append(f"printed |{name}| = {params[name + '_mag_tex']!r} differs from {correct[name + '_mag']}")
    freqs = [0.0, _f(inp["fw"]), _f(inp["f2"]), -_f(inp["f2"])]
    for row, f in zip(params["ev_rows"], freqs):
        segs = row["line"].split(" = ")
        z = Hd(spec, np.pi * f)[0]
        for seg in segs[1:]:
            if any(k in seg for k in (r"\sin", r"\cos", "|", r"\omega")):
                continue                     # a substitution step, not a number
            try:
                got = tex_eval(seg)
            except Exception as exc:
                probs.append(f"cannot evaluate a step of {row['line']!r}: {exc}")
                continue
            if abs(got - z) > 1e-9 * max(1, abs(z)):
                probs.append(f"step {seg!r} of {row['line']!r} is {got}, H_d is {z}")


def check(params, correct):
    probs = []
    spec = params["spec"]
    _printed(params, correct, probs)
    D, B, Cp, Cm, *_ = coeffs(params)
    real = is_real(spec)
    marked = correct["real"]["html"].strip()
    if marked != ("Yes" if real else "No"):
        probs.append(f"h real? numpy says {real}, marked {marked}")
    if abs(D.imag) > 1e-9 or abs(D.real - correct["D"]) > 1e-9 * max(1, abs(D)):
        probs.append(f"D = {D}, answer {correct['D']}")
    for name, z in (("B", B), ("Cp", Cp), ("Cm", Cm)):
        if abs(z) < 1e-6:
            probs.append(f"{name} is zero (angle undefined)")
            continue
        if abs(abs(z) - correct[f"{name}_mag"]) > 1e-9 * max(1, abs(z)):
            probs.append(f"|{name}| = {abs(z)}, answer {correct[f'{name}_mag']}")
        if not same_angle(ang(z), correct[f"{name}_ang"]) or not (-1 < correct[f"{name}_ang"] <= 1):
            probs.append(f"angle {name} = {ang(z)}, answer {correct[f'{name}_ang']}")
        q = Fraction(correct[f"{name}_ang"]).limit_denominator(24)
        if abs(float(q) - correct[f"{name}_ang"]) > 1e-12:
            probs.append(f"angle {name} = {correct[f'{name}_ang']} is not a simple multiple of pi")
    w2 = _f(params["inp"]["f2"])
    if not 0 < w2 < 1:
        probs.append("sinusoid frequency outside (0, pi)")
    f1 = _f(params["inp"]["f1"])
    if abs(f1) < 1e-12 or abs(abs(f1) - w2) < 1e-12:
        probs.append("exponential shares a frequency with another term")
    if spec["kind"] in ("jw", "wcos") and abs(abs(f1) - 1) < 1e-12:
        probs.append("exponential at pi for an H_d that jumps at +-pi")
    # time domain for FIR and LCCDE systems
    if spec["kind"] in ("fir", "lccde"):
        n_out = np.arange(0, 12)
        model = (correct["D"] + correct["B_mag"] * np.exp(1j * np.pi * correct["B_ang"]) * E_of_n(params, n_out)
                 + correct["Cp_mag"] * np.exp(1j * np.pi * correct["Cp_ang"]) * np.exp(1j * np.pi * w2 * n_out)
                 + correct["Cm_mag"] * np.exp(1j * np.pi * correct["Cm_ang"]) * np.exp(-1j * np.pi * w2 * n_out))
        if spec["kind"] == "fir":
            h, s0 = spec["h"], spec["start"]
            y = np.array([sum(h[k] * x_of_n(params, n - (s0 + k)) for k in range(len(h))) for n in n_out])
        else:
            g, q, a = _f(spec["g"]), spec["q"], _f(spec["a"])
            n_all = np.arange(-400, 12)
            yy = signal.lfilter([g, g * q], [1, -a], x_of_n(params, n_all))
            y = yy[-12:]
        if np.max(np.abs(y - model)) > 1e-8:
            probs.append(f"time-domain output differs from the answer by {np.max(np.abs(y - model))}")
    return probs


def _num(x):
    q = Fraction(x).limit_denominator(1000)
    if abs(float(q) - x) < 1e-12:
        return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return f"{x:.6g}"


def _differs(a, b):
    return abs(a - b) > 1e-2 * max(abs(b), 1e-3) + 1e-6


def _wrong(name, z, correct):
    """Overrides for a wrong coefficient z: the inputs (mag / angle) that would be marked wrong."""
    o, e = {}, {}
    if _differs(abs(z), correct[f"{name}_mag"]):
        o[f"{name}_mag"], e[f"{name}_mag"] = f"{abs(z):.6g}", 0
    a = ang(z)
    if abs(z) > 1e-9 and not same_angle(a, correct[f"{name}_ang"], 1e-2):
        o[f"{name}_ang"], e[f"{name}_ang"] = _num(a), 0
    return o, e


def submissions(params, correct):
    spec, inp = params["spec"], params["inp"]
    D, B, Cp, Cm, a1, cp, cm = coeffs(params)
    cases = []
    for name in ("B", "Cp", "Cm"):
        m, a = correct[f"{name}_mag"], correct[f"{name}_ang"]
        cases.append(({f"{name}_ang": _num(a)}, {f"{name}_ang": 1}))               # fraction
        cases.append(({f"{name}_ang": f"{a:.4f}"}, {f"{name}_ang": 1}))           # decimal
        cases.append(({f"{name}_mag": f"{m:.4g}"}, {f"{name}_mag": 1}))           # 4 significant digits
        if abs(abs(a) - 1) < 1e-12:
            cases.append(({f"{name}_ang": _num(-a)}, {f"{name}_ang": 1}))         # pi = -pi
        other = a - 2 if a > 0 else a + 2
        if not same_angle(other, a):
            cases.append(({f"{name}_ang": _num(other)}, {f"{name}_ang": 0}))      # not the principal value
    cases.append(({"D": _num(correct["D"])}, {"D": 1}))
    if _differs(float(inp["c"]), correct["D"]):
        cases.append(({"D": str(inp["c"])}, {"D": 0}))                          # forgot H_d(0)
    key = {o["html"].strip(): o["key"] for o in params["real"]}
    real = correct["real"]["html"].strip() == "Yes"
    cases.append(({"real": key["No" if real else "Yes"]}, {"real": 0}))
    w2 = np.pi * _f(inp["f2"])
    # C- with H_d(+w2) for both halves
    o, e = _wrong("Cm", cm * Hd(spec, w2)[0], correct)
    if o:
        cases.append((o, e))
    # the real-h shortcut C- = conj(C+) (wrong when h is complex)
    o, e = _wrong("Cm", np.conj(Cp), correct)
    if o:
        cases.append((o, e))
    # input coefficient of the -w2 half not conjugated
    o, e = _wrong("Cm", cp * Hd(spec, -w2)[0], correct)
    if o:
        cases.append((o, e))
    # ignoring periodicity: the formula evaluated at the written frequency outside [-pi, pi]
    fw = _f(inp["fw"])
    if spec["kind"] in ("absw", "jw", "wcos", "halfcos") and abs(fw) > 1:
        o, e = _wrong("B", a1 * Hd(spec, np.pi * fw, reduce=False)[0], correct)
        if o:
            cases.append((o, e))
    return cases
