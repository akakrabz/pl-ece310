"""Independent check of questions/freqresp/phase-and-group-delay.

H_d(w) is re-implemented with numpy from the family parameters (scipy.signal.freqz for the IIR).
Checked: |H_d| and np.angle at w_a, w_b (never +-pi; at w_b of a real-factor family the factor is
negative, so the no-jump answer is off by pi); the group delay from a centred difference of the
phase (and scipy.signal.group_delay for the IIR); the pi jumps found by scanning (0, pi) on a grid of
pi/240 for zeros of |H_d| across which the phase changes by about pi, against the checkbox's
correct options (every jump in (0, pi) must be listed; "none" is correct iff there is none)."""

import cmath
import math
import re
from fractions import Fraction

import numpy as np
from scipy import signal

NONE_TEXT = "none of the listed frequencies"


def H(params, w):
    s = params["sys"]
    w = np.atleast_1d(np.asarray(w, float))
    fam = s["fam"]
    if fam == "sin":
        return s["g"] * np.exp(-1j * s["M"] * w) * np.sin(s["K"] * w)
    if fam == "cos":
        return s["g"] * np.exp(-1j * s["M"] * w) * np.cos(s["K"] * w)
    if fam == "even":
        return s["g"] * (float(Fraction(s["b"])) + np.cos(w)) * np.exp(1j * s["s"] * np.abs(w))
    _, h = signal.freqz([float(Fraction(s["g"]))], [1, -float(Fraction(s["a"]))], worN=w)
    return h


def unit_factor(params, w):
    s = params["sys"]
    if s["fam"] in ("sin", "cos"):
        return np.exp(-1j * s["M"] * w)
    if s["fam"] == "even":
        return np.exp(1j * s["s"] * abs(w))
    return None


def tau_num(params, w, eps=1e-6):
    return float(-np.angle(H(params, w + eps)[0] / H(params, w - eps)[0]) / (2 * eps))


def parse_pi(html):
    """'$\\omega = \\tfrac{3\\pi}{4}$' -> 0.75."""
    t = html.split("=", 1)[1].replace("$", "").strip()
    m = re.fullmatch(r"(-?)\\tfrac\{(\d*)\\pi\}\{(\d+)\}", t)
    if m:
        return (-1 if m.group(1) else 1) * int(m.group(2) or 1) / int(m.group(3))
    m = re.fullmatch(r"(-?)(\d*)\\pi", t)
    if m:
        return (-1 if m.group(1) else 1) * int(m.group(2) or 1)
    raise ValueError(t)


def true_jumps(params):
    grid = np.arange(1, 240) * np.pi / 240
    out = []
    for w in grid:
        if abs(H(params, w)[0]) < 1e-9:
            d = np.angle(H(params, w + 1e-6)[0] / H(params, w - 1e-6)[0])
            if abs(abs(d) - np.pi) < 1e-3:
                out.append(w / np.pi)
    return out


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


def check(params, correct):
    probs = []
    negs = []
    for tag in ("a", "b"):
        f = float(Fraction(params[f"f{tag}"]))
        Hv = H(params, np.pi * f)[0]
        if abs(abs(Hv) - correct[f"mag_{tag}"]) > 1e-9 * max(1, abs(Hv)) or abs(Hv) < 1e-9:
            probs.append(f"|H_d({f} pi)| = {abs(Hv)}, answer {correct[f'mag_{tag}']}")
        ph = np.angle(Hv) / np.pi
        if abs(abs(ph) - 1) < 1e-9:
            probs.append(f"phase at {f} pi is +-pi")
        elif abs(ph - correct[f"ph_{tag}"]) > 1e-9:
            probs.append(f"phase at {f} pi is {ph}, answer {correct[f'ph_{tag}']}")
        u = unit_factor(params, np.pi * f)
        if u is not None:
            R = (Hv / u[0]).real if np.ndim(u) else (Hv / u).real      # the real factor at w_tag
            negs.append(R < 0)
            row = params["rows"][0 if tag == "a" else 1]
            try:
                if abs(tex_eval(row["R_val"]).real - R) > 1e-9 * max(1, abs(R)):
                    probs.append(f"printed R = {row['R_val']!r} at w_{tag}, true {R}")
            except Exception as exc:
                probs.append(f"cannot evaluate R = {row['R_val']!r}: {exc}")
        row = params["rows"][0 if tag == "a" else 1]
        if abs(tex_eval(row["mag_tex"]).real - abs(Hv)) > 1e-9 * max(1, abs(Hv)):
            probs.append(f"printed |H_d| = {row['mag_tex']!r} at w_{tag}, true {abs(Hv)}")
        if row.get("line"):                       # IIR: every step after the label is a value of H_d
            for seg in row["line"].split(" = ")[1:]:
                got = tex_eval(seg)
                if abs(got - Hv) > 1e-9 * max(1, abs(Hv)):
                    probs.append(f"step {seg!r} of {row['line']!r} is {got}, H_d is {Hv}")
    if params["sys"]["fam"] != "iir" and sorted(negs) != [False, True]:
        probs.append(f"exactly one of w_a, w_b must have a negative real factor: {negs}")
    fc = np.pi * float(Fraction(params["fc"]))
    lo, hi = np.angle(H(params, fc - 1e-6)[0]), np.angle(H(params, fc + 1e-6)[0])
    if abs(hi - lo) > 1e-3:
        probs.append(f"the principal phase breaks at w_c = {params['fc']} pi (jump or 2 pi wrap)")
    t = tau_num(params, fc)
    if abs(t - correct["tau"]) > 1e-5:
        probs.append(f"group delay at {params['fc']} pi is {t}, answer {correct['tau']}")
    if params["sys"]["fam"] == "iir":
        s = params["sys"]
        _, gd = signal.group_delay(([float(Fraction(s["g"]))], [1, -float(Fraction(s["a"]))]), w=[fc])
        if abs(gd[0] - correct["tau"]) > 1e-7:
            probs.append(f"scipy group_delay {gd[0]} != {correct['tau']}")
    jumps = true_jumps(params)
    opts = params["jumps"]
    marked = [o["html"] for o in correct["jumps"]]
    listed = [o["html"] for o in opts if o["html"].strip() != NONE_TEXT]
    lf = [parse_pi(x) for x in listed]
    for j in jumps:
        if not any(abs(j - x) < 1e-6 for x in lf):
            probs.append(f"jump at {j} pi is not listed")
    for x, html in zip(lf, listed):
        is_jump = any(abs(j - x) < 1e-6 for j in jumps)
        if is_jump != (html in marked):
            probs.append(f"option {html}: jump={is_jump} but marked {html in marked}")
    none_marked = NONE_TEXT in [m.strip() for m in marked]
    if none_marked != (not jumps):
        probs.append(f"'none' marked {none_marked} but jumps {jumps}")
    if len(listed) != 4:
        probs.append(f"{len(listed)} frequencies listed, expected 4")
    return probs


def _num(x):
    q = Fraction(x).limit_denominator(1000)
    if abs(float(q) - x) < 1e-12:
        return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return f"{x:.6g}"


def _differs(a, b):
    return abs(a - b) > 1e-2 * max(abs(b), 1e-3) + 1e-6


def submissions(params, correct):
    cases = []
    for tag in ("a", "b"):
        m, ph = correct[f"mag_{tag}"], correct[f"ph_{tag}"]
        cases.append(({f"ph_{tag}": f"{ph:.4f}"}, {f"ph_{tag}": 1}))          # decimal multiple of pi
        cases.append(({f"ph_{tag}": _num(ph)}, {f"ph_{tag}": 1}))             # fraction
        cases.append(({f"mag_{tag}": f"{m:.5g}"}, {f"mag_{tag}": 1}))
        other = ph - 2 if ph > 0 else ph + 2
        cases.append(({f"ph_{tag}": _num(other)}, {f"ph_{tag}": 0}))          # not the principal value
        if _differs(-ph, ph):
            cases.append(({f"ph_{tag}": _num(-ph)}, {f"ph_{tag}": 0}))        # sign slip
    if params["sys"]["fam"] != "iir":                                         # forgot the +-pi of R < 0
        for tag in ("a", "b"):
            w = np.pi * float(Fraction(params[f"f{tag}"]))
            u = unit_factor(params, w)
            if (H(params, w)[0] / u).real < 0:
                nojump = float(np.angle(u) / np.pi)
                if abs(nojump + 1) < 1e-12:
                    nojump = 1.0
                cases.append(({f"ph_{tag}": _num(nojump)}, {f"ph_{tag}": 0}))
    tau = correct["tau"]
    cases.append(({"tau": _num(tau)}, {"tau": 1}))
    if _differs(-tau, tau):
        cases.append(({"tau": _num(-tau)}, {"tau": 0}))                        # wrong sign
    # checkbox (net-correct: max(0, (right - wrong) / number correct))
    opts = params["jumps"]
    key = {o["html"].strip(): o["key"] for o in opts}
    right = [o["key"] for o in correct["jumps"]]
    cases.append(({"jumps": right}, {"jumps": 1}))
    c = len(right)
    wraps = [o["key"] for o in opts if o["html"].strip() != NONE_TEXT and o["key"] not in right]
    if wraps:
        extra = right + wraps[:1]
        cases.append(({"jumps": extra}, {"jumps": max(0.0, (c - 1) / c)}))
    none_key = key[NONE_TEXT]
    if none_key not in right:
        cases.append(({"jumps": [none_key]}, {"jumps": 0}))
        if c >= 2:
            cases.append(({"jumps": right[:1]}, {"jumps": 1 / c}))
    else:
        cases.append(({"jumps": [k for k in key.values() if k != none_key][:2]}, {"jumps": 0}))
    return cases
