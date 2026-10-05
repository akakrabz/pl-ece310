"""Independent check of questions/signals/complex-numbers.

The displayed LaTeX (w_tex, c_tex) is converted to a Python expression and evaluated with cmath, so
the check covers exactly what the student sees; roots come from numpy.roots."""

import cmath
import math
import re

import numpy as np


def _match_brace(s, i):
    """s[i] == '{' -> index of the matching '}'."""
    depth = 0
    for k in range(i, len(s)):
        if s[k] == "{":
            depth += 1
        elif s[k] == "}":
            depth -= 1
            if depth == 0:
                return k
    raise ValueError(f"unbalanced braces in {s!r}")


def tex_to_py(s):
    s = s.replace(r"\left(", "(").replace(r"\right)", ")").replace(r"\cdot", "*").replace(" ", "")
    s = s.replace(r"\tfrac{", r"\dfrac{")
    while r"\dfrac{" in s:                                   # \dfrac{A}{B} -> ((A)/(B))
        i = s.index(r"\dfrac{")
        j = _match_brace(s, i + 6)
        k = _match_brace(s, j + 1)
        s = s[:i] + "((" + s[i + 7:j] + ")/(" + s[j + 2:k] + "))" + s[k + 1:]
    s = re.sub(r"e\^\{(-?)j(\d*)\\pi(?:/(\d+))?\}",
               lambda m: f"exp({m.group(1)}1j*{m.group(2) or 1}*pi/{m.group(3) or 1})", s)
    s = s.replace(r"\pi", "pi")
    s = re.sub(r"\\sqrt\[(\d+)\]\{(\d+)\}", r"(\2**(1/\1))", s)
    s = re.sub(r"\\sqrt\{(\d+)\}", r"sqrt(\1)", s)
    s = re.sub(r"\^\{(-?\d+)\}", r"**(\1)", s)
    s = re.sub(r"(?<![\w.])j(?=\d|s|\()", "1j*", s)       # j4, j sqrt(3)
    s = re.sub(r"(?<![\w.])j", "1j", s)                     # lone j
    s = re.sub(r"(\d)(?=[sep(])", r"\1*", s)                # 2exp(, 4sqrt(, 3pi
    s = re.sub(r"\)(?=[(\w])", ")*", s)                     # )( , )exp
    return s


def ev(s):
    return complex(eval(tex_to_py(s), {"__builtins__": {}}, {"exp": cmath.exp, "sqrt": math.sqrt, "pi": math.pi}))


def principal_over_pi(z):
    a = cmath.phase(z) / math.pi
    return 1.0 if a <= -1 + 1e-12 else a


def _close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(b))


def check(params, correct):
    probs = []
    kind = params["kind"]
    if kind == "polar":
        w = ev(params["w_tex"])
        if not _close(abs(w), correct["mag"]):
            probs.append(f"|w| {correct['mag']} != {abs(w)} (from {tex_to_py(params['w_tex'])})")
        if not _close(principal_over_pi(w), correct["ang"]):
            probs.append(f"angle {correct['ang']} != {principal_over_pi(w)}")
        if abs(abs(w) - 1) < 1e-9 or abs(correct["ang"]) < 1e-9:
            probs.append("degenerate variant (|w| = 1 or angle 0)")
        if not (-1 < correct["ang"] <= 1):
            probs.append("angle not principal")
        if not _close(ev(params["mag_tex"]).real, correct["mag"]):
            probs.append(f"mag_tex {params['mag_tex']} != {correct['mag']}")
        if not _close(ev(params["ang_tex"]).real / math.pi, correct["ang"]):
            probs.append(f"ang_tex {params['ang_tex']} != {correct['ang']}")
    elif kind == "roots":
        c = ev(params["c_tex"])
        N = params["N"]
        roots = np.roots([1] + [0] * (N - 1) + [-c])
        mags = np.abs(roots)
        if not np.allclose(mags, mags[0]) or not _close(float(mags[0]), correct["mag"]):
            probs.append(f"root magnitudes {mags} vs {correct['mag']}")
        angs = sorted(principal_over_pi(complex(r)) for r in roots)
        sel = params["sel"]
        if sel == "minpos":
            want = min(a for a in angs if a > 1e-12)
        elif sel == "max":
            want = max(angs)
        else:
            want = max(a for a in angs if a < -1e-12)
        if not _close(want, correct["ang"], 1e-7):
            probs.append(f"selected root angle {correct['ang']} != {want} (all {angs})")
        polar = ev(params["c_polar_tex"]) if params["c_polar_tex"] else c
        if abs(polar - c) > 1e-9 * max(1, abs(c)):
            probs.append(f"c_tex {params['c_tex']} != c_polar_tex {params['c_polar_tex']}")
        if not _close(abs(c) ** (1 / N), correct["mag"]):
            probs.append("magnitude is not R^(1/N)")
    else:
        w = ev(params["w_tex"])
        if not _close(w.real, correct["re"]) or not _close(w.imag, correct["im"]):
            probs.append(f"Re/Im ({correct['re']}, {correct['im']}) != {w}")
        if not _close(ev(params["re_tex"]).real, correct["re"]) or not _close(ev(params["im_tex"]).real, correct["im"]):
            probs.append(f"re_tex/im_tex {params['re_tex']}, {params['im_tex']} disagree with the answers")
        for k in ("re", "im"):
            if abs(correct[k]) < 1e-9:
                probs.append(f"{k} is zero (degenerate)")
    return probs


def _num(x):
    return f"{x:.6g}"


def submissions(params, correct):
    kind = params["kind"]
    cases = []
    if kind in ("polar", "roots"):
        mag, ang = correct["mag"], correct["ang"]
        cases += [
            ({"ang": params["ang_frac"]}, {"ang": 1}),                          # exact fraction
            ({"ang": f"{ang:.5f}"}, {"ang": 1}),                                # decimal
            ({"mag": f"{mag:.5g}"}, {"mag": 1}),                                # 5 significant digits
            ({"ang": f"{ang:.4g}"}, {"ang": 1}),                                # 4 significant digits
            ({"ang": _num(ang * 180)}, {"ang": 0}),                             # degrees instead of multiples of pi
            ({"ang": _num(ang - 2 if ang > 0 else ang + 2)}, {"ang": 0}),       # off by 2 pi (not principal)
            ({"ang": _num(ang * math.pi)}, {"ang": 0}),                         # radians, not multiples of pi
        ]
        if abs(mag - round(mag)) > 1e-9:
            cases.append(({"mag": f"{mag:.4g}"}, {"mag": 1}))                    # 4 significant digits
        if abs(mag * mag - mag) > 1e-2 * mag:
            cases.append(({"mag": _num(mag * mag)}, {"mag": 0}))                 # squared magnitude
        if kind == "polar":
            raw = params["raw_ang"][0] / params["raw_ang"][1]
            if abs(raw - ang) > 1e-9:
                cases.append(({"ang": _num(raw)}, {"ang": 0}))                  # unreduced angle
        else:
            th = params["theta"][0] / params["theta"][1]
            if abs(th - ang) > 1e-6:
                cases.append(({"ang": _num(th)}, {"ang": 0}))                   # forgot to divide by N
    else:
        re_, im_ = correct["re"], correct["im"]
        cases += [
            ({"re": _num(re_), "im": _num(im_)}, {"re": 1, "im": 1}),            # 6-digit decimals
            ({"re": f"{re_:.5g}", "im": f"{im_:.5g}"}, {"re": 1, "im": 1}),      # 5 significant digits
            ({"re": f"{re_:.4g}", "im": f"{im_:.4g}"}, {"re": 1, "im": 1}),      # 4 significant digits
            ({"re": _num(im_), "im": _num(re_)}, {"re": 0, "im": 0}),            # swapped
        ]
        frac = re.fullmatch(r"(-?)\\tfrac\{(\d+)\}\{(\d+)\}|(-?\d+)", params["re_tex"])
        if frac:
            txt = frac.group(4) if frac.group(4) else f"{frac.group(1)}{frac.group(2)}/{frac.group(3)}"
            cases.append(({"re": txt}, {"re": 1}))                              # exact fraction
        # conjugate slip: e^{-j th} read as e^{+j th}
        flipped = sum(A * cmath.exp(1j * math.pi * abs(p / q)) for A, p, q in params["terms"])
        if abs(flipped.imag - im_) > 1e-3 * max(1, abs(im_)):
            cases.append(({"im": _num(flipped.imag)}, {"im": 0}))
    return cases
