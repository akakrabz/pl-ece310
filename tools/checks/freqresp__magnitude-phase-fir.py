"""Independent check of questions/freqresp/magnitude-phase-fir.

H_d(w) is evaluated with numpy straight from h[n] (DTFT sum over the true indices). Checked:
H_d(w) e^{j alpha w} (divided by j for antisymmetric h) is real on a grid and equals the A(w) of the
generated SymPy strings; |H_d| at the two frequencies; the principal phase (np.angle) at w1, where
the amplitude must be negative (so the no-jump answer is wrong) and the phase not +-pi; the group
delay from a centred finite difference of the phase at several frequencies; the first zero of |H_d|
on (0, pi] from a dense grid."""

import cmath
import math
import re
from fractions import Fraction

import numpy as np
from checklib import sympy

W = sympy.Symbol("w")


def _H(params, w):
    h = np.asarray(params["h"], float)
    n = np.arange(params["start"], params["start"] + len(h))
    w = np.atleast_1d(np.asarray(w, float))
    return np.exp(-1j * np.outer(w, n)) @ h


def _alpha(correct):
    return correct["alpha"]


def _eval(s, w):
    expr = sympy.sympify(s.replace("^", "**"), locals={"j": sympy.I, "w": W, "e": sympy.E})
    f = sympy.lambdify(W, expr, "numpy")
    return np.asarray(f(np.asarray(w, float)) * np.ones_like(w), complex)


def _amp(params, correct, w):
    """H_d e^{j alpha w} (/ j): must be real."""
    v = _H(params, w) * np.exp(1j * correct["alpha"] * np.asarray(w, float))
    return v / 1j if params["kind"] == "anti" else v


def _f(s):
    return float(Fraction(s))


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
    h = params["h"]
    kind = "sym" if h == h[::-1] else ("anti" if h == [-v for v in h[::-1]] else "none")
    if kind != params["kind"]:
        probs.append(f"h = {h} is {kind}, params say {params['kind']}")
    alpha = params["start"] + (len(h) - 1) / 2
    if abs(alpha - correct["alpha"]) > 1e-12 or abs(correct["tau"] - alpha) > 1e-12:
        probs.append(f"alpha/tau {correct['alpha']}/{correct['tau']} but the centre of h is {alpha}")
    grid = np.linspace(-np.pi, np.pi, 721)
    amp = _amp(params, correct, grid)
    if np.max(np.abs(amp.imag)) > 1e-9:
        probs.append("H_d e^{j alpha w} (/j) is not real: wrong alpha or form")
    for key in ("A_sym", "A_exp_sym"):
        got = _eval(params[key], grid)
        if np.max(np.abs(got - amp.real)) > 1e-9:
            probs.append(f"{key} = {params[key]!r} differs from the true amplitude")
    if not np.min(amp.real) < -1e-6:
        probs.append("A(w) is never negative: no pi to add")
    for name, f in (("mag_a", params["fa"]), ("mag_b", params["fb"])):
        m = float(np.abs(_H(params, math.pi * _f(f)))[0])
        if abs(m - correct[name]) > 1e-9 * max(1, m) or m < 1e-9:
            probs.append(f"{name}: |H_d({f} pi)| = {m}, answer {correct[name]}")
    f1 = math.pi * _f(params["f1"])
    H1 = _H(params, f1)[0]
    ph = np.angle(H1) / np.pi
    if abs(abs(ph) - 1) < 1e-9:
        probs.append("phase at w1 is +-pi (two spellings)")
    if abs(ph - correct["phase"]) > 1e-9:
        probs.append(f"phase at {params['f1']} pi is {ph}, answer {correct['phase']}")
    A1 = _amp(params, correct, [f1])[0].real
    if abs(A1) < 1e-9 or (A1 < 0) != bool(params["A1_neg"]):
        probs.append(f"A(w1) = {A1} but A1_neg = {params['A1_neg']}")
    # printed values in the solution
    for row, f in zip(params["mag_rows"], (params["fa"], params["fb"])):
        Av = _amp(params, correct, [math.pi * _f(f)])[0].real
        for key, want in (("A_val_tex", Av), ("mag_tex", abs(Av))):
            t = row[key]
            got = tex_eval(t[2:]).real if t.startswith("= ") else float(t.replace(r"\approx", "").strip())
            tol = 1e-9 if t.startswith("= ") else 6e-4 * max(1, abs(want))
            if abs(got - want) > tol * max(1, abs(want)):
                probs.append(f"printed {key} {t!r} at {f} pi, true {want}")
    t = params["A1_tex"]
    got = tex_eval(t[2:]).real if t.startswith("= ") else float(t.replace(r"\approx", "").strip())
    if abs(got - A1) > 6e-4 * max(1, abs(A1)):
        probs.append(f"printed A(w1) {t!r}, true {A1}")
    # group delay by a centred difference of the phase, away from zeros
    for w0 in (0.31, 1.07, 2.39, -0.77):
        if abs(_H(params, w0)[0]) < 1e-3:
            continue
        eps = 1e-6
        tau = -np.angle(_H(params, w0 + eps)[0] / _H(params, w0 - eps)[0]) / (2 * eps)
        if abs(tau - correct["tau"]) > 1e-5:
            probs.append(f"group delay at {w0} is {tau}, answer {correct['tau']}")
    # first zero on (0, pi]
    fine = np.linspace(1e-4, np.pi, 24001)
    mag = np.abs(_H(params, fine))
    idx = [i for i in range(1, len(fine) - 1) if mag[i] <= mag[i - 1] and mag[i] <= mag[i + 1] and mag[i] < 1e-3]
    if mag[-1] < 1e-9:
        idx.append(len(fine) - 1)
    if not idx:
        probs.append("no zero of |H_d| on (0, pi]")
    else:
        w_first = fine[idx[0]] / np.pi
        if abs(w_first - correct["wzero"]) > 1e-3:
            probs.append(f"first zero at {w_first} pi, answer {correct['wzero']}")
        if abs(_H(params, math.pi * correct["wzero"])[0]) > 1e-9:
            probs.append("|H_d| is not zero at the answer of (f)")
    return probs


def _differs(a, b):
    return abs(a - b) > 1e-2 * max(abs(b), 1e-3) + 1e-6


def _num(x):
    q = Fraction(x).limit_denominator(1000)
    if abs(float(q) - x) < 1e-12:
        return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return f"{x:.6g}"


def submissions(params, correct):
    alpha, ph, wz = correct["alpha"], correct["phase"], correct["wzero"]
    terms = [(Fraction(v), Fraction(k)) for v, k in params["terms"]]
    fn = "cos" if params["kind"] == "sym" else "sin"
    rev = " + ".join(f"{fn}({k}*w)*({v})" if k else f"({v})" for v, k in reversed(terms))
    cases = [
        ({"A": params["A_exp_sym"]}, {"A": 1}),                  # exponential (Euler) form
        ({"A": rev}, {"A": 1}),                                  # reordered, coefficient last
        ({"A": "-(" + params["A_sym"] + ")"}, {"A": 0}),         # sign of A flipped
        ({"A": "x + 1"}, {"A": "invalid"}),
        ({"alpha": f"{alpha:.4f}"}, {"alpha": 1}),
        ({"phase": f"{ph:.4f}"}, {"phase": 1}),                  # decimal multiple of pi
        ({"phase": _num(ph)}, {"phase": 1}),                     # exact fraction
        ({"wzero": f"{wz:.4f}"}, {"wzero": 1}),
        ({"tau": _num(alpha)}, {"tau": 1}),
    ]
    if any(k != 0 for _, k in terms):                            # forgot the 2 of 2cos / 2j sin
        half = " + ".join(f"({v / 2})*{fn}({k}*w)" if k else f"({v})" for v, k in terms)
        cases.append(({"A": half}, {"A": 0}))
    if params["kind"] == "anti":                                 # put the j inside A as well
        cases.append(({"A": "j*(" + params["A_sym"] + ")"}, {"A": 0}))
    if alpha != 0:
        cases.append(({"alpha": _num(-alpha)}, {"alpha": 0}))    # wrong sign of the delay
        cases.append(({"tau": _num(-alpha)}, {"tau": 0}))
        cases.append(({"A": f"exp(-j*({_num(alpha)})*w)*(" + params["A_sym"] + ")"}, {"A": 0}))   # kept e^{-j alpha w}
    nojump = Fraction(ph).limit_denominator(100) - 1
    nojump = nojump + 2 if nojump <= -1 else nojump
    if _differs(float(nojump), ph):
        cases.append(({"phase": _num(float(nojump))}, {"phase": 0}))   # the other decision about +-pi (forgot it / added it)
    other = ph - 2 if ph > 0 else ph + 2
    cases.append(({"phase": _num(other)}, {"phase": 0}))          # not the principal value
    for name in ("mag_a", "mag_b"):
        m = correct[name]
        if abs(m - round(m)) > 1e-9:
            cases.append(({name: f"{m:.4g}"}, {name: 1}))
    zs = [float(Fraction(z)) for z in params["zeros"]]
    if len(zs) > 1:
        cases.append(({"wzero": _num(zs[1])}, {"wzero": 0}))     # not the smallest zero
    if _differs(1.0, wz):
        cases.append(({"wzero": "1"}, {"wzero": 0}))
    return cases
