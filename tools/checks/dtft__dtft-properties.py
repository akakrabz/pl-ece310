"""Independent check of questions/dtft/dtft-properties.

y[n] is built sample by sample from the base sequence and the transformation (geometric tails are
truncated at n = 400, |a| <= 3/4, so the tail is below 1e-45), Y_d is computed from the samples by
the DTFT sum, and every multiple-choice option is evaluated from its DISPLAYED LaTeX (a small
converter for the LaTeX subset the generator emits -> SymPy). The correct option must match Y_d at
random frequencies, every distractor must differ from it and from the others. The numeric part is
recomputed from the samples (the Parseval integral by a periodic rectangle rule on an FFT grid)."""

import html
import re
from fractions import Fraction

import numpy as np
from checklib import sympy
from sympy.parsing.sympy_parser import implicit_multiplication_application, parse_expr, standard_transformations

W = sympy.Symbol("w")
TRANS = standard_transformations + (implicit_multiplication_application,)
NMAX = 400
RNG_W = [0.43, 1.21, 2.05, 2.77, -0.66, -1.9, -2.9, 0.08]


# ----------------------------------------------------------------------------- LaTeX -> SymPy
def _group(s, i):
    depth = 0
    for k in range(i, len(s)):
        if s[k] == "{":
            depth += 1
        elif s[k] == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:k], k + 1
    raise ValueError(f"unbalanced braces in {s!r}")


def _conv(s):
    out, i = [], 0
    while i < len(s):
        if s.startswith((r"\dfrac", r"\tfrac", r"\frac"), i):
            a, j = _group(s, s.index("{", i))
            b, j = _group(s, j)
            out.append(f"(({_conv(a)})/({_conv(b)}))")
            i = j
        elif s.startswith("e^{", i):
            a, j = _group(s, i + 2)
            out.append(f" exp({_conv(a)}) ")
            i = j
        elif s.startswith("^{", i):
            a, j = _group(s, i + 1)
            out.append(f"**({_conv(a)})")
            i = j
        elif s.startswith(r"\left(", i):
            out.append("(")
            i += 6
        elif s.startswith(r"\right)", i):
            out.append(")")
            i += 7
        elif s.startswith(r"\omega", i):
            out.append(" w ")
            i += 6
        elif s.startswith(r"\pi", i):
            out.append(" pi ")
            i += 3
        elif s.startswith(r"\cdot", i):
            out.append("*")
            i += 5
        elif s.startswith((r"\,", r"\;", r"\!"), i):
            out.append(" ")
            i += 2
        elif s[i] == "j":
            out.append(" I*")
            i += 1
        elif s[i] == "\\":
            raise ValueError(f"unknown LaTeX macro at {s[i:i + 15]!r}")
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


def latex_value(tex):
    """Function w -> complex for a LaTeX expression of the generator's subset."""
    tex = html.unescape(tex).strip()
    if tex.startswith("$") and tex.endswith("$"):
        tex = tex[1:-1]
    expr = parse_expr(_conv(tex), local_dict={"w": W, "I": sympy.I, "pi": sympy.pi, "exp": sympy.exp},
                      transformations=TRANS)
    f = sympy.lambdify(W, expr, "numpy")
    return lambda w: complex(f(w))


# ----------------------------------------------------------------------------- samples
def _x(params):
    """(values, start) of x[n]."""
    if params["base_kind"] == "geo":
        a, b = float(Fraction(params["a"])), float(Fraction(params["b"]))
        return b * a ** np.arange(NMAX + 1), 0
    return np.array(params["c"], dtype=float), 0


def _y(params):
    x, s = _x(params)
    n = np.arange(s, s + len(x))
    kind = params["kind"]
    if kind == "shift":
        return x.astype(complex), s + params["k"]
    if kind == "mod":
        w0 = float(Fraction(params["w0"])) * np.pi
        return x * np.exp(1j * w0 * n), s
    if kind == "alt":
        return x * (-1.0) ** n, s
    if kind == "cos":
        return x * np.cos(np.pi * n / 2), s
    if kind == "rev":
        return x[::-1].astype(complex), -(s + len(x) - 1)
    if kind == "conv":
        return np.convolve(x, x).astype(complex), 2 * s
    raise ValueError(kind)


def _dtft(vals, start, w):
    n = np.arange(start, start + len(vals))
    return complex(np.sum(vals * np.exp(-1j * w * n)))


def _X(params, w):
    x, s = _x(params)
    return _dtft(x, s, w)


def _numeric(params):
    at = params["at"]
    yv, ys = _y(params)
    if at == "energy":
        N = 8192
        n = np.arange(ys, ys + len(yv))
        grid = 2 * np.pi * np.arange(N) / N
        Y = np.array([np.sum(yv * np.exp(-1j * g * n)) for g in grid[::8]])       # coarse check below
        Yf = np.fft.fft(np.concatenate([yv, np.zeros(N - len(yv))])) * np.exp(-1j * grid * ys)
        assert np.allclose(Yf[::8], Y)
        return (2 / N) * np.sum(np.abs(Yf) ** 2)                                  # (1/pi) * integral
    w = {"pi": np.pi, "0": 0.0, "pi/2": np.pi / 2, "w0": float(Fraction(params["w0"])) * np.pi}[at]
    return _dtft(yv, ys, w)


def _options(params):
    return [(o["key"], o["html"]) for o in params["yd"]]


def check(params, correct):
    probs = []
    yv, ys = _y(params)
    Ytrue = [_dtft(yv, ys, w) for w in RNG_W]
    fX = latex_value(params["X_tex"])
    xv, xs = _x(params)
    if max(abs(fX(w) - _dtft(xv, xs, w)) for w in RNG_W) > 1e-8:
        probs.append("the stated X_d(w) does not match the DTFT of x[n]")
    opts = _options(params)
    if len(opts) != 4:
        probs.append(f"{len(opts)} options instead of 4")
    vals = {}
    for key, h in opts:
        try:
            f = latex_value(h)
            vals[key] = [f(w) for w in RNG_W]
        except Exception as exc:
            probs.append(f"option {h!r} cannot be evaluated: {exc}")
            return probs
    ckey = correct["yd"]["key"]
    if max(abs(a - b) for a, b in zip(vals[ckey], Ytrue)) > 1e-8 * max(1, max(abs(b) for b in Ytrue)):
        probs.append(f"the option marked correct ({correct['yd']['html']}) is not the DTFT of y")
    for key, h in opts:
        if key != ckey and max(abs(a - b) for a, b in zip(vals[key], Ytrue)) < 1e-6:
            probs.append(f"distractor {h} equals the DTFT of y")
    keys = list(vals)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if max(abs(a - b) for a, b in zip(vals[keys[i]], vals[keys[j]])) < 1e-6:
                probs.append(f"options {keys[i]} and {keys[j]} are the same function")
    if f"${params['Y_tex']}$" != html.unescape(correct["yd"]["html"]).strip():
        probs.append("the answer panel's Y_d differs from the correct option")
    num = _numeric(params)
    if abs(num.imag) > 1e-9 or abs(num.real - correct["num"]) > 1e-8 * max(1, abs(num.real)):
        probs.append(f"numeric part {correct['num']} != recomputed {num}")
    f = Fraction(correct["num"]).limit_denominator(1000)
    if abs(float(f) - correct["num"]) > 1e-12 or f.denominator > 64:
        probs.append(f"numeric answer {correct['num']} is not a nice fraction")
    if len(_slips(params, correct)) < 1:
        probs.append("no classic numeric slip differs from the answer")
    return probs


def _slips(params, correct):
    at, kind = params["at"], params["kind"]
    X0, Xpi, X2 = _X(params, 0.0).real, _X(params, np.pi).real, _X(params, np.pi / 2)
    c = []
    if kind == "shift":
        c.append(Xpi)                                   # dropped e^{-jk pi}
    elif kind == "mod":
        c.append(Xpi)                                   # X_d(w + w0) at w0: X_d(2 w0) = X_d(pi)
    elif kind == "alt":
        c.append(X0 if at == "0" else Xpi)              # forgot the shift by pi
    elif kind == "cos":
        c += [X0 + Xpi, X0] if at == "pi/2" else [2 * X2.real]      # missing 1/2, one term only
    elif kind == "rev":
        xv, _ = _x(params)
        c += [float(np.sum(np.abs(xv) ** 2)), 2 * float(np.sum(xv)) ** 2]   # forgot 2 pi; squared the sum
    elif kind == "conv":
        Xa = X0 if at == "0" else Xpi
        c += [2 * Xa, Xa]                               # 2 X_d instead of X_d^2; forgot to square
    ref = correct["num"]
    out = []
    for v in c:
        if abs(v - ref) > 1e-2 * max(abs(ref), 1e-3) and all(abs(v - u) > 1e-9 for u in out):
            out.append(v)
    return out


def _plain(v):
    f = Fraction(v).limit_denominator(1000)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def submissions(params, correct):
    ckey = correct["yd"]["key"]
    cases = [({"yd": ckey}, {"yd": 1})]
    cases += [({"yd": k}, {"yd": 0}) for k, _ in _options(params) if k != ckey]   # every distractor
    v = correct["num"]
    cases.append(({"num": _plain(v)}, {"num": 1}))
    f = Fraction(v).limit_denominator(1000)
    cases.append(({"num": f"{2 * f.numerator}/{2 * f.denominator}"}, {"num": 1}))   # unreduced fraction
    if abs(v - round(v)) > 1e-9:
        cases.append(({"num": f"{v:.4g}"}, {"num": 1}))                          # 4 significant digits
    else:
        cases.append(({"num": f"{v:.1f}"}, {"num": 1}))
    for s in _slips(params, correct):
        cases.append(({"num": _plain(s)}, {"num": 0}))
    cases.append(({"num": "1/0"}, {"num": "invalid"}))
    return cases
