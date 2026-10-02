"""Independent check of questions/dtft/dtft-of-finite-sequence.

The DTFT is evaluated with numpy straight from the definition; the generated SymPy strings (the
exp form used by the Test button and the closed cos / sin / Euler form shown in the solution) are
parsed with SymPy and compared at random frequencies; |X_d(w0)| is recomputed with numpy."""

import numpy as np
from checklib import sympy

W = sympy.Symbol("w")
W0 = {"0": 0.0, "pi/2": np.pi / 2, "pi": np.pi}


def _dtft(x, start, w):
    n = np.arange(start, start + len(x))
    return complex(np.sum(np.asarray(x, float) * np.exp(-1j * w * n)))


def _eval_sym(s, w):
    expr = sympy.sympify(s.replace("^", "**"), locals={"j": sympy.I, "w": W, "e": sympy.E})
    return complex(expr.subs(W, w).evalf(30))


def _num(v):
    return str(int(v)) if float(v).is_integer() else repr(float(v))


def check(params, correct):
    probs = []
    x, start, w0 = params["x"], params["start"], params["w0"]
    rng = np.random.default_rng(len(x) * 31 + start)
    for key in ("X_sym", "closed_sym"):
        for w in rng.uniform(-np.pi, np.pi, 5):
            got, want = _eval_sym(params[key], w), _dtft(x, start, w)
            if abs(got - want) > 1e-9 * max(1, abs(want)):
                probs.append(f"{key} = {params[key]!r} gives {got} at w={w:.3f}, DTFT is {want}")
                break
    pts = rng.uniform(-np.pi, np.pi, 4)          # the Test button's 'incorrect' answer must be wrong
    if max(abs(_eval_sym(params["X_sym_conj"], w) - _dtft(x, start, w)) for w in pts) < 1e-6:
        probs.append("X_sym_conj equals the DTFT (sequence even about n = 0?)")
    mag = abs(_dtft(x, start, W0[w0]))
    if abs(mag - correct["mag"]) > 1e-9 * max(1, mag):
        probs.append(f"|X_d({w0})| = {mag} but the correct answer is {correct['mag']}")
    if abs(mag - round(mag)) > 1e-9 or round(mag) <= 0:
        probs.append(f"|X_d({w0})| = {mag} is not a positive integer")
    if "X" in correct:
        probs.append("correct_answers['X'] is set; X must be graded by server.grade only")
    if len(_mag_mistakes(params, correct)) < 2:
        probs.append("fewer than two distinct classic mistakes for |X_d(w0)|")
    xs = list(x)
    symmetric, antisym = xs == xs[::-1], xs == [-v for v in xs[::-1]]
    if params["family"] == "gen" and (symmetric or antisym):
        probs.append(f"'gen' sequence {x} is (anti)symmetric: the solution would miss the closed form")
    if params["family"] == "sym" and not symmetric or params["family"] == "anti" and not antisym:
        probs.append(f"family {params['family']} does not match the sequence {x}")
    return probs


def _mag_mistakes(params, correct):
    x, start, w0 = params["x"], params["start"], params["w0"]
    X = _dtft(x, start, W0[w0])
    re, im = round(X.real), round(X.imag)
    mag = round(correct["mag"])
    cands = []
    if abs(X.imag) < 1e-9 and X.real < 0:
        cands.append(re)                                              # forgot the absolute value
    if w0 == "pi":
        cands.append(sum(x))                                          # forgot (-1)^n
    if w0 == "pi/2":
        cands.append(re * re + im * im)                               # forgot the square root
        if re and im:
            cands.append(abs(re) + abs(im))
    if w0 == "0":
        cands.append(sum(abs(v) for v in x))
    half = (len(x) - 1) // 2
    if params["family"] == "sym":                                     # forgot the 2 in 2a cos(kw)
        cands.append(round(abs(x[half] + sum(x[half - k] * np.cos(k * W0[w0]) for k in range(1, half + 1)))))
    if params["family"] == "anti":                                    # forgot the 2 in 2ja sin(kw)
        cands.append(round(abs(sum(x[half - k] * np.sin(k * W0[w0]) for k in range(1, half + 1)))))
    out = []
    for c in cands:
        if c != mag and c not in out:
            out.append(c)
    return out


def submissions(params, correct):
    x, start = params["x"], params["start"]
    ns = range(start, start + len(x))
    mag = round(correct["mag"])
    caret = " + ".join(f"({v})*e^(-j*({n})*w)" for n, v in zip(ns, x) if v)      # e^(...) syntax
    first = start
    factored = (f"e^(-j*({first})*w)*(" + " + ".join(f"({v})*e^(-j*{i}*w)" for i, v in enumerate(x) if v) + ")")
    sign_err = " + ".join(f"({v})*exp(j*({n})*w)" for n, v in zip(ns, x) if v)     # e^{+j w n}
    cases = [
        ({"X": caret}, {"X": 1}),
        ({"X": params["closed_sym"]}, {"X": 1}),          # cosine / sine / Euler form
        ({"X": factored}, {"X": 1}),
        ({"X": sign_err}, {"X": 0}),
        ({"mag": f"{2 * mag}/2"}, {"mag": 1}),
        ({"mag": f"{mag}.0"}, {"mag": 1}),
        ({"X": "exp(-j*w"}, {"X": "invalid"}),
        ({"X": "x + 1"}, {"X": "invalid"}),
    ]
    if start != 0:                                         # indexed from the first listed sample
        cases.append(({"X": " + ".join(f"({v})*exp(-j*{i}*w)" for i, v in enumerate(x) if v)}, {"X": 0}))
    if params["family"] == "sym":                          # dropped the linear-phase factor
        half = (len(x) - 1) // 2
        amp = f"({x[half]})" + "".join(f" + 2*({x[half - k]})*cos({k}*w)" for k in range(1, half + 1))
        cases.append(({"X": amp}, {"X": 0}))
    for m in _mag_mistakes(params, correct):
        cases.append(({"mag": _num(m)}, {"mag": 0}))
    return cases
