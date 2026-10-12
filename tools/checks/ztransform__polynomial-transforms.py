"""Independent check of questions/ztransform/polynomial-transforms.

The displayed LaTeX is parsed back (the sequence x_tex with its arrow, the Laurent polynomial H_tex),
so the check covers exactly what the student sees. X(z) is compared with SymPy, h[n] is read off
H(z), y = x * h comes from numpy.convolve, and Y(+-1) from a direct sum over the samples of y."""

import re

import numpy as np
import sympy
from checklib import S, Z

TERM = re.compile(r"([+-]?)(\d+)?(z)?(?:\^\{(-?\d+)\})?")


def parse_seq(tex):
    """'\\{1,\\ \\underset{\\uparrow}{2},\\ 3\\}' -> (values, start), trimmed to nonzero ends."""
    body = tex.strip()
    assert body.startswith(r"\{") and body.endswith(r"\}"), tex
    items = body[2:-2].split(r",\ ")
    vals, zero_at = [], None
    for i, it in enumerate(items):
        m = re.fullmatch(r"\\underset\{\\uparrow\}\{(-?\d+)\}", it)
        if m:
            zero_at = i
            vals.append(int(m.group(1)))
        else:
            vals.append(int(it))
    assert zero_at is not None, tex
    start = -zero_at
    while vals and vals[0] == 0:
        vals.pop(0)
        start += 1
    while vals and vals[-1] == 0:
        vals.pop()
    return vals, start


def parse_laurent(tex):
    """'2z^{2} - 1 + 3z^{-3}' -> {n: coefficient} with z^{-n} <-> n."""
    s = tex.replace(" ", "")
    out = {}
    for m in TERM.finditer(s):
        if not m.group(0) or m.group(0) in "+-":
            continue
        sign = -1 if m.group(1) == "-" else 1
        mag = int(m.group(2)) if m.group(2) else 1
        power = (int(m.group(4)) if m.group(4) else 1) if m.group(3) else 0
        out[-power] = out.get(-power, 0) + sign * mag
    return {n: c for n, c in out.items() if c != 0}


def as_vec(d):
    lo, hi = min(d), max(d)
    return [d.get(n, 0) for n in range(lo, hi + 1)], lo


def zpoly(vals, start):
    return sum(v * Z ** (-(start + i)) for i, v in enumerate(vals))


def sgn(n):
    return 1 if n % 2 == 0 else -1


def data(params):
    x, nx = parse_seq(params["x_tex"])
    h, nh = as_vec(parse_laurent(params["H_tex"]))
    y = np.convolve(x, h).astype(int).tolist()
    ny = nx + nh
    return x, nx, h, nh, y, ny


def at(vals, start, n):
    i = n - start
    return vals[i] if 0 <= i < len(vals) else 0


def yk_mistakes(params):
    x, nx, h, nh, y, ny = data(params)
    k = params["k"]
    right = at(y, ny, k)
    cands = [at(y, ny, -k)]                                                        # coefficient of z^{+k}
    cands.append(sum(at(x, nx, m) * at(h, nh, m - k) for m in range(nx - 10, nx + 20)))   # h not flipped
    cands.append(at(x, nx, k) * at(h, nh, k))                                      # pointwise product
    return [c for i, c in enumerate(cands) if c != right and c not in cands[:i]]


def Y_mistakes(params):
    x, nx, h, nh, y, ny = data(params)
    z0 = params["z0"]
    Xv = sum(v * (sgn(nx + i) if z0 == -1 else 1) for i, v in enumerate(x))
    Hv = sum(v * (sgn(nh + i) if z0 == -1 else 1) for i, v in enumerate(h))
    right = Xv * Hv
    cands = [Xv + Hv, Xv, Hv] if z0 == 1 else [sum(x) * sum(h), Xv + Hv, Xv]
    return [c for i, c in enumerate(cands) if c != right and c not in cands[:i]]


def check(params, correct):
    probs = []
    x, nx, h, nh, y, ny = data(params)
    if x != params["x_vals"] or nx != params["x_start"]:
        probs.append(f"x_tex parses to {x}@{nx}, params say {params['x_vals']}@{params['x_start']}")
    if sympy.expand(S(correct["X"]) - zpoly(x, nx)) != 0:
        probs.append(f"X = {correct['X']} != {zpoly(x, nx)}")
    d = {nx + i: v for i, v in enumerate(x)}
    if all(d.get(n, 0) == d.get(-n, 0) for n in d):
        probs.append("x is even about n = 0 (z^{+n} probe meaningless)")
    has_pos = any(v != 0 and nx + i > 0 for i, v in enumerate(x))
    has_neg = any(v != 0 and nx + i < 0 for i, v in enumerate(x))
    expect = {(True, True): r"i.e. $0 \lt \lvert z\rvert \lt \infty$", (True, False): r"i.e. $\lvert z\rvert \gt 0$",
              (False, True): r"i.e. $\lvert z\rvert \lt \infty$", (False, False): None}[(has_pos, has_neg)]
    html = correct["roc"]["html"].strip()
    if not (html.endswith(expect) if expect else html == "all $z$"):
        probs.append(f"ROC marked {html!r}, expected the option ending {expect!r}")
    if len(params["roc"]) != 4:
        probs.append("expected 4 ROC options")
    if correct["nh"] != nh:
        probs.append(f"first index of h: {correct['nh']} != {nh}")
    got_h = [int(round(v)) for v in correct["h"]["_value"][0]]
    if got_h != h:
        probs.append(f"h vector {got_h} != {h}")
    if nh == -(nh + len(h) - 1):
        probs.append("support of h symmetric about 0 (mirrored-index probe meaningless)")
    k = params["k"]
    if correct["yk"] != at(y, ny, k):
        probs.append(f"y[{k}] = {correct['yk']} but numpy gives {at(y, ny, k)}")
    nprod = sum(1 for m in range(nx, nx + len(x)) if at(x, nx, m) * at(h, nh, k - m) != 0)
    if nprod < 2 or k == 0 or at(y, ny, k) == 0:
        probs.append(f"y[{k}] is too easy ({nprod} products, value {at(y, ny, k)})")
    z0 = params["z0"]
    Yv = sum(v * (sgn(ny + i) if z0 == -1 else 1) for i, v in enumerate(y))
    if correct["Yz"] != Yv:
        probs.append(f"Y({z0}) = {correct['Yz']} but the samples of y give {Yv}")
    if len(yk_mistakes(params)) < 2 or len(Y_mistakes(params)) < 2:
        probs.append("fewer than two classic mistakes for y[k] or Y(z0)")
    return probs


def _vec(v):
    return "[" + ", ".join(map(str, v)) + "]"


def submissions(params, correct):
    x, nx, h, nh, y, ny = data(params)
    E = max(nx + len(x) - 1, 0)
    over = " + ".join(f"({v})*z^{E - (nx + i)}" for i, v in enumerate(x) if v)
    slash = " + ".join((f"({v})/z^{nx + i}" if nx + i > 0 else f"({v})*z^{-(nx + i)}") for i, v in enumerate(x) if v)
    plus = " + ".join(f"({v})*z^({nx + i})" for i, v in enumerate(x) if v)            # z^{+n}: sign mix-up
    cases = [
        ({"X": f"({over})/z^{E}"}, {"X": 1}),
        ({"X": slash}, {"X": 1}),
        ({"X": plus}, {"X": 0}),
        ({"X": "x + 1"}, {"X": "invalid"}),
    ]
    if nx != 0:                                                                       # indexed from the first sample
        cases.append(({"X": " + ".join(f"({v})*z^(-{i})" for i, v in enumerate(x) if v)}, {"X": 0}))
    for o in params["roc"]:
        cases.append(({"roc": o["key"]}, {"roc": 1 if o["key"] == correct["roc"]["key"] else 0}))
    eh = nh + len(h) - 1
    cases += [
        ({"h": "[" + " ".join(map(str, h)) + "]"}, {"h": 1}),                           # MATLAB style
        ({"h": _vec([f"{v}.0" for v in h])}, {"h": 1}),
        ({"nh": str(-eh)}, {"nh": 0}),                                                  # z^{2} read as n = 2
        ({"nh": "1.5"}, {"nh": "invalid"}),
    ]
    if h != h[::-1]:
        cases.append(({"h": _vec(h[::-1])}, {"h": 0}))
    if 0 in h:
        cases.append(({"h": _vec([v for v in h if v != 0])}, {"h": 0}))                # gaps dropped
    yk = correct["yk"]
    cases += [({"yk": f"+{yk}" if yk >= 0 else f"-{abs(yk)}"}, {"yk": 1}), ({"yk": f" {yk} "}, {"yk": 1})]
    cases += [({"yk": str(m)}, {"yk": 0}) for m in yk_mistakes(params)]
    Yv = correct["Yz"]
    cases += [({"Yz": f"+{Yv}" if Yv >= 0 else f"-{abs(Yv)}"}, {"Yz": 1}), ({"Yz": f" {Yv} "}, {"Yz": 1})]
    cases += [({"Yz": str(m)}, {"Yz": 0}) for m in Y_mistakes(params)]
    return cases
