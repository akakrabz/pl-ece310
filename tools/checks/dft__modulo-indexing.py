"""Independent check of questions/dft/modulo-indexing.

x is parsed from the displayed sequence. (a)-(c): the index expressions shown in the prompt are parsed
and evaluated with numpy (Python's % is n mod N), values compared with the answers. (d): the displayed
list of +-X_i is parsed, Y is built from X = np.fft.fft(x) and y = np.fft.ifft(Y); a real y must match
the matrix answer, a complex y must be reproduced by exactly one multiple-choice formula (each option
is parsed and evaluated numerically) and that one must be marked correct.
Probes: (a) plain flip, identity; (b) shift in the wrong direction, linear (zero-filled) advance;
(c) no reversal, reversal shifted the wrong way, plain flip; (d) the (-1)^k dropped, (-1)^k read as
(-1)^n, the negative, the reversal of x."""

import re

import checklib  # noqa: F401
import numpy as np
import prairielearn as pl


def _x(params):
    body = params["x_tex"].replace(r"\underset{\uparrow}", "")
    body = body.split(r"\{", 1)[1].rsplit(r"\}", 1)[0]
    return np.array([int(v.strip().strip("{}")) for v in body.split(r",\ ")], float)


def _index(tex, N):
    """'\\langle 2 - n\\rangle_{5}' -> integer index array."""
    m = re.fullmatch(r"\\langle (.*)\\rangle_\{(\d+)\}", tex.strip())
    if not m:
        assert tex.strip() == "n", tex
        return np.arange(N)
    assert int(m.group(2)) == N
    expr = m.group(1)
    assert re.fullmatch(r"[\dn +\-]+", expr), expr
    n = np.arange(N)
    return eval(expr, {"__builtins__": {}}, {"n": n}) % N


def _row(correct, name):
    return np.array(pl.from_json(correct[name]), dtype=complex).reshape(-1)


def _Y(params, x):
    items = params["Y_tex"].split(r"\{", 1)[1].rsplit(r"\}", 1)[0].split(r",\ ")
    X = np.fft.fft(x)
    out = []
    for it in items:
        m = re.fullmatch(r"\s*(-?)X_\{(\d+)\}\s*", it)
        out.append((-1 if m.group(1) else 1) * X[int(m.group(2))])
    return np.array(out)


def _opt_y(html, x):
    """'$y[n] = -e^{j\\frac{2\\pi}{3}n}\\,x[\\langle n - 3\\rangle_{6}]$' -> numpy array."""
    N = len(x)
    n = np.arange(N)
    s = html.strip().strip("$")
    m = re.fullmatch(r"y\[n\] = (-?)(?:e\^\{(-?)j(.*?)\}\\,)?x\[(.*)\]", s)
    assert m, s
    c = -1 if m.group(1) else 1
    ph = np.ones(N, complex)
    if m.group(3):
        body = m.group(3)
        q = re.fullmatch(r"\\frac\{(\d*)\\pi\}\{(\d+)\}n|(\d*)\\pi n", body)
        assert q, body
        if q.group(2):
            t = int(q.group(1) or 1) / int(q.group(2))
        else:
            t = int(q.group(3) or 1)
        ph = np.exp((-1 if m.group(2) else 1) * 1j * np.pi * t * n)
    return c * ph * x[_index(m.group(4), N)]


def _differs(a, b):
    return a.shape != b.shape or np.max(np.abs(a - b)) > 1e-3


def _probes(params, x):
    N = len(x)
    n = np.arange(N)
    mb, mc = params["mb"], params["mc"]
    lin = np.array([x[i + mb] if i + mb < N else 0 for i in range(N)])
    return {
        "ya": [x[::-1], x.copy()],
        "yb": [x[(n - mb) % N], lin],
        "yc": [x[(n - mc) % N], x[(-n - mc) % N], x[::-1]],
    }


def check(params, correct):
    probs = []
    x = _x(params)
    N = len(x)
    n = np.arange(N)
    if N not in (5, 6) or list(x) != params["x"]:
        probs.append(f"x = {x}")
    truth = {"ya": x[(-n) % N]}
    for name in ("ya", "yb", "yc"):
        idx = _index(params[f"{name}_idx_tex"], N)
        got = _row(correct, name)
        want = x[idx]
        if name == "ya" and not np.allclose(want, truth["ya"]):
            probs.append("(a) is not x[<-n>]")
        if got.shape != want.shape or not np.allclose(got, want):
            probs.append(f"{name}: answer {got.real} != {want}")
        wrong = [v for v in _probes(params, x)[name] if _differs(np.asarray(v, complex), want)]
        if len(wrong) < 2:
            probs.append(f"{name}: fewer than two distinct probes")
    answers = [tuple(_row(correct, k).real) for k in ("ya", "yb", "yc")]
    y = np.fft.ifft(_Y(params, x))
    is_real = np.max(np.abs(y.imag)) < 1e-9
    if is_real != params["d_real"]:
        probs.append(f"(d) d_real = {params['d_real']} but y = {np.round(y, 4)}")
    if params["d_real"]:
        got = _row(correct, "yd")
        if not np.allclose(got, y, atol=1e-9):
            probs.append(f"(d) answer {got.real} != ifft {np.round(y.real, 6)}")
        if not np.allclose(y.real, np.round(y.real)):
            probs.append("(d) y not integer")
        if tuple(np.round(y.real, 6)) in [tuple(np.round(a, 6)) for a in answers]:
            probs.append("(d) repeats an answer of (a)-(c)")
    else:
        opts = params["op"]
        hits = [o["key"] for o in opts if np.allclose(_opt_y(o["html"], x), y, atol=1e-9)]
        if len(hits) != 1:
            probs.append(f"(d) {len(hits)} options reproduce y")
        elif hits[0] != correct["op"]["key"]:
            probs.append("(d) the matching option is not marked correct")
        if len(opts) < 4:
            probs.append(f"(d) only {len(opts)} options")
    if len(set(answers)) < 3:
        probs.append("(a)-(c) answers repeat")
    return probs


def _fmt(v):
    return "[" + ", ".join(str(int(round(t.real))) for t in v) + "]"


def submissions(params, correct):
    x = _x(params)
    N = len(x)
    n = np.arange(N)
    cases = []
    for name in ("ya", "yb", "yc"):
        v = _row(correct, name)
        cases.append(({name: _fmt(v)}, {name: 1}))
        cases.append(({name: _fmt(v).replace(",", "")}, {name: 1}))                 # MATLAB style
        cases.append(({name: "[" + _fmt(v) + "]"}, {name: 1}))                       # Python style
        for w in _probes(params, x)[name]:
            if _differs(np.asarray(w, complex), v):
                cases.append(({name: _fmt(w)}, {name: 0}))
    cases.append(({"ya": "[1, 2,, 3]"}, {"ya": "invalid"}))
    if params["d_real"]:
        y = _row(correct, "yd")
        cases.append(({"yd": _fmt(y)}, {"yd": 1}))
        cases.append(({"yd": "[" + ", ".join(f"{t.real:.1f}" for t in y) + "]"}, {"yd": 1}))
        Yl = params["Y_list"]
        X = np.fft.fft(x)
        nosign = np.fft.ifft(np.array([X[i] for i in Yl["idx"]])).real          # (-1)^k dropped
        asn = (-1.0) ** n * nosign                                               # (-1)^k read as (-1)^n
        for w in (nosign, asn, -y.real, x[(-n) % N], x):
            if _differs(np.asarray(w, complex), y):
                cases.append(({"yd": _fmt(w)}, {"yd": 0}))
    else:
        for o in params["op"]:
            cases.append(({"op": o["key"]}, {"op": 1 if o["key"] == correct["op"]["key"] else 0}))
    return cases
