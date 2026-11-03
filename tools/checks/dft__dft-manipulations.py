"""Independent check of questions/dft/dft-manipulations.

Types dup, inter, dec, hw9: the displayed list {Y[k]} is parsed (entries 0, X_{i}, X_{i} + X_{j}) and built
from numpy.fft.fft(x); y = numpy.fft.ifft(Y) must be the answer. The rule itself is also rebuilt from the
type name for random real and complex x of several lengths and compared with the general statement
marked correct in (b); every other option of (b) must fail on random signals.
Type fa21: y[n] = x[<m - n>_8] c^n is built sample by sample, Y = numpy.fft.fft(y) at k = 0, 2, 4, 6 must be
the answer, and the formula options of (b) are evaluated for random complex x of length 8."""

import re

import checklib  # noqa: F401  (puts PrairieLearn's python library on sys.path)
import numpy as np
import prairielearn as pl

RNG = np.random.default_rng(2025)
C = {2: 1j, 6: -1j, 4: -1.0}


def _norm(s):
    return re.sub(r"\s+", "", s)


def _rand(N, cplx=True):
    v = RNG.standard_normal(N)
    return v + 1j * RNG.standard_normal(N) if cplx else v


# ----------------------------------------------------------------------------- rules and statements
def _rule(kind, x):
    """Y from x by the rule of the type (numpy), returned with its inverse DFT."""
    x = np.asarray(x, dtype=complex)
    N = len(x)
    X = np.fft.fft(x)
    if kind == "dup":
        Y = np.concatenate([X, X])
    elif kind == "inter":
        Y = np.zeros(2 * N, complex)
        Y[0::2] = X
    elif kind == "dec":
        Y = X[0::2]
    else:
        Y = X + X[(np.arange(N) - N // 2) % N]
    return Y, np.fft.ifft(Y)


def _gen(gid, x):
    """The general statement `gid` as a signal (independent numpy version of each option)."""
    x = np.asarray(x, dtype=complex)
    N = len(x)
    n = np.arange(N)
    h = N // 2
    z = np.zeros(N, complex)
    inter = np.zeros(2 * N, complex)
    inter[0::2] = x
    odd = np.zeros(2 * N, complex)
    odd[1::2] = x
    return {
        "dup": inter, "dup2": 2 * inter, "dupodd": odd,
        "half": np.tile(x, 2) / 2, "per": np.tile(x, 2), "zpad": np.concatenate([x, z]),
        "alias": x[:h] + x[h:2 * h], "alias_half": (x[:h] + x[h:2 * h]) / 2, "alias_diff": x[:h] - x[h:2 * h],
        "decim": x[0::2][:h], "trunc": x[:h],
        "mod1": (1 + (-1.0) ** n) * x, "mod1m": (1 - (-1.0) ** n) * x, "tsum": x + np.roll(x, h),
        "two": 2 * x, "modshift": (1 + (-1.0) ** n) * np.roll(x, h),
    }[gid]


def _fa21_formula(fid, X, m, ell):
    k = np.arange(8)
    ph = lambda arg, s=-1: np.exp(s * 2j * np.pi * m * arg / 8)       # noqa: E731
    if fid == "ok":
        return ph(k - ell) * X[(ell - k) % 8]
    if fid == "dir":
        return ph(k + ell) * X[(-ell - k) % 8]
    if fid == "norev":
        return ph(k - ell) * X[(k - ell) % 8]
    if fid == "psign":
        return ph(k - ell, 1) * X[(ell - k) % 8]
    if fid == "noshift":
        return ph(k) * X[(ell - k) % 8]
    return ph(k) * X[(-k) % 8]                                           # nomod


def _fa21_y(x, m, ell):
    n = np.arange(8)
    return np.asarray(x, complex)[(m - n) % 8] * C[ell] ** n


def _option_ids(params):
    """element key -> statement id, matched through the option text."""
    by_text = {_norm(o["text"]): gid for o, gid in zip(params["gen_choices"], params["gen_ids"])}
    return {o["key"]: by_text[_norm(o["html"])] for o in params["gen"]}


def _parse_list(tex, X):
    body = tex.split(r"\{", 1)[1].rsplit(r"\}", 1)[0]
    out = []
    for item in body.split(r",\ "):
        item = item.strip()
        if item == "0":
            out.append(0)
            continue
        idx = [int(i) for i in re.findall(r"X_\{(\d+)\}", item)]
        if not idx or _norm(item) != "+".join(f"X_{{{i}}}" for i in idx):
            raise ValueError(item)
        out.append(sum(X[i] for i in idx))
    return np.array(out, complex)


# ----------------------------------------------------------------------------- check
def check(params, correct):
    probs = []
    kind = params["kind"]
    ids = _option_ids(params)
    marked = ids[correct["gen"]["key"]]
    if kind == "fa21":
        x = np.array(params["x"], float)
        m, ell = params["m"], params["ell"]
        if not np.allclose(x[:4], x[4:]):
            probs.append("fa21: x is not 4-periodic")
        Y = np.fft.fft(_fa21_y(x, m, ell))
        got = np.array(pl.from_json(correct["Y"]), complex).reshape(-1)
        if not np.allclose(got, Y[[0, 2, 4, 6]], atol=1e-9):
            probs.append(f"fa21: answer {got} but fft gives {np.round(Y[[0, 2, 4, 6]], 6)}")
        if not np.allclose(Y[1::2], 0, atol=1e-9):
            probs.append("fa21: Y has nonzero odd bins")
        for _ in range(3):
            xr = _rand(8)
            Yr = np.fft.fft(_fa21_y(xr, m, ell))
            for key, fid in ids.items():
                ok = np.allclose(_fa21_formula(fid, np.fft.fft(xr), m, ell), Yr, atol=1e-9)
                if ok != (fid == marked):
                    probs.append(f"fa21 (b): option {fid} matches = {ok}, marked {marked}")
                    break
        if marked != "ok":
            probs.append(f"fa21 (b): marked {marked}")
        if len(_mistakes(params)) < 2:
            probs.append("fa21: fewer than two distinct wrong rows")
        return probs

    x = np.array(params["x"], float)
    N = len(x)
    Y = _parse_list(params["Y_list_tex"], np.fft.fft(x))
    Yr, yr = _rule(kind, x)
    if len(Y) != len(Yr) or not np.allclose(Y, Yr, atol=1e-9):
        probs.append(f"{kind}: displayed list is not the rule")
    y = np.fft.ifft(Y)
    got = np.array(pl.from_json(correct["y"]), float).reshape(-1)
    if got.shape != y.shape or not np.allclose(got, y.real, atol=1e-9) or not np.allclose(y.imag, 0, atol=1e-9):
        probs.append(f"{kind}: answer {got} but ifft gives {np.round(y, 6)}")
    if not np.allclose(got, np.round(got)):
        probs.append(f"{kind}: non-integer answer {got}")
    if params["M"] != len(Y):
        probs.append("M is not the length of the list")
    # the rule in general: random real and complex x of two lengths (even for dec, hw9)
    lengths = (6, 8) if kind in ("dec", "hw9") else (5, 8)
    for Nr in lengths:
        for cplx in (False, True):
            xr = _rand(Nr, cplx)
            _, yy = _rule(kind, xr)
            for key, gid in ids.items():
                ok = np.allclose(_gen(gid, xr), yy, atol=1e-9)
                if ok != (gid == marked):
                    probs.append(f"(b) N = {Nr}: option {gid} matches = {ok}, marked {marked}")
    if not np.allclose(_gen(marked, x), y, atol=1e-9):
        probs.append("(b) marked statement does not reproduce (a)")
    if len(_mistakes(params)) < 2:
        probs.append(f"{kind}: fewer than two distinct wrong rows")
    return probs


# ----------------------------------------------------------------------------- probes
def _mistakes(params):
    kind = params["kind"]
    x = np.array(params["x"], float)
    if kind == "fa21":
        m, ell = params["m"], params["ell"]
        X = np.fft.fft(x)
        right = np.fft.fft(_fa21_y(x, m, ell))[[0, 2, 4, 6]]
        cands = [_fa21_formula(f, X, m, ell)[[0, 2, 4, 6]] for f in ("dir", "norev", "psign", "noshift", "nomod")]
        cands += [X[[0, 2, 4, 6]], np.conj(right)]
    else:
        right = np.fft.ifft(_rule(kind, x)[0])
        others = {"dup": ["per", "half", "zpad", "dupodd", "dup2"], "inter": ["per", "dup", "zpad", "dupodd"],
                  "dec": ["decim", "alias_half", "alias_diff", "trunc"], "hw9": ["tsum", "two", "mod1m", "modshift"]}[kind]
        cands = [_gen(g, x) for g in others]
    out = []
    for v in cands:
        v = np.asarray(v, complex)
        if v.shape == right.shape and np.max(np.abs(v - right)) < 1e-3:
            continue
        if any(v.shape == u.shape and np.max(np.abs(v - u)) < 1e-9 for u in out):
            continue
        out.append(v)
    return out


def _fmt(z):
    z = complex(z)
    re_, im = round(z.real, 6), round(z.imag, 6)
    re_ = 0.0 if re_ == 0 else re_
    im = 0.0 if im == 0 else im

    def num(v):
        return str(int(v)) if float(v).is_integer() else f"{v:g}"

    if im == 0:
        return num(re_)
    jt = ("" if abs(im) == 1 else num(abs(im))) + "j"
    if re_ == 0:
        return ("-" if im < 0 else "") + jt
    return f"{num(re_)}{'-' if im < 0 else '+'}{jt}"


def _row(v, sep=", "):
    return "[" + sep.join(_fmt(z) for z in v) + "]"


def submissions(params, correct):
    name = "Y" if params["kind"] == "fa21" else "y"
    right = np.array(pl.from_json(correct[name]), complex).reshape(-1)
    cases = [
        ({name: _row(right, " ")}, {name: 1}),                                   # MATLAB style
        ({name: "[[" + ", ".join(_fmt(z) for z in right) + "]]"}, {name: 1}),    # Python style
        ({name: "[1, 2,, 3]"}, {name: "invalid"}),
    ]
    if name == "y":
        cases.append(({name: "[" + ", ".join(f"{z.real:.1f}" for z in right) + "]"}, {name: 1}))
        cases.append(({name: _row(np.array(params["x"], float))}, {name: 0}))   # x itself (wrong length or wrong values)
    else:
        cases.append(({name: _row(right).replace("j", "i")}, {name: 1}))
    for v in _mistakes(params):
        cases.append(({name: _row(v)}, {name: 0}))
    for o in params["gen"]:
        cases.append(({"gen": o["key"]}, {"gen": 1 if o["key"] == correct["gen"]["key"] else 0}))
    return cases
