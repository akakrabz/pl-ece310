"""Independent check of questions/dft/dft-values.

X[k] from numpy.fft.fft (not the exact Q(sqrt2)/Q(sqrt3) code path of server.py), Parseval from the
squared magnitudes of all N bins, the inverse direction from numpy.fft.ifft of the given V[k]."""

import numpy as np


def _V(params):
    return np.array([complex(a, b) for a, b in params["V"]])


def differs(wrong, right):
    """pl-number-input (relabs, rtol 1e-3, atol 1e-6) marks `wrong` incorrect, with margin."""
    return abs(float(wrong) - float(right)) > 2 * (1e-6 + 1e-3 * abs(float(right)))


def check(params, correct):
    probs = []
    x = np.array(params["x"], float)
    N, k0 = params["N"], params["k0"]
    X = np.fft.fft(x)
    if len(x) != N or N not in (4, 6, 8):
        probs.append(f"length {len(x)} vs N = {N}")
    for name, want in (("X0", X[0].real), ("XH", X[N // 2].real), ("re", X[k0].real), ("im", X[k0].imag),
                       ("P", np.sum(np.abs(X) ** 2))):
        if abs(correct[name] - want) > 1e-9 * max(1, abs(want)):
            probs.append(f"{name} = {correct[name]} but numpy gives {want}")
    if abs(X[k0].imag) < 1e-6:
        probs.append("Im X[k0] = 0: the sign-of-exponent probe is meaningless")
    if k0 in (0, N // 2) or not 0 < k0 < N:
        probs.append(f"k0 = {k0} is not an 'other' bin")
    if not differs(np.sum(np.abs(x)), X[0].real) or not differs(X[N // 2].real, X[0].real):
        probs.append("X[0] mistakes coincide with the answer")
    V = _V(params)
    v = np.fft.ifft(V)
    if params["inv_kind"] == "vec":
        got = np.array(correct["v"]["_value"][0], float)
        if len(V) != 4 or got.shape != (4,) or not np.allclose(got, v.real) or np.max(np.abs(v.imag)) > 1e-9:
            probs.append(f"v = {got} but ifft gives {v}")
        if not np.allclose(v.real, np.round(v.real)):
            probs.append(f"v = {v} is not integer")
        if np.allclose(np.roll(v.real[::-1], 1), v.real):
            probs.append("v is circularly even: the forward-transform probe is meaningless")
    else:
        if abs(correct["v0"] - v[0].real) > 1e-9 or abs(v[0].imag) > 1e-9:
            probs.append(f"v0 = {correct['v0']} but ifft gives {v[0]}")
        if not np.allclose(V, np.conj(V[(-np.arange(len(V))) % len(V)])):
            probs.append("the given V is not Hermitian (v would not be real)")
        if not differs(np.sum(V).real, v[0].real) or not differs(V[0].real / len(V), v[0].real):
            probs.append("v0 mistakes coincide with the answer")
    return probs


def _num(v):
    v = float(v)
    return str(int(round(v))) if abs(v - round(v)) < 1e-12 else f"{v:.6g}"


def _sig4(v):
    return f"{float(v):.4g}"


def _row(vals, sep=", "):
    return "[" + sep.join(_num(t) for t in vals) + "]"


def submissions(params, correct):
    x = np.array(params["x"], float)
    N, k0 = params["N"], params["k0"]
    X = np.fft.fft(x)
    X0, XH = round(X[0].real), round(X[N // 2].real)
    re, im = X[k0].real, X[k0].imag
    E = float(np.sum(x ** 2))
    cases = [
        ({"X0": f"{2 * X0}/2", "XH": f"{XH}.0"}, {"X0": 1, "XH": 1}),
        ({"X0": f"{X0}.00", "XH": f"{3 * XH}/3"}, {"X0": 1, "XH": 1}),
        ({"re": _sig4(re), "im": _sig4(im)}, {"re": 1, "im": 1}),                 # 4 significant digits
        ({"re": f"{re:.7g}", "im": f"{im:.7g}"}, {"re": 1, "im": 1}),
        ({"im": _num(-im)}, {"im": 0}),                                           # e^{+j...}: conjugate
        ({"P": f"{2 * N * E:g}/2"}, {"P": 1}),
        ({"P": f"{N * E:.1f}"}, {"P": 1}),
        ({"P": _num(E)}, {"P": 0}),                                               # forgot the factor N
        ({"P": _num(E / N)}, {"P": 0}),                                           # 1/N on the wrong side
        ({"X0": "1/0"}, {"X0": "invalid"}),
    ]
    if differs(np.sum(np.abs(x)), X0):
        cases.append(({"X0": _num(np.sum(np.abs(x)))}, {"X0": 0}))               # added magnitudes
    cases.append(({"X0": _num(XH)}, {"X0": 0}))                                  # X[N/2] for X[0]
    cases.append(({"XH": _num(X0)}, {"XH": 0}))                                  # forgot (-1)^n
    if XH != 0:
        cases.append(({"XH": _num(-XH)}, {"XH": 0}))                             # alternation started with -1
    if differs(im, re) and differs(re, im):
        cases.append(({"re": _sig4(im), "im": _sig4(re)}, {"re": 0, "im": 0}))  # swapped
    half = np.sum(x * np.exp(-1j * np.pi * k0 * np.arange(N) / N))               # e^{-j pi k n/N}: forgot the 2
    if differs(half.real, re) or differs(half.imag, im):
        cases.append(({"re": _sig4(half.real), "im": _sig4(half.imag)},
                      {"re": 0 if differs(half.real, re) else 1, "im": 0 if differs(half.imag, im) else 1}))
    V = _V(params)
    if params["inv_kind"] == "vec":
        v = np.round(np.fft.ifft(V).real).astype(int)
        fwd = np.round((np.fft.fft(V) / 4).real).astype(int)                     # e^{-j}: circularly reversed v
        cases += [
            ({"v": _row(v, " ")}, {"v": 1}),                                      # MATLAB style
            ({"v": "[[" + ", ".join(f"{t}.0" for t in v) + "]]"}, {"v": 1}),      # Python style, decimals
            ({"v": _row(fwd)}, {"v": 0}),
            ({"v": _row(4 * v)}, {"v": 0}),                                       # forgot the 1/N
            ({"v": _row(v[:3])}, {"v": 0}),                                       # wrong length
            ({"v": "[1, 2,, 3]"}, {"v": "invalid"}),
        ]
    else:
        v0 = np.fft.ifft(V)[0].real
        Nv = len(V)
        tot = round(float(np.sum(V).real))
        cases += [
            ({"v0": f"{tot}/{Nv}"}, {"v0": 1}),
            ({"v0": f"{2 * tot}/{2 * Nv}"}, {"v0": 1}),
            ({"v0": _sig4(v0)}, {"v0": 1}),
            ({"v0": _num(tot)}, {"v0": 0}),                                       # forgot the 1/N
            ({"v0": _num(V[0].real / Nv)}, {"v0": 0}),                            # V[0]/N
        ]
    return cases
