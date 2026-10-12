"""Independent check of questions/freqresp/pole-zero-and-magnitude.

Each plotted system is classified from its poles and zeros alone: |H_d| on a grid of pi/2400 (which
contains every pole/zero angle used) from the product of distances, then by shape: constant ->
all-pass; monotone with a zero at pi / at 0 -> low-pass / high-pass; zeros at 0 and pi with one
interior peak -> band-pass; a single interior zero with large gain at 0 and pi -> notch; three or more
zeros on [0, pi] -> comb. The class must equal the matching's correct option. The numeric part is
recomputed from K and the labelled real poles/zeros (which must be the plotted ones)."""

import math
from fractions import Fraction

import numpy as np

TEXT = {"lowpass": "low-pass", "highpass": "high-pass", "bandpass": "band-pass", "bandstop": "band-stop (notch)"}


def _mag(sysd, w):
    e = np.exp(1j * np.asarray(w, float))
    m = np.ones_like(e, dtype=float)
    for z in sysd["zeros"]:
        m = m * np.abs(e - complex(*z))
    for p in sysd["poles"]:
        m = m / np.abs(e - complex(*p))
    return m


def classify(sysd):
    w = np.linspace(0, np.pi, 2401)
    m = _mag(sysd, w)
    m = m / m.max()
    if m.max() - m.min() < 1e-9:
        return "allpass"
    zero = m < 1e-9
    nz = int(zero.sum())
    d = np.diff(m)
    if nz >= 3:
        return "comb"
    if np.all(d <= 1e-12) and zero[-1]:
        return "lowpass"
    if np.all(d >= -1e-12) and zero[0]:
        return "highpass"
    if zero[0] and zero[-1] and nz == 2:
        i = int(np.argmax(m))
        if np.all(d[:i] >= -1e-12) and np.all(d[i:] <= 1e-12):
            return "bandpass"
    if nz == 1 and not zero[0] and not zero[-1] and m[0] > 0.5 and m[-1] > 0.5:
        return "bandstop"
    return "other"


def _num_value(params, at_pi):
    K = float(Fraction(params["K"]))
    zinv = -1.0 if at_pi else 1.0
    num = np.prod([1 - float(Fraction(z)) * zinv for z in params["num_zeros"]])
    den = np.prod([1 - float(Fraction(q)) * zinv for q in params["num_poles"]])
    return abs(K * num / den)


def check(params, correct):
    probs = []
    stmts, opts = params["ftype"]
    if len(stmts) != 4:
        probs.append(f"{len(stmts)} statements, expected 4")
    texts = [o["html"].strip() for o in opts]
    kinds = []
    for i, sysd in enumerate(params["systems"]):
        k = classify(sysd)
        kinds.append(k)
        if k == "other":
            probs.append(f"system {i + 1} ({sysd['kind']}) has no clear shape")
            continue
        if k != sysd["kind"]:
            probs.append(f"system {i + 1}: generated as {sysd['kind']}, looks {k}")
        want = texts.index(TEXT[k])
        if correct["ftype"][i] != want:
            probs.append(f"system {i + 1}: correct option {texts[correct['ftype'][i]]}, numpy says {TEXT[k]}")
        if any(abs(complex(*p)) >= 1 - 1e-9 for p in sysd["poles"]):
            probs.append(f"system {i + 1} has a pole on/outside the unit circle (not stable)")
        if len(sysd["poles"]) < len(sysd["zeros"]):
            probs.append(f"system {i + 1}: more zeros than poles drawn (not causal as drawn)")
    if len(set(kinds)) < 2:
        probs.append(f"fewer than two distinct shapes {kinds}")
    if len(texts) != 4 or sorted(texts) != sorted(TEXT.values()):
        probs.append(f"options {texts}")
    if params["systems"][params["num_label"] - 1]["kind"] not in ("lowpass", "highpass"):
        probs.append("numeric part not on a low-pass / high-pass system")
    if not ({"lowpass", "highpass"} & set(kinds)):
        probs.append("no low-pass or high-pass system")
    # numeric part
    i0 = params["num_label"] - 1
    sysd = params["systems"][i0]
    plotted_z = sorted(round(complex(*z).real, 9) for z in sysd["zeros"] if abs(complex(*z)) > 1e-12)
    plotted_p = sorted(round(complex(*p).real, 9) for p in sysd["poles"] if abs(complex(*p)) > 1e-12)
    n_z = len(sysd["zeros"])
    n_p = len(sysd["poles"])
    if n_p < n_z:
        probs.append(f"system {i0 + 1}: {n_z} zeros but {n_p} poles drawn (H(z) would not be causal)")
    if plotted_z != sorted(round(float(Fraction(z)), 9) for z in params["num_zeros"]) or \
            plotted_p != sorted(round(float(Fraction(q)), 9) for q in params["num_poles"]):
        probs.append("the numeric part's poles/zeros are not the plotted ones")
    if any(abs(complex(*z).imag) > 1e-12 for z in sysd["zeros"] + sysd["poles"]):
        probs.append("numeric part on a system with complex poles/zeros")
    fig = params["stmts"][i0]["fig_html"]
    for v in params["num_zeros"] + params["num_poles"]:
        q = Fraction(v)
        if abs(q) != 1:
            lab = ("\u2212" if q < 0 else "") + (str(abs(q.numerator)) if q.denominator == 1
                                                 else f"{abs(q.numerator)}/{q.denominator}")
            if lab not in fig:
                probs.append(f"value {v} is not labelled on the plot")
    val = _num_value(params, params["at_pi"])
    if abs(val - correct["mag"]) > 1e-9 * max(1, val):
        probs.append(f"|H_d| = {val}, answer {correct['mag']}")
    # the magnitude from the plotted geometry (distances) agrees up to the gain
    w0 = math.pi if params["at_pi"] else 0.0
    geo = float(Fraction(params["K"])) * _mag(sysd, [w0])[0]
    if abs(geo - correct["mag"]) > 1e-9 * max(1, geo):
        probs.append(f"distance product {geo} != answer {correct['mag']}")
    return probs


def _num(x):
    q = Fraction(x).limit_denominator(1000)
    if abs(float(q) - x) < 1e-12:
        return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
    return f"{x:.6g}"


def _differs(a, b):
    return abs(a - b) > 1e-2 * max(abs(b), 1e-3) + 1e-6


def submissions(params, correct):
    ans = correct["ftype"]
    full = {f"ftype-dropdown-{i}": str(a) for i, a in enumerate(ans)}
    i, j = next((i, j) for i in range(len(ans)) for j in range(i + 1, len(ans)) if ans[i] != ans[j])
    swapped = dict(full)
    swapped[f"ftype-dropdown-{i}"], swapped[f"ftype-dropdown-{j}"] = str(ans[j]), str(ans[i])
    n_opts = len(params["ftype"][1])
    shifted = {f"ftype-dropdown-{i}": str((a + 2) % n_opts) for i, a in enumerate(ans)}
    cases = [
        (full, {"ftype": 1}),
        (swapped, {"ftype": 0.5}),                         # two systems swapped: 2 of 4 right
        (shifted, {"ftype": 0}),
    ]
    m = correct["mag"]
    cases.append(({"mag": _num(m)}, {"mag": 1}))
    cases.append(({"mag": f"{m:.4f}"}, {"mag": 1}))
    K = float(Fraction(params["K"]))
    if K != 1 and _differs(m / K, m):
        cases.append(({"mag": _num(m / K)}, {"mag": 0}))   # forgot the gain
    other = _num_value(params, not params["at_pi"])
    if _differs(other, m):
        cases.append(({"mag": _num(other)}, {"mag": 0}))   # evaluated at the other end (z = 1 vs z = -1)
    return cases
