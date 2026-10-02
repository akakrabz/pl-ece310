"""Independent check of questions/lccde/two-sided-recursion.

h: scipy.signal.residuez on (b, a) + our own side assignment (|p| < 1 right-sided, |p| > 1
left-sided), confirmed by a truncated z-sum at a point of the annulus ROC. The answers alpha, beta
are tested by actually running the backward recursion (and lfilter forward for the causal term) on
a long window and comparing y_c + y_a with the direct two-sided convolution x * h."""

import html
from fractions import Fraction

import numpy as np
from checklib import zsum
from scipy.signal import lfilter, residuez

W = 80                      # window n = -W .. W


def _setup(params):
    b = [float(Fraction(v)) for v in params["b"]]
    a = [float(Fraction(v)) for v in params["a"]]
    r, p, k = residuez(b, a)
    terms = [(complex(ri).real, complex(pi).real) for ri, pi in zip(r, p)]
    inside = [t for t in terms if abs(t[1]) < 1]
    outside = [t for t in terms if abs(t[1]) > 1]
    x0, x1 = (float(Fraction(v)) for v in params["x"])
    return b, a, inside, outside, list(np.atleast_1d(k)), x0, x1


def _h_two_sided(inside, outside, direct, n):
    v = sum(A * p ** n for A, p in inside) if n >= 0 else -sum(A * p ** n for A, p in outside)
    if 0 <= n < len(direct):
        v += float(np.real(direct[n]))
    return v


def _h_causal(inside, outside, n):
    return sum(A * p ** n for A, p in inside + outside) if n >= 0 else 0.0


def _h_swapped(inside, outside, n):
    return sum(A * p ** n for A, p in outside) if n >= 0 else -sum(A * p ** n for A, p in inside)


def _y_from(hfun, x0, x1, n):
    return x0 * hfun(n) + x1 * hfun(n - 1)


def check(params, correct):
    probs = []
    b, a, inside, outside, direct, x0, x1 = _setup(params)
    if len(inside) != 1 or len(outside) != 1:
        return [f"expected one pole inside and one outside the unit circle, got {inside}, {outside}"]
    (A_in, p_in), (A_out, p_out) = inside[0], outside[0]
    if any(abs(d) > 1e-12 for d in direct):
        probs.append(f"unexpected direct terms {direct}")
    h = lambda n: _h_two_sided(inside, outside, direct, n)                  # noqa: E731
    # truncated z-sum at a point of the ROC |p_in| < |z| < |p_out|
    z0 = np.sqrt(abs(p_in) * abs(p_out)) * np.exp(0.9j)
    Hz = np.polyval(b[::-1], 1 / z0) / np.polyval(a[::-1], 1 / z0)
    zs = zsum(h, z0, -400, 400)
    if abs(zs - Hz) > 1e-8 * max(1, abs(Hz)):
        probs.append(f"z-sum {zs} != H(z0) {Hz}")
    # (b): the backward recursion with the submitted-correct alpha, beta reproduces the anti-causal part
    alpha, beta = correct["alpha"], correct["beta"]
    if abs(alpha) >= 1:
        probs.append(f"backward recursion gain |alpha| = {abs(alpha)} >= 1 (not stable backward)")
    n = np.arange(-W, W + 1)
    x = np.where(n == 0, x0, 0.0) + np.where(n == 1, x1, 0.0)
    yc = lfilter([A_in], [1.0, -p_in], x)                                    # forward, from rest
    ya = np.zeros_like(x)                                                    # backward, from rest at n = W
    for i in range(len(n) - 1, 0, -1):
        ya[i - 1] = alpha * ya[i] + beta * x[i]
    hv = np.array([h(m) for m in range(-2 * W, 2 * W + 1)])
    ydirect = np.convolve(x, hv)[2 * W: 2 * W + len(n)]                      # same n axis as x
    if not np.allclose(yc + ya, ydirect, rtol=1e-9, atol=1e-9):
        probs.append("forward + backward recursions do not reproduce x * h")
    for name, m in (("ym1", -1), ("y0", 0), ("y1", 1)):
        want = ydirect[m + W]
        if abs(correct[name] - want) > 1e-9 * max(1, abs(want)):
            probs.append(f"{name} = {correct[name]} but x * h gives {want}")
    # (a): the marked option runs the inside pole forward
    meta = {html.unescape(k).strip(): v for k, v in params["dir_meta"].items()}
    marked = meta.get(html.unescape(correct["dir"]["html"]).strip())
    if marked != ["fwd_in"]:
        probs.append(f"marked direction option {correct['dir']['html']!r} -> {marked}")
    return probs


def _key(params, tag):
    meta = {html.unescape(k).strip(): v for k, v in params["dir_meta"].items()}
    for o in params["dir"]:
        if meta.get(html.unescape(o["html"]).strip()) == [tag]:
            return o["key"]
    raise KeyError(tag)


def _plain(v):
    f = Fraction(v).limit_denominator(10000)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def submissions(params, correct):
    b, a, inside, outside, direct, x0, x1 = _setup(params)
    (A_in, p_in), (A_out, p_out) = inside[0], outside[0]
    names = ["alpha", "beta", "ym1", "y0", "y1"]
    vals = {k: correct[k] for k in names}
    cases = [
        ({k: _plain(vals[k]) for k in names}, {k: 1 for k in names}),           # exact fractions
        ({k: f"{vals[k]:.4g}" for k in names}, {k: 1 for k in names}),           # 4 significant digits
    ]

    def same(u, v):
        return abs(u - v) <= 1e-3 * max(abs(v), 1e-3)

    def ys(hfun):
        return {"ym1": _y_from(hfun, x0, x1, -1), "y0": _y_from(hfun, x0, x1, 0), "y1": _y_from(hfun, x0, x1, 1)}

    # the causal (unstable) system: both terms right-sided
    yc = ys(lambda n: _h_causal(inside, outside, n))
    cases.append(({k: _plain(v) for k, v in yc.items()}, {k: (1 if same(v, vals[k]) else 0) for k, v in yc.items()}))
    # swapped sides: inside pole anti-causal, outside pole causal
    ysw = ys(lambda n: _h_swapped(inside, outside, n))
    cases.append(({k: _plain(v) for k, v in ysw.items()}, {k: (1 if same(v, vals[k]) else 0) for k, v in ysw.items()}))
    cases.append(({"dir": _key(params, "fwd_out")}, {"dir": 0}))
    cases.append(({"alpha": _plain(1 / p_in)}, {"alpha": 0}))                     # backward recursion on the wrong pole
    # sign slip in the backward recursion, and dividing only half the equation (beta = -A_out)
    cases.append(({"beta": _plain(-vals["beta"])}, {"beta": 0}))
    if not same(-A_out, vals["beta"]):
        cases.append(({"beta": _plain(-A_out)}, {"beta": 0}))
    # running the whole difference equation in one direction
    cases.append(({"dir": _key(params, "both_f")}, {"dir": 0}))
    cases.append(({"dir": _key(params, "both_b")}, {"dir": 0}))
    cases.append(({"dir": _key(params, "fwd_in")}, {"dir": 1}))
    cases.append(({"ym1": "y[-1]"}, {"ym1": "invalid"}))
    return cases
