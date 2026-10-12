"""Independent check of questions/sampling/nyquist-rate.

Sums and products of sinusoids: the band edge is measured with an FFT of one second of the signal
(all tones are integer Hz, so every tone sits exactly on a bin). Sinc and abstract families: the
CTFT rules (box, modulation, product = convolution of spectra, convolution = product, time scaling)
recomputed here from the stored numbers. The displayed formula is parsed back to its frequencies,
and the checkbox options are re-derived from f_s > 2B (strict, Lecture 17)."""

import math
import re
from fractions import Fraction as F

import numpy as np


def _fft_edge(freqs_fn, fmax):
    fsn = int(10 * fmax) + 10
    t = np.arange(fsn) / fsn
    x = freqs_fn(t)
    X = np.abs(np.fft.rfft(x)) / fsn
    nz = np.nonzero(X > 1e-6 * X.max())[0]
    return F(int(nz.max()))


def band_edge(d):
    fam = d["fam"]
    f = [F(v) for v in d["f"]]
    if fam == "sum":
        return _fft_edge(lambda t: sum(np.cos(2 * np.pi * float(v) * t + 0.3) for v in f), float(max(f)))
    if fam == "product":
        f3 = F(d["f3"]) if "f3" in d else None
        fn = (lambda t: np.cos(2 * np.pi * float(f[0]) * t) * np.sin(2 * np.pi * float(f[1]) * t)
              + (np.cos(2 * np.pi * float(f3) * t) if f3 else 0))
        return _fft_edge(fn, float(f[0] + f[1] + (f3 or 0)))
    if fam == "sinc":
        return f[0]
    if fam == "sinc2":
        return f[0] + f[0]                       # box * box
    if fam == "sincmod":
        return f[0] + F(d["f0"])                 # box shifted to +-Omega_0
    if fam == "sincprod":
        return f[0] + f[1]
    if fam == "sincconv":
        return min(f)
    if fam == "scale":
        return F(d["c"]) * f[0]
    if fam == "square":
        return f[0] + f[0]
    if fam == "modulate":
        return f[0] + F(d["f0"])
    if fam == "absprod":
        return f[0] + F(d["f1"])
    if fam == "absconv":
        return min(f[0], F(d["f1"]))
    raise KeyError(fam)


def shown_freqs(tex):
    """Frequencies (Hz) written in a formula: '300\\pi' -> 150, '2\\pi\\cdot 150' -> 150."""
    out = [F(m) for m in re.findall(r"2\\pi\\cdot (\d+)", tex)]
    rest = re.sub(r"2\\pi\\cdot \d+", "", tex)
    out += [F(m) / 2 for m in re.findall(r"(\d+)\\pi", rest)]
    out += [F(m) for m in re.findall(r"= (\d+)\$ Hz", tex)]
    return sorted(set(out))


def _rate(html):
    return F(re.search(r"f_s = ([\d.]+)", html).group(1))


def check(params, correct):
    probs = []
    d = params["data"]
    B = band_edge(d)
    if abs(correct["B"] - float(B)) > 1e-9 * float(B):
        probs.append(f"B = {correct['B']} but recomputed {B} ({d})")
    if abs(correct["fN"] - 2 * float(B)) > 1e-9 * float(B):
        probs.append("Nyquist rate is not 2B")
    if abs(correct["TN"] - 1 / (2 * float(B))) > 1e-12:
        probs.append("T_N is not 1/(2B)")
    expect = set()
    for k in ("f", "f0", "f1", "f3"):
        if k in d:
            expect |= {F(v) for v in (d[k] if isinstance(d[k], list) else [d[k]])}
    shown = set(shown_freqs(params["sig_tex"] + " " + params["given_html"]))
    if not expect <= shown:
        probs.append(f"frequencies {sorted(expect)} not all shown in the formula ({sorted(shown)})")
    opts = params["rates"]
    if len(opts) != 5:
        probs.append("expected 5 rate options")
    marked = {o["html"] for o in correct["rates"]}
    for o in opts:
        r = _rate(o["html"])
        if (r > 2 * B) != (o["html"] in marked):
            probs.append(f"option {o['html']} marked {o['html'] in marked} but f_s > 2B is {r > 2 * B}")
    if not any(_rate(o["html"]) == 2 * B for o in opts):
        probs.append("no option exactly at the Nyquist rate")
    return probs


def _wrong_B(d, B):
    """Family-specific wrong bandwidths."""
    fam = d["fam"]
    f = [F(v) for v in d["f"]]
    if fam == "sum":
        w = [min(f)]
    elif fam == "product":
        w = [max(f + ([F(d["f3"])] if "f3" in d else []))]
    elif fam == "sinc":
        w = [2 * f[0]]                               # read the coefficient of pi as hertz
    elif fam in ("sinc2", "square"):
        w = [f[0]]                                   # forgot that squaring doubles the band
    elif fam in ("sincmod", "modulate"):
        w = [F(d["f0"]), f[0]]                       # only the carrier / only the baseband edge
    elif fam in ("sincprod", "sincconv"):
        w = [max(f)] if fam == "sincprod" else [max(f), f[0] + f[1]]
    elif fam == "scale":
        w = [f[0] / F(d["c"])]                       # divided instead of multiplied
    elif fam == "absprod":
        w = [max(f[0], F(d["f1"]))]
    else:
        w = [max(f[0], F(d["f1"])), f[0] + F(d["f1"])]
    return [v for v in w if abs(v - B) > F(1, 50) * B]


def _s(q):
    q = F(q)
    return str(q.numerator) if q.denominator == 1 else repr(float(q))


def submissions(params, correct):
    d = params["data"]
    B = band_edge(d)
    cases = [
        ({"B": f"{float(B)}0", "fN": f"{2 * float(B):.6g}"}, {"B": 1, "fN": 1}),            # decimals
        ({"B": f"{2 * B.numerator}/{2 * B.denominator}"}, {"B": 1}),                          # unreduced fraction
        ({"B": repr(2 * math.pi * float(B))}, {"B": 0}),                                      # rad/s, not Hz
        ({"fN": _s(B)}, {"fN": 0}),                                                           # B instead of 2B
        ({"fN": repr(4 * math.pi * float(B))}, {"fN": 0}),                                    # 2 Omega_max in rad/s
        ({"TN": f"1/{_s(2 * B)}"}, {"TN": 1}),                                                # exact fraction
        ({"TN": f"{1 / (2 * float(B)):.4g}"}, {"TN": 1}),                                     # 4 significant digits
        ({"TN": f"{1 / (2 * float(B)):.6e}"}, {"TN": 1}),                                     # scientific notation
        ({"TN": f"1/{_s(B)}"}, {"TN": 0}),                                                    # 1/B
        ({"TN": _s(2 * B)}, {"TN": 0}),                                                       # the rate, not the period
        ({"TN": "1/0"}, {"TN": "invalid"}),
    ]
    for wb in _wrong_B(d, B):
        cases.append(({"B": _s(wb)}, {"B": 0}))
    opts = params["rates"]
    good = [o["key"] for o in opts if _rate(o["html"]) > 2 * B]
    at = [o["key"] for o in opts if _rate(o["html"]) == 2 * B]
    below = [o["key"] for o in opts if _rate(o["html"]) < 2 * B]
    c = len(good)
    cases.append(({"rates": good}, {"rates": 1}))
    cases.append(({"rates": list(reversed(good))}, {"rates": 1}))
    cases.append(({"rates": good + at}, {"rates": (c - 1) / c}))                              # counted f_s = 2B
    cases.append(({"rates": good + at + below}, {"rates": max(0.0, (c - 1 - len(below)) / c)}))
    cases.append(({"rates": at + below}, {"rates": 0}))
    return cases
