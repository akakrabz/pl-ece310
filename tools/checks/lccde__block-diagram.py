"""Independent check of questions/lccde/block-diagram.

The diagram itself is parsed: gain labels (accent-coloured text) and their triangles' positions
are read from the SVG, and every z^-1 box is counted, so the check covers what the student sees.
Direct form I: triangles pointing right are feed-forward taps (row = delay index), triangles
pointing left are feedback taps g_i. From the taps: b = feed-forward gains, a = [1, -feedback gains];
h by scipy.signal.lfilter; FIR/IIR by running the recursion exactly (Fractions) and checking that
h vanishes from some index on (FIR) or keeps nonzero samples (IIR)."""

import re
from fractions import Fraction

import numpy as np
from scipy.signal import lfilter

DY, Y0 = 78, 46


def _parse(svg):
    tris = []
    for m in re.finditer(r'<polygon points="([\d.]+),([\d.]+) ([\d.]+),([\d.]+) ([\d.]+),([\d.]+)"', svg):
        x1, y1, x2, y2, x3, y3 = map(float, m.groups())
        tris.append((x3 > x1, (y1 + y2) / 2, (x1 + x3) / 2))     # (points right?, y centre, x centre)
    labels = []
    for m in re.finditer(r'<text x="([\d.]+)" y="([\d.]+)"[^>]*fill="#d1501f"[^>]*>([^<]+)</text>', svg):
        labels.append((float(m.group(1)), float(m.group(2)), m.group(3)))
    ndelay = len(re.findall(r"<rect ", svg))
    taps = {"ff": {}, "fb": {}}
    for right, yc, xc in tris:
        lab = [lbl for lx, ly, lbl in labels if abs(lx - xc) < 1 and abs(ly - (yc - 19)) < 1]
        assert len(lab) == 1, (xc, yc, labels)
        val = Fraction(lab[0].replace("−", "-"))
        row = round((yc - Y0) / DY)
        taps["ff" if right else "fb"][row] = val
    return taps, ndelay


def _vectors(params):
    taps, ndelay = _parse(params["fig_html"])
    nff = max(taps["ff"]) if taps["ff"] else 0
    b = [taps["ff"].get(j, Fraction(0)) for j in range(nff + 1)]
    K = max(taps["fb"]) if taps["fb"] else 0
    a = [Fraction(1)] + [-taps["fb"].get(i, Fraction(0)) for i in range(1, K + 1)]
    return b, a, ndelay


def _exact_h(b, a, N):
    h = []
    for n in range(N):
        v = (b[n] if n < len(b) else 0) - sum(a[k] * h[n - k] for k in range(1, len(a)) if n - k >= 0)
        h.append(Fraction(v))
    return h


def check(params, correct):
    probs = []
    b, a, ndelay = _vectors(params)
    gb = np.array(correct["bvec"]["_value"][0], float)
    ga = np.array(correct["avec"]["_value"][0], float)
    if gb.shape != (len(b),) or not np.allclose(gb, [float(v) for v in b]):
        probs.append(f"b {gb} != diagram {b}")
    if ga.shape != (len(a),) or not np.allclose(ga, [float(v) for v in a]):
        probs.append(f"a {ga} != diagram {a}")
    if params["form"] != "I" or ndelay != len(b) - 1 + len(a) - 1:   # direct form I: one box per delayed x and y
        probs.append(f"SVG has {ndelay} delays; direct form I needs {len(b) - 1 + len(a) - 1}")
    x = np.zeros(20)
    x[0] = 1
    h = lfilter([float(v) for v in b], [float(v) for v in a], x)
    for k, n in (("h0", 0), ("h1", 1), ("h2", 2)):
        if abs(h[n] - correct[k]) > 1e-9:
            probs.append(f"{k} = {correct[k]} but lfilter gives {h[n]}")
    he = _exact_h(b, a, 40)
    fir = all(v == 0 for v in he[20:])
    if not fir and all(v == 0 for v in he[10:40]):
        probs.append("impulse response vanishes on 10..39 but not later?")
    marked = correct["firiir"]["html"].strip()
    if marked.startswith("FIR") != fir:
        probs.append(f"marked {marked} but exact recursion says FIR={fir}")
    if params["family"] == "cancel" and len(a) == 1:
        probs.append("cancel family without feedback")
    for v in list(b) + list(a):
        if Fraction(v).denominator not in (1, 2, 4):
            probs.append(f"coefficient {v} is not a short decimal")
    return probs


def _row(v, sep=", "):
    return "[" + sep.join(f"{float(t):g}" for t in v) + "]"


def _far(u, v):
    return abs(u - v) > 0.02 * max(1.0, abs(v))


def _h2_mistakes(b, a, h2):
    """h[2] from the recursion with the feedback signs flipped (a_i used as drawn gains), and with
    the feedback ignored (h[2] = b_2)."""
    x = np.zeros(4)
    x[0] = 1
    out = []
    if len(a) > 1:
        flipped = [1.0] + [-float(v) for v in a[1:]]
        out.append(float(lfilter([float(v) for v in b], flipped, x)[2]))
        out.append(float(b[2]) if len(b) > 2 else 0.0)
    else:
        out.append(float(b[1]) if len(b) > 1 else 0.0)                 # shifted by one sample
        out.append(float(b[0]))                                        # gains read in reverse order
    out.append(-h2)                                                    # sign slip
    res = []
    for w in out:
        if _far(w, h2) and all(_far(w, r) for r in res):
            res.append(w)
    return res


def submissions(params, correct):
    b, a, ndelay = _vectors(params)
    h0, h1, h2 = correct["h0"], correct["h1"], correct["h2"]
    cases = [
        ({"bvec": _row(b, " "), "avec": _row(a, " ")}, {"bvec": 1, "avec": 1}),             # MATLAB style
        ({"avec": "[" + ", ".join(f"{float(t):.2f}" for t in a) + "]"}, {"avec": 1}),      # 0.50 instead of 0.5
        ({"h0": f"{2 * Fraction(h0).limit_denominator(64).numerator}/{2 * Fraction(h0).limit_denominator(64).denominator}"},
         {"h0": 1}),
        ({"h2": f"{Fraction(h2).limit_denominator(256)}"}, {"h2": 1}),                       # exact fraction
        ({"h2": f"{h2:.4g}"}, {"h2": 1}),                                                    # 4 significant digits
        ({"h1": "1/0.5x"}, {"h1": "invalid"}),
    ]
    if len(a) > 1:
        flipped = [a[0]] + [-v for v in a[1:]]
        cases.append(({"avec": _row(flipped)}, {"avec": 0}))                               # feedback signs not flipped
        cases.append(({"avec": _row(a[1:])}, {"avec": 0}))                                 # dropped a0 = 1
    if [float(v) for v in a] != [float(v) for v in b]:
        cases.append(({"bvec": _row(a), "avec": _row(b)}, {"bvec": 0, "avec": 0}))         # b and a swapped
    if params["family"] == "cancel":                                                        # cancelled B/A: not as drawn
        q, r = np.polydiv([float(v) for v in b][::-1], [float(v) for v in a][::-1])
        qz = np.trim_zeros(q[::-1], "b")
        cases.append(({"bvec": _row(qz), "avec": "[1]"}, {"bvec": 0, "avec": 0}))
    wrong_key = [o["key"] for o in params["firiir"] if o["key"] != correct["firiir"]["key"]][0]
    cases.append(({"firiir": wrong_key}, {"firiir": 0}))
    b1 = float(b[1]) if len(b) > 1 else 0.0
    g1 = -float(a[1]) if len(a) > 1 else 0.0
    for wrong in (b1 - g1 * h0, b1):                                                        # sign error / no feedback
        if _far(wrong, h1):
            cases.append(({"h1": repr(wrong)}, {"h1": 0}))
            break
    for w in _h2_mistakes(b, a, h2):
        cases.append(({"h2": repr(w)}, {"h2": 0}))
    return cases
