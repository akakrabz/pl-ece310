"""Independent check of questions/dft/periodic-sinc-reading.

|X_d(w)| of x[n] = A cos(w0 n), 0 <= n <= N-1, is summed directly with numpy (no closed form) from the
ANSWERS (read) or the given A, N, w0 (predict):
  read     the plotted polyline must match the direct sum at every vertex; the tick labels must be the
           two zeros of |X_d| next to the main lobe (numerically zero, a sign of the local minimum on both
           sides), the y label must be the main-lobe maximum to within 2 %, the maximum must sit at w0,
           and the labels must not collide;
  predict  the height must be the main-lobe maximum to within 2 %, the zeros must be zeros of |X_d|
           bracketing the maximum with no other zero between them, and the marked MC option must be
           the one with twice the samples (checked by recomputing the main-lobe width of every option)."""

import html
import re

import numpy as np

W = 640
LEFT, RIGHT, TOP, BOT, H = 48, 34, 30, 30, 260


def xd(A, w0, N, w):
    n = np.arange(N)
    w = np.atleast_1d(np.asarray(w, float))
    return np.abs((A * np.cos(w0 * n))[None, :] @ np.exp(-1j * np.outer(n, w)))[0]


def lobe(A, w0, N):
    """(max, argmax) of |X_d| within the main lobe, and the nearest zeros left/right of the maximum."""
    w = np.linspace(w0 - 3 * np.pi / N, w0 + 3 * np.pi / N, 6001)
    v = xd(A, w0, N, w)
    i = int(np.argmax(v))
    j = i
    while j > 0 and v[j - 1] < v[j]:
        j -= 1
    m = i
    while m < len(v) - 1 and v[m + 1] < v[m]:
        m += 1
    return v[i], w[i], w[j], w[m]


def parse_pi(lab):
    lab = lab.replace("−", "-")
    m = re.fullmatch(r"(-?)(\d*)π(?:/(\d+))?", lab)
    if not m:
        return None
    s = -1 if m.group(1) else 1
    num = int(m.group(2)) if m.group(2) else 1
    den = int(m.group(3)) if m.group(3) else 1
    return s * num / den


def check(params, correct):
    probs = []
    N, A = params["N"], params["A"]
    w0 = np.pi * params["r"][0] / params["r"][1]
    Hh = A * N / 2
    if params["read"]:
        w0, N, A = np.pi * correct["w0"], correct["N"], correct["A"]
        if N != params["N"] or A != params["A"]:
            probs.append("answers disagree with params")
    peak, at, zl, zr = lobe(A, w0, N)
    if abs(peak - Hh) > 0.02 * Hh:
        probs.append(f"main-lobe max {peak:.3f} vs AN/2 = {Hh}")
    if abs(at - w0) > 0.02 * 4 * np.pi / N:
        probs.append(f"max at {at / np.pi:.4f} pi, w0 = {w0 / np.pi:.4f} pi")
    if not (w0 - 4 * np.pi / N >= -1e-9 and w0 + 4 * np.pi / N <= np.pi + 1e-9):
        probs.append("main lobe too close to 0 or pi")
    lo, hi = (np.pi * a / b for a, b in params["lohi"])
    for z, name in ((lo, "lo"), (hi, "hi")):
        if xd(A, w0, N, z)[0] > 1e-9 * N * A:
            probs.append(f"{name} = {z / np.pi:.4f} pi is not a zero: {xd(A, w0, N, z)[0]:.3g}")
    if abs(zl - lo) > 2 * np.pi / N * 1e-3 * 6 or abs(zr - hi) > 2 * np.pi / N * 1e-3 * 6:
        probs.append(f"nearest minima {zl / np.pi:.4f}, {zr / np.pi:.4f} pi vs labels {lo / np.pi:.4f}, {hi / np.pi:.4f}")
    if params["read"]:
        fig = params["fig_html"]
        xlabels = re.findall(r'<text x="([\d.]+)" y="[\d.]+" text-anchor="middle">([^<]+)</text>', fig)
        vals = [(float(x), parse_pi(t) if t != "0" else 0.0, t) for x, t in xlabels]
        got = sorted(v for _, v, _ in vals if v is not None)
        want = sorted([0.0, lo / np.pi, hi / np.pi, 1.0])
        if len(got) != 4 or not np.allclose(got, want, atol=1e-9):
            probs.append(f"x tick labels {got} != {want}")
        vals.sort()
        for (x0, _, a), (x1, _, b) in zip(vals, vals[1:]):
            if x1 - x0 < 3.75 * (len(a) + len(b)) + 6:
                probs.append(f"x labels {a} and {b} only {x1 - x0:.1f} px apart")
        if not re.search(rf'text-anchor="end">{int(Hh)}</text>', fig):
            probs.append(f"height {Hh} not written on the y axis")
        pts = re.search(r'<polyline points="([^"]+)"', fig).group(1).split()
        ymax = params["ymax"]
        bad = 0
        for s in pts:
            px, py = map(float, s.split(","))
            w = np.pi * (px - LEFT) / (W - LEFT - RIGHT)
            y = (H - BOT - py) * ymax / (H - TOP - BOT)
            if abs(y - xd(A, w0, N, w)[0]) > 0.6 * ymax / (H - TOP - BOT) + 0.002 * ymax:
                bad += 1
        if bad:
            probs.append(f"{bad} polyline vertices off the direct-sum |X_d|")
        if len(pts) < 600:
            probs.append("polyline too coarse")
    else:
        if abs(correct["Hp"] - Hh) > 1e-9 or abs(correct["lo"] - lo / np.pi) > 1e-12 or abs(correct["hi"] - hi / np.pi) > 1e-12:
            probs.append("predict answers wrong")
        marked = html.unescape(correct["width"]["html"])
        if f"N' = {2 * N}$" not in marked:
            probs.append(f"marked option is not 2N samples: {marked}")
        for o in params["width"]:
            t = html.unescape(o["html"])
            m = re.search(r"N' = (\d+)\$", t)
            if m and o["key"] != correct["width"]["key"]:
                Np = int(m.group(1))
                if abs(4 * np.pi / Np - 2 * np.pi / N) < 1e-12:
                    probs.append(f"distractor {t} also halves the width")
    return probs


def _num(v):
    v = float(v)
    return str(int(round(v))) if abs(v - round(v)) < 1e-12 else f"{v:.6g}"


def submissions(params, correct):
    N, A = params["N"], params["A"]
    Hh = A * N // 2
    r = params["r"]
    (la, lb), (ha, hb) = params["lohi"]
    if params["read"]:
        w0 = correct["w0"]
        return [
            ({"w0": f"{r[0]}/{r[1]}"}, {"w0": 1}),
            ({"w0": f"{2 * r[0]}/{2 * r[1]}", "A": f"{2 * Hh}/{N}"}, {"w0": 1, "A": 1}),
            ({"w0": f"{w0:.4g}", "A": f"{A:.1f}"}, {"w0": 1, "A": 1}),
            ({"N": f"+{N}"}, {"N": 1}),
            ({"w0": f"{la}/{lb}"}, {"w0": 0}),                         # a zero, not the centre
            ({"w0": f"{ha}/{hb}"}, {"w0": 0}),
            ({"w0": _num(np.pi * w0)}, {"w0": 0}),                     # radians in the pi box
            ({"w0": _num(w0 / 2)}, {"w0": 0}),
            ({"N": str(N // 2)}, {"N": 0}) if N % 2 == 0 else ({"N": str(2 * N)}, {"N": 0}),   # width 2pi/N
            ({"N": str(2 * N)}, {"N": 0}),
            ({"A": _num(Hh / N)}, {"A": 0}),                           # height AN
            ({"A": _num(Hh)}, {"A": 0}),                               # raw height
            ({"A": _num(4 * Hh / N)}, {"A": 0}),
            ({"N": f"{N}.0"}, {"N": "invalid"}),
        ]
    H2 = A * N
    cases = [
        ({"Hp": f"{H2}/2"}, {"Hp": 1}),
        ({"Hp": f"{Hh}.0", "lo": f"{2 * la}/{2 * lb}"}, {"Hp": 1, "lo": 1}),
        ({"lo": f"{correct['lo']:.4g}", "hi": f"{correct['hi']:.4g}"}, {"lo": 1, "hi": 1}),
        ({"Hp": _num(H2)}, {"Hp": 0}),                                 # AN (exponential)
        ({"Hp": _num(N)}, {"Hp": 0 if A != 2 else 1}),                 # height N of C alone
        ({"Hp": _num(A / 2)}, {"Hp": 0}),
        ({"lo": _num(correct["lo"] + 1 / N), "hi": _num(correct["hi"] - 1 / N)}, {"lo": 0, "hi": 0}),   # +-pi/N
        ({"lo": _num(correct["lo"] - 2 / N), "hi": _num(correct["hi"] + 2 / N)}, {"lo": 0, "hi": 0}),   # +-4pi/N
        ({"lo": _num(np.pi * correct["lo"])}, {"lo": 0}),
    ]
    for o in params["width"]:
        cases.append(({"width": o["key"]}, {"width": 1 if o["key"] == correct["width"]["key"] else 0}))
    return cases
