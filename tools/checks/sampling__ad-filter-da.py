"""Independent check of questions/sampling/ad-filter-da.

analysis  x_a(t) is parsed from the displayed LaTeX, sampled for one second (N = f_s samples, so DFT
          bin k is k Hz), the ideal H_d is applied bin by bin at w_k = 2 pi k / N reduced to (-pi, pi],
          and the ideal D/A output amplitude at f < f_s/2 is 2|Y[f]|/N (nothing at or above f_s/2).
design    the largest T is found by brute force: a sampling rate f_s works iff every input frequency
          on a fine grid of (0, B) lands where it should (kept band: unaliased, below f_s/2 and inside
          the digital passband w_c = 2 pi F_c / f_s; removed band: principal |w| outside the passband);
          every rate scanned below the answer must fail and the answer must work. The analog cutoff
          of part (b) is found by scanning tones below f_s/2 through H_d(2 pi f T)."""

import math
import re
from fractions import Fraction

import numpy as np


def _frac_pi(s):
    s = s.strip()
    neg = s.startswith("-")
    s = s.lstrip("-").strip()
    m = re.fullmatch(r"\\tfrac\{(\d*)\\pi\}\{(\d+)\}", s)
    if m:
        v = Fraction(int(m.group(1) or 1), int(m.group(2)))
    else:
        m = re.fullmatch(r"(\d*)\\pi", s)
        assert m, s
        v = Fraction(int(m.group(1) or 1))
    return -v if neg else v


def _tones(x_tex):
    out = []
    for amp, om in re.findall(r"(\d*)\\cos\((\d+)\\pi t\)", x_tex):
        out.append((int(om) / 2, int(amp or 1)))        # (Hz, amplitude)
    return out


def _H(kind, cuts, w):
    """Ideal filter at |w| (radians, array)."""
    w = np.abs(w)
    c = [float(Fraction(v)) * math.pi for v in cuts]
    if kind == "LPF":
        return (w <= c[0]).astype(float)
    if kind == "HPF":
        return (w >= c[0]).astype(float)
    return ((w >= c[0]) & (w <= c[1])).astype(float)


def _check_analysis(an, correct):
    probs = []
    fs = an["fs"]
    tones = _tones(an["x_tex"])
    if sorted(tones) != sorted((t["f"], t["A"]) for t in an["tones"]):
        probs.append(f"displayed x_a {an['x_tex']} != stored tones {an['tones']}")
    N = fs
    n = np.arange(N)
    x = sum(A * np.cos(2 * math.pi * f * n / fs) for f, A in tones)
    X = np.fft.fft(x)
    k = np.arange(N)
    w = 2 * math.pi * k / N
    w = np.where(w > math.pi, w - 2 * math.pi, w)
    Y = X * _H(an["kind"], an["cuts"], w)
    for i, f in enumerate(an["listed"]):
        want = 2 * abs(Y[int(f)]) / N if f < fs / 2 else 0.0
        got = correct[f"a{i + 1}"]
        if abs(want - got) > 1e-6:
            probs.append(f"amplitude at {f} Hz: simulated {want:.6f}, answer {got}")
    # no tone on w = 0 or pi, none on a cutoff edge
    for f, A in tones:
        r = (f / fs) % 1
        if min(abs(r), abs(r - 0.5), abs(r - 1)) < 1e-9:
            probs.append(f"tone {f} Hz lands on w = 0 or pi")
        wp = abs(2 * math.pi * r - (2 * math.pi if r > 0.5 else 0))
        for c in an["cuts"]:
            if abs(wp - float(Fraction(c)) * math.pi) < 1e-9:
                probs.append(f"tone {f} Hz lands on the cutoff")
    nonzero = sum(1 for i in range(5) if correct[f"a{i + 1}"])
    if nonzero < 3:
        probs.append(f"only {nonzero} of the five amplitudes are nonzero (an all-zero guess scores too much)")
    if sum(1 for f in an["listed"] if f >= fs / 2) > 1:
        probs.append("more than one listed frequency at or above f_s/2")
    # every blocked tone's would-be output frequency is listed (so ignoring H_d costs points)
    for f, A in tones:
        r = (f / fs) % 1
        fp = fs * min(r, 1 - r)
        wp = 2 * math.pi * min(r, 1 - r)
        if _H(an["kind"], an["cuts"], np.array([wp]))[0] == 0 and not any(abs(fp - fl) < 1e-9 for fl in an["listed"]):
            probs.append(f"blocked tone {f} Hz: its frequency after sampling, {fp} Hz, is not listed")
    if len(set(an["listed"])) != 5:
        probs.append("listed frequencies not distinct")
    return probs


def _num(tex):
    m = re.search(r"([\d.]+)\\ \\text\{(k?Hz)\}", tex)
    v = float(m.group(1))
    return v * 1000 if m.group(2) == "kHz" else v


def _works(kind, B, Fc, fs, grid):
    """Can a sampling rate fs realize the ideal analog LPF/HPF with cutoff Fc on inputs below B?"""
    wc = 2 * math.pi * Fc / fs
    if Fc >= fs / 2:
        return False
    wraw = 2 * math.pi * grid / fs
    wp = np.abs((wraw + math.pi) % (2 * math.pi) - math.pi)          # principal |w|
    unaliased = grid < fs / 2
    if kind == "LPF":
        keep, drop = grid < Fc, grid > Fc
        return bool(np.all(unaliased[keep] & (wp[keep] < wc)) and np.all(wp[drop] > wc))
    keep, drop = grid > Fc, grid < Fc
    return bool(np.all(unaliased[keep] & (wp[keep] > wc)) and np.all(wp[drop] < wc))


def _check_design(de, correct):
    probs = []
    kind = de["kind"]
    B, Fc = _num(de["B_tex"]), _num(de["Fc_tex"])
    if (B, Fc) != (de["B"], de["Fc"]):
        probs.append("displayed B / F_c differ from the stored values")
    grid = np.linspace(0, B, 40001)[1:-1]
    grid = grid[np.abs(grid - Fc) > B / 20000]
    fs_ans = 1 / correct["Tmax"]
    if not _works(kind, B, Fc, fs_ans * (1 + 1e-9), grid):
        probs.append(f"f_s = {fs_ans} (T = {correct['Tmax']}) does not realize the {kind}")
    for fs in np.linspace(2 * Fc * 1.0005, fs_ans * 0.998, 300):
        if fs < fs_ans * 0.998 and _works(kind, B, Fc, fs, grid):
            probs.append(f"a slower rate f_s = {fs:.2f} Hz also works: T is not the largest")
            break
    wc = 2 * math.pi * Fc / fs_ans / math.pi
    if abs(wc - correct["wc"]) > 1e-9:
        probs.append(f"w_c/pi should be {wc}, answer {correct['wc']}")
    # (b) analog cutoff of a given digital filter, no aliasing
    fs2 = de["fs2"]
    wc2 = float(_frac_pi(de["wc2_tex"])) * math.pi
    f = np.linspace(0, fs2 / 2, 200001)[1:-1]
    passed = (2 * math.pi * f / fs2 <= wc2) if kind == "LPF" else (2 * math.pi * f / fs2 >= wc2)
    edge = f[passed].max() if kind == "LPF" else f[passed].min()
    if abs(edge - correct["fc"]) > fs2 / 2 / 100000:
        probs.append(f"(b) analog cutoff by scan {edge}, answer {correct['fc']}")
    return probs


def check(params, correct):
    if params["is_analysis"]:
        return _check_analysis(params["an"], correct)
    return _check_design(params["de"], correct)


def _fmt(v):
    v = float(v)
    return str(int(v)) if v.is_integer() else repr(v)


def _differs(w, r, atol=1e-6):
    return abs(float(w) - float(r)) > 2 * (atol + 1e-3 * abs(float(r)))


def submissions(params, correct):
    cases = []
    if params["is_analysis"]:
        an = params["an"]
        fs, listed = an["fs"], an["listed"]
        names = [f"a{i + 1}" for i in range(5)]
        right = [correct[nm] for nm in names]
        cases.append(({nm: _fmt(v) + ".0" for nm, v in zip(names, right)}, {nm: 1 for nm in names}))
        cases.append(({nm: f"{2 * int(v)}/2" for nm, v in zip(names, right)}, {nm: 1 for nm in names}))
        tones = _tones(an["x_tex"])
        cuts = [float(Fraction(c)) * fs / 2 for c in an["cuts"]]      # analog cutoffs in Hz
        kind = an["kind"]

        def naive_pass(f):
            if kind == "LPF":
                return f < cuts[0]
            if kind == "HPF":
                return f > cuts[0]
            return cuts[0] < f < cuts[1]

        naive = [sum(A for f, A in tones if f == fl and naive_pass(f)) for fl in listed]     # forgot aliasing
        cases.append(({nm: _fmt(v) for nm, v in zip(names, naive)},
                      {nm: (1 if v == r else 0) for nm, v, r in zip(names, naive, right)}))
        # filtered correctly, but the output put back at the input frequency
        H = []
        for f, A in tones:
            r = (f / fs) % 1
            wp = 2 * math.pi * min(r, 1 - r)
            H.append(float(_H(kind, an["cuts"], np.array([wp]))[0]))
        back = [sum(A * h for (f, A), h in zip(tones, H) if f == fl) for fl in listed]
        cases.append(({nm: _fmt(v) for nm, v in zip(names, back)},
                      {nm: (1 if v == r else 0) for nm, v, r in zip(names, back, right)}))
        allp = [sum(A for (f, A) in tones if abs(fs * min((f / fs) % 1, 1 - (f / fs) % 1) - fl) < 1e-9)
                for fl in listed]                                                         # ignored H_d
        cases.append(({nm: _fmt(v) for nm, v in zip(names, allp)},
                      {nm: (1 if v == r else 0) for nm, v, r in zip(names, allp, right)}))
        cases.append(({nm: "0" for nm in names}, {nm: (1 if r == 0 else 0) for nm, r in zip(names, right)}))
        cases.append(({names[0]: "-"}, {names[0]: "invalid"}))
        return cases
    de = params["de"]
    T, wc, fc = correct["Tmax"], correct["wc"], correct["fc"]
    B, Fc = de["B"], de["Fc"]
    cases.append(({"Tmax": f"1/{round(1 / T)}", "wc": str(Fraction(wc).limit_denominator(50)),
                   "fc": _fmt(fc) + ".0"}, {"Tmax": 1, "wc": 1, "fc": 1}))
    cases.append(({"Tmax": f"{T:.4g}", "wc": f"{wc:.4g}", "fc": f"{int(2 * fc)}/2"}, {"Tmax": 1, "wc": 1, "fc": 1}))
    if de["kind"] == "LPF":
        nyq = 1 / (2 * B)
        cases.append(({"Tmax": f"1/{2 * B}"}, {"Tmax": 0}))                            # insisted on Nyquist
        if _differs(Fc / B, wc):
            cases.append(({"wc": str(Fraction(Fc, B))}, {"wc": 0}))                     # w_c at the Nyquist T
        if _differs(1 / (2 * Fc), T, 1e-12):
            cases.append(({"Tmax": f"1/{2 * Fc}"}, {"Tmax": 0}))                        # Nyquist of the passband only
        del nyq
        cases.append(({"wc": str(Fraction(Fc, B + Fc))}, {"wc": 0}))                    # dropped the 2 of 2 pi F_c T
    else:
        cases.append(({"Tmax": f"1/{B + Fc}"}, {"Tmax": 0}))                            # the low-pass trick
        cases.append(({"Tmax": f"1/{B}"}, {"Tmax": 0}))                                 # Nyquist without the 2
        cases.append(({"wc": str(Fraction(Fc, 2 * B))}, {"wc": 0}))                     # dropped the 2 of 2 pi F_c T
        if _differs(Fraction(2 * Fc, B + Fc), wc):
            cases.append(({"wc": str(Fraction(2 * Fc, B + Fc))}, {"wc": 0}))
    wc2 = float(_frac_pi(de["wc2_tex"]))
    cases.append(({"fc": f"{wc2 / 2:.6g}"}, {"fc": 0}))                                 # w_c/(2 pi) without T
    cases.append(({"fc": f"{wc2 * math.pi * de['fs2']:.6g}"}, {"fc": 0}))               # rad/s instead of Hz
    cases.append(({"fc": _fmt(wc2 * de["fs2"])}, {"fc": 0}))                             # forgot the factor 1/2
    return cases
