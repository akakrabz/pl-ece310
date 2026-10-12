"""Digital filtering of analog signals: ideal A/D (period T) -> ideal digital filter H_d(w) -> ideal
D/A (same T). Midterm 2 family "A/D -> H_d -> D/A" (6 of 7 past exams), taught after Lecture 17
(tag `ahead`). Two variant families:

analysis  x_a(t) = sum of 4-5 zero-phase cosines (written in rad/s, integer-Hz frequencies) and an
          ideal LPF / HPF / BPF with cutoffs at multiples of pi that never coincide with a tone.
          Each tone lands at w = 2 pi f T, reduced to its principal value |w| in (0, pi); it survives
          iff H_d(|w|) = 1 and comes out at |w| / (2 pi T) Hz; tones landing on the same |w| add.
          Every variant has three passed tones at three distinct output frequencies (at least one of
          them aliased), one or two blocked tones, and in a little over half of the variants a
          fourth passed tone that folds onto one of the outputs (amplitudes add). The student gives
          the amplitude in y_a(t) at five listed frequencies (Hz): the three outputs (nonzero), the
          would-be output frequency of every blocked tone (0), and, if one slot is left, one
          original frequency above f_s/2 (0). An all-zero answer scores 2/5; ignoring H_d, ignoring
          aliasing, or putting tones back at their input frequency scores at most 4/5. Variants where
          one of these, or the elimination guess (0 in the stopband and at or above f_s/2, a direct
          tone's amplitude at its own frequency, a most common amplitude of the tones above f_s/2
          elsewhere), would be fully right are rejected.
          Models: SP2025 #5(a), FA2024 #5(b), SP2023 #7, FA2021 #4, FA2019 #9(b).
design    x_a(t) band-limited to B Hz (X_a = 0 for |Omega| >= 2 pi B). For an overall ideal analog
          low-pass with cutoff F_c (aliasing allowed outside the passband) the largest period is
          T = 1/(B + F_c) and w_c = 2 pi F_c T (SP2025 #5(b), FA2024 #5(c), FA2023 #5(a): w_c = 2pi/3 when
          F_c = B/2); for a high-pass no aliasing is allowed, T = 1/(2B) and w_c = 2 pi F_c T
          (FA2023 #5(b)). Also the analog cutoff w_c / (2 pi T) of a given digital filter
          (FA2019 #9(b), 2019 notes section 5.3: H_eff(Omega) = H_d(Omega T)).
"""

import random
from fractions import Fraction as F

from ece310 import fmt, svg
from ece310.w7_analog import chain_svg, plain, principal, rad_tex, tex_pi

# ----------------------------------------------------------------------------- analysis family
# Tones sit on the grid w = 2 pi j / D (j = 1 .. D/2 - 1); every filter below has at least three grid
# points in its passband and at least one in its stopband, and no cutoff on a grid point.
FS_FOR = {10: [500, 1000, 2000, 3000, 5000, 10000],
          12: [600, 1200, 2400, 3600, 6000, 12000],
          16: [800, 1600, 3200, 4800, 8000, 16000]}
SPECS = {
    10: [("LPF", (F(2, 3),)), ("LPF", (F(3, 4),)), ("HPF", (F(1, 3),)), ("HPF", (F(1, 4),))],
    12: [("LPF", (F(3, 4),)), ("HPF", (F(1, 4),)), ("BPF", (F(1, 4), F(3, 4)))],
    16: [("LPF", (F(3, 5),)), ("LPF", (F(2, 3),)), ("LPF", (F(4, 5),)), ("HPF", (F(1, 5),)), ("HPF", (F(1, 3),)),
         ("HPF", (F(2, 5),)), ("BPF", (F(1, 5), F(3, 5))), ("BPF", (F(1, 3), F(2, 3))), ("BPF", (F(1, 6), F(5, 6))),
         ("BPF", (F(2, 5), F(4, 5)))],
}


def _hd(kind, cuts, w):
    """Ideal filter value at |w| (multiple of pi, 0 < w < 1); cutoffs never coincide with tones."""
    if kind == "LPF":
        return 1 if w < cuts[0] else 0
    if kind == "HPF":
        return 1 if w > cuts[0] else 0
    return 1 if cuts[0] < w < cuts[1] else 0


def _hd_tex(kind, cuts):
    a = tex_pi(cuts[0])
    if kind == "LPF":
        return (rf"H_d(\omega) = \begin{{cases}} 1, & \lvert\omega\rvert \le {a} \\ 0, & {a} \lt \lvert\omega\rvert \le \pi"
                r"\end{cases}")
    if kind == "HPF":
        return (rf"H_d(\omega) = \begin{{cases}} 0, & \lvert\omega\rvert \lt {a} \\ 1, & {a} \le \lvert\omega\rvert \le \pi"
                r"\end{cases}")
    b = tex_pi(cuts[1])
    return (rf"H_d(\omega) = \begin{{cases}} 1, & {a} \le \lvert\omega\rvert \le {b} \\ 0, & \text{{otherwise in }} [-\pi, \pi]"
            r"\end{cases}")


def _hd_svg(kind, cuts):
    c = [float(v) for v in cuts]
    if kind == "LPF":
        pts = [(-1, 0), (-c[0], 0), (-c[0], 1), (c[0], 1), (c[0], 0), (1, 0)]
    elif kind == "HPF":
        pts = [(-1, 1), (-c[0], 1), (-c[0], 0), (c[0], 0), (c[0], 1), (1, 1)]
    else:
        pts = [(-1, 0), (-c[1], 0), (-c[1], 1), (-c[0], 1), (-c[0], 0), (c[0], 0), (c[0], 1), (c[1], 1), (c[1], 0),
               (1, 0)]
    ticks = [(-1, "−π"), (1, "π"), (0, "0")]
    for v in cuts:
        ticks += [(float(v), svg.fmt_pi(v)), (-float(v), svg.fmt_pi(-v))]
    fig = svg.plot([pts], xlim=(-1.15, 1.15), ylim=(0, 1.35), xticks=ticks, yticks=[(1, "1")],
                   xlabel="ω", ylabel="Hd(ω)", width=520, height=200, max_width=520,
                   title="ideal digital filter, one period")
    # subscript d in the axis label (svg.plot escapes its labels, so patch the text node)
    return fig.replace(">Hd(ω)<", '>H<tspan font-size="10" dy="4">d</tspan><tspan dy="-4">(ω)</tspan><', 1)


def _tone(j, g, fs, alias):
    """An analog frequency landing on grid point j: j g itself, or k f_s +- j g (k = 1, 2) when aliased."""
    if not alias:
        return j * g
    return random.choice([1, 1, 2]) * fs + random.choice([1, -1]) * j * g


def _analysis():
    while True:
        D = random.choice([10, 12, 16, 16])
        fs = random.choice(FS_FOR[D])
        g = fs // D
        kind, cuts = random.choice(SPECS[D])
        cuts = list(cuts)
        pas = [j for j in range(1, D // 2) if _hd(kind, cuts, F(2 * j, D))]
        stop = [j for j in range(1, D // 2) if not _hd(kind, cuts, F(2 * j, D))]
        extra = random.choices(["plain", "fold", "two_blocked"], weights=[30, 30, 40])[0]
        if extra == "two_blocked" and len(stop) < 2:
            extra = "fold"
        jp = random.sample(pas, 3)                         # three distinct output frequencies
        jb = random.sample(stop, 2 if extra == "two_blocked" else 1)
        specs = [(j, True) for j in jp] + [(j, False) for j in jb]
        if extra == "fold":
            specs.append((random.choice(jp), True))        # a second tone folding onto an output frequency
        tones = [(_tone(j, g, fs, random.random() < 0.55), random.choice([1, 2, 3, 4]), j, ok) for j, ok in specs]
        fr = [t[0] for t in tones]
        if len(set(fr)) < len(tones):
            continue
        if not any(f > fs / 2 and ok for f, A, j, ok in tones) or sum(1 for f in fr if f > fs / 2) < 2:
            continue                                       # aliasing must matter for a nonzero answer
        if extra == "fold":
            pair = [t for t in tones if t[2] == specs[-1][0]]
            if all(f > fs / 2 for f, A, j, ok in pair) and random.random() < 0.5:
                continue                                   # prefer folds of a direct tone and an alias
        out = {}                                           # output frequency (Hz) -> amplitude
        for f, A, j, ok in tones:
            if ok:
                out[j * g] = out.get(j * g, 0) + A
        listed = sorted(out) + sorted({j * g for j in jb})
        if len(listed) < 5:                                # one original frequency above f_s/2
            big = [f for f, A, j, ok in tones if f > fs / 2 and ok] or [f for f, A, j, ok in tones if f > fs / 2]
            listed.append(random.choice(big))
        listed = sorted(listed)
        amps = [out.get(f, 0) for f in listed]
        # the classic wrong answers must really be wrong somewhere
        fc_an = [c * fs / 2 for c in cuts]                 # analog cutoffs in Hz
        naive = {}
        for f, A, j, ok in tones:                          # forgot aliasing: filter the analog frequency itself
            okn = (f < fc_an[0]) if kind == "LPF" else (f > fc_an[0]) if kind == "HPF" else (fc_an[0] < f < fc_an[1])
            if okn:
                naive[f] = naive.get(f, 0) + A
        orig = {}                                          # filtered correctly, "the D/A returns the input frequency"
        allp = {}                                          # ignored H_d
        for f, A, j, ok in tones:
            if ok:
                orig[f] = orig.get(f, 0) + A
            allp[j * g] = allp.get(j * g, 0) + A
        if any([naive.get(f, 0) for f in listed] == amps or [d.get(f, 0) for f in listed] == amps for d in (orig, allp)):
            continue
        # elimination guess without computing any alias: 0 in the stopband and at or above f_s/2, a direct
        # tone's amplitude at its own frequency, and a most common amplitude of the tones above f_s/2 in the
        # remaining slots (every tied mode counts)
        direct = {}
        for f, A, j, ok in tones:
            if f < fs / 2:
                direct[f] = direct.get(f, 0) + A
        high = [A for f, A, j, ok in tones if f > fs / 2]
        top = max(high.count(a) for a in high)

        def _elim(m):
            return [0 if (fl >= fs / 2 or not _hd(kind, cuts, F(2 * fl, fs))) else direct.get(fl, m) for fl in listed]

        if any(_elim(m) == amps for m in {a for a in high if high.count(a) == top}):
            continue
        break

    order = list(range(len(tones)))
    random.shuffle(order)
    x_tex = " + ".join(
        (f"{tones[i][1]}" if tones[i][1] != 1 else "") + rf"\cos({rad_tex(tones[i][0])} t)" for i in order)
    rows = []
    for i in order:
        f, A, j, ok = tones[i]
        w = F(2 * f, fs)
        wp = principal(w)
        rows.append({
            "f_tex": str(f), "w_tex": tex_pi(w), "wp_tex": tex_pi(abs(wp)), "A": A,
            "alias": "yes" if f > fs / 2 else "no", "H": "1" if ok else "0",
            "out_tex": (str(j * g) + r"\ \text{Hz}") if ok else r"\text{removed}",
        })
    y_terms = [(out[f], rf"\cos({rad_tex(f)} t)") for f in sorted(out)]
    fcs = [fmt.tex_num(c * fs / 2) for c in cuts]
    folds = [f for f in out if sum(1 for t in tones if t[3] and t[2] * g == f) > 1]
    p = {
        "is_analysis": True, "kind": kind, "cuts": [plain(c) for c in cuts], "fs": fs, "D": D,
        "tones": [{"f": f, "A": A} for f, A, j, ok in tones],
        "listed": listed,
        "x_tex": x_tex, "fs_tex": str(fs), "T_tex": rf"T = \tfrac{{1}}{{{fs}}}\ \text{{s}}",
        "hd_tex": _hd_tex(kind, cuts), "hd_svg_html": _hd_svg(kind, cuts),
        "kind_word": {"LPF": "low-pass", "HPF": "high-pass", "BPF": "band-pass"}[kind],
        "half_tex": str(fs // 2),
        "rows": rows,
        "y_tex": fmt.tex_sum(y_terms),
        "fc_an_phrase_html": (rf"the analog cutoffs ${fcs[0]}$ and ${fcs[1]}$ Hz" if kind == "BPF"
                              else rf"the analog cutoff ${fcs[0]}$ Hz"),
        "fold_html": (rf"Two tones land on the same digital frequency and come out together at ${folds[0]}$ Hz, "
                      r"so their amplitudes add. " if folds else ""),
        "labels": [{"name": f"a{i + 1}", "f": f} for i, f in enumerate(listed)],
        "ans_list": [{"f": f, "A": a} for f, a in zip(listed, amps)],
    }
    ans = {f"a{i + 1}": float(a) for i, a in enumerate(amps)}
    return p, ans


# ----------------------------------------------------------------------------- design family
RATIOS = [F(1, 7), F(1, 5), F(1, 4), F(1, 3), F(1, 2), F(3, 5), F(2, 3), F(3, 4)]
BS = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 2100, 2400, 3000, 4000, 5000, 6000,
      8000, 10000, 12000]
FS2 = [400, 600, 800, 1000, 1200, 1600, 2000, 3000, 4000, 6000, 8000, 10000]
WC2 = [F(1, 6), F(1, 5), F(1, 4), F(1, 3), F(2, 5), F(1, 2), F(3, 5), F(2, 3), F(3, 4), F(4, 5), F(5, 6)]


def _hz_tex(f):
    f = F(f)
    if f.denominator == 1 and f >= 1000 and f % 1000 == 0:
        return rf"{f // 1000}\ \text{{kHz}}"
    return rf"{fmt.tex_num(f)}\ \text{{Hz}}"


def _design():
    kind = random.choice(["LPF", "LPF", "HPF"])
    while True:
        r = random.choice(RATIOS)
        B = random.choice(BS)
        Fc = r * B
        if Fc.denominator == 1:
            break
    if kind == "LPF":
        fs_min = B + Fc
        wc = 2 * Fc / fs_min                       # multiple of pi
    else:
        fs_min = 2 * B
        wc = Fc / B
    while True:
        fs2 = random.choice(FS2)
        wc2 = random.choice(WC2)
        fc2 = wc2 * fs2 / 2
        if fc2.denominator == 1 and fc2 != Fc:
            break
    kind_word = "low-pass" if kind == "LPF" else "high-pass"
    p = {
        "is_design": True, "kind": kind, "kind_word": kind_word, "is_lpf": kind == "LPF", "is_hpf": kind == "HPF",
        "B": int(B), "Fc": int(Fc),
        "B_tex": _hz_tex(B), "Fc_tex": _hz_tex(Fc), "OmB_tex": rad_tex(B),
        "fsmin_tex": str(int(fs_min)), "Tmax_tex": rf"\tfrac{{1}}{{{int(fs_min)}}}",
        "wc_tex": tex_pi(wc), "wc_plain": plain(wc),
        "BplusFc_tex": str(int(B + Fc)), "twoB_tex": str(int(2 * B)),
        "Tnyq_tex": rf"\tfrac{{1}}{{{int(2 * B)}}}",
        "fs2": fs2, "wc2": plain(wc2), "fs2_tex": str(fs2), "T2_tex": rf"T = \tfrac{{1}}{{{fs2}}}\ \text{{s}}",
        "wc2_tex": tex_pi(wc2), "fc2_tex": str(int(fc2)), "Omc2_tex": rad_tex(fc2),
    }
    ans = {"Tmax": float(F(1) / fs_min), "wc": float(wc), "fc": float(fc2)}
    return p, ans


def generate(data):
    p = data["params"]
    p["chain_html"] = chain_svg(True)
    if random.random() < 0.5:
        pa, ans = _analysis()
        p["an"] = pa
        p["is_analysis"], p["is_design"] = True, False
    else:
        pd, ans = _design()
        p["de"] = pd
        p["is_analysis"], p["is_design"] = False, True
    data["correct_answers"].update(ans)
