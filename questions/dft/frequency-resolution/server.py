"""Frequency resolution of the DFT (Lecture 23 Sec. 2, eqs. 8-15; Lecture 22 Sec. 1.2).

Two variant types (Midterm 2 family "DFT spectral analysis and resolution"):
  l23   Lecture 23 Exercise 1 / slide Example 1 (guitar E2 = 82 Hz, F2 = 87 Hz at f_s = 8 kHz: full-lobe
        N > 4 pi / |w2 - w1| = 3200 samples, 0.4 s; doubling or halving f_s leaves the time unchanged),
        HW9 #5 (piano notes, f_s = 1.2 kHz) and FA2025 MT2 #8 / SP2025 MT2 #1(f) (Hamming: main lobe twice as
        wide). Two or three equal-amplitude tones in Hz (the closest pair decides), f_s given. Convention
        of the notes (eqs. 9-10, 15: "N > 3200"): the criteria are strict inequalities and the question asks
        for the BOUND: full-lobe N > N_full = 4 pi / dw = 2 f_s / df, half-lobe N > N_half = f_s / df,
        recording time NT > T_full = 2 / df, Hamming full-lobe N > 8 pi / dw = 4 f_s / df; MC: how N and the
        time change when f_s doubles or halves. f_s / df is an integer, so every bound is "exam nice".
  fa23  FA2023 MT2 #6 (a cosine's frequency from a leaky 32-point DFT: tallest bin k = 7, |X[8]| > |X[6]|
        -> lambda0 ~ 7pi/16 and 7pi/16 <= lambda0 <= 15pi/32; the friend's 64-point plot is the record
        zero-padded with 32 zeros), with the FA2025 MT2 #8 alternative "32 samples instead of 16" (taller,
        narrower lobes) and FA2021 MT2 #7(b) (padding gives a finer grid, not a better peak):
        x[n] = cos(lambda0 n), N in {16, 20, 24, 32}, lambda0 = 2 pi (k* + d) / N with |d| in
        [0.15, 0.35]; (a) the estimate 2 pi k*/N, (b) the bracket between bin k* and the half bin on the side
        of the taller neighbour, (c) MC: which processing produced a second, 2N-point plot (zero padding
        with N zeros, or 2N samples of the cosine; distractors: upsampling by 2, a Hamming window).
  rec   FA2021 MT2 #7 setting (cos(10 pi t) sampled at 100 Hz for 5 s, L = 500, zero-padded to N = 1024)
        asked the other way round, with FA2023 MT2 #6(c) and FA2019 MT2 #1(a) (padding does not change the
        DTFT): f_s and the duration D given; the smallest spacing (Hz) of two equal tones that full-lobe and
        half-lobe separation guarantee (2/D, 1/D), the full-lobe bound as a digital frequency (4 pi / L),
        the Hamming full-lobe spacing (4/D), and an MC on what zero padding to N does (nothing).
"""

import cmath
import math
import random
from fractions import Fraction as Fr

from ece310 import fmt, mc, svg, w8_dft as W

NOTES3 = [("C3", 131), ("C#3", 139), ("D3", 146), ("D#3", 156), ("E3", 165), ("F3", 175), ("F#3", 185),
          ("G3", 196), ("G#3", 208), ("A3", 220), ("A#3", 233), ("B3", 247), ("C4", 262), ("C#4", 277),
          ("D4", 294), ("D#4", 311), ("E4", 330), ("F4", 349), ("F#4", 370), ("G4", 392), ("G#4", 415),
          ("A4", 440), ("A#4", 466), ("B4", 494)]


def _khz(fs):
    return f"{fs // 1000}\\ \\text{{kHz}}" if fs % 1000 == 0 else (
        f"{fs / 1000:g}\\ \\text{{kHz}}" if fs % 100 == 0 else f"{fs}\\ \\text{{Hz}}")


def _tones(p):
    """(story, freqs, names or None, fs)."""
    story = random.choice(["generic", "generic", "guitar", "piano"])
    if story == "guitar":
        fs = random.choice([4000, 8000, 8000, 16000, 44100])
        return story, [82, 87, 110], ["E2", "F2", "A2"], fs
    if story == "piano":
        while True:
            i = random.randint(0, len(NOTES3) - 5)
            picks = sorted(random.sample(range(i, i + 5), 3))
            names = [NOTES3[j][0] for j in picks]
            fr = [NOTES3[j][1] for j in picks]
            g1, g2 = fr[1] - fr[0], fr[2] - fr[1]
            df = min(g1, g2)
            fs = random.choice([1200, 2400, 4800, 8000, 9600])
            if g1 != g2 and fs % df == 0 and fr[2] < fs / 2:
                return story, fr, names, fs
    while True:
        df = random.choice([2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 25, 40, 50])
        fs = random.choice([4000, 6000, 8000, 10000, 12000, 16000, 20000, 24000, 44100, 48000])
        if fs % df:
            continue
        f1 = random.randint(10, 200) * random.choice([1, 2, 5])
        if random.random() < 0.3:
            fr = [f1, f1 + df]
        else:
            g = df + random.choice([1, 2, 3, 4, 5, 6, 8, 10, 15, 20]) * max(1, df // 4)
            fr = [f1, f1 + df, f1 + df + g] if random.random() < 0.5 else [f1, f1 + g, f1 + g + df]
        if fr[-1] < fs / 4:
            return story, fr, None, fs


def _l23(p, c):
    story, fr, names, fs = _tones(p)
    gaps = [b - a for a, b in zip(fr, fr[1:])]
    df = min(gaps)
    i = gaps.index(df)
    f1, f2 = fr[i], fr[i + 1]
    for st in ("generic", "guitar", "piano"):
        p["story_" + st] = story == st
    p["n_tones"] = len(fr)
    p["fs"], p["fs_tex"] = fs, _khz(fs)
    if names:
        p["tones_tex"] = ", ".join(f"\\text{{{n}}} \\approx {f}\\ \\text{{Hz}}".replace("#", "\\#") for n, f in zip(names, fr))
    else:
        p["tones_tex"] = ", ".join(f"{f}\\ \\text{{Hz}}" for f in fr)
    p["two_tones"] = len(fr) == 2
    p["f1"], p["f2"], p["df"] = f1, f2, df
    p["gaps_tex"] = ",\\ ".join(f"{b} - {a} = {b - a}" for a, b in zip(fr, fr[1:]))
    Nf, Nh, Tf, Nham = Fr(2 * fs, df), Fr(fs, df), Fr(2, df), Fr(4 * fs, df)
    p["Nf"], p["Nh"], p["Nham"] = int(Nf), int(Nh), int(Nham)
    p["Tf_tex"] = fmt.tex_num(Tf)
    p["Tf_dec"] = f"{float(Tf):.4g}"
    p["Th_tex"] = fmt.tex_num(Fr(1, df))
    p["dw_tex"] = rf"\frac{{2\pi\cdot {df}}}{{ {fs} }}"
    r = random.choice([2, Fr(1, 2)])
    fs2 = fs * r
    p["fs2_tex"] = _khz(int(fs2)) if fs2.denominator == 1 else f"{float(fs2):g}\\ \\text{{Hz}}"
    p["fs_change"] = "doubled" if r == 2 else "halved"
    p["Nf2_tex"] = fmt.tex_num(Fr(2 * fs2, df))
    opts = {
        (2, 1): r"$N_{\text{full}}$ doubles; the minimum recording time is unchanged.",
        (Fr(1, 2), 1): r"$N_{\text{full}}$ halves; the minimum recording time is unchanged.",
        (1, 2): r"$N_{\text{full}}$ is unchanged; the minimum recording time doubles.",
        (1, Fr(1, 2)): r"$N_{\text{full}}$ is unchanged; the minimum recording time halves.",
        (2, 2): r"$N_{\text{full}}$ and the minimum recording time both double.",
        (Fr(1, 2), Fr(1, 2)): r"$N_{\text{full}}$ and the minimum recording time both halve.",
        (1, 1): r"Neither changes: both depend only on the spacing of the tones in hertz.",
    }
    right = opts[(r, 1)]
    wrong = [v for k, v in opts.items() if k != (r, 1)]
    keep = [opts[(1, 1 / r)], opts[(r, r)]]                      # the classic slips always shown
    others = [w for w in wrong if w not in keep]
    p["fs_choices"] = mc.choices(right, keep + random.sample(others, 1), n=4)
    p["Nf2_word"] = "doubles" if r == 2 else "halves"
    c["Nf"], c["Nh"], c["Tf"], c["Nham"] = float(Nf), float(Nh), float(Tf), float(Nham)


# ----------------------------------------------------------------------------- fa23
def _C(w, L):
    s = math.sin(w / 2)
    ratio = L * math.cos(L * w / 2) / math.cos(w / 2) if abs(s) < 1e-12 else math.sin(L * w / 2) / s
    return cmath.exp(-0.5j * w * (L - 1)) * ratio


def _mag(lam, L, M, k):
    """|X[k]| of the M-point DFT of cos(lam n), 0 <= n <= L-1."""
    w = 2 * math.pi * k / M
    return abs(0.5 * (_C(w - lam, L) + _C(w + lam, L)))


WIDTH, LEFT, RIGHT, TOP, BOT, HEIGHT = 680, 48, 34, 30, 30, 230


def _stem(vals, labels, ticks, title):
    M = len(vals)
    ymax = 1.18 * max(vals)
    fig = svg.stem(list(range(M)), [round(v, 4) for v in vals], xlim=(-1, M), ylim=(0, ymax),
                   xticks=[(k, str(k)) for k in ticks], xlabel="k", ylabel="|X[k]|",
                   width=WIDTH, max_width=WIDTH, height=HEIGHT, title=title)

    def X(k):
        return LEFT + (k + 1) * (WIDTH - LEFT - RIGHT) / (M + 1)

    def Y(v):
        return HEIGHT - BOT - v * (HEIGHT - TOP - BOT) / ymax

    extra = "".join(f'<text x="{X(k):.1f}" y="{Y(vals[k]) - 8:.1f}" text-anchor="middle" font-size="11">{vals[k]:.1f}</text>'
                    for k in labels)
    return fig.replace("</svg>", extra + "</svg>"), ymax


def _ticks(M, important):
    """Index labels: the important bins, plus 0 and M-1 when at least 3 bins from every important one."""
    out = set(important)
    for t in (0, M - 1):
        if all(abs(t - u) >= 3 for u in important):
            out.add(t)
    return sorted(out)


def _bracket_safe(N, ks, d):
    """Every frequency within 1.5 bins of bin ks whose N-point DFT has the same tallest bin (ks in
    0..N/2) and the same taller neighbour lies in the bracket [ks, ks + 1/2] or [ks - 1/2, ks] (in bins),
    up to 0.0075 bin: the key's "sure" interval (FA2023 #6(b)) is then really sure, mirror lobe included."""
    side = d > 0
    lo_b, hi_b = (ks, ks + 0.5) if side else (ks - 0.5, ks)
    for i in range(-300, 301):
        kk = ks + i * 1.5 / 300
        lam = 2 * math.pi * kk / N
        v = [_mag(lam, N, N, k) for k in range(N // 2 + 1)]
        if max(range(len(v)), key=lambda k: v[k]) == ks and (v[ks + 1] > v[ks - 1]) == side:
            if not lo_b - 0.0075 <= kk <= hi_b + 0.0075:
                return False
    return True


def _fa23(p, c):
    while True:
        N = random.choice([16, 20, 24, 32])
        ks = random.randint(3, N // 2 - 3)
        d = random.choice([1, -1]) * random.choice([Fr(3, 20), Fr(1, 5), Fr(3, 10), Fr(7, 20)])
        lam = 2 * math.pi * float(ks + d) / N
        vals = [_mag(lam, N, N, k) for k in range(N)]
        half = vals[: N // 2 + 1]
        order = sorted(range(len(half)), key=lambda k: -half[k])
        nb_hi, nb_lo = half[ks + 1], half[ks - 1]
        side_ok = (nb_hi > 1.12 * nb_lo) if d > 0 else (nb_lo > 1.12 * nb_hi)
        if not (order[0] == ks and half[ks] > 1.5 * half[order[1]] and side_ok and _bracket_safe(N, ks, d)):
            continue
        proc = random.choice(["pad", "long"])
        M = 2 * N
        vals2 = [_mag(lam, N if proc == "pad" else M, M, k) for k in range(M)]
        ratio = max(vals2[: N + 1]) / half[ks]
        if (proc == "pad" and ratio <= 1.2) or (proc == "long" and ratio >= 1.55):
            break
    k2 = max(range(N + 1), key=lambda k: vals2[k])
    fig1, ymax1 = _stem(vals, [ks - 1, ks, ks + 1], _ticks(N, [ks, N - ks]), f"DFT magnitude |X[k]|, N = {N}")
    fig2, ymax2 = _stem(vals2, [k2], _ticks(M, [k2, M - k2]), f"The friend's DFT: |X[k]|, k = 0, …, {M - 1}")
    p["fig1_html"], p["fig2_html"] = fig1, fig2
    p["vals1"] = [round(v, 4) for v in vals]
    p["vals2"] = [round(v, 4) for v in vals2]
    p["N"], p["Nm1"], p["M"], p["Mm1"], p["ks"], p["Nmks"] = N, N - 1, M, M - 1, ks, N - ks
    p["d"] = [d.numerator, d.denominator]
    p["proc_kind"] = proc
    p["right"] = d > 0
    p["v_ks"], p["v_lo"], p["v_hi"] = f"{vals[ks]:.1f}", f"{vals[ks - 1]:.1f}", f"{vals[ks + 1]:.1f}"
    p["ksm1"], p["ksp1"] = ks - 1, ks + 1
    p["k2"], p["v_k2"] = k2, f"{vals2[k2]:.1f}"
    p["two_ks"] = 2 * ks
    p["v_2ks"] = f"{vals2[2 * ks]:.1f}"
    est = Fr(2 * ks, N)
    a, b = (est, Fr(2 * ks + 1, N)) if d > 0 else (Fr(2 * ks - 1, N), est)
    p["est_tex"] = W.tex_pi_frac(est)
    p["est_raw_tex"] = rf"\frac{{2\pi\cdot {ks}}}{{ {N} }}"
    p["a_tex"], p["b_tex"] = W.tex_pi_frac(a), W.tex_pi_frac(b)
    p["half_tex"] = rf"\frac{{2\pi\,({ks} {'+' if d > 0 else '-'} \frac12)}}{{ {N} }}"
    p["lam_true_tex"] = W.tex_pi_frac(Fr(2, N) * (ks + d))
    p["lam_true_dec"] = f"{float(Fr(2, N) * (ks + d)):.4g}"
    p["kfrac_tex"] = f"{float(ks + d):g}"
    opts = {
        "pad": rf"She zero-padded the ${N}$ samples with ${N}$ zeros and computed the ${M}$-point DFT.",
        "long": rf"She used ${M}$ samples of the cosine, $\{{x[n]\}}_{{n=0}}^{{ {M - 1} }}$, and computed the ${M}$-point DFT.",
        "up": rf"She upsampled $x[n]$ by $2$ (a zero between consecutive samples) and computed the ${M}$-point DFT.",
        "ham": rf"She multiplied the ${N}$ samples by a Hamming window and computed the ${M}$-point DFT of the result, zero-padded.",
    }
    p["proc_choices"] = mc.choices(opts[proc], [v for k, v in opts.items() if k != proc], n=4)
    p["proc_pad"], p["proc_long"] = proc == "pad", proc == "long"
    c["lam"], c["a"], c["b"] = float(est), float(a), float(b)


def _rec(p, c):
    """FA2021 MT2 #7 setting (x_c(t) sampled at 100 Hz for 5 s, L = 500, zero-padded to N = 1024) asked the
    other way round: which spacing of two equal tones does this record resolve (Lecture 23 eq. 9: "as we
    increase the length of our signal, we may distinguish frequencies that are closer"), and does the
    zero padding help (FA2023 #6(c), FA2019 #1(a): it does not change the DTFT)."""
    while True:
        fs = random.choice([100, 200, 250, 400, 500, 800, 1000, 8000])
        D = random.choice([Fr(1, 4), Fr(1, 2), Fr(1), Fr(2), Fr(5, 2), Fr(4), Fr(5)])
        L = fs * D
        if L.denominator == 1 and 64 <= L <= 4000:
            break
    L = int(L)
    Np = 1 << (L - 1).bit_length()
    Np = Np * 2 if Np == L or random.random() < 0.4 else Np
    p["fs"], p["L"], p["Np"] = fs, L, Np
    p["D_tex"] = fmt.tex_num(D) if D.denominator == 1 else f"{float(D):g}"
    p["D_unit"] = "second" if D == 1 else "seconds"
    p["dff_tex"], p["dfh_tex"], p["dfham_tex"] = fmt.tex_num(2 / D), fmt.tex_num(1 / D), fmt.tex_num(4 / D)
    p["dwf_tex"] = W.tex_pi_frac(Fr(4, L))
    p["bin_tex"] = fmt.tex_num(Fr(fs, Np))
    p["pad_shrunk_tex"] = fmt.tex_num(Fr(2 * fs, Np))
    opts = {
        "same": r"The same as without padding: zero padding does not change the DTFT, so the main lobes keep their width.",
        "lobe": r"It shrinks by the factor $\frac{L}{N}$, because the main lobes become $\frac{4\pi}{N}$ wide.",
        "bin": r"It shrinks to the bin spacing $\frac{f_s}{N}$ of the padded DFT.",
        "grow": r"It grows by the factor $\frac{N}{L}$, because the appended zeros dilute the signal.",
    }
    p["pad_choices"] = mc.choices(opts["same"], [opts["lobe"], opts["bin"], opts["grow"]], n=4)
    c["dff"], c["dfh"], c["dwf"], c["dfham"] = float(2 / D), float(1 / D), float(Fr(4, L)), float(4 / D)


def generate(data):
    p, c = data["params"], data["correct_answers"]
    kind = random.choices(["l23", "fa23", "rec"], weights=[4, 35 / 10, 25 / 10])[0]
    p["l23"], p["fa23"], p["rec"] = kind == "l23", kind == "fa23", kind == "rec"
    {"l23": _l23, "fa23": _fa23, "rec": _rec}[kind](p, c)


def grade(data):
    """Bounds N_full, N_half, N_Ham: the notes write "N > 3200" (enter the bound); the smallest integer N
    (bound + 1) is also accepted, so neither convention loses points."""
    if not data["params"].get("l23"):
        return
    changed = False
    for name in ("Nf", "Nh", "Nham"):
        ps = data["partial_scores"].get(name)
        sub = data["submitted_answers"].get(name)
        if not ps or (ps.get("score") or 0) >= 1 or name in data["format_errors"]:
            continue
        try:
            got = float(sub)
        except (TypeError, ValueError):
            continue
        want = float(data["params"][name]) + 1
        # written so that NaN / inf fail: every comparison with NaN is False
        if abs(got - want) <= 1e-9 * want:
            ps["score"] = 1.0
            changed = True
    if changed:
        import prairielearn as pl
        pl.set_weighted_score_data(data)
