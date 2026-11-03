"""Four DFT plots, mixed up (Midterm 2 family "DFT spectral analysis"; Lectures 22-23).

Same shape as FA2025 MT2 #8: x_a(t) = A1 cos(Ω1 t) + A2 cos(Ω2 t) is sampled at T; with N = 16,
  DFT I   = 16-point DFT of {x[n]}_{n=0}^{15},
  DFT II  = those 16 samples zero-padded to 32, 32-point DFT,
  DFT III = 32-point DFT of {x[n]}_{n=0}^{31},
  DFT IV  = Hamming window h[n] = 0.54 - 0.46 cos(2πn/31) times the 32 samples, 32-point DFT.
The four |X[k]| stem plots (ece310.svg.stem, one common y scale with grid lines, 16 or 32 stems) are
shown in random order as plots (a)-(d). (a) one dropdown per plot (DFT I-IV); (b) A1, A2, Ω1, Ω2
(Ω as multiples of π rad/s) read off the plot matched to DFT III, graded with rtol 5 %.
FA2025 key: (a) I, IV, II, III for the exam's order; (b) A1 ≈ 1, A2 ≈ 3/4, Ω1 ≈ 175π/4, Ω2 ≈ 125π/2.

Both tones sit on the 32-point grid (bins k1 < k2 in [3, 13], at least one odd, so DFT I leaks), which
makes the four plots unmistakable and the readings exact: DFT I is the only 16-stem plot; DFT II's even
bins repeat DFT I (same DTFT, evaluated twice as densely); DFT III has single stems of height 16 A_i
(twice DFT II's 8 A_i) and exact zeros elsewhere; DFT IV has three-stem main lobes of about 8.4 A_i and
side lobes below 3 % of its peak.

Variant types (wording and number ranges of each source):
  fa25    FA2025 MT2 #8 verbatim: T = 1/100 s, 0 < Ω_i < 100π, A_i in {1/2, 3/4, 1, 5/4, 3/2, 2},
          |k2 - k1| >= 3; 20 % of them are the key's own reading (k = 7, 10; A = 1, 3/4).
  fa24    FA2024 MT2 #9 sampling and reading wording ("sampled at a rate of f_s Hz for D seconds to
          produce {x[n]}_{n=0}^{31}", "A1 and A2 are integers and A1 > A2"): f_s in {32, 64, 80, 128,
          160, 320} Hz, A in {1..5}; 20 % are the exam's numbers (64 Hz, 0.5 s, 3 at k = 3, 2 at k = 8).
  sp23    SP2023 MT2 #10 wording ("composed of spectral components", T = 1/12,000 s, "sampled above
          the Nyquist rate"): T = 1/F with F in {4000, 6000, 8000, 12000, 16000, 24000}, A in
          {1/2, 3/4, 1, 3/2} (the exam's 1/2, 1, 3/4).
  lecture Lecture 22 eqs. 16-17 (s[n] = A1 cos(ω1 n) + A2 cos(ω2 n), A1 = 1, A2 = 1/2 there) as a
          digital signal: ω_i asked as multiples of π (rad/sample), A in {1/2, 1, 3/2, 2}.
In fa25, sp23 and lecture the tones are numbered by frequency (Ω1 < Ω2); in fa24 by amplitude."""

import math
import random
from fractions import Fraction as F

import numpy as np

from ece310 import fmt, svg

DFTS = ["I", "II", "III", "IV"]
LETTERS = ["a", "b", "c", "d"]


def _bins(min_gap):
    while True:
        k1, k2 = sorted(random.sample(range(3, 14), 2))
        if k2 - k1 >= min_gap and (k1 % 2 or k2 % 2):
            return k1, k2


def _amps(pool, ratio_max):
    while True:
        a, b = random.sample(pool, 2)
        if max(a, b) / min(a, b) <= ratio_max:
            return F(a), F(b)


def _ystep(peaks):
    """Largest tick step in {16, 8, 4} dividing every DFT III peak, with at most 8 ticks."""
    for step in (16, 8, 4):
        if all(p % step == 0 for p in peaks) and max(peaks) / step <= 8:
            return step
    return None


def _setup(p):
    u = random.random()
    kind = "fa25" if u < 0.35 else "fa24" if u < 0.6 else "sp23" if u < 0.8 else "lecture"
    p["kind"] = kind
    for k in ("fa25", "fa24", "sp23", "lecture"):
        p["is_" + k] = kind == k
    p["analog"] = kind != "lecture"
    exam_case = False
    while True:
        if kind == "fa25":
            inv_T = 100
            if random.random() < 0.2:
                (k1, k2), (A1, A2), exam_case = (7, 10), (F(1), F(3, 4)), True
            else:
                k1, k2 = _bins(3)
                A1, A2 = _amps([F(1, 2), F(3, 4), F(1), F(5, 4), F(3, 2), F(2)], F(8, 3))
        elif kind == "fa24":
            inv_T = random.choice([32, 64, 64, 80, 128, 160, 320])
            if random.random() < 0.2:
                inv_T, (k1, k2), (A1, A2), exam_case = 64, (3, 8), (F(3), F(2)), True
            else:
                k1, k2 = _bins(4)
                A1, A2 = _amps([1, 2, 3, 4, 5], 3)
        elif kind == "sp23":
            inv_T = random.choice([4000, 6000, 8000, 12000, 12000, 16000, 24000])
            k1, k2 = _bins(4)
            A1, A2 = _amps([F(1, 2), F(3, 4), F(1), F(3, 2)], 2)
        else:
            inv_T = None
            k1, k2 = _bins(4)
            A1, A2 = _amps([F(1, 2), F(1), F(3, 2), F(2)], 3)
        if _ystep([16 * A1, 16 * A2]):
            break
    return kind, inv_T, k1, k2, A1, A2, exam_case


def _dfts(k1, k2, A1, A2):
    n = np.arange(32)
    x = float(A1) * np.cos(2 * np.pi * k1 * n / 32) + float(A2) * np.cos(2 * np.pi * k2 * n / 32)
    ham = 0.54 - 0.46 * np.cos(2 * np.pi * n / 31)
    raw = {"I": np.fft.fft(x[:16]), "II": np.fft.fft(np.concatenate([x[:16], np.zeros(16)])),
           "III": np.fft.fft(x), "IV": np.fft.fft(ham * x)}
    out = {}
    for k, v in raw.items():
        m = np.round(np.abs(v), 3)
        m[m < 5e-3] = 0.0
        out[k] = [float(a) for a in m]
    return out, float(ham.sum())


def _figure(vals, ytop, step, letter):
    M = len(vals)
    if M == 16:
        xt = [(k, str(k)) for k in range(16)]
    else:
        xt = [(k, str(k) if k % 2 == 0 else "") for k in range(32)]
    yt = [(v, str(v)) for v in range(step, ytop + 1, step)]
    return svg.stem(list(range(M)), vals, ylim=(0, ytop * 1.06), xticks=xt, yticks=yt, xlabel="k", ylabel="|X[k]|",
                    width=600, height=200, max_width=600, grid=True, title=f"Plot ({letter})")


def _pi_tex(q):
    """Multiple of π as LaTeX: \\tfrac{175}{4}\\pi, 12\\pi."""
    q = F(q)
    if q == 1:
        return r"\pi"
    return fmt.tex_num(q) + r"\pi"


def _num(v, d=2):
    s = f"{v:.{d}f}".rstrip("0").rstrip(".")
    return s


def generate(data):
    p, c = data["params"], data["correct_answers"]
    kind, inv_T, k1, k2, A1, A2, exam_case = _setup(p)
    p["exam_case"] = exam_case
    D, sum_h = _dfts(k1, k2, A1, A2)
    order = random.sample(DFTS, 4)                   # order[i] = DFT shown as plot LETTERS[i]
    peaks = [16 * A1, 16 * A2]
    step = _ystep(peaks)
    ytop = int(max(peaks))
    plot_of = {}
    for letter, name in zip(LETTERS, order):
        p["fig_" + letter + "_html"] = _figure(D[name], ytop, step, letter)
        p["vals_" + letter] = D[name]
        p["ch_" + letter] = [{"text": f"DFT {d}", "correct": d == name} for d in DFTS]
        plot_of[name] = letter
    p["order"] = order
    p["plot_of"] = plot_of
    p["ystep"] = step
    # ---- the signal and its description
    p["k1"], p["k2"], p["mk1"], p["mk2"] = k1, k2, 32 - k1, 32 - k2
    p["kA"] = {"k1": k1, "k2": k2, "A1": float(A1), "A2": float(A2)}
    if kind == "lecture":
        W1, W2 = F(k1, 16), F(k2, 16)                # ω / π
    else:
        W1, W2 = F(k1 * inv_T, 16), F(k2 * inv_T, 16)    # Ω / π = 2 k / (32 T) / π
        p["inv_T"] = inv_T
        p["T_tex"] = rf"\tfrac{{1}}{{{inv_T}}}" if kind != "sp23" else rf"\tfrac{{1}}{{{inv_T:,}}}".replace(",", "{,}")
        p["Om_bin_tex"] = _pi_tex(F(inv_T, 16))       # Ω per bin of the 32-point DFT
        p["nyq_tex"] = _pi_tex(F(inv_T))              # π / T
    if kind == "fa24":
        p["fs"] = inv_T
        D_s = F(32, inv_T)
        p["D_tex"] = fmt.tex_num(D_s) if D_s.denominator == 1 else f"{float(D_s):g}"
        p["D_unit"] = "second" if D_s == 1 else "seconds"
        p["fs_half"] = _num(inv_T / 2)
    # numbering: by frequency (fa25, sp23, lecture) or by amplitude, A1 > A2 (fa24)
    tones = [(k1, A1, W1), (k2, A2, W2)]
    if kind == "fa24":
        tones.sort(key=lambda t: -t[1])
    (ka, Aa, Wa), (kb, Ab, Wb) = tones
    c["A1"], c["A2"], c["W1"], c["W2"] = float(Aa), float(Ab), float(Wa), float(Wb)
    p["ka"], p["kb"], p["mka"], p["mkb"] = ka, kb, 32 - ka, 32 - kb
    p["ha"], p["hb"] = fmt.tex_num(16 * Aa), fmt.tex_num(16 * Ab)
    p["Aa_tex"], p["Ab_tex"] = fmt.tex_num(Aa), fmt.tex_num(Ab)
    p["Wa_tex"], p["Wb_tex"] = _pi_tex(Wa), _pi_tex(Wb)
    if kind != "lecture":
        p["fa_hz"], p["fb_hz"] = _num(float(Wa) / 2), _num(float(Wb) / 2)
    # ---- worked solution: one line per plot
    I, II, III, IV = D["I"], D["II"], D["III"], D["IV"]
    m_show = max(range(1, 8), key=lambda m: I[m])
    odd = [k for k in (k1, k2) if k % 2]
    lo, hi = sorted([odd[0] // 2, odd[0] // 2 + 1])
    why = {
        "I": (rf"the only plot with $16$ bins ($k = 0, \dots, 15$): $16$ samples, $16$-point DFT. "
              rf"The tone at $k = {odd[0]}$ of the 32-point grid falls between its bins ${lo}$ and ${hi}$ "
              rf"($\omega = \tfrac{{2\pi\cdot{odd[0]}}}{{32}} = \tfrac{{2\pi\cdot{odd[0] / 2:g}}}{{16}}$) and leaks into every bin."),
        "II": (rf"$32$ bins whose even bins repeat plot ({plot_of['I']}): $\lvert X[2m]\rvert$ here equals $\lvert X[m]\rvert$ there "
               rf"(e.g. ${_num(II[2 * m_show])}$ at $k = {2 * m_show}$ vs ${_num(I[m_show])}$ at $k = {m_show}$). Zero padding evaluates the "
               rf"same DTFT of the same $16$ samples twice as densely; the peaks stay near $\tfrac{{16A}}{{2}} = 8A$ "
               rf"(${_num(II[k1])}$ at $k = {k1}$, ${_num(II[k2])}$ at $k = {k2}$) and the side lobes stay."),
        "III": (rf"$32$ bins, single stems of heights ${fmt.tex_num(16 * A1)}$ and ${fmt.tex_num(16 * A2)}$ at $k = {k1}, {k2}$ "
                rf"(mirrors ${32 - k1}, {32 - k2}$) and exact zeros elsewhere: $32$ samples hold whole periods of both tones, and the "
                rf"peaks $\tfrac{{32A}}{{2}} = 16A$ are twice those of DFT II: more data, taller and narrower lobes."),
        "IV": (rf"$32$ bins, three-stem main lobes (wider: $\tfrac{{8\pi}}{{32}}$) with almost nothing beyond them (side lobes near "
               rf"$-42$ dB); the peaks are scaled by $\sum_n h[n] / 32 = {sum_h / 32:.2f}$: ${_num(IV[k1])}$ at $k = {k1}$ "
               rf"against ${fmt.tex_num(16 * A1)}$ in DFT III."),
    }
    p["sol_items"] = [{"letter": L, "dft": name, "why": why[name]} for L, name in zip(LETTERS, order)]
    p["III_letter"], p["II_letter"], p["IV_letter"] = plot_of["III"], plot_of["II"], plot_of["IV"]
    p["sum_h"] = f"{sum_h:.2f}"
    p["IV_peak_a"] = _num(IV[ka])
    p["II_peak_a"] = _num(II[ka])
    p["dk_even"] = (k2 - k1) % 2 == 0
    p["dk_odd"] = not p["dk_even"]
