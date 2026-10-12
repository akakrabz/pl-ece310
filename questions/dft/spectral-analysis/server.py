"""DFT spectral analysis (Midterm 2 family "DFT spectral analysis and resolution", 6 of 7 past exams).

(a) Bins and hertz (FA2019 MT2 #2: F_s = 40 kHz, N = 1000, bins 300 and 800 -> 12 kHz and -8 kHz,
    spacing 40 Hz; FA2024 #6: which analog frequencies land in bin k when aliasing is allowed):
    f_s and N given; asked the bin spacing f_s/N, the frequency in (-f_s/2, f_s/2] of a bin k_hi > N/2
    (a negative frequency, -(N - k_hi) f_s/N), and the smallest frequency above f_s/2 that also puts a
    peak in bin k_lo (the alias f_s - k_lo f_s/N, through its negative-frequency half).
(b) Reading a DFT magnitude plot (FA2024 #9: 64 Hz for 0.5 s, peaks 48 and 32 -> A = 3, 2 and
    Omega = 12 pi, 32 pi rad/s; SP2023 #10; SP2025 #9): x_c(t) = A1 cos(Omega1 t + phi1) +
    A2 cos(Omega2 t + phi2) sampled at f_s for D seconds (N = 16, 24 or 32 samples) with both
    frequencies on DFT bins and below f_s/2; the plot (ece310.svg.stem) shows |X[k]|, with
    |X[k_i]| = |X[N - k_i]| = N A_i / 2. Asked A1, A2 and Omega1 < Omega2 in rad/s (entered as
    multiples of pi). 2019 course notes, Sec. 2.6 and 2.6.1 (|X[k]| peaks of N A / 2)."""

import cmath
import math
import random

from ece310 import fmt, svg


def _part_a(p):
    while True:
        N = random.choice([200, 250, 400, 500, 800, 1000])
        df = random.choice([5, 8, 10, 16, 20, 25, 40, 50])
        fs = N * df
        if fs <= 50000:
            break
    k_lo = random.randint(N // 20, N // 2 - N // 20)
    k_hi = random.randint(N // 2 + N // 20, N - N // 20)
    p["Na"], p["dfa"], p["fsa"] = N, df, fs
    p["fsa_tex"] = f"{fs // 1000}\\ \\text{{kHz}}" if fs % 1000 == 0 else f"{fs}\\ \\text{{Hz}}"
    p["k_lo"], p["k_hi"] = k_lo, k_hi
    p["Nm_hi"] = N - k_hi
    p["f_hi"] = -(N - k_hi) * df
    p["f_lo"] = k_lo * df
    p["f_alias"] = fs - k_lo * df
    p["fs_half"] = fs // 2
    p["f_hi_abs"] = (N - k_hi) * df
    return df, -(N - k_hi) * df, fs - k_lo * df


def _part_b(p):
    N = random.choice([16, 24, 32])
    df = random.choice([1, 2, 4, 5, 10, 20, 25, 50, 100])
    fs = N * df
    D = fmt.Q(1) / df                                   # seconds: N / fs
    while True:
        k1, k2 = sorted(random.sample(range(2, N // 2), 2))   # k1 >= 2: the labels 0 and k1 never touch
        A1, A2 = random.sample([1, 2, 3, 4, 5], 2)
        if k2 - k1 >= 2:
            break
    ph = [random.choice([0, 1, 2, 3, 4, 5, 6, 7]) / 4 for _ in range(2)]   # phases (multiples of pi), irrelevant to |X|
    mag = []
    for k in range(N):
        X = 0j
        for n in range(N):
            xn = (A1 * math.cos(2 * math.pi * k1 * n / N + math.pi * ph[0])
                  + A2 * math.cos(2 * math.pi * k2 * n / N + math.pi * ph[1]))
            X += xn * cmath.exp(-2j * math.pi * k * n / N)
        mag.append(abs(X))
    h = {k1: N * A1 // 2, k2: N * A2 // 2}
    exact = [h.get(k, h.get(N - k, 0)) if k else 0 for k in range(N)]
    if max(abs(a - b) for a, b in zip(mag, exact)) > 1e-6:
        raise RuntimeError("on-bin cosines leaked")       # cannot happen: whole numbers of periods
    ticks = {0, k1, k2, N - k2, N - k1}                  # pairwise >= 2 bins apart (k1 >= 2, k2 - k1 >= 2, k2 < N/2)
    for extra in (N // 2, N - 1):                        # only where the label cannot fuse with a neighbour
        if min(abs(extra - t) for t in ticks) >= 2:
            ticks.add(extra)
    ticks = sorted(ticks)
    width = max(560, 60 + 18 * N)
    p["xticks"] = ticks
    p["fig_html"] = svg.stem(list(range(N)), exact, xticks=[(k, str(k)) for k in ticks], width=width, max_width=width,
                             yticks=[(v, str(v)) for v in sorted(set(exact) - {0})],
                             values=True, xlabel="k", ylabel="|X[k]|", title=f"DFT magnitude, N = {N}")
    p["Nb"], p["Nbm1"], p["Nb_half"], p["dfb"], p["fsb"] = N, N - 1, N // 2, df, fs
    p["Db_tex"] = fmt.tex_num(D) if D.denominator == 1 else f"{float(D):g}"
    p["Db_unit"] = "second" if D == 1 else "seconds"
    p["k1"], p["k2"], p["h1"], p["h2"] = k1, k2, N * A1 // 2, N * A2 // 2
    p["mk1"], p["mk2"] = N - k1, N - k2
    p["mag"] = exact
    p["W1v"], p["W2v"] = 2 * k1 * df, 2 * k2 * df         # Omega / pi
    p["A1v"], p["A2v"] = A1, A2
    p["Om_per_bin"] = 2 * df
    p["Om_per_bin_tex"] = f"{2 * df}\\pi"
    p["f1"], p["f2"] = k1 * df, k2 * df
    p["halfN"] = N // 2
    return A1, A2, 2 * k1 * df, 2 * k2 * df


def generate(data):
    p, c = data["params"], data["correct_answers"]
    df, f_hi, f_alias = _part_a(p)
    c["df"], c["f_hi"], c["f_alias"] = float(df), float(f_hi), float(f_alias)
    A1, A2, W1, W2 = _part_b(p)
    c["A1"], c["A2"], c["W1"], c["W2"] = float(A1), float(A2), float(W1), float(W2)
