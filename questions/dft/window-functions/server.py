"""Window functions (Midterm 2 family "DFT spectral analysis"; Lecture 23 §2-3.1, eqs. 8-21, slides 7-13).

Facts used (L23 eqs. 20-21, slide 13): rectangular window main lobe 4π/N wide, tallest side lobe -13 dB
(about 0.22 of the main-lobe peak); Hamming h[n] = 0.54 - 0.46 cos(2πn/(N-1)), 0 <= n <= N-1: main lobe
8π/N, tallest side lobe -42 dB (about 0.008).

Past-exam shapes and lecture examples reproduced (one variant type each):
(a) dB <-> amplitude ratio (L23 lecture of Oct 28: "-13 dB is something like between 0.1 and 0.2",
    "-42 dB is under 1/100", the 1/50 tone). (i) a side-lobe level of -L dB as a fraction of the main
    lobe, 10^(-L/20): type "lecture" = the rectangular (13 dB) or Hamming (42 dB) window, type "other" =
    another level; (ii) an amplitude ratio 1/m in dB, 20 log10(1/m) (m = 50 is the lecture's tone).
(b) A weak tone beside a strong one: type "lecture" = L23 eq. 16 / Figures 2 and 4 exactly,
    cos(0.5πn) + (1/10)cos(0.625πn), N = 64 (rectangular masks it, Hamming shows it); type "random":
    A2/A1 in {1/10, 1/20, 1/50, 1/2}, |ω2 - ω1| >= 0.25π and > 2·8π/N. Checkbox "with which window does
    the weaker tone's main lobe rise above the stronger tone's tallest side lobe": the window-functions
    page's rule A_weak/A_strong > 10^(-13/20) ≈ 0.22 (rectangular), > 10^(-42/20) ≈ 0.008 (Hamming); every
    ratio is at least a factor 2.2 away from both thresholds.
(c) (i) Hamming main-lobe width 8π/N (as a multiple of π); (ii) smallest N for full-lobe separation,
    rectangular N > 4π/Δω, Hamming N > 8π/Δω (L23 eqs. 8-10, slide 13): type "digital" = two equal tones
    at ω1, ω2 = ω1 + Δω (site quick check 0.30π / 0.34π shape); type "notes" = L23 Exercise 1 (guitar
    chord at 8 kHz: N > 2 f_s/Δf) with two adjacent notes from the HW9 piano table. The bounds 4π/Δω and
    8π/Δω are never integers (Δf or 100·Δω/π has a factor 3, 7, 11 or 13), so "smallest N" has no
    strict/non-strict ambiguity.
(d) type "sp25" = SP2025 MT2 #1(f) True/False, "Applying the Hamming window for spectral analysis leads
    to lower sidelobes and narrower main lobes in spectral components." (False), weighted 40 %, and its
    three siblings (lower/wider True, higher/narrower False, higher/wider False); type "fa19" = FA2019
    MT2 #1(c) True/False on the DTFT of a finite piece of a cosine (the exam's impulse formula with
    ω0 = π/32, N = 16 and random ω0, N: False; the shifted-window form ½R_d(ω ∓ ω0): True); type
    "choice" = L23 lecture of Oct 28, "which window would you choose" for a stated two-tone signal
    (close tones of similar amplitude, 4π/N < Δω < 8π/N: rectangular; well-separated tones with a
    weak one, 1/20 ... 1/50: Hamming).
Window statements on past exams: SP2025 #1(f) only (T/F bank Q40); FA2023 #1(g) and FA2019 #1(f) are
DFT-leakage statements, FA2021 #7 is peak location after padding (other questions)."""

import math
import random
from decimal import Decimal
from fractions import Fraction as F

from ece310 import fmt, mc

RECT_DB, HAM_DB = 13, 42

# HW9 piano table (Hz), octaves 3 and 4
NOTES = [("C3", 131), ("C♯3", 139), ("D3", 146), ("D♯3", 156), ("E3", 165), ("F3", 175), ("F♯3", 185),
         ("G3", 196), ("G♯3", 208), ("A3", 220), ("A♯3", 233), ("B3", 247), ("C4", 262), ("C♯4", 277),
         ("D4", 294), ("D♯4", 311), ("E4", 330), ("F4", 349), ("F♯4", 370), ("G4", 392), ("G♯4", 415),
         ("A4", 440), ("A♯4", 466), ("B4", 494)]


def dec(q) -> str:
    """Exact decimal of a Fraction whose denominator is 2^a 5^b: 0.625, 0.08, 3."""
    q = F(q)
    if q.denominator == 1:
        return str(q.numerator)
    s = format(Decimal(q.numerator) / Decimal(q.denominator), "f")
    return s.rstrip("0").rstrip(".") if "." in s else s


def pi_dec(q) -> str:
    """0.5\\pi, \\pi (q = 1)."""
    return r"\pi" if F(q) == 1 else dec(q) + r"\pi"


def pi_frac(q) -> str:
    r"""\tfrac{\pi}{32}, \tfrac{3\pi}{8}, \pi."""
    q = F(q)
    num = r"\pi" if q.numerator == 1 else f"{q.numerator}\\pi"
    return num if q.denominator == 1 else rf"\tfrac{{{num}}}{{{q.denominator}}}"


def coef(a) -> str:
    """Amplitude in front of \\cos: '' for 1, else the decimal."""
    return "" if F(a) == 1 else dec(a)


def sig3(x: float) -> str:
    """3 significant digits, plain: 0.224, 0.00794, -27.96 -> -28.0."""
    if x == 0:
        return "0"
    d = 2 - int(math.floor(math.log10(abs(x))))
    v = round(x, d)
    return f"{v:.{max(d, 0)}f}"


# ----------------------------------------------------------------------------- (a)
def part_a(p, c):
    if random.random() < 0.35:
        win, L = random.choice([("rectangular", RECT_DB), ("Hamming", HAM_DB)])
        p["a_lect"], p["a_win"] = True, win
    else:
        L = random.choice([15, 18, 22, 25, 26, 28, 31, 32, 35, 38, 45, 46, 48, 52, 55, 58])
        p["a_lect"] = False
    p["a_other"] = not p["a_lect"]
    r = 10 ** (-L / 20)
    p["aL"], p["a_ratio_s"] = L, sig3(r)
    c["ratio"] = r
    m = random.choice([4, 5, 6, 8, 12, 15, 25, 30, 40, 50, 50, 60, 80, 150, 200, 250, 300, 400, 500, 800])
    db = 20 * math.log10(1 / m)
    p["am"], p["a_db_s"] = m, sig3(db)
    p["a_log_s"] = f"{math.log10(m):.4f}"
    p["a_lecture_m"] = m == 50
    c["db"] = db


# ----------------------------------------------------------------------------- (b)
def part_b(p):
    if random.random() < 0.2:                       # L23 eq. 16, Figures 2 and 4
        A_s, ratio, w_s, w_w, N = F(1), F(1, 10), F(1, 2), F(5, 8), 64
        p["b_lect"] = True
        sig = r"\cos(0.5\pi n) + \tfrac{1}{10}\cos(0.625\pi n)"
    else:
        p["b_lect"] = False
        N = random.choice([64, 100, 128, 200, 256])
        A_s = F(random.choice([1, 2, 4, 5, 10]))
        ratio = random.choice([F(1, 10), F(1, 20), F(1, 50), F(1, 2)])
        grid = [F(k, 20) for k in range(3, 18)]     # 0.15 ... 0.85
        while True:
            w_s, w_w = random.sample(grid, 2)
            if abs(w_w - w_s) >= F(1, 4) and abs(w_w - w_s) * N > 16:   # more than two Hamming lobes apart
                break
        A_w = A_s * ratio
        terms = sorted([(w_s, A_s), (w_w, A_w)])
        sig = " + ".join(f"{coef(a)}\\cos({pi_dec(w)} n)" for w, a in terms)
    A_w = A_s * ratio
    p["b_sig"], p["bN"], p["bNm1"] = sig, N, N - 1
    p["b_As"], p["b_Aw"] = dec(A_s), (r"\tfrac{1}{10}" if p["b_lect"] else dec(A_w))
    p["b_ws"], p["b_ww"] = pi_dec(w_s), pi_dec(w_w)
    p["b_ratio_tex"] = fmt.tex_num(ratio)
    p["b_ratio"] = float(ratio)
    p["b_ratio_db"] = f"{20 * math.log10(float(ratio)):.0f}"
    p["b_dw"] = pi_dec(abs(w_w - w_s))
    rect_ok = float(ratio) > 10 ** (-RECT_DB / 20)
    ham_ok = float(ratio) > 10 ** (-HAM_DB / 20)
    p["b_rect_ok"], p["b_rect_no"] = rect_ok, not rect_ok
    p["b_choices"] = [{"text": "rectangular window", "correct": rect_ok},
                      {"text": "Hamming window", "correct": ham_ok},
                      {"text": "neither window", "correct": not (rect_ok or ham_ok)}]
    p["b_spec"] = {"A_s": float(A_s), "A_w": float(A_w), "w_s": float(w_s), "w_w": float(w_w), "N": N}


# ----------------------------------------------------------------------------- (c)
def part_c(p, c):
    Nc = random.choice([16, 20, 25, 32, 40, 50, 64, 80, 100, 128, 160, 200, 256])
    p["cN"] = Nc
    p["c_wham_tex"] = fmt.tex_num(F(8, Nc))
    p["c_wrect_tex"] = fmt.tex_num(F(4, Nc))
    c["W_ham"] = float(F(8, Nc))
    if random.random() < 0.5:
        p["c_dig"], p["c_notes"] = True, False
        d = random.choice([F(3, 100), F(7, 100), F(9, 100), F(11, 100), F(3, 50), F(3, 25), F(3, 20), F(7, 200),
                           F(9, 200), F(13, 200), F(3, 40), F(9, 50), F(7, 50)])
        w1 = F(random.randint(4, 12), 20)
        p["c_w1"], p["c_w2"], p["c_dw"] = pi_dec(w1), pi_dec(w1 + d), pi_dec(d)
        p["c_dw_frac"] = fmt.tex_num(d)
        b_rect, b_ham = 4 / d, 8 / d                       # in units where Δω = d·π
        p["c_b_rect"] = fmt.tex_num(b_rect)
        p["c_b_ham"] = fmt.tex_num(b_ham)
        p["c_spec"] = {"kind": "digital", "w1": float(w1), "w2": float(w1 + d)}
    else:
        p["c_dig"], p["c_notes"] = False, True
        while True:
            i = random.randrange(len(NOTES) - 1)
            (n1, f1), (n2, f2) = NOTES[i], NOTES[i + 1]
            df = f2 - f1
            if any(df % q == 0 for q in (3, 7, 11, 13, 17, 19, 23)):
                break
        fs = random.choice([4000, 8000, 8000, 16000])
        p["c_n1"], p["c_n2"], p["c_f1"], p["c_f2"], p["c_df"] = n1, n2, f1, f2, df
        p["c_fs_tex"] = f"{fs // 1000}\\ \\text{{kHz}}"
        p["c_fs"] = fs
        b_rect, b_ham = F(2 * fs, df), F(4 * fs, df)
        p["c_b_rect"] = fmt.tex_num(b_rect)
        p["c_b_ham"] = fmt.tex_num(b_ham)
        p["c_spec"] = {"kind": "notes", "f1": f1, "f2": f2, "fs": fs}
        p["c_T_rect"] = f"{float((math.floor(b_rect) + 1) / fs):.3g}"
    N_rect = math.floor(b_rect) + 1
    N_ham = math.floor(b_ham) + 1
    assert b_rect.denominator != 1 and b_ham.denominator != 1
    p["c_b_rect_dec"] = f"{float(b_rect):.2f}"
    p["c_b_ham_dec"] = f"{float(b_ham):.2f}"
    p["c_Nrect"], p["c_Nham"] = N_rect, N_ham
    c["N_rect"], c["N_ham"] = N_rect, N_ham


# ----------------------------------------------------------------------------- (d)
TF_SP25 = [("lower", "narrower", 4), ("lower", "wider", 2), ("higher", "narrower", 2), ("higher", "wider", 2)]
CHOICE_TEXT = {
    "rect_main": "Rectangular window: its narrower main lobe is needed to separate the two tones.",
    "ham_side": "Hamming window: its lower side lobes are needed to reveal the weaker tone.",
    "ham_main": "Hamming window: its narrower main lobe is needed to separate the two tones.",
    "rect_side": "Rectangular window: its lower side lobes are needed to reveal the weaker tone.",
}


def part_d(p):
    u = random.random()
    p["d_sp25"] = p["d_fa19"] = p["d_choice"] = False
    if u < 0.4:
        p["d_sp25"] = True
        side, main = random.choices([(s, m) for s, m, _ in TF_SP25], weights=[w for *_, w in TF_SP25])[0]
        truth = side == "lower" and main == "wider"
        p["d_side"], p["d_main"] = side, main
        p["d_exam_stmt"] = side == "lower" and main == "narrower"
        p["d_statement"] = (f"Applying the Hamming window for spectral analysis leads to {side} sidelobes and "
                            f"{main} main lobes in spectral components.")
        p["d_choices"] = mc.true_false(truth)
        p["d_spec"] = {"kind": "sp25", "side": side, "main": main}
        p["d_truth_word"] = "True" if truth else "False"
    elif u < 0.6:
        p["d_fa19"] = True
        if random.random() < 0.4:
            w0, N = F(1, 32), 16
        else:
            w0 = random.choice([F(1, 16), F(1, 8), F(1, 5), F(1, 4), F(1, 3), F(3, 8), F(1, 2), F(2, 3)])
            N = random.choice([8, 12, 16, 20, 24, 32, 50, 64])
        impulses = random.random() < 0.6
        w0t = pi_frac(w0)
        p["d_w0"], p["d_N"], p["d_Nm1"] = w0t, N, N - 1
        p["d_exam_stmt"] = impulses and w0 == F(1, 32) and N == 16
        if impulses:
            p["d_formula"] = (rf"X_d(\omega) = \pi\left[\delta\left(\omega - {w0t}\right) + "
                              rf"\delta\left(\omega + {w0t}\right)\right]")
        else:
            p["d_formula"] = (rf"X_d(\omega) = \tfrac12 R_d\left(\omega - {w0t}\right) + "
                              rf"\tfrac12 R_d\left(\omega + {w0t}\right)")
        p["d_impulses"], p["d_shifted"] = impulses, not impulses
        p["d_choices"] = mc.true_false(not impulses)
        p["d_spec"] = {"kind": "fa19", "w0": float(w0), "N": N, "impulses": impulses}
        p["d_truth_word"] = "False" if impulses else "True"
        p["d_main_w"] = fmt.tex_num(F(4, N))
    else:
        p["d_choice"] = True
        if random.random() < 0.5:                    # close tones, similar amplitudes -> rectangular
            while True:
                N = random.choice([40, 50, 64, 80, 100, 128, 200])
                d = F(random.randint(2, 20), 100)
                if 5 <= d * N <= 7:
                    break
            A1, A2 = random.choice([(3, 2), (4, 3), (5, 4), (5, 3), (2, F(3, 2)), (5, F(7, 2))])
            w1 = F(random.randint(4, 12), 20)
            w2 = w1 + d
            key = "rect_main"
            p["d_close"], p["d_weak"] = True, False
        else:                                        # well-separated, weak tone -> Hamming
            while True:
                N = random.choice([64, 100, 128, 200, 256])
                grid = [F(k, 20) for k in range(3, 18)]
                w1, w2 = sorted(random.sample(grid, 2))
                if w2 - w1 >= F(1, 4) and (w2 - w1) * N >= 16:
                    break
            A1 = F(random.choice([1, 2, 4, 5, 10]))
            A2 = A1 * random.choice([F(1, 20), F(1, 25), F(1, 40), F(1, 50)])
            if random.random() < 0.5:
                A1, A2 = A2, A1                      # weak tone at the lower frequency
            key = "ham_side"
            p["d_close"], p["d_weak"] = False, True
            d = w2 - w1
        A1, A2 = F(A1), F(A2)
        p["d_sig"] = f"{coef(A1)}\\cos({pi_dec(w1)} n) + {coef(A2)}\\cos({pi_dec(w2)} n)"
        p["d_N"], p["d_Nm1"] = N, N - 1
        p["d_dw"] = pi_dec(d)
        p["d_rect_w"] = fmt.tex_num(F(4, N))
        p["d_ham_w"] = fmt.tex_num(F(8, N))
        p["d_rect_w_dec"] = dec(F(4, N))             # N = 2^a 5^b: finite decimals
        p["d_ham_w_dec"] = dec(F(8, N))
        small, big = sorted([A1, A2])
        p["d_ratio_tex"] = fmt.tex_num(small / big)
        p["d_ratio_db"] = f"{20 * math.log10(float(small / big)):.1f}"
        p["d_choices"] = [{"text": CHOICE_TEXT[k], "correct": k == key} for k in CHOICE_TEXT]
        random.shuffle(p["d_choices"])
        p["d_spec"] = {"kind": "choice", "A1": float(A1), "A2": float(A2), "w1": float(w1), "w2": float(w2), "N": N}
        p["d_exam_stmt"] = False
    p["d_tf"] = p["d_sp25"] or p["d_fa19"]


def generate(data):
    p, c = data["params"], data["correct_answers"]
    part_a(p, c)
    part_b(p)
    part_c(p, c)
    part_d(p)
