"""Ideal D/A reconstruction (sinc interpolation), the zero-order hold, and which analog cosines an
ideal A/D -> ideal D/A chain maps onto a given output (Midterm 2 family "sampling"; taught after
Lecture 17, so tagged `ahead`).

Models: FA2019 MT2 #8 (D/A of delta[n-3]: ZOH and ideal D/A), SP2025 MT2 #1(c) (x_a(3) = 0 for
delta[n-2] + 2 delta[n-4], T = 1 s), FA2023 #1(e) / FA2024 #1(c) (one sample -> one sinc),
SP2021 MT2 #9 (which inputs give y_a(t) = 2cos(300 pi t) at F_s = 400 Hz), SP2023 #8 (a phase fixes
the sign of the alias). Ideal D/A as in the 2019 notes eq. (3.3) and the keys:
y_a(t) = sum_n x[n] sin(pi(t - nT)/T) / (pi(t - nT)/T); ZOH g(t) = u(t) - u(t - T).

(a) x[n] = A cos(w1 n + theta) with w1 = w0 (principal), w0 + 2pi or 2pi k - w0 (0 < w0 < pi):
    y_a(t) = A cos(2 pi f0 t + phi), f0 = w0 / (2 pi T), phi = theta (or -theta for 2pi k - w0).
(b) a short x[n] (c delta[n-k], two deltas, or three samples): the ZOH value at t1 = (m1 + alpha)T,
    the ideal D/A value at a sample instant m2 T (= x[m2]) and halfway, at (m3 + 1/2)T, where
    y_a = (2/pi) sum_n x[n] (-1)^(m3-n) / (2(m3-n) + 1)  (exact; entered as a decimal).
(c) checkbox: inputs B cos(2 pi F t + phi) sampled at f_s whose ideal D/A output equals
    B cos(2 pi f0 t + psi): exactly F = k f_s + f0 with phi = psi, or F = k f_s - f0 with phi = -psi.
"""

import math
import random
from fractions import Fraction as F

from ece310 import fmt
from ece310.w7_analog import plain, principal, rad_tex, tex_pi

# ----------------------------------------------------------------------------- (a)
W0 = [F(1, 6), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(5, 6), F(1, 5), F(2, 5), F(3, 5), F(4, 5),
      F(1, 8), F(3, 8), F(5, 8), F(7, 8)]
FS = [200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400, 3000, 4000, 6000, 8000]
THETA = [s * a for a in (F(1, 6), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(5, 6)) for s in (1, -1)]


def _phase_tex(th):
    """' + \\tfrac{\\pi}{4}' / ' - \\tfrac{\\pi}{3}' / '' for 0."""
    if th == 0:
        return ""
    return (" - " if th < 0 else " + ") + tex_pi(abs(th))


def _part_a():
    c = random.choice(W0)
    fs = random.choice([f for f in FS if (c * f / 2).denominator == 1])
    A = random.choice([2, 3, 4, 5, 6])
    th = random.choice(THETA)
    rep = random.choices(["principal", "plus", "minus", "minus2"], weights=[22, 26, 36, 16])[0]
    if rep == "principal":
        w1, phi = c, th
    elif rep == "plus":
        w1, phi = c + 2, th
    elif rep == "minus":
        w1, phi = 2 - c, -th
    else:
        w1, phi = 4 - c, -th
    f0 = c * fs / 2
    p = {
        "A": A, "fs": fs, "w1": plain(w1), "theta": plain(th), "rep": rep,
        "x_tex": rf"{A}\cos\left({tex_pi(w1)} n{_phase_tex(th)}\right)",
        "T_tex": rf"T = \tfrac{{1}}{{{fs}}}\ \text{{s}}",
        "w0_tex": tex_pi(c), "w1_tex": tex_pi(w1), "th_tex": tex_pi(th), "phi_tex": tex_pi(phi),
        "f0_tex": str(f0), "fs_tex": str(fs), "Om0_tex": rad_tex(f0),
        "y_tex": rf"{A}\cos\left({rad_tex(f0)} t{_phase_tex(phi)}\right)",
        "is_principal": rep == "principal",
        "is_alias": rep != "principal",
    }
    if rep == "plus":
        p["reduce_tex"] = (rf"{tex_pi(w1)} - 2\pi = {tex_pi(c)}, \text{{ so }} "
                           rf"x[n] = {A}\cos\left({tex_pi(c)} n{_phase_tex(th)}\right)")
        p["reduce_why_html"] = r"Adding $2\pi$ to a digital frequency does not change the samples: $\cos\big((\omega_0 + 2\pi)n + \theta\big) = \cos(\omega_0 n + \theta)$ for every integer $n$."
    elif rep in ("minus", "minus2"):
        k2 = 2 if rep == "minus" else 4
        p["reduce_tex"] = (rf"{tex_pi(w1)} - {tex_pi(k2)} = -{tex_pi(c)}, \text{{ so }} "
                           rf"x[n] = {A}\cos\left(-{tex_pi(c)} n{_phase_tex(th)}\right) = "
                           rf"{A}\cos\left({tex_pi(c)} n{_phase_tex(-th)}\right)")
        p["reduce_why_html"] = (r"The principal value of the digital frequency is negative, and a cosine is even: "
                           r"$\cos(-\omega_0 n + \theta) = \cos(\omega_0 n - \theta)$. The phase changes sign along with the frequency.")
    ans = {"f0": float(f0), "amp": float(A), "phi": float(phi)}
    return p, ans


# ----------------------------------------------------------------------------- (b)
TS = [(F(1), r"\text{s}", "s"), (F(2), r"\text{s}", "s"), (F(1, 2), r"\text{s}", "s"), (F(1), r"\text{ms}", "ms"),
      (F(2), r"\text{ms}", "ms"), (F(4), r"\text{ms}", "ms"), (F(5), r"\text{ms}", "ms"), (F(10), r"\text{ms}", "ms")]


def _dec(q):
    """Terminating decimal of a Fraction (denominators 2^a 5^b only): 6.5, 1.625, 3."""
    q = F(q)
    if q.denominator == 1:
        return str(q.numerator)
    s = f"{float(q):.6f}".rstrip("0").rstrip(".")
    assert F(s) == q, (q, s)
    return s


def _x_terms_tex(xs):
    """2\\delta[n-3] - \\delta[n-5]."""
    terms = []
    for k in sorted(xs):
        d = r"\delta[n]" if k == 0 else rf"\delta[n-{k}]"
        terms.append((xs[k], d))
    return fmt.tex_sum(terms)


def _kT(k):
    """0, T, 2T, ..."""
    return "0" if k == 0 else ("T" if k == 1 else f"{k}T")


def _xv(xs, n):
    return xs.get(n, 0)


def _half_R(xs, m3):
    """y_a((m3 + 1/2)T) = (2/pi) * R with R = sum_n x[n] (-1)^(m3-n) / (2(m3-n)+1)."""
    return sum((F(v) * (-1) ** abs(m3 - n) / (2 * (m3 - n) + 1) for n, v in xs.items()), F(0))


def _part_b():
    while True:
        fam = random.choice(["delta", "pair", "seq", "seq"])
        if fam == "delta":
            xs = {random.randint(1, 5): random.choice([1, 2, 3, 4, -1, -2, -3])}
        elif fam == "pair":
            k1, d = random.randint(0, 3), random.choice([2, 3])
            xs = {k1: random.choice([1, 2, 3, -1, -2]), k1 + d: random.choice([1, 2, 3, -1, -2, -3])}
        else:
            n0 = random.randint(0, 2)
            vals = [random.choice([1, 2, 3, -1, -2, -3]), random.choice([-2, -1, 0, 1, 2, 3]),
                    random.choice([1, 2, 3, -1, -2, -3])]
            xs = {n0 + i: v for i, v in enumerate(vals) if v != 0}
        lo, hi = min(xs), max(xs)
        # ZOH instant: x[m1] != x[m1 + 1] so that "the wrong sample" is a real mistake
        c1 = [m for m in range(max(0, lo - 1), hi + 1) if _xv(xs, m) != _xv(xs, m + 1)]
        m1 = random.choice(c1)
        alpha = random.choice([F(1, 4), F(1, 2), F(3, 4)])
        # ideal D/A at a sample instant
        if fam == "pair" and random.random() < 0.65:
            m2 = random.randint(lo + 1, hi - 1)                       # between the two samples: 0 (SP2025 #1(c))
        elif fam == "delta":
            k = lo
            m2 = random.choice([k, k, k - 1, k + 1, k + 2])
        else:
            m2 = random.randint(lo, hi)
        if m2 < 0:
            continue
        # ideal D/A halfway between two sample instants
        m3 = random.choice([m for m in range(max(0, lo - 2), hi + 2)])
        R = _half_R(xs, m3)
        if R == 0 or abs(2 * float(R) / math.pi) < 0.05:
            continue
        lin = F(_xv(xs, m3) + _xv(xs, m3 + 1), 2)                  # linear interpolation (wrong)
        if abs(2 * float(R) / math.pi - float(lin)) < 0.03 * max(0.1, abs(2 * float(R) / math.pi)):
            continue
        break
    Tv, Tunit, Tword = random.choice(TS)
    t1, t2, t3 = (m1 + alpha) * Tv, m2 * Tv, (m3 + F(1, 2)) * Tv
    if fam == "seq":
        x_tex = fmt.tex_seq([_xv(xs, n) for n in range(lo, hi + 1)], lo)
        x_intro = r"the sequence $x[n] = " + x_tex + r"$ (arrow at $n = 0$; all samples not shown are zero)"
    else:
        x_tex = _x_terms_tex(xs)
        x_intro = r"$x[n] = " + x_tex + r"$"
    # worked halfway sum
    hterms = []
    for n in sorted(xs):
        j = m3 - n
        coef = F((-1) ** abs(j), 2 * j + 1)
        hterms.append((F(xs[n]) * coef, ""))
    half_terms_tex = r"\tfrac{2}{\pi}\left(" + " + ".join(
        rf"{fmt.tex_num(xs[n], paren=True)}\cdot{fmt.tex_num(F((-1) ** abs(m3 - n), 2 * (m3 - n) + 1), paren=True)}"
        for n in sorted(xs)) + r"\right)"
    val3 = 2 * float(R) / math.pi
    p = {
        "family": fam, "xs": {str(k): v for k, v in sorted(xs.items())},
        "x_intro_html": x_intro, "x_tex": x_tex,
        "Tval": plain(Tv), "Tunit": Tword,
        "Tb_tex": rf"T = {_dec(Tv)}\ {Tunit}",
        "t1_tex": rf"t_1 = {_dec(t1)}\ {Tunit}", "t2_tex": rf"t_2 = {_dec(t2)}\ {Tunit}",
        "t3_tex": rf"t_3 = {_dec(t3)}\ {Tunit}",
        "m1": m1, "alpha": plain(alpha), "m2": m2, "m3": m3,
        "t1T_tex": fmt.tex_num(m1 + alpha), "t2T_tex": _kT(m2), "t3T_tex": fmt.tex_num(m3 + F(1, 2)),
        "ratio2_tex": (r"\dfrac{\sin(\pi n)}{\pi n}" if m2 == 0 else
                       rf"\dfrac{{\sin\big(\pi({m2} - n)\big)}}{{\pi({m2} - n)}}"),
        "sign3_tex": "(-1)^{n}" if m3 == 0 else rf"(-1)^{{{m3} - n}}",
        "den3_tex": rf"{2 * m3 + 1} - 2n",
        "m1p1": m1 + 1, "x_m1_tex": fmt.tex_num(_xv(xs, m1)), "x_m2_tex": fmt.tex_num(_xv(xs, m2)),
        "half_terms_tex": half_terms_tex, "R_tex": fmt.tex_num(R, paren=True), "R": plain(R),
        "m1T_tex": _kT(m1), "m1p1T_tex": _kT(m1 + 1),
        "val3_tex": f"{val3:.4f}",
        "half_pair_tex": rf"t_3 = {fmt.tex_num(m3 + F(1, 2))}\,T",
    }
    ans = {"zoh": float(_xv(xs, m1)), "ideal_s": float(_xv(xs, m2)), "ideal_h": val3}
    return p, ans


# ----------------------------------------------------------------------------- (c)
FSC = [200, 300, 400, 500, 600, 800, 1000, 1200]
PSI = [F(0), F(0), F(0), F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(1, 2), F(-1, 2)]


def _out_of(F_hz, phi, fs):
    """(output frequency in Hz, output phase) of the ideal D/A for input cos(2 pi F t + phi) sampled at fs,
    or None when the digital frequency lands on 0 or pi (degenerate)."""
    r = F(F_hz, fs) % 1
    if r == 0 or r == F(1, 2):
        return None
    if r < F(1, 2):
        return r * fs, principal(phi)
    return (1 - r) * fs, principal(-phi)


def _amp(B):
    return "" if B == 1 else str(B)


def _cos_tex(B, F_hz, phi):
    return rf"{_amp(B)}\cos\left({rad_tex(F_hz)} t{_phase_tex(phi)}\right)"


def _part_c():
    while True:
        fs = random.choice(FSC)
        D = random.choice([8, 10, 12])
        if fs % D:
            continue
        g = fs // D
        f0 = g * random.randint(1, D // 2 - 1)
        B = random.choice([1, 2, 3, 4])
        psi = random.choice(PSI)
        target = (F(f0), principal(psi))
        good = [(f0, psi), (fs - f0, -psi), (fs + f0, psi), (2 * fs - f0, -psi), (2 * fs + f0, psi)]
        bad = [(fs - f0, psi), (fs + f0, -psi), (f0 + F(fs, 2), psi), (F(fs, 2) - f0, psi), (2 * f0, psi),
               (f0, -psi), (3 * f0, psi), (fs - 2 * f0, psi), (F(fs, 2) + 2 * f0, psi)]
        good = [(F(a), principal(b)) for a, b in good if a > 0]
        bad = [(F(a), principal(b)) for a, b in bad if a > 0]
        good = [o for o in good if _out_of(o[0], o[1], fs) == target]
        bad = [o for o in bad if _out_of(o[0], o[1], fs) not in (target, None) and o not in good
               and (o[0] * 2).denominator == 1]
        bad = list(dict.fromkeys(bad))
        n_good = random.choice([2, 2, 3])
        if len(good) < n_good or len(bad) < 5 - n_good:
            continue
        chosen_good = random.sample(good, n_good)
        trap = (F(fs - f0), principal(psi))                     # folded alias with the phase NOT flipped
        pool = list(bad)
        chosen_bad = []
        if psi != 0 and trap in pool and random.random() < 0.8:
            chosen_bad.append(trap)
            pool.remove(trap)
        chosen_bad += random.sample(pool, 5 - n_good - len(chosen_bad))
        break
    opts = [(o, True) for o in chosen_good] + [(o, False) for o in chosen_bad]
    random.shuffle(opts)
    choices = [{"text": f"$x_a(t) = {_cos_tex(B, o[0], o[1])}$", "correct": ok} for o, ok in opts]
    rows = []
    for (Fh, ph), ok in sorted(opts, key=lambda t: t[0][0]):
        r = F(Fh, fs)
        w = 2 * r                                     # digital frequency / pi
        wp = principal(w)
        out = _out_of(Fh, ph, fs)
        samp = rf"{_amp(B)}\cos\left({tex_pi(wp)} n{_phase_tex(ph)}\right)" if wp > 0 else \
            rf"{_amp(B)}\cos\left({tex_pi(-wp)} n{_phase_tex(-ph)}\right)"
        rows.append({
            "in_tex": _cos_tex(B, Fh, ph), "F_tex": fmt.tex_num(Fh), "w_tex": tex_pi(w), "wp_tex": tex_pi(wp),
            "samp_tex": samp, "out_tex": _cos_tex(B, out[0], out[1]), "verdict": "yes" if ok else "no",
        })
    p = {
        "fsc": fs, "f0c": plain(f0), "B": B, "psi": plain(psi),
        "yc_tex": _cos_tex(B, f0, psi),
        "w0c_tex": tex_pi(F(2 * f0, fs)),
        "f0c_tex": str(f0), "fsc_tex": str(fs),
        "options": [{"F": plain(o[0]), "phi": plain(o[1]), "correct": ok} for o, ok in opts],
        "inputs_choices": choices,
        "crows": rows,
    }
    return p


def generate(data):
    pa, ans_a = _part_a()
    pb, ans_b = _part_b()
    pc = _part_c()
    p = data["params"]
    p.update({"a": pa, "b": pb, "c": pc})
    p["inputs_choices"] = pc.pop("inputs_choices")
    c = data["correct_answers"]
    c.update(ans_a)
    c.update(ans_b)
