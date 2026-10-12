"""The DFT samples the DTFT: X[k] = X_d(2 pi k / N) for a signal of length at most N, and zero
padding samples the same DTFT more densely (Midterm 2 family "DFT computation and properties").

Models: FA2024 MT2 #7 (smallest N whose DFT contains X_d(2pi/5) and X_d(pi/3): N = lcm(5, 6) = 30),
FA2023 #1(h) (X[11] of a length-14 DFT is X_d(-3pi/7)), SP2025 #7(c) and SP2023 #9(c) (zero padding:
X[2] = Z[15], X[6] = Z^*[15]), FA2024 #1(d) (X[k] = Y[2k] after padding to 2N), the zero-padding
True/False statements FA2019 #1(a), SP2023 #1(c), FA2023 #1(f), SP2025 #1(d), and FA2023 #6(c)
(padding adds bins, not resolution). 2019 course notes, Sec. 2.5.1.1 (sampling property) and 2.6.2
(zero padding).

Parts
  (a) which bin k of an N-point DFT holds X_d(w0) (w0 = 2 pi k0/N, sometimes written as the negative
      frequency w0 - 2 pi);
  (b) the smallest N whose DFT grid contains two given frequencies: lcm of the denominators of w/(2 pi);
  (c) x of length L in {6, 8} zero-padded to M = (L/2) q (a multiple of L or not): which of six relations
      X[a] = Y[b], X[a] = Y^*[b] are true (checkbox; direct, mirrored and conjugate partners, and the traps
      "no index scaling", "ratio of the appended zeros", "both last bins", a bin with no partner, a
      superfluous conjugate);
  (d) one multiple-choice statement on what zero padding does (and does not do)."""

import math
import random
from fractions import Fraction as Fr

import numpy as np
from ece310 import fmt, mc, w8_dft as W

TRUE_STMTS = [
    r"The DTFT is unchanged, $Y_d(\omega) = X_d(\omega)$; the ${M}$-point DFT samples that same DTFT on a grid ${r}$ times finer.",
    r"It gives more closely spaced samples of the same DTFT (useful for plotting, or to reach a length the FFT likes), but no better resolution.",
    r"$X_d(0) = \sum_n x[n]$ and the energy $\sum_n \lvert x[n]\rvert^2$ are unchanged by the zero padding.",
    r"The appended zeros only add terms $0\cdot e^{-j\omega n}$ to the DTFT sum, so $Y_d(\omega) = X_d(\omega)$ at every $\omega$.",
    r"The width of a sinusoid's main lobe in $\lvert X_d(\omega)\rvert$ is set by the ${L}$ original samples; the padding leaves it unchanged.",
]
FALSE_STMTS = [
    r"It improves the frequency resolution: two sinusoids too close to be told apart from the ${L}$ samples become separable.",
    r"It changes the DTFT: $Y_d(\omega) = X_d({r}\omega)$.",
    r"It removes the leakage of a sinusoid whose frequency falls between the DFT frequencies $\frac{2\pi k}{L}$.",
    r"It makes the main lobe of each sinusoid's peak in $\lvert X_d(\omega)\rvert$ narrower.",
    r"It changes the DTFT by adding new frequency content between the old DFT frequencies.",
]


def _omega_tex(t):
    """w = t*pi in LaTeX."""
    return W.tex_pi_frac(t)


def _part_a(p):
    N = random.choice([8, 10, 12, 16, 18, 20, 24])
    while True:
        k = random.randint(1, N - 1)
        if 2 * k != N and math.gcd(k, N) < N // 2:
            break
    t = Fr(2 * k, N)                                   # w0 / pi in (0, 2)
    neg = t > 1 and random.random() < 0.6
    shown = t - 2 if neg else t
    p["Na"], p["ka"], p["Nam1"] = N, k, N - 1
    p["wa"] = [shown.numerator, shown.denominator]
    p["wa_tex"] = _omega_tex(shown)
    p["wa_neg"] = neg
    p["wa_2pi_tex"] = _omega_tex(t)
    p["ka_frac_tex"] = (rf"\frac{{N\omega_0}}{{2\pi}} = \frac{{{N}}}{{2\pi}}\cdot {_omega_tex(t)} = {k}")
    return k


def _part_b(p):
    fracs = [Fr(a, q) for q in (3, 4, 5, 6, 8, 9, 10, 12) for a in range(1, q) if math.gcd(a, q) == 1 and 2 * a < q]
    while True:
        f1, f2 = sorted(random.sample(fracs, 2))
        q1, q2 = f1.denominator, f2.denominator
        if q1 == q2:
            continue
        N = q1 * q2 // math.gcd(q1, q2)
        if N > 60:
            continue
        t1, t2 = 2 * f1, 2 * f2                         # w / pi, in (0, 1)
        d1, d2 = t1.denominator, t2.denominator
        trap_pi = d1 * d2 // math.gcd(d1, d2)
        traps = {trap_pi, q1 * q2, max(q1, q2)}
        traps.discard(N)
        if len(traps) >= 2 and trap_pi != N:
            break
    if random.random() < 0.5:
        f1, f2, t1, t2, q1, q2 = f2, f1, t2, t1, q2, q1
    p["wb"] = [[t1.numerator, t1.denominator], [t2.numerator, t2.denominator]]
    p["wb1_tex"], p["wb2_tex"] = _omega_tex(t1), _omega_tex(t2)
    p["qb"] = [q1, q2]
    p["fb1_tex"], p["fb2_tex"] = fmt.tex_num(f1), fmt.tex_num(f2)
    p["q1"], p["q2"] = q1, q2
    p["kb1"], p["kb2"] = int(f1 * N), int(f2 * N)
    i = 0 if t1.denominator != q1 else 1
    wt, ft, qt, dt = (p["wb1_tex"], f1, q1, t1.denominator) if i == 0 else (p["wb2_tex"], f2, q2, t2.denominator)
    p["trap_note_html"] = (rf"Work with $\frac{{\omega}}{{2\pi}}$, not $\frac{{\omega}}{{\pi}}$: ${wt} = 2\pi\cdot{fmt.tex_num(ft)}$ needs a "
                          rf"multiple of ${qt}$, not just of ${dt}$ (that slip gives $N = {trap_pi}$).")
    p["Nb"] = N
    return N


def _wfrac(num, den):
    """2*pi*num/den as a multiple of pi, reduced to [0, 2)."""
    return Fr(2 * num, den) % 2


def _rel_tex(a, b, conj):
    return rf"X[{a}] = Y{'^*' if conj else ''}[{b}]"


def _part_c(p):
    """SP2025 #7(c) / SP2023 #9(c): x of length L in {6, 8} zero-padded to M = (L/2) q; a checkbox of
    relations X[a] = Y[b] and X[a] = Y^*[b]. M is not always a multiple of L (q odd: X[a] with odd a has
    no partner in Y), and the conjugate relations use X_d(-w) = X_d^*(w) for the real x."""
    while True:
        L = random.choice([6, 8])
        x = [random.randint(-3, 4) for _ in range(L)]
        if x[0] == 0 or x[-1] == 0 or sum(1 for v in x if v == 0) > 2:
            continue
        h = L // 2
        if random.random() < 0.5:
            q = random.choice([3, 5, 7, 9] if L == 6 else [3, 5, 7, 9, 11, 13, 15])
        else:
            q = random.choice([4, 6, 8, 10, 12, 14, 16, 18, 20] if L == 6 else [4, 6, 8, 10, 12, 14])
        M = h * q
        X, Y = np.fft.fft(x), np.fft.fft(x, M)
        rho = Fr(M, L)

        def partner(k):
            b = k * rho
            return int(b) % M if b.denominator == 1 else None

        paired = [k for k in range(1, L) if k != h and partner(k) is not None]
        unpaired = [k for k in range(1, L) if partner(k) is None]
        TP, TC, FP, FC = [], [], [], []                 # true/false x plain/conjugate: (a, b, conj, kind)
        for k in paired:
            bk = partner(k)
            TP.append((k, bk, False, "direct"))
            TC += [(k, (M - bk) % M, True, "mirror"), (L - k, bk, True, "mirrorX")]
            FP += [(k, k, False, "no_scaling"), (k, (M - bk) % M, False, "plain_mirror")]      # Y[M - k rho] = X^*[k]
            FC += [(k, bk, True, "conj_wrong"),                                                  # superfluous conjugate
                   (L - k, (M - bk) % M, True, "conj_mirrorX_wrong"),                             # = X[k], not X[L - k]
                   (k, (M - k) % M, True, "conj_no_scaling")]
            z = k * (rho - 1)                           # ratio of appended zeros instead of the total length
            if z.denominator == 1 and int(z) % M != bk:
                FP.append((k, int(z) % M, False, "zeros_ratio"))
                FC.append((k, (M - int(z)) % M, True, "conj_zeros_ratio"))                       # SP2023 #9(c) v
        if partner(h) is not None:
            TC.append((h, partner(h), True, "half_conj"))
        for k in unpaired:
            bk = k * rho
            FP += [(k, math.floor(bk), False, "no_partner"), (k, math.ceil(bk) % M, False, "no_partner")]
            FC.append((k, (M - math.floor(bk)) % M, True, "no_partner_conj"))
        FP.append((L - 1, M - 1, False, "last_bins"))

        def gap(r):
            a_, b_, conj, _ = r
            v = np.conj(Y[b_]) if conj else Y[b_]
            return abs(X[a_] - v)

        TP = [r for r in TP if gap(r) < 1e-6]
        TC = [r for r in TC if gap(r) < 1e-6]
        FP = [r for r in FP if gap(r) > 0.5]
        FC = [r for r in FC if gap(r) > 0.5]
        nt = random.choice([1, 2, 3, 4])               # true relations
        nc = random.choice([2, 3, 4])                  # conjugate relations, about half of them true
        ntc = min(nt, len(TC), random.choice([nc // 2, (nc + 1) // 2]))
        nfc, ntp = nc - ntc, nt - ntc
        nfp = 6 - nt - nfc
        if min(nfc, ntp, nfp) < 0 or ntp > len(TP) or nfc > len(FC) or nfp > len(FP):
            continue
        pick_t = random.sample(TC, ntc) + random.sample(TP, ntp)
        pick_f = random.sample(FC, nfc) + random.sample(FP, nfp)
        rels = pick_t + pick_f
        if len({(r[0], r[1], r[2]) for r in rels}) != 6:
            continue
        if unpaired and not any(r[3].startswith("no_partner") for r in rels):
            continue                                    # M not a multiple of L: show a bin without a partner

        def shortcut(r):                                # "plain iff b = a M/L, conjugate iff b != a M/L"
            on_grid = (r[0] * rho).denominator == 1 and int(r[0] * rho) % M == r[1]
            return (not on_grid) if r[2] else on_grid

        if all(shortcut(r) == (r in pick_t) for r in rels):
            continue                                    # the index-only shortcut must fail
        break
    random.shuffle(rels)
    p["xc"], p["Lc"], p["Lcm1"], p["M"], p["Mm1"], p["Mz"] = x, L, L - 1, M, M - 1, M - L
    p["xc_tex"] = W.tex_list([str(v) for v in x])
    p["ratio_tex"] = fmt.tex_num(rho)
    p["ratio_int"] = rho.denominator == 1
    p["rels"] = [{"a": r[0], "b": r[1], "conj": r[2], "kind": r[3]} for r in rels]
    p["rel_choices"] = [{"text": "$" + _rel_tex(*r[:3]) + "$", "correct": r in pick_t} for r in rels]
    lines = []
    for a, b, conj, kind in rels:
        wa = _wfrac(a, L)
        wb = _wfrac(b, M)
        if conj:
            wneg = (-wb) % 2
            why = (rf"$Y^*[{b}] = X_d^*\!\left({W.tex_pi_frac(wb)}\right) = X_d\!\left({W.tex_pi_frac(wneg)}\right)$ ($x$ real, $2\pi$-periodic) "
                   rf"and $X[{a}] = X_d\!\left({W.tex_pi_frac(wa)}\right)$")
            ok = wneg == wa
        else:
            why = rf"$X[{a}] = X_d\!\left({W.tex_pi_frac(wa)}\right)$ and $Y[{b}] = X_d\!\left({W.tex_pi_frac(wb)}\right)$"
            ok = wb == wa
        if ok:
            why += ": the same frequency."
        elif kind.startswith("no_partner"):
            why += rf": different frequencies; no bin of $Y$ sits at ${W.tex_pi_frac(wa)}$, since ${a}\cdot\frac{{{M}}}{{{L}}}$ is not an integer."
        else:
            why += ": different frequencies."
        lines.append({"rel_tex": _rel_tex(a, b, conj), "verdict": "true" if ok else "false", "why_html": why})
    p["rel_lines"] = lines
    return L


def generate(data):
    p, c = data["params"], data["correct_answers"]
    c["ka"] = _part_a(p)
    c["Nb"] = _part_b(p)
    L = _part_c(p)
    M = p["M"]

    def fill(t):
        return t.replace("{M}", str(M)).replace("{r}", p["ratio_tex"]).replace("{L}", str(L))

    p["zp_choices"] = mc.choices(fill(random.choice(TRUE_STMTS)), [fill(t) for t in random.sample(FALSE_STMTS, 3)], n=4)
