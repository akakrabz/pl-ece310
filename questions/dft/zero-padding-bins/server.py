"""Which bins of a DFT hold which DTFT values, and how many zeros to pad so that the DFT evaluates
X_d at two given frequencies (Midterm 2 family "DFT values and zero padding").

Same shape as FA2025 MT2 #6 (length-32 real x: (a) five True/False relations X[0] = X_d(2pi),
X[2] = X_d(pi/8), X[4] = X_d(pi/8), X[6] = X_d^*(-3pi/8), X[8] = X_d^*(-pi/2); (b) which of
L = 4, 16, 22, 40, 48 zeros make Y[k] evaluate X_d(pi/3) and X_d(4pi/9): 32 + L divisible by 18,
so 4, 22, 40). Part (c) is FA2024 MT2 #7 (smallest N capturing X_d(2pi/5) and X_d(pi/3): lcm(5, 6) = 30)
asked for the padding instead of the length; part (d) adds the Hermitian partner of a bin, as in
FA2025 #6(a) iv-v and FA2023 #1(h) (X[11] = X_d(-3pi/7) for N = 14). The "hz" variant type states
the two frequencies in hertz for samples taken at f_s, in the style of FA2019 MT2 #2 (F_s = 40 kHz,
N = 1000: bins 300 and 800 are 12 kHz and -8 kHz) and FA2024 #6.

Parts
  (a) checkbox, five relations X[k] = X_d(theta) / X[k] = X_d^*(theta) for the length-N real x
      (true kinds: X[0] = X_d(2pi), the plain bin, the bin written as a negative frequency, a
      conjugate at -theta, the mirror bin N - k; false kinds: the pi k/N grid, a conjugate without
      the sign flip, the negative frequency without a conjugate, the 4 pi k/N grid);
  (b) checkbox, five candidate L: N + L must be a multiple of q = lcm(q1, q2), where w_i/(2pi) = a_i/q_i.
      Distractors: a multiple of the lcm of the denominators of w_i/pi only (the "pi" slip), a length
      that fits one frequency only, an L that is itself a multiple of q;
  (c) the smallest L >= 0 (integer);
  (d) for that L (M = N + L): k with Y[k] = X_d(w_b) and k' in 0..M-1 with Y[k'] = X_d^*(w_b): k' = M - k.
N in {16, 20, 24, 32}; w_i/(2pi) = a/q with q in {4, 5, 6, 8, 9, 10, 12, 16, 18, 20, 24}, q1 and q2 not
dividing each other, lcm at most 40, N not a multiple of it (so L_min > 0)."""

import math
import random
from fractions import Fraction as Fr

from ece310 import w8_dft as W

QS = (4, 5, 6, 8, 9, 10, 12, 16, 18, 20, 24)
POOL = [Fr(a, q) for q in QS for a in range(1, q) if math.gcd(a, q) == 1 and 2 * a < q]
LMAX = 100


def _lcm(a, b):
    return a * b // math.gcd(a, b)


def _pi(t):
    return W.tex_pi_frac(t)


def _num(v):
    """Exact number in LaTeX (integers plain, fractions \\tfrac)."""
    v = Fr(v)
    if v.denominator == 1:
        return str(v.numerator)
    return rf"\tfrac{{{v.numerator}}}{{{v.denominator}}}"


# ----------------------------------------------------------------------------- (a) relations
def _rels(N):
    """Five relations (k, conj, theta, truth, kind); theta in units of pi."""
    h = N // 2
    inner = [k for k in range(1, N) if k != h]
    while True:
        nt = random.choice([2, 3, 3, 4])
        tkinds = random.sample(["plain", "wrap", "conjneg", "mirror"], nt - 1)
        tkinds.append("zero2pi" if random.random() < 0.4 else random.choice(["plain", "conjneg", "wrap", "mirror"]))
        fkinds = random.sample(["grid", "conjsame", "neg", "double"], 5 - nt)
        out = []
        for kind in tkinds + fkinds:
            if kind == "zero2pi":
                out.append((0, False, Fr(2), True, kind))
                continue
            k = random.choice(inner)
            if kind == "plain":
                out.append((k, False, Fr(2 * k, N), True, kind))
            elif kind == "wrap":
                k = random.choice([v for v in inner if v > h])
                out.append((k, False, Fr(2 * k, N) - 2, True, kind))
            elif kind == "conjneg":
                out.append((k, True, Fr(-2 * k, N), True, kind))
            elif kind == "mirror":
                k = random.choice([v for v in inner if v < h])
                out.append((N - k, True, Fr(2 * k, N), True, kind))
            elif kind == "grid":
                k = random.choice([v for v in inner if v % 2 == 0] or inner)
                out.append((k, False, Fr(k, N), False, kind))
            elif kind == "conjsame":
                k = random.choice([v for v in inner if v < h])
                out.append((k, True, Fr(2 * k, N), False, kind))
            elif kind == "neg":
                k = random.choice([v for v in inner if v < h])
                out.append((k, False, Fr(-2 * k, N), False, kind))
            else:                                               # double: the 4 pi k / N grid
                k = random.choice([v for v in inner if 4 * v < N])
                out.append((k, False, Fr(4 * k, N), False, kind))
        if len({r[0] for r in out}) != 5:
            continue
        if sum(1 for r in out if r[1]) < 1 or sum(1 for r in out if not r[1]) < 2:
            continue
        return sorted(out, key=lambda r: r[0])


def _rel_tex(k, conj, theta):
    return rf"X[{k}] = X_d{'^*' if conj else ''}\left({_pi(theta)}\right)"


def _rel_why(N, k, conj, theta, kind):
    wk = Fr(2 * k, N)
    base = rf"$\omega_{{{k}}} = \frac{{2\pi\cdot {k}}}{{{N}}} = {_pi(wk)}$"
    if kind == "zero2pi":
        return r"$X[0] = X_d(0)$, and $X_d$ is $2\pi$-periodic: $X_d(0) = X_d(2\pi)$."
    if kind == "plain":
        return base + "."
    if kind == "wrap":
        return base + rf", and ${_pi(wk)} - 2\pi = {_pi(theta)}$: $X_d$ is $2\pi$-periodic."
    if kind == "conjneg":
        return base + rf"; $x$ is real, so $X_d^*\left({_pi(theta)}\right) = X_d\left({_pi(-theta)}\right)$, the same frequency."
    if kind == "mirror":
        return (rf"$\omega_{{{k}}} = \frac{{2\pi\cdot {k}}}{{{N}}} = {_pi(wk)} \equiv {_pi(wk - 2)}$ ($2\pi$-periodicity); "
                rf"$x$ is real, so $X_d^*\left({_pi(theta)}\right) = X_d\left({_pi(-theta)}\right)$, the same frequency.")
    if kind == "grid":
        return base + rf", not ${_pi(theta)}$: the bins are $\frac{{2\pi}}{{{N}}}$ apart, so ${_pi(theta)}$ is bin ${Fr(k, 2)}$."
    if kind == "conjsame":
        return (base + rf", but $X_d^*\left({_pi(theta)}\right) = X_d\left({_pi(-theta)}\right) = X^*[{k}]$: a conjugate needs the "
                "frequency's sign flipped too.")
    if kind == "neg":
        return (base + rf", but $X_d\left({_pi(theta)}\right) = X_d^*\left({_pi(-theta)}\right) = X^*[{k}] \ne X[{k}]$ in general "
                "(a negative frequency needs a conjugate).")
    return base + rf", not ${_pi(theta)}$ (that is bin ${2 * k}$)."


# ----------------------------------------------------------------------------- (b), (c), (d)
def _pick_freqs(N):
    while True:
        f1, f2 = random.sample(POOL, 2)
        q1, q2 = f1.denominator, f2.denominator
        if q1 % q2 == 0 or q2 % q1 == 0:
            continue
        q = _lcm(q1, q2)
        if q > 40 or N % q == 0:
            continue
        d = _lcm((2 * f1).denominator, (2 * f2).denominator)
        if d == q and random.random() < 0.85:
            continue                                        # prefer pairs where the w/pi slip gives a wrong length
        return f1, f2, q1, q2, q, d


def _candidates(N, q1, q2, q, d, Lmin):
    def works(L):
        return (N + L) % q == 0

    good = list(range(Lmin, LMAX + 1, q))
    pi_trap = [L for L in range(1, LMAX + 1) if (N + L) % d == 0 and not works(L)]
    one_only = [L for L in range(1, LMAX + 1) if ((N + L) % q1 == 0) != ((N + L) % q2 == 0)]
    l_trap = [L for L in range(q, LMAX + 1, q) if not works(L)]
    while True:
        nt = min(random.choice([1, 2, 2, 3]), len(good))
        if random.random() < 0.5 or len(good) - 1 < nt:
            tr = [Lmin] + random.sample(good[1:], nt - 1)
        else:
            tr = random.sample(good[1:], nt)
        fl = []
        for src in random.sample([pi_trap, one_only, l_trap], 3):
            opts = [L for L in src if L not in fl]
            if opts and len(fl) < 5 - nt:
                fl.append(random.choice(opts))
        while len(fl) < 5 - nt:
            fl.append(random.choice([L for L in one_only + pi_trap if L not in fl] or
                                    [L for L in range(1, LMAX + 1) if not works(L) and L not in fl]))
        cands = sorted(set(tr + fl))
        if len(cands) == 5:
            return cands


def _cand_row(N, L, fs):
    M = N + L
    vals = [M * f for f in fs]
    ok = all(v.denominator == 1 for v in vals)
    cells = "".join(f"<td>${_num(v)}$</td>" for v in vals)
    return (f"<tr><td>${L}$</td><td>${M}$</td>{cells}<td><b>{'works' if ok else 'fails'}</b></td></tr>")


def generate(data):
    p, c = data["params"], data["correct_answers"]
    N = random.choice([16, 20, 24, 32])
    hz = random.random() < 0.35
    f1, f2, q1, q2, q, d = _pick_freqs(N)
    Lmin = (-N) % q
    M = N + Lmin
    p["N"], p["Nm1"], p["hz"], p["dig"] = N, N - 1, hz, not hz

    # (a)
    rels = _rels(N)
    p["rels"] = [{"k": r[0], "conj": r[1], "theta": [r[2].numerator, r[2].denominator], "kind": r[4]} for r in rels]
    p["rel_choices"] = [{"text": "$" + _rel_tex(*r[:3]) + "$", "correct": r[3]} for r in rels]
    p["rel_lines"] = [{"rel_tex": _rel_tex(*r[:3]), "verdict": "True" if r[3] else "False",
                       "why_html": _rel_why(N, r[0], r[1], r[2], r[4])} for r in rels]
    p["Nh"], p["dw_tex"] = N // 2, _pi(Fr(2, N))

    # (b)
    p["f"] = [[f1.numerator, f1.denominator], [f2.numerator, f2.denominator]]
    p["q1"], p["q2"], p["q"], p["d"] = q1, q2, q, d
    p["w1_tex"], p["w2_tex"] = _pi(2 * f1), _pi(2 * f2)
    p["f1_tex"], p["f2_tex"] = _num(f1), _num(f2)
    if hz:
        m = random.choice([v for v in (10, 20, 25, 50, 100, 200) if 200 <= q * v <= 8000])
        fs = q * m
        p["fs"] = fs
        p["F1"], p["F2"] = int(f1 * fs), int(f2 * fs)
        p["spacing_tex"] = rf"\frac{{f_s}}{{M}} = \frac{{{fs}}}{{M}}"
    cands = _candidates(N, q1, q2, q, d, Lmin)
    p["L_choices"] = [{"text": f"$L = {L}$", "correct": (N + L) % q == 0} for L in cands]
    p["cands"] = cands
    hdr = ("<tr><th>$L$</th><th>$M = N + L$</th>"
           + (r"<th>$\frac{M f_1}{f_s}$</th><th>$\frac{M f_2}{f_s}$</th>" if hz else
              r"<th>$M\cdot\frac{\omega_1}{2\pi}$</th><th>$M\cdot\frac{\omega_2}{2\pi}$</th>")
           + "<th></th></tr>")
    p["cand_table_html"] = ('<table class="table table-sm table-bordered text-center" style="width:auto">' + hdr
                            + "".join(_cand_row(N, L, (f1, f2)) for L in cands) + "</table>")
    pi_bad = [L for L in range(1, LMAX + 1) if (N + L) % d == 0 and (N + L) % q]
    p["pi_note"] = d != q
    p["pi_trap_L"] = pi_bad[0] if pi_bad else 0
    p["pi_trap_M"] = N + p["pi_trap_L"]
    c["Lmin"] = Lmin
    p["Lmin"] = Lmin

    # (c), (d)
    p["Lmin_M"] = M
    b = random.choice([0, 1])
    fb = (f1, f2)[b]
    k = int(M * fb)
    p["b"] = b + 1
    p["wb_tex"] = _pi(2 * fb)
    p["wbneg_tex"] = _pi(-2 * fb)
    p["wbc_tex"] = _pi(2 - 2 * fb)
    p["fb_tex"] = _num(fb)
    if hz:
        p["Fb"] = int(fb * p["fs"])
    p["kb"], p["kbc"], p["Mm1"] = k, M - k, M - 1
    c["k"], c["kc"] = k, M - k
