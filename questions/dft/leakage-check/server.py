"""DFT of a truncated cosine: exactly two nonzero bins, or leakage (Lecture 22 Sec. 1.2, eqs. 11-14).

X[k] = (A/2) C(2pi k/N - w0) + (A/2) C(2pi k/N + w0). With k0 = w0 N / (2 pi) (the number of periods in
the N samples) an integer, the DFT evaluates the two lobes at their peaks and at their zeros: X[k0] =
X[N - k0] = A N / 2 and every other bin is 0. Otherwise every bin is (in general) nonzero, the tallest at
the bin nearest k0.

(a) checkbox, four statements in the wording of the past True/False items FA2019 MT2 #1(f) ("The 64-point
    DFT of x[n] = cos(pi n / 16), 0 <= n <= 63 has only two nonzero elements": true), FA2023 MT2 #1(g)
    (cos(pi n / 8) on 24 samples: 1.5 periods, false), SP2025 MT2 #1(e) (cos(50 pi t) sampled at 50 Hz
    for 1 s: cos(pi n), a single bin) and Lecture 22 slide Example 2 (sin(pi n / 11), N = 33: k0 = 1.5;
    sin(3 pi n / 5), N = 20: k0 = 6; cos(pi n / 4) + cos(pi n / 3), N = 48: four bins). Statement kinds:
    on-bin sinusoid / sum of two on-bin sinusoids / on-bin sampled analog cosine (true); half-period
    record, a record one or two samples off a whole number of periods, a sum with one off-bin term,
    the f0 = f_s / 2 cosine with one bin, a sampled analog cosine with a fractional number of
    periods (false). 1-3 true statements.
(b), (c) two variant types:
  l22   Lecture 22 eqs. 13-14 (cos(pi n / 3) kept for 24 and for 22 samples) with the (b) wording of
        SP2025 MT2 #9 ("Determine the DFT of the truncated sequence {x[n]}_{n=0}^{15}" of cos(pi n / 2):
        X[4] = X[12] = 8): (b) an on-bin record of A cos(w0 n): the two nonzero bins and |X[k]| = AN/2;
        (c) the same cosine on a record a few samples shorter or longer: the number of periods in the
        record (a fraction such as 11/3) and the tallest bin in 0 <= k <= N/2.
  fa21  FA2021 MT2 #7 (cos(10 pi t) sampled at 100 Hz for 5 s, L = 500: (a) peaks at k = 25 and 475,
        |X| = 250; (b) zero-padded to N = 1024: k = 51.2, peaks at 51 and 973): (b) the L-point DFT of
        A cos(Omega0 t) sampled at f_s for D seconds (a whole number of periods): the two nonzero bins
        and their magnitude; (c) after zero padding to a power of two: the fractional index
        kappa = N w0 / (2 pi) and the tallest bin in 0 <= k <= N/2.
Off-bin variants keep the fractional part of k0 away from 1/2 (and are checked: the tallest bin is the
nearest one, at least 5 % taller than every other bin in 0..N/2)."""

import cmath
import math
import random
from fractions import Fraction as Fr

from ece310 import fmt, w1_signals as S, w8_dft as W

FRACS = [Fr(p, q) for q in (2, 3, 4, 5, 6, 8, 10, 11, 12, 16) for p in range(1, q) if math.gcd(p, q) == 1]


def period(t):
    """Fundamental period of cos(t pi n), t = p/q in lowest terms."""
    return 2 * t.denominator // math.gcd(t.numerator, 2 * t.denominator)


def _C(w, L):
    s = math.sin(w / 2)
    ratio = L * math.cos(L * w / 2) / math.cos(w / 2) if abs(s) < 1e-12 else math.sin(L * w / 2) / s
    return cmath.exp(-0.5j * w * (L - 1)) * ratio


def mag(A, w0, L, Npts, k):
    """|X[k]| of the Npts-point DFT of A cos(w0 n), 0 <= n <= L-1 (zero-padded when Npts > L)."""
    w = 2 * math.pi * k / Npts
    return abs(A / 2 * (_C(w - w0, L) + _C(w + w0, L)))


def tallest_ok(A, w0, L, Npts, kap):
    """Tallest bin in 0..Npts/2 is round(kap), at least 5 % above every other bin there."""
    vals = [mag(A, w0, L, Npts, k) for k in range(Npts // 2 + 1)]
    order = sorted(range(len(vals)), key=lambda k: -vals[k])
    return order[0] == round(kap) and vals[order[0]] >= 1.05 * vals[order[1]]


# ----------------------------------------------------------------------------- (a) statements
def _amp_tex(a):
    return "" if a == 1 else str(a)


def _sin_tex(a, fn, t):
    return f"{_amp_tex(a)}\\{fn}\\!\\left({S.tex_pi_n(t)}\\right)"


def _disc_stmt(terms, N):
    sig = " + ".join(_sin_tex(*tm) for tm in terms)
    cnt = "two" if len(terms) == 1 else "four"
    if random.random() < 0.5:
        text = rf"The ${N}$-point DFT of $x[n] = {sig}$, $0 \le n \le {N - 1}$, has exactly {cnt} nonzero values."
    else:
        text = rf"The DFT of $\{{x[n]\}}_{{n=0}}^{{ {N - 1} }}$ with $x[n] = {sig}$ has exactly {cnt} nonzero values."
    lines = []
    ok = True
    for a, fn, t in terms:
        k0 = t * N / 2
        ok &= k0.denominator == 1
        lines.append(rf"${S.tex_pi(t)}$: $k_0 = \frac{{\omega_0 N}}{{2\pi}} = {fmt.tex_num(k0)}$"
                     + (rf" (bins ${k0}$ and ${N - k0}$)" if k0.denominator == 1 else " (not an integer: leakage)"))
    return {"text": text, "N": N, "terms": [[a, fn, t.numerator, t.denominator] for a, fn, t in terms],
            "claim": 2 * len(terms), "analog": None, "why": "; ".join(lines), "truth": ok}


def _analog_stmt(f0, fs, D):
    N = int(fs * D)
    W0 = 2 * f0
    om = f"{fmt.tex_num(W0)}\\pi t" if W0 != 1 else r"\pi t"
    Dt = fmt.tex_num(D) if Fr(D).denominator == 1 else f"{float(D):g}"
    text = (rf"An analog signal $x_a(t) = \cos({om})$ is sampled at $f_s = {fs}$ Hz for ${Dt}$ s to obtain "
            rf"$\{{x[n]\}}_{{n=0}}^{{ {N - 1} }}$. Its DFT $X[k]$ is nonzero for exactly two values of $k$, $0 \le k \le {N - 1}$.")
    t = Fr(2 * f0, fs)                                         # w0 / pi
    k0 = Fr(f0) * D
    if t == 1:
        why = (rf"$\omega_0 = \frac{{\Omega_0}}{{f_s}} = \pi$: $x[n] = \cos(\pi n) = (-1)^n = e^{{j\frac{{2\pi\cdot {N // 2}}}{{ {N} }}n}}$, "
               rf"a single bin $k = {N // 2}$ (one value, not two)")
        ok = False
    else:
        ok = k0.denominator == 1
        why = (rf"$\omega_0 = \frac{{\Omega_0}}{{f_s}} = {S.tex_pi(t)}$, $k_0 = f_0 D = {fmt.tex_num(Fr(f0))}\cdot {Dt} = {fmt.tex_num(k0)}$"
               + (rf" (bins ${k0}$ and ${N - k0}$)" if ok else " periods (not an integer: leakage)"))
    return {"text": text, "N": N, "terms": [[1, "cos", t.numerator, t.denominator]], "claim": 2,
            "analog": [float(f0), fs, float(D)], "why": why, "truth": ok}


def _pick_t(maxP=32):
    while True:
        t = random.choice(FRACS)
        if period(t) <= maxP:
            return t


def _stmt(kind):
    for _ in range(500):
        a = random.choice([1, 1, 2, 3])
        fn = random.choice(["cos", "sin"])
        if kind in ("on1", "half1", "near1"):
            t = _pick_t()
            P = period(t)
            if kind == "on1":
                Ns = [c * P for c in range(1, 20) if 8 <= c * P <= 64]
            elif kind == "half1":
                if t.numerator % 2 == 0:
                    continue
                Ns = [(2 * c + 1) * t.denominator for c in range(1, 20) if 8 <= (2 * c + 1) * t.denominator <= 64]
            else:
                Ns = [c * P + d for c in range(1, 20) for d in (-2, -1, 1, 2) if 8 <= c * P + d <= 64 and c * P >= 8]
            Ns = [N for N in Ns if (t * N / 2).denominator != 1] if kind != "on1" else Ns
            if not Ns:
                continue
            return _disc_stmt([(a, fn, t)], random.choice(Ns))
        if kind in ("on2", "mixed2"):
            t1, t2 = _pick_t(16), _pick_t(16)
            if t1 == t2:
                continue
            P1, P2 = period(t1), period(t2)
            L = P1 * P2 // math.gcd(P1, P2)
            if kind == "on2":
                Ns = [c * L for c in range(1, 10) if 8 <= c * L <= 64]
            else:
                Ns = [c * P1 for c in range(1, 20) if 8 <= c * P1 <= 64 and (t2 * c * P1 / 2).denominator != 1]
            if not Ns:
                continue
            terms = [(a, fn, t1), (random.choice([1, 2]), random.choice(["cos", "sin"]), t2)]
            random.shuffle(terms)
            return _disc_stmt(terms, random.choice(Ns))
        fs = random.choice([20, 40, 50, 60, 80, 100])
        D = random.choice([Fr(1, 2), Fr(1), Fr(1), Fr(2)])
        N = fs * D
        if N.denominator != 1 or not 10 <= N <= 100:
            continue
        if kind == "pi1":
            return _analog_stmt(Fr(fs, 2), fs, D)
        f0 = Fr(random.randint(2, fs - 1), 2)                     # half-hertz grid below f_s / 2
        k0 = f0 * D
        if k0 < 1:
            continue
        if kind == "aon" and k0.denominator == 1:
            return _analog_stmt(f0, fs, D)
        if kind == "aoff" and k0.denominator != 1:
            return _analog_stmt(f0, fs, D)
    raise RuntimeError(kind)


def _part_a(p):
    while True:
        nt = random.choice([1, 2, 2, 3])
        tk = random.sample(["on1", "on2", "aon", "on1"], nt)
        fk = random.sample(["half1", "near1", "mixed2", "pi1", "aoff", "half1"], 4 - nt)
        st = [_stmt(k) for k in tk + fk]
        if len({s["text"] for s in st}) == 4 and all(s["truth"] == (i < nt) for i, s in enumerate(st)):
            break
    random.shuffle(st)
    p["stmts"] = st
    p["stmt_choices"] = [{"text": s["text"], "correct": s["truth"]} for s in st]
    p["stmt_lines"] = [{"text": s["text"], "verdict": "True" if s["truth"] else "False", "why": s["why"]} for s in st]


# ----------------------------------------------------------------------------- (b), (c)
def _l22(p, c):
    while True:
        t = Fr(2, random.choice([4, 5, 6, 8, 10, 12]))           # w0 = 2 pi / P: N / P periods
        P = period(t)
        A = random.randint(1, 5)
        N1 = random.choice([m * P for m in range(2, 9) if 12 <= m * P <= 64] or [0])
        if not N1:
            continue
        N2 = N1 + random.choice([-3, -2, -1, 1, 2, 3])
        k01, k02 = t * N1 / 2, t * N2 / 2
        fr = k02 - math.floor(k02)
        if N2 < 10 or k01.denominator != 1 or fr == 0 or abs(fr - Fr(1, 2)) < Fr(1, 8):
            continue
        if not 1 <= round(k02) < N2 / 2 or not tallest_ok(A, math.pi * t, N2, N2, float(k02)):
            continue
        break
    p["A"], p["N1"], p["N2"], p["N1m1"], p["N2m1"], p["P"] = A, N1, N2, N1 - 1, N2 - 1, P
    p["x_tex"] = (f"{A}" if A != 1 else "") + rf"\cos\!\left({S.tex_pi_n(t)}\right)"
    p["w0_tex"] = S.tex_pi(t)
    p["t"] = [t.numerator, t.denominator]
    p["k01"], p["N1mk01"], p["k02_tex"] = int(k01), int(N1 - k01), fmt.tex_num(k02)
    p["H_tex"] = fmt.tex_num(Fr(A * N1, 2))
    p["k02_floor"], p["k02_dec"] = math.floor(k02), f"{float(k02):.4g}"
    p["kmax"] = round(k02)
    p["N2half"] = N2 // 2
    pk = mag(A, math.pi * t, N2, N2, round(k02))
    p["peak2"] = f"{pk:.2f}"
    p["peak2_rel"] = "below" if pk < A * N2 / 2 - 0.05 else "close to"
    p["peak2_full_tex"] = fmt.tex_num(Fr(A * N2, 2))
    c["k1"], c["k2"], c["Xk"] = int(k01), int(N1 - k01), float(Fr(A * N1, 2))
    c["kap"], c["kmax"] = float(k02), round(k02)


def _fa21(p, c):
    while True:
        fs = random.choice([40, 50, 64, 80, 100, 128, 200, 250, 256])
        D = random.choice([Fr(1, 2), Fr(1), Fr(3, 2), Fr(2), Fr(5, 2), Fr(4), Fr(5)])
        L = fs * D
        if L.denominator != 1 or not 32 <= L <= 1000:
            continue
        L = int(L)
        f0 = Fr(random.randint(2, 2 * fs // 2 - 1), 2)        # half-hertz grid
        k0 = f0 * D
        if f0 > Fr(3 * fs, 10) or k0.denominator != 1 or not 2 <= k0 <= L // 2 - 2:
            continue
        Np = 1 << (L - 1).bit_length()
        if Np == L:
            Np *= 2
        Np = random.choice([Np, 2 * Np]) if Np <= 1024 else Np
        kap = f0 * Np / fs
        fr = kap - math.floor(kap)
        if kap > 80 or min(fr, 1 - fr) < Fr(3, 20) or abs(fr - Fr(1, 2)) < Fr(1, 8):
            continue
        A = random.randint(1, 4)
        w0 = 2 * math.pi * float(f0) / fs
        if not tallest_ok(A, w0, L, Np, float(kap)):
            continue
        break
    W0 = 2 * f0
    p["A"], p["fs"], p["L"], p["Lm1"], p["Np"], p["Np_half"] = A, fs, L, L - 1, Np, Np // 2
    p["D_tex"] = fmt.tex_num(D) if D.denominator == 1 else f"{float(D):g}"
    p["D_unit"] = "second" if D == 1 else "seconds"
    p["xc_tex"] = (f"{A}" if A != 1 else "") + rf"\cos({fmt.tex_num(W0) if W0 != 1 else ''}\pi t)"
    p["f0_tex"] = fmt.tex_num(f0)
    p["W0_tex"] = (fmt.tex_num(W0) if W0 != 1 else "") + r"\pi"
    t = Fr(2 * f0, fs)
    p["w0_tex"] = S.tex_pi(t)
    p["t"] = [t.numerator, t.denominator]
    p["k0"], p["Lmk0"] = int(k0), int(L - k0)
    p["H_tex"] = fmt.tex_num(Fr(A * L, 2))
    p["kap_tex"] = fmt.tex_num(kap) if kap.denominator <= 20 else f"{float(kap):.4g}"
    p["kap_dec"] = f"{float(kap):.4g}"
    p["kmax"] = round(kap)
    p["kmax_mirror"] = Np - round(kap)
    p["lobe_bins"] = fmt.tex_num(Fr(2 * Np, L)) if Fr(2 * Np, L).denominator <= 10 else f"{2 * Np / L:.3g}"
    pk = mag(A, w0, L, Np, round(kap))
    p["peak_pad"] = f"{pk:.1f}"
    p["peak_pad_rel"] = "below" if pk < A * L / 2 - 0.05 else "close to (not exactly)"
    c["k1"], c["k2"], c["Xk"] = int(k0), int(L - k0), float(Fr(A * L, 2))
    c["kap"], c["kmax"] = float(kap), round(kap)


def generate(data):
    p, c = data["params"], data["correct_answers"]
    _part_a(p)
    kind = random.choice(["l22", "fa21"])
    p["l22"], p["fa21"] = kind == "l22", kind == "fa21"
    (_l22 if kind == "l22" else _fa21)(p, c)
