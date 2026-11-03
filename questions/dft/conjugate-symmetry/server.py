"""DFT conjugate symmetry (Midterm 2 family "DFT computation and properties"). Three variant types, one
per past exam problem, each reproducing that problem's shape, wording and number ranges.

fa23 -- same shape as FA2023 MT2 #7 (prove the conjugation property, B^*[N] = B^*[0], two real DFTs
        from one complex DFT):
  (a) MC: the N-point DFT of b^*[n] is B^*[<-k>_N] (written <-k>_N or <N-k>_N); distractors B^*[k]
      (the DFT of b^*[<-n>_N]), B[<-k>_N] (of b[<-n>_N]), B^*[N+k] (an out-of-range index, = B^*[k]),
      -B[k] and the plain flip B^*[N-1-k] (four of the five shown);
  (b) x, y real, z = x + j y with a given N-point DFT Z (N in {4, 6}; Gaussian-integer entries from
      integer x, y for N = 4, from Hermitian Gaussian-integer X, Y for N = 6): Re and Im of X[k1] and
      Y[k1], k1 not in {0, N/2}, via X[k] = (Z[k] + Z^*[<-k>_N])/2, Y[k] = (Z[k] - Z^*[<-k>_N])/(2j);
  (c) X[k2], Y[k2] at a self-paired bin k2 in {0, N/2} (#7(b): the exam's Z^*[N] is Z^*[0]): Re Z, Im Z.
sp21 -- same shape as SP2021 MT2 #5 ({X[k]} = {-2, -4j, 3, 5, 3, 4j, -2j}, N = 7: is x real?), N in 5..8:
  (a) checkbox: five lists, 1-3 of them Hermitian (X[<-k>_N] = X^*[k]: X[0] real, X[N/2] real for
      even N, X[k] and X[N-k] conjugate; general, all-real and imaginary-pair styles). Near misses: the
      plain flip k <-> N-1-k (the exam's list, always shown), pairs equal instead of conjugate or a
      negated pair (at least one: |X[k]| = |X[N-k]| alone does not decide), X[0] complex, X[N/2]
      complex, all-real values that are not circularly even, a broken real part;
  (b) Re and Im of x[0] = (1/N) sum_k X[k] for one more list (the exam key's route);
  (c) MC (fixed order): is that x real? yes / no: x[0] is complex / no: x[0] is real but a pair is not
      conjugate (the key's "real X[0] proves nothing"). Cases: Hermitian; the exam's list (plain flip
      with X[N-1] complex), equal pairs, X[0] complex; full plain flip, broken real part, real values.
fa25 -- same shape as FA2025 MT2 #6(a) (N = 32: X[0] = X_d(2pi), X[2] = X_d(pi/8), X[4] = X_d(pi/8),
        X[6] = X_d^*(-3pi/8), X[8] = X_d^*(-pi/2)), real x of length N in {16, 20, 24, 32, 40, 48}:
  (a) checkbox, fixed order by k: five relations X[k] = X_d(w) or X_d^*(w), 2-4 true. True forms:
      w = w_k, w_k -+ 2pi (incl. X[0] = X_d(2pi)), X_d^*(-w_k), X_d^*(2pi - w_k), X[N-m] = X_d^*(w_m).
      False forms: the pi k/N grid (exam iii), a superfluous conjugate X_d^*(w_k), a missing conjugate
      X_d(-w_k), 4 pi k/N, X[N-m] = X_d(w_m), X_d^*(-pi k/N). Neither "conjugate <=> true" nor
      "ignore the conjugate" nor "only the literal w_k" answers all five;
  (b) the bins k with X[k] = X_d^*(w_a), w_a = 2 pi m_a/N in (0, pi) (k = N - m_a), and
      X[k] = X_d^*(w_b), w_b = -2 pi m_b/N or 2 pi - 2 pi m_b/N (k = m_b); bins and frequencies unused in (a).
Lecture 21, Table 1 (conjugation x^*[n] <-> X^*[<-k>_N]); Lecture 20 (X[k] = X_d(2 pi k/N))."""

import random
from fractions import Fraction as Fr

from ece310 import mc, w8_dft as W

KINDS = ["fa23", "sp21", "fa25"]


# ----------------------------------------------------------------------------- Gaussian rationals (re, im)
def _g(a, b=0):
    return (Fr(a), Fr(b))


def _add(u, v):
    return (u[0] + v[0], u[1] + v[1])


def _sub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def _cj(u):
    return (u[0], -u[1])


def _sc(u, c):
    c = Fr(c)
    return (u[0] * c, u[1] * c)


def _mulj(u):
    """j u"""
    return (-u[1], u[0])


def _div2j(u):
    """u / (2j) = -j u / 2"""
    return (u[1] / 2, -u[0] / 2)


def _tex(u):
    return W.tex_cplx(u[0], u[1])


def _par(u):
    """LaTeX of u, in parentheses unless it is a single non-negative term."""
    t = _tex(u)
    single = (u[1] == 0 and u[0] >= 0) or (u[0] == 0 and u[1] > 0)
    return t if single and "frac" not in t else r"\left(" + t + r"\right)"


def _list_tex(vals):
    return W.tex_list([_tex(v) for v in vals])


def _is_herm(L):
    N = len(L)
    return all(L[k] == _cj(L[(-k) % N]) for k in range(N))


def _first_break(L):
    N = len(L)
    return next(k for k in range(N) if L[k] != _cj(L[(-k) % N]))


_J4 = {0: (1, 0), 1: (0, -1), 2: (-1, 0), 3: (0, 1)}       # e^{-j pi m/2} = (-j)^m


def _dft4(v):
    """Exact 4-point DFT of an integer list (Gaussian integers)."""
    out = []
    for k in range(4):
        s = _g(0)
        for n, a in enumerate(v):
            c, d = _J4[(k * n) % 4]
            s = _add(s, _g(a * c, a * d))
        out.append(s)
    return out


NZ = [-4, -3, -2, -1, 1, 2, 3, 4]


def _herm(N, lo=-3, hi=5, style="mixed"):
    """A random Hermitian list of Gaussian integers (the DFT of a real signal), no zero entries.
    style: mixed (general pairs), real (all values real: circularly even), imag (pairs +-jb)."""
    npairs = (N - 1) // 2
    nzr = [v for v in range(lo, hi + 1) if v]
    while True:
        X0 = random.choice(nzr)
        if style == "mixed":
            pairs = [(random.randint(lo, hi), random.randint(-4, 4)) for _ in range(npairs)]
            if sum(1 for a, b in pairs if b != 0) < max(1, npairs - 1):
                continue
        elif style == "real":
            pairs = [(random.choice(nzr), 0) for _ in range(npairs)]
        else:
            pairs = [(0, random.choice(NZ)) for _ in range(npairs)]
        mid = [(random.choice(nzr), 0)] if N % 2 == 0 else []
        half = [(X0, 0)] + pairs + mid
        if (0, 0) in half or len(set(half)) < len(half):
            continue
        full = half + [(a, -b) for a, b in reversed(pairs)]
        return [_g(a, b) for a, b in full]


# ============================================================================= fa23
CONJ_RIGHT = [r"B^*[\langle -k\rangle_N]", r"B^*[\langle N-k\rangle_N]"]
CONJ_WRONG = [r"B^*[k]", r"B[\langle -k\rangle_N]", r"B^*[N+k]", r"-B[k]", r"B^*[N-1-k]"]


def _seq(t):
    return r"$\{" + t + r"\}_{k=0}^{N-1}$"


def _kind_fa23(p, c):
    N = random.choice([4, 6])
    while True:
        if N == 4:
            x = [random.randint(-3, 4) for _ in range(4)]
            y = [random.randint(-3, 4) for _ in range(4)]
            if sum(v == 0 for v in x + y) > 2:
                continue
            X, Y = _dft4(x), _dft4(y)
        else:
            X, Y = _herm(6, -3, 4), _herm(6, -3, 4)
        Z = [_add(X[k], _mulj(Y[k])) for k in range(N)]
        if any(abs(v) > 12 for u in Z for v in u) or len(set(Z)) < N or _g(0) in Z:
            continue
        h = N // 2
        k1s = [k for k in range(1, N) if k != h and all(v != 0 for v in X[k] + Y[k])
               and abs(X[k][0]) != abs(Y[k][0]) and Z[k] != _cj(Z[(-k) % N])]
        k2s = [k for k in (0, h) if Z[k][0] != 0 and Z[k][1] != 0 and abs(Z[k][0]) != abs(Z[k][1])]
        if k1s and k2s:
            break
    k1, k2 = random.choice(k1s), random.choice(k2s)
    k1m = (-k1) % N
    p["N"], p["Nm1"], p["k1"], p["k1m"], p["k2"] = N, N - 1, k1, k1m, k2
    p["Z"] = [[str(a), str(b)] for a, b in Z]
    p["Z_tex"] = _list_tex(Z)
    p["k2_zero"] = k2 == 0
    p["k2_mid"] = k2 != 0

    corr = random.choice(CONJ_RIGHT)
    p["conj_right_tex"] = corr
    p["conj_choices"] = mc.choices(_seq(corr), [_seq(t) for t in random.sample(CONJ_WRONG, len(CONJ_WRONG))], n=5)

    Zk, Zm = Z[k1], Z[k1m]
    s = _add(Zk, _cj(Zm))
    d = _sub(Zk, _cj(Zm))
    Xk, Yk = _sc(s, Fr(1, 2)), _div2j(d)
    assert Xk == X[k1] and Yk == Y[k1]
    p["X_line_tex"] = (rf"X[{k1}] = \tfrac12\left(Z[{k1}] + Z^*[{k1m}]\right) = \tfrac12\left[{_par(Zk)} + {_par(_cj(Zm))}\right]"
                       rf" = \tfrac12\left({_tex(s)}\right) = {_tex(Xk)}")
    p["Y_line_tex"] = (rf"Y[{k1}] = \tfrac{{1}}{{2j}}\left(Z[{k1}] - Z^*[{k1m}]\right) = \tfrac{{1}}{{2j}}\left[{_par(Zk)} - {_par(_cj(Zm))}\right]"
                       rf" = \tfrac{{1}}{{2j}}\left({_tex(d)}\right) = -\tfrac{{j}}{{2}}\left({_tex(d)}\right) = {_tex(Yk)}")
    p["Zk1_tex"], p["Zk1m_tex"], p["Zk1mc_tex"] = _tex(Zk), _tex(Zm), _tex(_cj(Zm))
    p["ReZk1_tex"], p["ImZk1_tex"] = W.tex_cplx(Zk[0], 0), W.tex_cplx(Zk[1], 0)
    Z2 = Z[k2]
    p["Zk2_tex"], p["Zk2c_tex"] = _tex(Z2), _tex(_cj(Z2))
    p["Xs_tex"], p["Ys_tex"] = W.tex_cplx(Z2[0], 0), W.tex_cplx(Z2[1], 0)
    c["reX"], c["imX"] = float(Xk[0]), float(Xk[1])
    c["reY"], c["imY"] = float(Yk[0]), float(Yk[1])
    c["Xs"], c["Ys"] = float(Z2[0]), float(Z2[1])


# ============================================================================= sp21
def _plain_flip(N, lo=-3, hi=5):
    """L[k] = L^*[N-1-k] for every k (the pairing of SP2021 #5), L[0] = L[N-1] real, and the self-paired
    entry of the plain flip ((N-1)/2 for odd N, or L[N/2-1] = L[N/2] for even N) real too."""
    nzr = [v for v in range(lo, hi + 1) if v]
    while True:
        L = [None] * N
        a = random.choice(nzr)
        L[0] = L[N - 1] = _g(a)
        if N % 2:
            L[(N - 1) // 2] = _g(random.choice(nzr))
        else:
            b = random.choice(nzr)
            L[N // 2 - 1] = L[N // 2] = _g(b)
        for k in range(1, N):
            if L[k] is None:
                u = _g(random.randint(lo, hi), random.choice(NZ))
                L[k], L[N - 1 - k] = u, _cj(u)
        if not _is_herm(L) and len(set(L)) >= N // 2 + 1:
            return L


def _bad_list(kind, N):
    h = N // 2
    while True:
        if kind == "plain_flip":
            L = _plain_flip(N)
        elif kind == "x0_complex":
            L = _herm(N)
            L[0] = (L[0][0], Fr(random.choice([-3, -2, -1, 1, 2, 3])))
        elif kind == "mid_complex":
            L = _herm(N)
            L[h] = (L[h][0], Fr(random.choice([-3, -2, -1, 1, 2, 3])))
        elif kind == "equal_pairs":
            L = _herm(N)
            for k in range(1, (N + 1) // 2):
                L[N - k] = L[k]
        elif kind == "neg_pair":
            L = _herm(N)
            ks = [k for k in range(1, (N + 1) // 2) if L[k][0] != 0 and L[k][1] != 0]
            if not ks:
                continue
            k = random.choice(ks)
            L[N - k] = _sc(_cj(L[k]), -1)
        elif kind == "re_broken":
            L = _herm(N)
            k = random.randint(1, (N - 1) // 2)
            L[N - k] = _add(L[N - k], _g(random.choice([-2, -1, 1, 2])))
        elif kind == "realent_nonsym":
            nzr = [v for v in range(-3, 6) if v]
            L = [_g(random.choice(nzr)) for _ in range(N)]
        elif kind == "sp21":                      # the exam's list: plain flip, X[N-1] complex
            L = _plain_flip(N)
            L[N - 1] = _g(random.randint(-3, 3), random.choice(NZ))
        else:
            raise KeyError(kind)
        if not _is_herm(L) and _g(0) not in L:
            return L


LEAD = {
    "good": r"Hermitian: $X[0]${mid} real and $X[{N} - k] = X^*[k]$ for every pair, so $x[n]$ is real.",
    "x0_complex": r"$X[0] = \sum_n x[n]$ must be real.",
    "mid_complex": r"$X[{h}] = \sum_n (-1)^n x[n]$ pairs with itself ($\langle -{h}\rangle_{{N}} = {h}$) and must be real.",
    "equal_pairs": r"The pairs are equal, $X[{N} - k] = X[k]$, not conjugate.",
    "plain_flip": r"Conjugates sit at $k$ and ${Nm1} - k$ (a plain flip, as in SP2021 #5), not at $k$ and $\langle -k\rangle_{{N}} = {N} - k$.",
    "sp21": r"Paired by a plain flip ($k \leftrightarrow {Nm1} - k$, as in SP2021 #5), and $X[0]$, $X[{Nm1}]$ do not even match that way.",
    "realent_nonsym": r"All values are real, but that does not make $x[n]$ real: the real parts must also be circularly even.",
    "neg_pair": r"One pair has $X[{N} - k] = -X^*[k]$: the sign is wrong.",
    "re_broken": r"The imaginary parts are odd, but the real parts are not even.",
}


def _fill(t, N):
    mid = rf", $X[{N // 2}]$" if N % 2 == 0 else ""
    return (t.replace("{mid}", mid).replace("{Nm1}", str(N - 1)).replace("{N}", str(N))
            .replace("{h}", str(N // 2)))


def _detail(L):
    N = len(L)
    k = _first_break(L)
    if k == 0 or 2 * k == N:
        return rf"Here $X[{k}] = {_tex(L[k])}$ is not real."
    return rf"Here $X[{k}] = {_tex(L[k])}$ but $X^*[{(-k) % N}] = {_tex(_cj(L[(-k) % N]))}$."


def _kind_sp21(p, c):
    N = random.choice([5, 6, 7, 7, 8])
    h = N // 2
    pool = ["x0_complex", "x0_complex", "equal_pairs", "realent_nonsym", "re_broken", "neg_pair"]
    if N % 2 == 0:
        pool += ["mid_complex", "mid_complex"]
    while True:
        n_true = random.choice([1, 2, 2, 3])
        styles = [random.choice(["mixed", "mixed", "real", "imag"]) for _ in range(n_true)]
        bad = [random.choice(["plain_flip", "plain_flip", "sp21"]), random.choice(["equal_pairs", "equal_pairs", "neg_pair"])]
        while len(bad) < 5 - n_true:
            k = random.choice(pool)
            if k not in bad:
                bad.append(k)
        bad = bad[:5 - n_true]
        opts = [("good", _herm(N, style=s)) for s in styles] + [(k, _bad_list(k, N)) for k in bad]
        if len({tuple(L) for _, L in opts}) < 5:
            continue
        # shortcuts that must fail: X[0] (and X[N/2]) real only; that plus |X[k]| = |X[N-k]|
        selfp = [L[0][1] == 0 and (N % 2 or L[h][1] == 0) for _, L in opts]
        mags = [s and all(L[k][0] ** 2 + L[k][1] ** 2 == L[-k % N][0] ** 2 + L[-k % N][1] ** 2 for k in range(N))
                for s, (_, L) in zip(selfp, opts)]
        if sum(selfp) > n_true and sum(mags) > n_true:
            break
    random.shuffle(opts)
    p["N"], p["Nm1"], p["h"], p["N_even"] = N, N - 1, h, N % 2 == 0
    p["real_choices"] = [{"text": "$" + _list_tex(L) + "$", "correct": k == "good"} for k, L in opts]
    p["real_kinds"] = [k for k, _ in opts]
    p["real_reasons"] = [{"list_tex": _list_tex(L), "verdict": "real" if k == "good" else "not real",
                          "why_html": _fill(LEAD[k], N) + ("" if k == "good" else " " + _detail(L))} for k, L in opts]

    # (b), (c): one more list, the exam's question
    case = random.choice(["R", "C", "P"])
    while True:
        if case == "R":
            sub, L = "good", _herm(N)
        elif case == "C":
            sub = random.choice(["sp21", "sp21", "equal_pairs", "x0_complex"])
            L = _bad_list(sub, N)
        else:
            sub = random.choice(["plain_flip", "plain_flip", "re_broken", "realent_nonsym"])
            L = _bad_list(sub, N)
        S = (sum(u[0] for u in L), sum(u[1] for u in L))
        x0 = _sc(S, Fr(1, N))
        if (x0[1] != 0) != (case == "C") or S[0] == 0 or x0[0] == L[0][0]:
            continue
        if any(tuple(L) == tuple(M) for _, M in opts):
            continue
        break
    p["case"], p["case_kind"] = case, sub
    p["XL"] = [[str(a), str(b)] for a, b in L]
    p["XL_tex"] = _list_tex(L)
    terms = " + ".join(_tex(u) if i == 0 else _par(u) for i, u in enumerate(L))
    p["x0_line_tex"] = (rf"x[0] = \frac{{1}}{{{N}}}\sum_{{k=0}}^{{{N - 1}}} X[k]\,e^{{j\frac{{2\pi}}{{{N}}}k\cdot 0}}"
                        rf" = \frac{{1}}{{{N}}}\left[{terms}\right] = \frac{{1}}{{{N}}}\left({_tex(S)}\right)"
                        + ("" if x0 == S else rf" = {_tex(x0)}"))
    p["x0_tex"] = _tex(x0)
    c["rex0"], c["imx0"] = float(x0[0]), float(x0[1])
    opts_c = [
        rf"Yes: $X[k] = X^*[\langle -k\rangle_{{{N}}}]$ for every $k$, so $x[n]$ is real.",
        r"No: $x[0]$ is not real, so $x[n]$ is not real.",
        rf"No: $x[0]$ is real, but $X[k] \ne X^*[\langle -k\rangle_{{{N}}}]$ for some $k$.",
    ]
    right = {"R": 0, "C": 1, "P": 2}[case]
    p["isreal_choices"] = [{"text": t, "correct": i == right} for i, t in enumerate(opts_c)]
    if case == "R":
        why = _fill(LEAD["good"], N)
    elif case == "C":
        why = (rf"$x[0] = {_tex(x0)}$ is complex, and a real signal has a real $x[0]$. (The symmetry test agrees: "
               + _detail(L)[0].lower() + _detail(L)[1:] + ")")
    else:
        why = (rf"$x[0] = {_tex(x0)}$ is real, but that proves nothing: the test is $X[k] = X^*[\langle -k\rangle_{{{N}}}]$ for every $k$. "
               + _fill(LEAD[sub], N) + " " + _detail(L))
    p["isreal_why_html"] = why
    p["isreal_verdict"] = ["real", "not real", "not real"][right]


# ============================================================================= fa25
def _kind_fa25(p, c):
    while True:
        if _try_fa25(p, c):
            return


def _try_fa25(p, c):
    N = random.choice([16, 20, 24, 32, 40, 48])
    h = N // 2

    def w(k):
        return Fr(2 * k, N)

    lo = list(range(1, h))

    def nice(k):                                          # exam-like frequencies: w_k/pi with denominator <= 8
        return Fr(2 * k, N).denominator <= 8
    T = {   # name: (true?, candidate k, k -> (shown k, conjugate?, w/pi))
        "direct": (True, [k for k in range(1, N) if k != h], lambda k: (k, False, w(k))),
        "x0_2pi": (True, [0], lambda k: (0, False, Fr(random.choice([2, 2, -2])))),
        "neg": (True, list(range(h + 1, N)), lambda k: (k, False, w(k) - 2)),
        "plus2pi": (True, lo, lambda k: (k, False, w(k) + 2)),
        "conj_neg": (True, lo, lambda k: (k, True, -w(k))),
        "conj_mirror": (True, lo, lambda k: (k, True, 2 - w(k))),
        "conj_upper": (True, lo, lambda m: (N - m, True, w(m))),
        "pi_grid": (False, lo, lambda k: (k, False, Fr(k, N))),
        "conj_same": (False, lo, lambda k: (k, True, w(k))),
        "noconj_neg": (False, lo, lambda k: (k, False, -w(k))),
        "double": (False, list(range(1, (h + 1) // 2)), lambda k: (k, False, 2 * w(k))),
        "upper_noconj": (False, lo, lambda m: (N - m, False, w(m))),
        "conj_pi_grid": (False, lo, lambda k: (k, True, -Fr(k, N))),
    }
    tnames = [n for n, v in T.items() if v[0]]
    fnames = [n for n, v in T.items() if not v[0]]
    n_true = random.choice([2, 3, 3, 4, 4])
    names = random.sample(tnames, n_true)
    fn = [n for n in fnames if n != "pi_grid"]
    if random.random() < 0.7:
        names += ["pi_grid"] + random.sample(fn, 4 - n_true)
    else:
        names += random.sample(fn, 5 - n_true)
    rels, used = [], set()
    for nm in names:
        truth, cands, f = T[nm]
        cands = [k for k in cands if k not in used and (N - k) not in used and f(k)[0] not in used]
        if not cands:
            return False
        pretty = [k for k in cands if nice(k)]
        k, conj, t = f(random.choice(pretty if pretty and random.random() < 0.75 else cands))
        used |= {k, (N - k) % N}
        eff = -t if conj else t
        ok = (eff - w(k)) % 2 == 0
        assert ok == truth, (nm, k, conj, t)
        rels.append((k, conj, t, ok, nm))
    shortcuts = [lambda r: r[1], lambda r: not r[1], lambda r: (r[2] - w(r[0])) % 2 == 0,
                 lambda r: r[2] == w(r[0]) and not r[1]]
    if any(all(s(r) == r[3] for r in rels) for s in shortcuts):
        return False
    rels.sort(key=lambda r: r[0])
    # (b): two conjugate look-ups at bins and frequencies that (a) does not use
    touched = set()
    for k, conj, t, ok, nm in rels:
        touched |= {k, (N - k) % N}
        q = t * N / 2
        if q.denominator == 1:
            touched |= {int(q) % N, int(-q) % N}
    free = [m for m in lo if m not in touched and (N - m) not in touched and 3 * m != N]
    if len(free) < 2:
        return False
    pretty = [m for m in free if nice(m)]
    ma, mb = random.sample(pretty if len(pretty) >= 2 and random.random() < 0.75 else free, 2)
    tb = random.choice([-w(mb), -w(mb), 2 - w(mb)])
    p["N"], p["Nm1"], p["h"] = N, N - 1, h
    p["wk_tex"] = W.tex_pi_frac(Fr(2, N))
    p["rel_choices"], p["rel_lines"], p["rels"] = [], [], []
    for k, conj, t, ok, nm in rels:
        star = "^*" if conj else ""
        rt = rf"X[{k}] = X_d{star}\!\left({W.tex_pi_frac(t)}\right)"
        p["rel_choices"].append({"text": "$" + rt + "$", "correct": ok})
        p["rels"].append({"k": k, "conj": conj, "t": [t.numerator, t.denominator], "kind": nm})
        wk = W.tex_pi_frac(w(k))
        line = rf"$\omega_{{{k}}} = \frac{{2\pi\cdot {k}}}{{{N}}} = {wk}$; " if k else r"$\omega_0 = 0$; "
        if conj:
            line += (rf"$x$ real, so $X_d^*\!\left({W.tex_pi_frac(t)}\right) = X_d\!\left({W.tex_pi_frac(-t)}\right)$"
                     rf" and ${W.tex_pi_frac(-t)} \equiv {W.tex_pi_frac((-t) % 2)}$ (mod $2\pi$): ")
        else:
            line += rf"${W.tex_pi_frac(t)} \equiv {W.tex_pi_frac(t % 2)}$ (mod $2\pi$): "
        line += (rf"the same frequency, $= \omega_{{{k}}}$." if ok else rf"a different frequency, $\ne \omega_{{{k}}}$.")
        if nm == "pi_grid":
            line += rf" (The bins are $\frac{{2\pi}}{{{N}}}$ apart, not $\frac{{\pi}}{{{N}}}$.)"
        elif nm == "conj_same":
            line += r" (The conjugate moves the frequency to $-\omega$; it is not free.)"
        elif nm in ("noconj_neg", "upper_noconj"):
            line += rf" (Without the conjugate this is $X[{(N - k) % N}]$, not $X[{k}]$.)"
        p["rel_lines"].append({"rel_tex": rt, "verdict": "True" if ok else "False", "why_html": line})
    p["wa_tex"], p["wb_tex"] = W.tex_pi_frac(w(ma)), W.tex_pi_frac(tb)
    p["wa_red_tex"] = W.tex_pi_frac((-w(ma)) % 2)
    p["wb_red_tex"] = W.tex_pi_frac((-tb) % 2)
    p["wa_chain_tex"] = rf"X_d^*\!\left({p['wa_tex']}\right) = X_d\!\left({W.tex_pi_frac(-w(ma))}\right) = X_d\!\left({p['wa_red_tex']}\right)"
    p["wb_chain_tex"] = (rf"X_d^*\!\left({p['wb_tex']}\right) = X_d\!\left({W.tex_pi_frac(-tb)}\right)"
                         + ("" if -tb == (-tb) % 2 else rf" = X_d\!\left({p['wb_red_tex']}\right)"))
    p["ma"], p["mb"] = ma, mb
    p["ta"], p["tb"] = [w(ma).numerator, w(ma).denominator], [tb.numerator, tb.denominator]
    c["ka"], c["kb"] = N - ma, mb
    p["ka_line_tex"] = rf"k = \frac{{{N}}}{{2\pi}}\cdot {p['wa_red_tex']} = {N - ma}"
    p["kb_line_tex"] = rf"k = \frac{{{N}}}{{2\pi}}\cdot {p['wb_red_tex']} = {mb}"
    return True


def generate(data):
    p, c = data["params"], data["correct_answers"]
    kind = random.choice(KINDS)
    p["kind"] = kind
    for k in KINDS:
        p["k_" + k] = kind == k
    {"fa23": _kind_fa23, "sp21": _kind_sp21, "fa25": _kind_fa25}[kind](p, c)
