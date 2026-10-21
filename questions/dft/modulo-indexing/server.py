"""Modulo indexing <n>_N and a DFT-domain permutation (Lecture 21, Example 1 and Example 3).

x: length N in {5, 6}, small integers, x[0] != 0, all 2N rotations/reflections distinct (so every
wrong index rule gives a different row).
(a) x[<-n>_N]   (b) x[<n + m>_N]   (c) x[<m2 - n>_N]   (1 x N rows, pl-matrix-input)
(d) Y[k] given as a list of +-X_i (Example 3 style): Y[k] = (-1)^{a k} X[<tau k + e>_N].
    X[<tau k + e>] <-> e^{j 2 pi l n/N} x[<tau n>] with l = -tau e; (-1)^k = e^{-j 2 pi k (N/2)/N} then
    delays by N/2:  y[n] = c e^{j 2 pi l n/N} x[<sigma n + d>_N],  sigma = tau,
    c = (-1)^l, d = -sigma N/2 if a = 1 (else c = 1, d = 0).
    Real kinds (N = 6, l in {0, 3}: y real integer) -> matrix row 'yd':
        alt     (-1)^k X[k]          -> x[<n - 3>]
        half    X[<k - 3>]           -> (-1)^n x[n]
        altrev  (-1)^k X[<-k>]       -> x[<3 - n>]
        althalf (-1)^k X[<k - 3>]    -> -(-1)^n x[<n - 3>]
    Complex kinds (frequency shifts, y complex) -> multiple choice 'op' among formulas
    c e^{j 2 pi l n/N} x[<sigma n + d>]:
        fsh X[<k - m>], fshm X[<k + m>], frev X[<m - k>]  (N = 5, 6)
        ex3 (-1)^k X[<k + m>], ex3b (-1)^k X[<k - m>]       (N = 6; Example 3 is ex3 with m = 2)
    (The input of (d) depends on the kind, as the coordinator asked: a real y is typed, a complex one chosen.)"""

import random
from fractions import Fraction as F

import numpy as np
import prairielearn as pl
from ece310 import fmt, mc
from ece310.w5_freqresp import tex_ej
from ece310.w8_dft import tex_list, tex_mod


def _rand_x(N):
    while True:
        x = [random.randint(-2, 6) for _ in range(N)]
        if x[0] == 0 or sum(v == 0 for v in x) > 1 or len(set(x)) < N - 1:
            continue
        orbit = {tuple(x[(s * n + d) % N] for n in range(N)) for s in (1, -1) for d in range(N)}
        if len(orbit) == 2 * N:
            return x


def _idx(sigma, d, var="n"):
    """sigma*var + d as an index expression (without the <>)."""
    if sigma == 1:
        return var if d == 0 else (f"{var} + {d}" if d > 0 else f"{var} - {-d}")
    if d == 0:
        return f"-{var}"
    return f"{d} - {var}" if d > 0 else f"-{var} - {-d}"


def _idx_tex(sigma, d, N, var="n"):
    if sigma == 1 and d == 0:
        return var
    return tex_mod(_idx(sigma, d, var), N)


def _apply(x, c, ell, sigma, d):
    """y[n] = c e^{j 2 pi ell n/N} x[<sigma n + d>_N] as complex numpy array."""
    N = len(x)
    n = np.arange(N)
    return c * np.exp(2j * np.pi * ell * n / N) * np.asarray(x, float)[(sigma * n + d) % N]


def _formula_tex(c, ell, sigma, d, N):
    ell = _norm(ell, N)
    e = tex_ej(F(2 * ell, N), "n") if ell else ""
    return ("-" if c < 0 else "") + (e + r"\," if e else "") + f"x[{_idx_tex(sigma, d, N)}]"


def _real_formula_tex(c, ell, sigma, d, N):
    """For ell in {0, N/2}: (-1)^n instead of e^{j pi n}."""
    pre = ("-" if c < 0 else "") + (r"(-1)^n\," if ell % N else "")
    return pre + f"x[{_idx_tex(sigma, d, N)}]"


def _norm(ell, N):
    ell %= N
    return ell - N if 2 * ell > N else ell


def _formula(a, tau, e, N):
    ell = _norm(-tau * e, N)
    sigma = tau
    if a:
        return (-1) ** (ell % 2), ell, sigma, -sigma * (N // 2)
    return 1, ell, sigma, 0


REAL_KINDS = {"alt": (1, 1, 0), "half": (0, 1, -3), "altrev": (1, -1, 0), "althalf": (1, 1, -3)}


def _row(v):
    return pl.to_json(np.array([[float(t) for t in v]]))


def _seq_steps(x, N, sigma, d):
    idx = [(sigma * n + d) % N for n in range(N)]
    return idx, [x[i] for i in idx]


def generate(data):
    p, c = data["params"], data["correct_answers"]
    N = random.choice([5, 6])
    x = _rand_x(N)
    p["N"], p["Nm1"], p["x"] = N, N - 1, x
    p["x_tex"] = fmt.tex_seq(x, 0)

    # (d) kind first: (b), (c) must not repeat its answer
    if N == 6 and random.random() < 0.55:
        kind = random.choice(list(REAL_KINDS))
        a, tau, e = REAL_KINDS[kind]
    else:
        kind = random.choice(["fsh", "fshm", "frev"] + (["ex3", "ex3b", "ex3"] if N == 6 else []))
        ms = [1, 2, 4, 5] if N == 6 else [1, 2, 3, 4]
        m = random.choice(ms)
        a, tau, e = {"fsh": (0, 1, -m), "fshm": (0, 1, m), "frev": (0, -1, m), "ex3": (1, 1, m), "ex3b": (1, 1, -m)}[kind]
    real = kind in REAL_KINDS
    ms_bc = [v for v in range(1, N) if not (N == 6 and v == 3)]
    mb = random.choice(ms_bc)
    mc_ = random.choice(ms_bc)

    rows = []
    for name, (sigma, d), lab in (("ya", (-1, 0), "a"), ("yb", (1, mb), "b"), ("yc", (-1, mc_), "c")):
        idx, vals = _seq_steps(x, N, sigma, d)
        c[name] = _row(vals)
        p[f"{name}_idx_tex"] = tex_mod(_idx(sigma, d), N)
        rows.append({"lab": lab, "expr_tex": f"x[{tex_mod(_idx(sigma, d), N)}]",
                     "idx_tex": tex_list([str(i) for i in idx]), "val_tex": tex_list([str(v) for v in vals])})
    p["rows"] = rows
    p["mb"], p["mc"] = mb, mc_

    # (d)
    Yidx = [(tau * k + e) % N for k in range(N)]
    Ysgn = [(-1) ** (a * k) for k in range(N)]
    p["Y_list"] = {"idx": Yidx, "sgn": Ysgn}
    p["Y_tex"] = tex_list([("-" if s < 0 else "") + f"X_{{{i}}}" for i, s in zip(Yidx, Ysgn)])
    p["X_tex"] = tex_list([f"X_{{{k}}}" for k in range(N)])
    cc, ell, sigma, d = _formula(a, tau, e, N)
    p["d_formula"] = {"c": cc, "ell": ell, "sigma": sigma, "d": d}
    p["d_real"], p["d_mc"], p["kind"] = real, not real, kind
    pat = ("(-1)^k\\," if a else "") + f"X[{_idx_tex(tau, e, N, 'k')}]"
    p["pattern_tex"] = pat
    p["check_tex"] = rf"Y[1] = {'-' if a else ''}X_{{{Yidx[1]}}}"
    steps = []
    if tau == 1 and e != 0:
        steps.append(rf"X[{_idx_tex(1, e, N, 'k')}] \leftrightarrow {tex_ej(F(2 * _norm(-e, N), N), 'n')}\,x[n]"
                     r"\quad(\text{frequency shift: } e^{j\frac{2\pi}{N}\ell n}x[n] \leftrightarrow X[\langle k - \ell\rangle_N])")
    elif tau == -1 and e == 0:
        steps.append(rf"X[{tex_mod('-k', N)}] \leftrightarrow x[{tex_mod('-n', N)}]\quad(\text{{circular reversal}})")
    elif tau == -1:
        steps.append(rf"X[{tex_mod(f'{e} - k', N)}] = Z[{tex_mod(f'k - {e}', N)}],\ Z[k] = X[{tex_mod('-k', N)}] \leftrightarrow x[{tex_mod('-n', N)}]"
                     rf"\ \Rightarrow\ {tex_ej(F(2 * _norm(e, N), N), 'n')}\,x[{tex_mod('-n', N)}]")
    if a:
        steps.append(rf"(-1)^k = e^{{-j\frac{{2\pi}}{{{N}}}k\cdot {N // 2}}}:\ \text{{circular delay by }} {N // 2}\ (n \to n - {N // 2})")
    shown = _real_formula_tex(cc, ell, sigma, d, N) if real else _formula_tex(cc, ell, sigma, d, N)
    if a and ell % N:
        steps.append(rf"{tex_ej(F(2 * ell, N), f'(n - {N // 2})')} = (-1)^{{{ell}}}\,{tex_ej(F(2 * ell, N), 'n')}")
    p["d_steps"] = [{"s": s} for s in steps]
    p["y_formula_tex"] = shown
    if real:
        y = _apply(x, cc, ell, sigma, d).real
        yi = [int(round(v)) for v in y]
        c["yd"] = _row(yi)
        p["yd_tex"] = tex_list([str(v) for v in yi])
    else:
        cands = [(cc, -ell, sigma, d),                      # wrong direction
                 (1, 0, 1, -e if tau == 1 else e),          # time shift instead of a frequency shift
                 (1, 0, 1, e if tau == 1 else -e)]
        if sigma == -1:
            cands += [(cc, ell, 1, d), (cc, -ell, -1, d)]
        if a:
            cands += [(cc, ell, sigma, 0), (-cc, ell, sigma, d), (1, 0, 1, -(N // 2))]
        else:
            cands += [(cc, ell, sigma, -(N // 2)) if N == 6 else (cc, ell, -sigma, 0), (1, 0, -1, 0)]
        right = _apply(x, cc, ell, sigma, d)
        uniq = []
        for t in cands:
            v = _apply(x, *t)
            if np.max(np.abs(v - right)) > 1e-6 and all(np.max(np.abs(v - _apply(x, *u))) > 1e-6 for u in uniq):
                uniq.append(t)
        random.shuffle(uniq)
        dis = uniq[:4]
        p["op_choices"] = mc.choices(f"$y[n] = {shown}$", [f"$y[n] = {_formula_tex(*t, N)}$" for t in dis], n=5)
