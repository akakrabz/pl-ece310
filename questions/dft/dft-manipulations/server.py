"""DFT manipulations: a new DFT is built from the values X[k] of a short real signal; find the signal
(Midterm 2 family "DFT properties"; the shapes not already in dft/dft-properties and dft/modulo-indexing).

Variant types, one per exam/homework problem:
  dup   FA2025 MT2 #7 made numeric: {Y[k]}_{k=0}^{2N-1} = {X_0..X_{N-1}, X_0..X_{N-1}} (the DFT written twice)
        -> y[2n] = x[n], y[2n+1] = 0 (N = 4, 5, 6).
  inter the dual of FA2025 #7 (zeros between the DFT values): Y[2k] = X[k], Y[2k+1] = 0 -> y[n] = x[<n>_N]/2
        (x has even entries so y is an integer row).
  dec   the dual of zero padding (FA2025 #6(b), FA2024 #1(d) "X[k] = Y[2k]" read backwards): keep every other
        DFT value, Y[k] = X[2k], k < N/2 -> time aliasing y[n] = x[n] + x[n + N/2] (N = 6, 8).
  hw9   HW9 #2(b) in general form: Y[k] = X[k] + X[<k - N/2>_N] (for X = {A, B, C, D, 0, 0, 0, 0} this is
        {A, B, C, D, A, B, C, D}) -> y[n] = (1 + (-1)^n) x[n] (N = 4, 6, 8).
  fa21  FA2021 MT2 #3(b): x of length 8 with x[n + 4] = x[n] (X[k] = 0 for odd k), y[n] = x[<m - n>_8] c^n
        with c in {j, -j, -1}; Y[k] for k = 0, 2, 4, 6 (circular reversal + shift as in FA2024 #8, then a
        modulation: Y[k] = e^{-j 2 pi m (k - l)/8} X[<l - k>_8], c = e^{j 2 pi l / 8}).
(a) the numbers: y[n] as a row (types dup..hw9; integers), or [Y[0], Y[2], Y[4], Y[6]] (fa21; Gaussian
    integers, pl-matrix-input with allow-complex).
(b) multiple choice: the general statement that holds for every x (the FA2025 #7 "show that" result for
    dup), with distractors taken from the other rules (the dual rule, zero padding vs. upsampling, decimation
    in time instead of aliasing, the wrong modulation direction, a missing reversal).
The phase-ramp / sign-flip / permuted-list shapes (FA2019 #4, SP2021 #4, SP2025 #7(b) and #8, SP2023 #9(b),
FA2024 #8 forwards) are in dft/dft-properties and dft/modulo-indexing; conjugation (FA2023 #7, SP2021 #5)
is in dft/conjugate-symmetry."""

import random
from fractions import Fraction as Fr

import numpy as np
import prairielearn as pl
from ece310 import fmt, mc, w8_dft as W

KINDS = ["dup", "dup", "dup", "inter", "inter", "dec", "dec", "hw9", "hw9", "fa21", "fa21"]

N2 = r"\tfrac{N}{2}"
GEN = {   # id -> general statement (LaTeX/HTML), for every length-N signal x[n]
    "dup": r"$y[2n] = x[n]$ and $y[2n+1] = 0$ for $0 \le n \le N-1$",
    "dup2": r"$y[2n] = 2x[n]$ and $y[2n+1] = 0$ for $0 \le n \le N-1$",
    "dupodd": r"$y[2n] = 0$ and $y[2n+1] = x[n]$ for $0 \le n \le N-1$",
    "half": r"$y[n] = \tfrac12\,x[\langle n\rangle_N]$ for $0 \le n \le 2N-1$",
    "per": r"$y[n] = x[\langle n\rangle_N]$ for $0 \le n \le 2N-1$ (two copies of $x$)",
    "zpad": r"$y[n] = x[n]$ for $0 \le n \le N-1$ and $y[n] = 0$ for $N \le n \le 2N-1$",
    "alias": rf"$y[n] = x[n] + x[n + {N2}]$ for $0 \le n \le {N2} - 1$",
    "alias_half": rf"$y[n] = \tfrac12\left(x[n] + x[n + {N2}]\right)$ for $0 \le n \le {N2} - 1$",
    "alias_diff": rf"$y[n] = x[n] - x[n + {N2}]$ for $0 \le n \le {N2} - 1$",
    "decim": rf"$y[n] = x[2n]$ for $0 \le n \le {N2} - 1$",
    "trunc": rf"$y[n] = x[n]$ for $0 \le n \le {N2} - 1$",
    "mod1": r"$y[n] = \left(1 + (-1)^n\right)x[n]$ for $0 \le n \le N-1$",
    "mod1m": r"$y[n] = \left(1 - (-1)^n\right)x[n]$ for $0 \le n \le N-1$",
    "tsum": rf"$y[n] = x[n] + x[\langle n - {N2}\rangle_N]$ for $0 \le n \le N-1$",
    "two": r"$y[n] = 2x[n]$ for $0 \le n \le N-1$",
    "modshift": rf"$y[n] = \left(1 + (-1)^n\right)x[\langle n - {N2}\rangle_N]$ for $0 \le n \le N-1$",
}
DISTRACT = {
    "dup": ["half", "per", "zpad", "dup2", "dupodd"],
    "inter": ["per", "dup", "zpad", "dup2"],
    "dec": ["decim", "alias_half", "alias_diff", "trunc"],
    "hw9": ["tsum", "two", "mod1m", "modshift"],
}
RIGHT = {"dup": "dup", "inter": "half", "dec": "alias", "hw9": "mod1"}
CPOW = {2: ("j^n", 1j), 6: ("(-j)^n", -1j), 4: ("(-1)^n", -1)}      # l -> (c^n LaTeX, c)


def _rand_x(N, vals):
    while True:
        x = [random.choice(vals) for _ in range(N)]
        if x[0] == 0 or x[-1] == 0 or sum(v == 0 for v in x) > 1 or len(set(x)) < N - 1:
            continue
        return x


def _items(names):
    return W.tex_list(names)


def _row(vals):
    return pl.to_json(np.array([[float(v) for v in vals]]))


# ----------------------------------------------------------------------------- types dup .. hw9
def _real_kind(p, c, kind):
    if kind == "dup":
        N = random.choice([4, 5, 6])
        x = _rand_x(N, [-3, -2, -1, 1, 2, 3, 4, 5])
        M = 2 * N
        Y = [f"X_{{{k % N}}}" for k in range(M)]
        y = [x[n // 2] if n % 2 == 0 else 0 for n in range(M)]
    elif kind == "inter":
        N = random.choice([4, 5, 6])
        x = _rand_x(N, [-6, -4, -2, 2, 4, 6, 8])
        M = 2 * N
        Y = [f"X_{{{k // 2}}}" if k % 2 == 0 else "0" for k in range(M)]
        y = [x[n % N] // 2 for n in range(M)]
    elif kind == "dec":
        N = random.choice([6, 8])
        while True:
            x = _rand_x(N, [-3, -2, -1, 1, 2, 3, 4, 5])
            h = N // 2
            y = [x[n] + x[n + h] for n in range(h)]
            if sum(v != 0 for v in y) >= h - 1 and y != [x[2 * n] for n in range(h)]:
                break
        M = N // 2
        Y = [f"X_{{{2 * k}}}" for k in range(M)]
    else:
        N = random.choice([4, 6, 8])
        x = _rand_x(N, [-3, -2, -1, 1, 2, 3, 4, 5])
        h = N // 2
        M = N
        Y = [f"X_{{{k}}} + X_{{{(k - h) % N}}}" for k in range(N)]
        y = [(1 + (-1) ** n) * x[n] for n in range(N)]
    p["N"], p["Nm1"], p["M"], p["Mm1"], p["h"] = N, N - 1, M, M - 1, N // 2
    p["x"] = x
    p["Y_list_tex"] = _items(Y)
    p["y"] = y
    p["y_tex"] = W.tex_list([str(v) for v in y])
    c["y"] = _row(y)

    right = RIGHT[kind]
    dis = random.sample(DISTRACT[kind], min(4, len(DISTRACT[kind])))
    opts = mc.choices(GEN[right], [GEN[d] for d in dis], n=5)
    inv = {v: k for k, v in GEN.items()}
    p["gen_choices"] = opts
    p["gen_ids"] = [inv[o["text"]] for o in opts]

    # worked numbers
    if kind == "dup":
        p["work_tex"] = W.tex_list([f"x[{n // 2}]" if n % 2 == 0 else "0" for n in range(M)])
    elif kind == "inter":
        p["work_tex"] = W.tex_list([rf"\tfrac12 x[{n % N}]" for n in range(M)])
    elif kind == "dec":
        h = N // 2
        p["work_tex"] = W.tex_list([f"x[{n}] + x[{n + h}]" for n in range(h)])
        p["work2_tex"] = W.tex_list([f"{x[n]} {'+' if x[n + h] >= 0 else '-'} {abs(x[n + h])}" for n in range(h)])
    else:
        p["work_tex"] = W.tex_list([f"2x[{n}]" if n % 2 == 0 else "0" for n in range(N)])


# ----------------------------------------------------------------------------- type fa21
def _ell_tex(ell):
    """(k - l) as written in an exponent, l in {2, 6, 4}: (k - 2), (k + 2), (k - 4)."""
    return {2: "(k - 2)", 6: "(k + 2)", 4: "(k - 4)"}[ell]


def _idx_tex(ell):
    """l - k inside <>_8: 2 - k, -2 - k, 4 - k."""
    return {2: "2 - k", 6: "-2 - k", 4: "4 - k"}[ell]


def _phase_tex(m, arg, sign="-"):
    t = W.tex_pi_frac(Fr(2 * m, 8))
    return rf"e^{{{sign}j{t}{arg}}}"


def _fa21_formula(kid, m, ell):
    """LaTeX of a candidate general formula for Y[k] (kid: ok, dir, norev, psign, noshift, nomod)."""
    if kid == "ok":
        return rf"$Y[k] = {_phase_tex(m, _ell_tex(ell))}\,X[\langle {_idx_tex(ell)}\rangle_8]$"
    if kid == "dir":
        e2 = (-ell) % 8
        return rf"$Y[k] = {_phase_tex(m, _ell_tex(e2))}\,X[\langle {_idx_tex(e2)}\rangle_8]$"
    if kid == "norev":
        idx = {2: "k - 2", 6: "k + 2", 4: "k - 4"}[ell]
        return rf"$Y[k] = {_phase_tex(m, _ell_tex(ell))}\,X[\langle {idx}\rangle_8]$"
    if kid == "psign":
        return rf"$Y[k] = {_phase_tex(m, _ell_tex(ell), sign='')}\,X[\langle {_idx_tex(ell)}\rangle_8]$"
    if kid == "noshift":
        return rf"$Y[k] = {_phase_tex(m, 'k')}\,X[\langle {_idx_tex(ell)}\rangle_8]$"
    return rf"$Y[k] = {_phase_tex(m, 'k')}\,X[\langle -k\rangle_8]$"


def _fa21(p, c):
    while True:
        base = [random.randint(-2, 4) for _ in range(4)]
        if base[0] == 0 or sum(v == 0 for v in base) > 1 or len(set(base)) < 3:
            continue
        x = base + base
        ell = random.choice([2, 2, 6, 4])
        m = random.choice([v for v in range(1, 8) if v != 4])
        X = np.fft.fft(np.array(x, float))        # exact for this small integer x (Gaussian integers at even k)
        Xe = [complex(round(X[k].real), round(X[k].imag)) for k in range(8)]
        Y = {}
        for k in (0, 2, 4, 6):
            ph = (-1j) ** ((m * (k - ell) // 2) % 4)
            Y[k] = ph * Xe[(ell - k) % 8]
        vals = [Y[k] for k in (0, 2, 4, 6)]
        if sum(abs(v) > 0 for v in vals) >= 3 and len({(v.real, v.imag) for v in vals}) >= 3:
            break
    p["N"], p["Nm1"] = 8, 7
    p["x"] = x
    p["m"], p["ell"] = m, ell
    p["cn_tex"] = CPOW[ell][0]
    p["c_tex"] = {2: "j", 6: "-j", 4: "-1"}[ell]
    p["y_def_tex"] = rf"y[n] = x[\langle {m} - n\rangle_8]\cdot {CPOW[ell][0]}"
    p["ell_signed"] = {2: 2, 6: -2, 4: 4}[ell]
    p["shift_tex"] = {2: r"\langle k - 2\rangle_8", 6: r"\langle k + 2\rangle_8", 4: r"\langle k - 4\rangle_8"}[ell]
    p["cexp_tex"] = {2: r"j^n = e^{j\frac{2\pi}{8}\cdot 2n}", 6: r"(-j)^n = e^{-j\frac{2\pi}{8}\cdot 2n}",
                     4: r"(-1)^n = e^{j\frac{2\pi}{8}\cdot 4n}"}[ell]
    p["Z_tex"] = rf"Z[k] = {_phase_tex(m, 'k')}\,X[\langle -k\rangle_8]"
    p["Xeven_tex"] = ", ".join(rf"X[{k}] = {W.tex_cplx(int(Xe[k].real), int(Xe[k].imag))}" for k in (0, 2, 4, 6))
    p["base_tex"] = W.tex_list([str(v) for v in base])
    rows = []
    for k in (0, 2, 4, 6):
        i = (ell - k) % 8
        ph = (-1j) ** ((m * (k - ell) // 2) % 4)
        rows.append(f"<tr><td>${k}$</td><td>${i}$</td><td>${W.tex_cplx(int(Xe[i].real), int(Xe[i].imag))}$</td>"
                    f"<td>${W.tex_cplx(int(ph.real), int(ph.imag))}$</td>"
                    f"<td>${W.tex_cplx(int(Y[k].real), int(Y[k].imag))}$</td></tr>")
    p["table_html"] = ('<table class="table table-sm table-bordered text-center" style="width:auto">'
                       rf"<tr><th>$k$</th><th>$\langle {_idx_tex(ell)}\rangle_8$</th><th>$X[\langle {_idx_tex(ell)}\rangle_8]$</th>"
                       rf"<th>${_phase_tex(m, _ell_tex(ell))}$</th><th>$Y[k]$</th></tr>" + "".join(rows) + "</table>")
    p["Yrow_tex"] = W.tex_list([W.tex_cplx(int(v.real), int(v.imag)) for v in vals])
    p["Yvals"] = [[int(v.real), int(v.imag)] for v in vals]
    c["Y"] = pl.to_json(np.array([[complex(v) for v in vals]]))
    kids = ["norev", "psign", "nomod"] + (["dir"] if ell != 4 else []) + (["noshift"] if (m * ell) % 8 else [])
    dis = random.sample(kids, min(4, len(kids)))
    texts = {k: _fa21_formula(k, m, ell) for k in ["ok"] + kids}
    opts = mc.choices(texts["ok"], [texts[k] for k in dis], n=5)
    inv = {v: k for k, v in texts.items()}
    p["gen_choices"] = opts
    p["gen_ids"] = [inv[o["text"]] for o in opts]
    p["ok_formula_tex"] = texts["ok"]


def generate(data):
    p, c = data["params"], data["correct_answers"]
    kind = random.choice(KINDS)
    p["kind"] = kind
    for k in ("dup", "inter", "dec", "hw9", "fa21"):
        p["is_" + k] = kind == k
    p["real_kind"] = kind != "fa21"
    if kind == "fa21":
        _fa21(p, c)
    else:
        _real_kind(p, c, kind)
    p["x_tex"] = W.tex_list([fmt.tex_num(v) for v in p["x"]])
    p["X_names_tex"] = W.tex_list([f"X_{{{k}}}" for k in range(p["N"])])
