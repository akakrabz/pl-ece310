"""Two sinusoids through an ideal A/D with period T1 and an ideal D/A with period T2 (Lectures 17-19;
Midterm 2 family "sampling"). Same shape as FA2025 MT2 #4: x_a(t) = A1 f1(Omega_1 t) + A2 f2(Omega_2 t)
(sines, cosines or mixed; integer-Hz frequencies f1 < f2, written in rad/s).

(a) T1 != T2, both tones below f_s/2 (no aliasing): the digital frequencies w_i = Omega_i T1 of x[n]
    (multiples of pi) and the output frequencies Omega_i' = w_i / T2 (rad/s, multiples of pi).
    FA2025 #4(a): sin(20 pi t) + sin(90 pi t), T1 = 1/100, T2 = 1/20 -> sin(4 pi t) + sin(18 pi t).
(b) T1 = T2 = T with 2 f1 < f_s < 2 f2 (below the Nyquist rate of the second tone only): y_a(t) written
    as b1 f1(Omega_1' t) + b2 f2(Omega_2' t), Omega' >= 0; asked: both Omega' and both signed b. The second
    tone folds (sine: sign flip), wraps (no flip) or wraps and folds (flip); a cosine never flips.
    FA2025 #4(b): T = 1/60 -> sin(20 pi t) - sin(30 pi t).
(c) T1 = T2 = T and a stated output y_a(t): every T that produces it; asked: the LARGEST such T and
    how many there are. One variant type per past exam problem of this shape:
      kill   FA2025 #4(c): y_a = the first tone only; the second (a sine) must be sampled at multiples
             of pi (sin(k pi n) = 0), the first must not alias: T = k/(2 f2), k < f2/f1.
             (FA2025: T = k/90, k = 1..4, largest 2/45 s.)
      dc     SP2023 #8(b): x_a = A1 cos(Omega_1 t + phi) + A2 cos(Omega_2 t), y_a = A1 cos(Omega_1 t + phi) + A2:
             the second cosine must land on DC (Omega_2 T = 2 pi k; odd multiples of pi give (-1)^n,
             rebuilt as cos(pi t/T)), the first must not alias: T = k/f2, k < f2/(2 f1).
             (SP2023: 2cos(100 pi t + pi/4) + cos(300 pi t) -> T = 1/150 s only.)
      alias  SP2025 #6(c): sin + sin -> y_a = sin(Omega_3 t) + sin(Omega_1 t): the second sine must alias
             onto +Omega_3 (not -Omega_3: a sine's sign matters), the first must not alias:
             f_s = (f2 - f3)/k > 2 max(f1, f3). (SP2025: sin(30 pi t) + sin(100 pi t) -> T = 1/45 s.)
    Each variant's solution set is computed twice: from the type's formula and by an exact search over
    every candidate rate (f_i +- g)/k that could put a tone on a target frequency, DC or pi; variants
    where the two disagree (e.g. equal amplitudes letting the tones swap roles) are rejected.
About 8 % of variants use the exam's own signal (FA2025 #4 with all its periods; SP2023 #8 and
SP2025 #6 signals with new periods in (a), (b)). Related shapes covered elsewhere: which T / which
Omega_0 give a stated x[n] (SP2025 #6(a), FA2023 #3, SP2023 #8(a)) in sampling/sinusoid-aliasing;
which inputs give a stated output (SP2021 #9) in sampling/ideal-reconstruction."""

import math
import random
from fractions import Fraction as F

from ece310.w7_analog import chain_svg, plain, principal, tex_pi

NICE = [n for n in range(10, 12001) if n % 5 == 0]
F1S = [5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100]


# ----------------------------------------------------------------------------- exact chain model
def image(comp, fs):
    """Output term of one tone A fn(2 pi f t + ph pi) through ideal A/D -> ideal D/A, both at f_s.
    Returns (coef, fn, W, ph) with Omega_out = W pi rad/s, or None if the samples vanish."""
    A, fn, f, ph = comp
    r = principal(2 * F(f) / F(fs))
    if r == 0 or r == 1:
        th = ph - F(1, 2) if fn == "sin" else ph            # A cos(r pi n + th pi)
        val = A * math.cos(math.pi * float(th))
        if abs(val) < 1e-12:
            return None
        return (val, "const", F(0), F(0)) if r == 0 else (val, "cos", F(fs), F(0))   # (-1)^n -> cos(pi t/T)
    if r > 0:
        return (A, fn, r * F(fs), F(ph))
    return (-A if fn == "sin" else A, fn, -r * F(fs), -F(ph))


def chain_out(comps, fs):
    return [t for t in (image(c, fs) for c in comps) if t is not None]


_TS = [0.0137 * k + 0.00191 * k * k for k in range(80)]


def _eval(terms, t):
    s = 0.0
    for c, fn, W, ph in terms:
        if fn == "const":
            s += float(c)
        else:
            arg = math.pi * (float(W) * t + float(ph))
            s += float(c) * (math.sin(arg) if fn == "sin" else math.cos(arg))
    return s


def same(a, b):
    return all(abs(_eval(a, t) - _eval(b, t)) < 1e-9 for t in _TS)


def search(comps, target):
    """Every sampling rate f_s at which A/D -> D/A (same T) maps x_a onto the target exactly."""
    gs = [F(W) / 2 for c, fn, W, ph in target]
    gmax = max(gs)
    fset = set(gs) | {F(0)} | {F(c[2]) for c in comps}
    cands = set()
    for A, fn, f, ph in comps:
        for g in fset:
            for num in (F(f) + g, F(f) - g):
                k = 1
                while num > 0 and num / k > 2 * gmax:
                    cands.add(num / k)
                    k += 1
    return sorted(fs for fs in cands if same(chain_out(comps, fs), target))


# ----------------------------------------------------------------------------- LaTeX
def coef_tex(c, first):
    c = int(c)
    mag = "" if abs(c) == 1 else str(abs(c))
    if first:
        return ("-" if c < 0 else "") + mag
    return (" - " if c < 0 else " + ") + mag


def arg_tex(W, var, ph=F(0)):
    """'20\\pi t', '\\tfrac{\\pi}{5} n', with an optional phase (multiple of pi)."""
    s = f"{tex_pi(W)} {var}"
    if ph:
        s += (" + " if ph > 0 else " - ") + tex_pi(abs(ph))
    return s


def term_tex(c, fn, W, ph, var, first):
    if fn == "const":
        c = int(c)
        return str(c) if first else (" - " if c < 0 else " + ") + str(abs(c))
    return coef_tex(c, first) + rf"\{fn}\!\left({arg_tex(W, var, ph)}\right)"


def sum_tex(terms, var):
    return "".join(term_tex(c, fn, W, ph, var, i == 0) for i, (c, fn, W, ph) in enumerate(terms)) or "0"


def sym_term(c, fn, ph, wname, var, first):
    """A form with a symbolic frequency: 2\\sin(\\omega_1 n + \\tfrac{\\pi}{4})."""
    a = wname + " " + var
    if ph:
        a += (" + " if ph > 0 else " - ") + tex_pi(abs(ph))
    return coef_tex(c, first) + rf"\{fn}\!\left({a}\right)"


def b_sym(bname, fn, ph, wname):
    """b_1\\sin(\\Omega_1' t + phase)."""
    a = wname + " t"
    if ph:
        a += (" + " if ph > 0 else " - ") + tex_pi(abs(ph))
    return rf"{bname}\{fn}\!\left({a}\right)"


def T_tex(fs):
    T = 1 / F(fs)
    return rf"\tfrac{{{T.numerator}}}{{{T.denominator}}}" if T.denominator != 1 else str(T.numerator)


def rad(f):
    return 2 * F(f)


# ----------------------------------------------------------------------------- signals per type
def _kill():
    f1 = random.choice(F1S)
    f2s = [m for m in range(math.floor(F(3, 2) * f1) + 1, math.ceil(F(13, 2) * f1))
           if m % f1 and (m % 5 == 0 or f1 < 15) and m <= 600]
    f2 = random.choice(f2s)
    fn1 = random.choices(["sin", "cos"], weights=[70, 30])[0]
    A1, A2 = random.choice([1, 1, 1, 2, 3]), random.choice([1, 1, 1, 2, 3])
    return [(A1, fn1, F(f1), F(0)), (A2, "sin", F(f2), F(0))]


def _dc():
    f1 = random.choice(F1S)
    f2s = [m for m in range(2 * f1 + 1, 9 * f1) if m % (2 * f1) and (m % f1 == 0 or m % 5 == 0) and m <= 900]
    f2 = random.choice(f2s)
    ph = random.choice([F(0), F(0), F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(1, 6), F(-1, 6)])
    A1, A2 = random.choice([1, 2, 2, 3]), random.choice([1, 1, 2, 3])
    return [(A1, "cos", F(f1), ph), (A2, "cos", F(f2), F(0))]


def _alias():
    while True:
        f1 = random.choice(F1S)
        fs = random.choice([n for n in NICE if 2 * f1 < n <= min(8 * f1, 600)])
        f3 = random.choice([m for m in range(5, (fs + 1) // 2, 5) if m != f1] or [0])
        if not f3 or 2 * f3 == fs:
            continue
        f2 = random.choice([1, 1, 2]) * fs + f3
        if f2 > 1000:
            continue
        fn1 = random.choices(["sin", "cos"], weights=[75, 25])[0]
        A1, A2 = random.choice([1, 1, 1, 2]), random.choice([1, 1, 1, 2])
        return [(A1, fn1, F(f1), F(0)), (A2, "sin", F(f2), F(0))], F(f3)


def _target(kind, comps, f3=None):
    (A1, fn1, f1, ph1), (A2, fn2, f2, ph2) = comps
    if kind == "kill":
        return [(A1, fn1, rad(f1), ph1)]
    if kind == "dc":
        return [(A1, fn1, rad(f1), ph1), (A2, "const", F(0), F(0))]
    return sorted([(A1, fn1, rad(f1), ph1), (A2, "sin", rad(f3), F(0))], key=lambda t: t[2])


def _formula(kind, comps, f3=None):
    """The type's closed-form solution set (sampling rates), or None at a boundary case."""
    f1, f2 = comps[0][2], comps[1][2]
    if kind == "kill":                          # f_s = 2 f2 / k > 2 f1
        if (f2 / f1).denominator == 1:
            return None
        return sorted(2 * f2 / k for k in range(1, math.floor(f2 / f1) + 1))
    if kind == "dc":                            # f_s = f2 / k > 2 f1
        if (f2 / (2 * f1)).denominator == 1:
            return None
        return sorted(f2 / k for k in range(1, math.floor(f2 / (2 * f1)) + 1))
    lo = 2 * max(f1, f3)                        # f_s = (f2 - f3)/k > 2 max(f1, f3)
    out = []
    k = 1
    while (f2 - f3) / k >= lo:
        if (f2 - f3) / k == lo:
            return None
        out.append((f2 - f3) / k)
        k += 1
    return sorted(out)


NAMED = {
    "kill": ([(1, "sin", F(10), F(0)), (1, "sin", F(45), F(0))], None, (100, 20, 60), "FA2025 Midterm 2 #4"),
    "dc": ([(2, "cos", F(50), F(1, 4)), (1, "cos", F(150), F(0))], None, None, "SP2023 Midterm 2 #8"),
    "alias": ([(1, "sin", F(15), F(0)), (1, "sin", F(50), F(0))], F(5), None, "SP2025 Midterm 2 #6(c)"),
}


def _den_ok(q, d=20):
    return F(q).denominator <= d


def _part_a(comps):
    f1, f2 = comps[0][2], comps[1][2]
    while True:
        c1 = [n for n in NICE if 2 * f2 < n <= 6 * f2 and _den_ok(2 * f1 / n) and _den_ok(2 * f2 / n)]
        if not c1:
            return None
        fs1 = random.choice(c1)
        r1, r2 = 2 * f1 / fs1, 2 * f2 / fs1
        c2 = [n for n in NICE if n != fs1 and (r1 * n).denominator == 1 and (r2 * n).denominator == 1 and n <= 2 * fs1
              and (n % 10 == 0 or n < 100)]
        slow = [n for n in c2 if n < fs1]
        fast = [n for n in c2 if n > fs1]
        pool = slow if slow and (random.random() < 0.7 or not fast) else fast
        if pool:
            return fs1, random.choice(pool)


def _part_b(comps, sols, fs1):
    f1, f2 = comps[0][2], comps[1][2]
    by_reg = {}
    for n in NICE:
        if not (2 * f1 < n < 2 * f2) or F(n) in sols or n == fs1:
            continue
        r = 2 * f2 / n
        p = principal(r)
        if p in (0, 1) or abs(p) * n / 2 == f1 or not (_den_ok(2 * f1 / n, 24) and _den_ok(r, 24)):
            continue
        reg = "fold" if r < 2 else "wrap" if r < 3 else "wrapfold" if r < 4 else "far"
        by_reg.setdefault(reg, []).append(n)
    regs = [g for g in ("fold", "wrap", "wrapfold") if g in by_reg]
    if not regs:
        return None
    reg = random.choices(regs, weights=[{"fold": 55, "wrap": 25, "wrapfold": 20}[g] for g in regs])[0]
    return random.choice(by_reg[reg])


def _variant():
    while True:
        kind = random.choices(["kill", "dc", "alias"], weights=[45, 30, 25])[0]
        named = random.random() < 0.08
        f3, fixed, src = None, None, None
        if named:
            comps, f3, fixed, src = NAMED[kind]
        elif kind == "kill":
            comps = _kill()
        elif kind == "dc":
            comps = _dc()
        else:
            comps, f3 = _alias()
        target = _target(kind, comps, f3)
        sols = search(comps, target)
        if not sols or sols != _formula(kind, comps, f3):
            continue
        if fixed:
            fs1, fs2, fsb = (F(v) for v in fixed)
        else:
            pa = _part_a(comps)
            if pa is None:
                continue
            fs1, fs2 = F(pa[0]), F(pa[1])
            fsb = _part_b(comps, sols, fs1)
            if fsb is None:
                continue
            fsb = F(fsb)
        return kind, comps, f3, target, sols, fs1, fs2, fsb, src


# ----------------------------------------------------------------------------- solution text
def _reduce_html(comp, fs, idx):
    """How tone idx lands after sampling at f_s (part b)."""
    A, fn, f, ph = comp
    r = 2 * f / fs
    p = principal(r)
    w = rf"\omega_{idx} = \Omega_{idx} T = {tex_pi(rad(f))}\cdot{T_tex(fs)} = {tex_pi(r)}"
    if r < 1:
        return rf"${w} \lt \pi$: no aliasing, the term passes unchanged."
    shift = r - p
    s = rf"${w} \gt \pi$: aliased. Subtract ${tex_pi(shift)}$ (whole turns change no sample): ${tex_pi(p)}$"
    if p > 0:
        return s + ", already in $(0, \\pi)$, so no sign change."
    if fn == "sin":
        return s + (r", a <b>negative</b> frequency, and $\sin(-\theta) = -\sin\theta$: the term changes sign and comes back at "
                    rf"${tex_pi(-p)}$.")
    return s + (r", a negative frequency, but $\cos(-\theta) = \cos\theta$: a cosine keeps its sign"
                + (" (its phase changes sign)" if ph else "") + rf" and comes back at ${tex_pi(-p)}$.")


def _c_html(kind, comps, f3, sols):
    (A1, fn1, f1, ph1), (A2, fn2, f2, ph2) = comps
    Ts = [1 / s for s in sols]
    lst = ",\\ ".join(T_tex(1 / T) for T in sorted(Ts))
    K = len(sols)
    if kind == "kill":
        nxt = 2 * f2 / (K + 1)
        r = 2 * f1 / nxt
        return (rf"Only the first term may come out, unchanged, so (i) it must not alias: $\Omega_1 T = {tex_pi(rad(f1))}\,T \lt \pi$, "
                rf"i.e. $T \lt \tfrac{{1}}{{{2 * f1}}}$ s; and (ii) the second term must vanish. A sine sampled at a multiple of "
                rf"$\pi$ gives $\sin(k\pi n) = 0$ for every $n$ (a different frequency would leave a tone behind), so "
                rf"$\Omega_2 T = {tex_pi(rad(f2))}\,T = k\pi$, $T = \tfrac{{k}}{{{2 * f2}}}$. Both: "
                rf"$\tfrac{{k}}{{{2 * f2}}} \lt \tfrac{{1}}{{{2 * f1}}} \iff k \lt \tfrac{{{2 * f2}}}{{{2 * f1}}} = "
                rf"{_num(f2 / f1)}$, so {_klist(K)}: "
                rf"$$T \in \left\{{{lst}\right\}}\ \text{{s}}.$$ The largest is $\mathbf{{{T_tex(sols[0])}}}$ s and there are "
                rf"$\mathbf{{{K}}}$. The next one, $k = {K + 1}$ ($T = {T_tex(nxt)}$ s), also kills the second term, but then "
                rf"$\Omega_1 T = {tex_pi(r)} \ge \pi$ aliases the first.")
    if kind == "dc":
        nxt = f2 / (K + 1)
        return (rf"The first term must come out unchanged: $\Omega_1 T = {tex_pi(rad(f1))}\,T \lt \pi$, i.e. "
                rf"$T \lt \tfrac{{1}}{{{2 * f1}}}$ s. The constant ${A2}$ can only come from the second cosine landing on "
                rf"$\omega = 0$: $\Omega_2 T = {tex_pi(rad(f2))}\,T = 2\pi k$, $T = \tfrac{{k}}{{{f2}}}$, and then "
                rf"$\cos(2\pi k n) = 1$ (an odd multiple of $\pi$ would give ${A2}(-1)^n$ instead, which the D/A rebuilds as a cosine "
                rf"at $\tfrac{{\pi}}{{T}}$, not a constant). Both: $\tfrac{{k}}{{{f2}}} \lt \tfrac{{1}}{{{2 * f1}}} \iff "
                rf"k \lt \tfrac{{{f2}}}{{{2 * f1}}} = {_num(f2 / (2 * f1))}$, so {_klist(K)}: "
                rf"$$T \in \left\{{{lst}\right\}}\ \text{{s}}.$$ The largest is $\mathbf{{{T_tex(sols[0])}}}$ s and there "
                + (r"is only $\mathbf{1}$" if K == 1 else rf"are $\mathbf{{{K}}}$")
                + rf". $T = {T_tex(nxt)}$ s also sends the second cosine to DC but aliases the first.")
    lo = 2 * max(f1, f3)
    return (rf"Both output frequencies must lie below $\tfrac{{f_s}}{{2}}$, so $f_s \gt {_num(lo)}$ Hz, and the first term then "
            rf"passes unchanged. The term ${term_tex(A2, 'sin', rad(f3), F(0), 't', True)}$ must be the alias of "
            rf"${term_tex(A2, 'sin', rad(f2), F(0), 't', True)}$ with a <b>plus</b> sign: "
            rf"${tex_pi(rad(f2))}\,T = {tex_pi(rad(f3))}\,T + 2\pi k$, i.e. $f_s = \tfrac{{{_num(f2)} - {_num(f3)}}}{{k}} = "
            rf"\tfrac{{{_num(f2 - f3)}}}{{k}}$ Hz. (The folded solutions $f_s = \tfrac{{{_num(f2 + f3)}}}{{k}}$ give "
            rf"$\omega_2 \equiv -{tex_pi(2 * f3)}\,T$ and $\sin(-\theta) = -\sin\theta$: the wrong sign.) With $f_s \gt {_num(lo)}$: "
            + (rf"only $k = 1$, $f_s = {_num(sols[0])}$ Hz," if K == 1 else rf"$k = 1, \dots, {K}$,")
            + rf" $$T \in \left\{{{lst}\right\}}\ \text{{s}}.$$ The largest is $\mathbf{{{T_tex(sols[0])}}}$ s and there "
            + (r"is only $\mathbf{1}$." if K == 1 else rf"are $\mathbf{{{K}}}$."))


def _klist(K):
    return {1: "only $k = 1$", 2: "$k = 1, 2$"}.get(K, rf"$k = 1, \dots, {K}$")


def _num(q):
    q = F(q)
    return str(q.numerator) if q.denominator == 1 else rf"\tfrac{{{q.numerator}}}{{{q.denominator}}}"


# ----------------------------------------------------------------------------- generate
def generate(data):
    kind, comps, f3, target, sols, fs1, fs2, fsb, src = _variant()
    (A1, fn1, f1, ph1), (A2, fn2, f2, ph2) = comps
    p = data["params"]
    c = data["correct_answers"]
    p["kind"] = kind
    p["comps"] = [{"A": A, "fn": fn, "f": plain(f), "ph": plain(ph)} for A, fn, f, ph in comps]
    p["target"] = [{"c": int(cc), "fn": fn, "W": plain(W), "ph": plain(ph)} for cc, fn, W, ph in target]
    p["fs"] = {"a1": plain(fs1), "a2": plain(fs2), "b": plain(fsb)}
    p["src_html"] = src or ""
    p["named"] = src is not None
    p["chain_html"] = chain_svg(False, title="x_a(t) into an ideal A/D with period T1, x[n] into an ideal D/A with period T2, "
                                             "output y_a(t)").replace(
        'font-style="italic">T</text>', 'font-style="italic">T<tspan font-size="10" dy="4">1</tspan></text>', 1).replace(
        'font-style="italic">T</text>', 'font-style="italic">T<tspan font-size="10" dy="4">2</tspan></text>', 1)
    W1, W2 = rad(f1), rad(f2)
    xa = [(A1, fn1, W1, ph1), (A2, fn2, W2, ph2)]
    p["xa_tex"] = "x_a(t) = " + sum_tex(xa, "t")
    p["comp1_tex"] = term_tex(A1, fn1, W1, ph1, "t", True)
    p["comp2_tex"] = term_tex(A2, fn2, W2, ph2, "t", True)
    p["T1_tex"], p["T2_tex"], p["Tb_tex"] = T_tex(fs1), T_tex(fs2), T_tex(fsb)

    # (a)
    r1, r2 = W1 / fs1, W2 / fs1
    Wa1, Wa2 = r1 * fs2, r2 * fs2
    p["ax_form_tex"] = "x[n] = " + sym_term(A1, fn1, ph1, r"\omega_1", "n", True) + sym_term(A2, fn2, ph2, r"\omega_2", "n", False)
    p["ay_form_tex"] = "y_a(t) = " + sym_term(A1, fn1, ph1, r"\Omega_1'", "t", True) + sym_term(A2, fn2, ph2, r"\Omega_2'", "t", False)
    p["a_w_tex"] = (rf"\omega_1 = \Omega_1 T_1 = {tex_pi(W1)}\cdot{T_tex(fs1)} = {tex_pi(r1)},\qquad "
                    rf"\omega_2 = \Omega_2 T_1 = {tex_pi(W2)}\cdot{T_tex(fs1)} = {tex_pi(r2)}")
    p["a_x_tex"] = "x[n] = " + sum_tex([(A1, fn1, r1, ph1), (A2, fn2, r2, ph2)], "n")
    p["a_W_tex"] = (rf"\Omega_1' = \frac{{\omega_1}}{{T_2}} = {tex_pi(r1)}\cdot {_num(fs2)} = {tex_pi(Wa1)},\qquad "
                    rf"\Omega_2' = \frac{{\omega_2}}{{T_2}} = {tex_pi(r2)}\cdot {_num(fs2)} = {tex_pi(Wa2)}")
    p["a_y_tex"] = "y_a(t) = " + sum_tex([(A1, fn1, Wa1, ph1), (A2, fn2, Wa2, ph2)], "t")
    p["a_ratio_html"] = (rf"The samples taken at {_num(fs1)} Hz are played back at {_num(fs2)} Hz, so every frequency is "
                         + (rf"divided by ${_num(fs1 / fs2)}$" if fs1 > fs2 else rf"multiplied by ${_num(fs2 / fs1)}$")
                         + r"; using $T_1$ in the D/A would just give $x_a(t)$ back.")
    c["w1"], c["w2"] = float(r1), float(r2)
    c["Wa1"], c["Wa2"] = float(Wa1), float(Wa2)

    # (b)
    out = [image(cp, fsb) for cp in comps]
    p["b_form_tex"] = "y_a(t) = " + b_sym("b_1", fn1, ph1, r"\Omega_1'") + " + " + b_sym("b_2", fn2, ph2, r"\Omega_2'")
    p["b1_html"] = _reduce_html(comps[0], fsb, 1)
    p["b2_html"] = _reduce_html(comps[1], fsb, 2)
    p["b_x_tex"] = "x[n] = " + sum_tex([(o[0], o[1], o[2] / fsb, o[3]) for o in out], "n")
    p["b_y_tex"] = "y_a(t) = " + sum_tex(out, "t")
    p["b_W_tex"] = rf"\Omega' = \frac{{\omega}}{{T}} = {_num(fsb)}\,\omega"
    p["b_flip"] = out[1][0] < 0
    c["Wb1"], c["b1"] = float(out[0][2]), float(out[0][0])
    c["Wb2"], c["b2"] = float(out[1][2]), float(out[1][0])

    # (c)
    p["yc_tex"] = "y_a(t) = " + sum_tex(target, "t")
    p["c_html"] = _c_html(kind, comps, f3, sols)
    p["sols"] = [plain(1 / s) for s in sols]
    p["is_kill"], p["is_dc"], p["is_alias"] = kind == "kill", kind == "dc", kind == "alias"
    c["Tmax"] = float(1 / sols[0])
    c["count"] = len(sols)
