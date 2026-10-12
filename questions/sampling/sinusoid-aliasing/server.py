"""Sampling sinusoids (Lecture 17: omega = Omega T; a sampled sinusoid knows its frequency only
modulo 2 pi and, for a cosine, up to sign).

(a) x_a(t) = A cos(2 pi f0 t + theta) or A sin(2 pi f0 t + theta) sampled at f_s: write
    x[n] = A cos(omega n + phi) with 0 < omega < pi, phi principal (never +-pi). Folding
    (omega_raw in (pi, 2pi) mod 2pi) flips the sign of the phase; a sine first becomes
    cos(. - pi/2). Regimes: no aliasing (f0 < f_s/2), folding, wrap-around, wrap + fold.
(c) checkbox (FA2023 MT2 #3 / SP2025 MT2 #6(a) shapes): which listed analog frequencies (or which
    sampling periods) give exactly a stated x[n]; a sine (or a phased cosine) accepts only
    omega_raw = +omega0 + 2 pi k, a plain cosine also -omega0 + 2 pi k.
(b) the lowest analog frequency consistent with the samples of (a): omega f_s / (2 pi).
Models: Lecture 17 notes section 1.1 (eqs. 12-14, omega = Omega T, from sampling a sinusoid) and the
study site's Lecture 17 page section 8; SP2025 MT2 #6(a), SP2023 MT2 #8(a), FA2023 MT2 #3, SP2021 MT2
#1(a) and #9(b)."""

import random
from fractions import Fraction as F

from ece310 import fmt

FS = [100, 120, 150, 200, 240, 300, 400, 500, 600, 800, 1000, 1200, 2000]
DENS = [3, 4, 5, 6, 8, 10, 12]
THETAS = [F(1, 6), F(-1, 6), F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(1, 2), F(-1, 2), F(2, 3), F(-2, 3), F(3, 4), F(-3, 4)]


def principal(t):
    t = F(t)
    while t <= -1:
        t += 2
    while t > 1:
        t -= 2
    return t


def tex_pi(q, paren=False):
    """q*pi: \\tfrac{3\\pi}{4}, -\\tfrac{\\pi}{2}, \\pi, 0."""
    q = F(q)
    if q == 0:
        return "0"
    s = "-" if q < 0 else ""
    p, d = abs(q.numerator), q.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    out = s + (num if d == 1 else r"\tfrac{" + num + "}{" + str(d) + "}")
    return rf"\left({out}\right)" if paren and q < 0 else out


def arg_tex(w, phase):
    """'\\tfrac{3\\pi}{4}n + \\tfrac{\\pi}{3}' for omega n + phase (multiples of pi)."""
    s = tex_pi(w) + "n" if w != 1 else r"\pi n"
    if phase:
        s += (" + " if phase > 0 else " - ") + tex_pi(abs(phase))
    return s


def dec(q):
    q = F(q)
    return fmt.plain_num(q) if q.denominator == 1 else f"{float(q):g}"


def reduce_cos(r, th):
    """cos(r pi n + th pi) = cos(w pi n + phi pi), 0 < w < 1. Returns (w, phi, kind of step)."""
    m = r % 2
    if m < 1:
        return m, principal(th), ("none" if r < 1 else "wrap")
    return 2 - m, principal(-th), ("fold" if r < 2 else "wrapfold")


def _part_a():
    while True:
        fs = F(random.choice(FS))
        reg = random.choices(["none", "fold", "wrap", "wrapfold"], weights=[25, 40, 15, 20])[0]
        lo = {"none": 0, "fold": 1, "wrap": 2, "wrapfold": 3}[reg]
        den = random.choice(DENS)
        r = F(lo * den + random.randint(1, den - 1), den)
        f0 = r * fs / 2
        if f0.denominator != 1:
            continue
        kind = random.choice(["cos", "cos", "sin"])
        th = random.choice(THETAS + ([F(0)] if kind == "sin" else []))
        thc = th - F(1, 2) if kind == "sin" else th
        w, phi, step = reduce_cos(r, thc)
        if phi in (F(0), F(1)) or principal(thc) == phi and step in ("fold", "wrapfold"):
            continue                     # phi = 0 or pi: nothing to flip / ambiguous; or a flip that changes nothing
        A = random.randint(1, 5)
        return dict(fs=fs, f0=f0, r=r, kind=kind, th=th, thc=thc, w=w, phi=phi, step=step, A=A)


def _ok_digital(d, w0, kind):
    """Does a digital frequency d*pi reproduce x[n] = A kind(w0 pi n + ph pi)? (sin / phased cos: only
    d = w0 mod 2; plain cos: d = +-w0 mod 2.)"""
    return (d - w0) % 2 == 0 or (kind == "cos" and (d + w0) % 2 == 0)


def _choose(cands, w0, kind, n_opts):
    """cands: {value: digital frequency d (multiple of pi)}. Pick 1-3 correct options and distractors:
    the sign flip (for a sine or a phased cosine), a genuinely different frequency (shifted by pi,
    or doubled), and a multiple of pi (all samples 0 for a sine)."""
    good = [v for v, d in cands.items() if _ok_digital(d, w0, kind)]
    bad = [v for v, d in cands.items() if not _ok_digital(d, w0, kind)]
    flip = [v for v in bad if (cands[v] + w0) % 2 == 0]
    other = [v for v in bad if cands[v] % 1 != 0 and v not in flip]
    halves = [v for v in bad if cands[v] % 1 == 0]
    n_good = random.choice([1, 2, 2, 3])
    if len(good) < n_good:
        return None
    pick = random.sample(good, n_good)
    need = n_opts - n_good
    pool = []
    if flip:
        pool += random.sample(flip, min(len(flip), 2 if kind != "cos" else 0))
    if other:
        pool += random.sample(other, min(len(other), 2))
    if halves:
        pool += random.sample(halves, 1)
    rest = [v for v in bad if v not in pool]
    random.shuffle(rest)
    pool = (pool + rest)[:need]
    if len(pool) < need or not any(v in other for v in pool):
        return None
    return sorted(pick + pool), set(pick)


def _part_b():
    """Checkbox: which analog frequencies (or periods) give exactly x[n] = A kind(w0 pi n + ph pi)."""
    while True:
        shape = random.choice(["freq", "freq", "period"])
        kind = random.choice(["sin", "sin", "cosph", "cos"])
        den = random.choice([3, 4, 5, 6, 8])
        w0 = F(random.randint(1, den - 1), den)             # 0 < w0 < 1
        ph = random.choice([F(1, 4), F(-1, 3), F(1, 6)]) if kind == "cosph" else F(0)
        if shape == "freq":
            fs = F(random.choice(FS))
            base = w0 * fs / 2                               # Hz for w0 pi
            if base.denominator != 1:
                continue
            cands = {}
            for k in range(-2, 3):
                for f in (base + k * fs, -base + k * fs, base + fs / 2 + k * fs, 2 * base + k * fs, k * fs + fs / 2):
                    if f != 0 and abs(f) <= 2 * fs and f.denominator == 1:
                        cands[f] = 2 * f / fs
            got = _choose(cands, w0, kind, 6)
            if got is None:
                continue
            opts, good = got
            return dict(shape="freq", kind=kind, w0=w0, ph=ph, fs=fs, opts=opts, good=good, A=random.randint(1, 4))
        f0 = F(random.choice([10, 20, 25, 40, 50, 60, 75, 100, 150, 200]))
        # 2 pi f0 T = d pi  ->  T = d/(2 f0) for digital frequencies d = +-w0 + 2k, w0 + 1 + 2k, 2 w0, m
        cands = {}
        for k in range(0, 3):
            for d in (w0 + 2 * k, 2 * k - w0, w0 + 1 + 2 * k, 2 * w0 + 2 * k, F(k + 1)):
                T = d / (2 * f0)
                if T > 0 and T.denominator <= 2000:
                    cands[T] = d
        got = _choose(cands, w0, kind, 5)
        if got is None:
            continue
        opts, good = got
        return dict(shape="period", kind=kind, w0=w0, ph=ph, f0=f0, opts=opts, good=good, A=random.randint(1, 4))


def _sig_b(b):
    A = "" if b["A"] == 1 else str(b["A"])
    fn = "sin" if b["kind"] == "sin" else "cos"
    xn = rf"x[n] = {A}\{fn}\!\left({arg_tex(b['w0'], b['ph'])}\right)"
    return fn, A, xn


def generate(data):
    a = _part_a()
    p = data["params"]
    A = "" if a["A"] == 1 else str(a["A"])
    arg = rf"2\pi\cdot {dec(a['f0'])}\,t" + ((" + " if a["th"] > 0 else " - ") + tex_pi(abs(a["th"])) if a["th"] else "")
    p["xa_tex"] = rf"x_a(t) = {A}\{a['kind']}\!\left({arg}\right)"
    p["fs_tex"] = dec(a["fs"])
    p["A_tex"] = A
    p["a"] = {"fs": str(a["fs"]), "f0": str(a["f0"]), "kind": a["kind"], "th": str(a["th"]), "A": a["A"]}
    # solution (a)
    r = a["r"]
    p["raw_tex"] = rf"\omega_0 = \dfrac{{2\pi f_0}}{{f_s}} = \dfrac{{2\pi\cdot {dec(a['f0'])}}}{{{dec(a['fs'])}}} = {tex_pi(r)}"
    p["is_sin"] = a["kind"] == "sin"
    p["sin_tex"] = (rf"{A}\sin\!\left({arg_tex(r, a['th'])}\right) = {A}\cos\!\left({arg_tex(r, a['thc'])}\right)"
                    if a["kind"] == "sin" else "")
    steps = {
        "none": rf"Since $f_0 = {dec(a['f0'])} \lt \tfrac{{f_s}}{{2}} = {dec(a['fs'] / 2)}$ Hz, $\omega_0 = {tex_pi(r)}$ already lies in "
                r"$(0, \pi)$: no aliasing, and the phase is unchanged.",
        "fold": rf"$\omega_0 = {tex_pi(r)}$ lies between $\pi$ and $2\pi$. Subtract $2\pi$: ${tex_pi(r - 2)}$, a <b>negative</b> "
                r"frequency; the cosine is even, $\cos(-\alpha) = \cos\alpha$, so negate the whole argument: the frequency becomes "
                rf"${tex_pi(2 - r)}$ and the phase changes sign.",
        "wrap": rf"$\omega_0 = {tex_pi(r)}$ is more than $2\pi$: subtracting $2\pi$ (one full turn per sample changes nothing) "
                rf"gives ${tex_pi(r - 2)}$, already in $(0, \pi)$; the phase is unchanged.",
        "wrapfold": rf"$\omega_0 = {tex_pi(r)}$: subtracting $4\pi$ gives ${tex_pi(r - 4)}$, negative; since the cosine is even, "
                    rf"negate the argument: frequency ${tex_pi(4 - r)}$, and the phase changes sign.",
    }
    p["step_html"] = steps[a["step"]]
    p["xn_tex"] = rf"x[n] = {A}\cos\!\left({arg_tex(a['w'], a['phi'])}\right)"
    p["w_frac"] = fmt.plain_num(a["w"])
    p["phi_frac"] = fmt.plain_num(a["phi"])
    fmin = a["w"] * a["fs"] / 2
    p["fmin_tex"] = rf"\dfrac{{\omega}}{{2\pi}}f_s = \dfrac{{{tex_pi(a['w'])}}}{{2\pi}}\cdot {dec(a['fs'])} = {dec(fmin)}"
    p["aliased"] = a["step"] != "none"

    b = _part_b()
    fn, Ab, xn = _sig_b(b)
    p["b"] = {"shape": b["shape"], "kind": b["kind"], "w0": str(b["w0"]), "ph": str(b["ph"]), "A": b["A"],
              "fs": str(b.get("fs", "")), "f0": str(b.get("f0", ""))}
    p["xn_b_tex"] = xn
    p["b_freq"] = b["shape"] == "freq"
    p["b_period"] = b["shape"] == "period"
    if b["shape"] == "freq":
        p["b_fs_tex"] = dec(b["fs"])
        p["b_sig_tex"] = rf"{Ab}\{fn}\!\left(2\pi f t" + (("+" if b["ph"] > 0 else "-") + tex_pi(abs(b["ph"])) if b["ph"] else "") + r"\right)"
        p["b_choices"] = [{"text": rf"$f = {dec(f)}$ Hz", "correct": f in b["good"]} for f in b["opts"]]
        rule = (rf"\dfrac{{2\pi f}}{{{dec(b['fs'])}}} = {tex_pi(b['w0'])} + 2\pi k" +
                (r"\ \text{ or }\ -" + tex_pi(b["w0"]) + r" + 2\pi k" if b["kind"] == "cos" else ""))
        p["b_rule_tex"] = rule
        p["b_lines"] = [{"line_html": rf"$f = {dec(f)}$ Hz: $\omega = \tfrac{{2\pi f}}{{f_s}} = {tex_pi(2 * f / b['fs'])}$ &rarr; "
                         + ("gives exactly $x[n]$" if f in b["good"] else _why_not(2 * f / b["fs"], b))} for f in b["opts"]]
    else:
        f0 = b["f0"]
        p["b_sig_tex"] = rf"{Ab}\{fn}\!\left(2\pi\cdot {dec(f0)}\,t" + (("+" if b["ph"] > 0 else "-") + tex_pi(abs(b["ph"])) if b["ph"] else "") + r"\right)"
        p["b_choices"] = [{"text": rf"$T = {_frac_tex(T)}$ s", "correct": T in b["good"]} for T in b["opts"]]
        p["b_rule_tex"] = (rf"2\pi\cdot {dec(f0)}\,T = {tex_pi(b['w0'])} + 2\pi k" +
                           (r"\ \text{ or }\ -" + tex_pi(b["w0"]) + r" + 2\pi k" if b["kind"] == "cos" else ""))
        p["b_lines"] = [{"line_html": rf"$T = {_frac_tex(T)}$ s: $\omega = 2\pi f_0 T = {tex_pi(2 * f0 * T)}$ &rarr; "
                         + ("gives exactly $x[n]$" if T in b["good"] else _why_not(2 * f0 * T, b))} for T in b["opts"]]
    p["b_kind_html"] = {"sin": "a sine is odd, so the sign of the digital frequency matters: only $+\\omega_0 + 2\\pi k$ works",
                        "cosph": "the phase fixes the sign of the frequency: only $+\\omega_0 + 2\\pi k$ works "
                                 "($-\\omega_0$ would flip the phase)",
                        "cos": "a plain cosine is even, so both $\\pm\\omega_0 + 2\\pi k$ work"}[b["kind"]]

    c = data["correct_answers"]
    c["w"] = float(a["w"])
    c["phi"] = float(a["phi"])
    c["fmin"] = float(fmin)


def _frac_tex(T):
    T = F(T)
    return fmt.tex_num(T) if T.denominator != 1 else str(T.numerator)


def _why_not(r, b):
    m = principal(r)
    if m in (0, 1):
        return r"every sample is $\sin(m\pi n) = 0$" if b["kind"] == "sin" else r"a multiple of $\pi$: not the right sequence"
    if m == -b["w0"]:
        return {"sin": r"$\sin(-\omega_0 n) = -\sin(\omega_0 n)$: the wrong sign",
                "cosph": r"the phase would come out with the wrong sign",
                "cos": "fine for a plain cosine"}[b["kind"]]
    if m == r:
        return "a different frequency"
    return rf"$\equiv {tex_pi(m)} \pmod{{2\pi}}$: a different frequency"
