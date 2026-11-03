"""A triangle spectrum through ideal A/D -> ideal digital filter H_d -> ideal D/A, shared period T
(Lectures 17-19; Midterm 2 family "A/D -> H_d -> D/A"). X_a(Omega) is a triangle of height 1 that reaches
0 at |Omega| = 2 pi B (B in Hz). Exact values from ece310.w6_sampling (X_d = (1/T) sum_k X_a((w - 2 pi k)/T)).

(a) for a given B (or T): X_d(0) = 1/T, the band edge w_B = 2 pi B T of the copy centred at 0 (a
    multiple of pi, above pi when the copies overlap) and X_d at one more frequency (inside the overlap
    when there is one, so the copies must be added).
(b) Y_d(w_c) (the cutoff belongs to the passband), the cutoff Omega_c = w_c/T of Y_a in rad/s, Y_a(Omega_c)
    = T Y_d(w_c), and Yes/No: is y_a(t) x_a(t) passed through the ideal analog filter of cutoff w_c/T?
    (Yes iff no copy reaches the passband.)
(c) the largest bandwidth / period. One variant type per past exam problem of this shape:
      fa25  FA2025 MT2 #5 (and SP2025 #5(a) for part (a)-(b)): T and an ideal LPF w_c fixed. Largest B
            with no aliasing in Y_d: 2 pi - 2 pi B T >= w_c, B = (1 - w_c/(2 pi))/T (the key writes
            "B < 50 Hz" and answers 50; the question says to enter the bound), and largest B with no
            aliasing in X_d: B = 1/(2T). FA2025: T = 1/60, w_c = pi/3, B = 20 -> X_d height 60, edge 2pi/3,
            Y_d(pi/3) = 30, Y_a cut at 20 pi rad/s with value 1/2; 50 Hz and 30 Hz. SP2025 #5(a): T = 1/300,
            B = 100 Hz.
      fa19  FA2019 MT2 #9: x_a band-limited to Omega_max = 2 pi B rad/s, LPF w_c fixed, T the unknown:
            the Nyquist period T_NYQ = pi/Omega_max and the largest T for which the overall system is
            still LTI (no aliasing in the passband): T = (2 pi - w_c)/Omega_max. FA2019: Omega_max =
            1000 pi, w_c = pi/3 -> 1 ms and 1/600 s.
      hpf   the high-pass lesson of FA2023 #5(b) and FA2021 #4 on the same setup: an ideal HPF keeps
            w_c <= |w| <= pi, where the copies meet first, so any overlap is aliasing in Y_d: both
            answers are 1/(2T) (the filter buys nothing).
Regimes of the given B (or T), b = 2BT: copies apart (w_c < b pi <= pi), overlap outside the
passband (LPF only, pi < b pi < 2 pi - w_c), overlap inside the passband (b pi > 2 pi - w_c, or > pi
for the HPF); b < 2 always (so X_d(0) = 1/T) and b != 2 - w_c/pi (no copy edge on the cutoff).
About 8 % of variants use the exam's numbers (FA2025 #5, SP2025 #5(a), FA2019 #9 with a new T).
Not duplicated here: the largest T for a cutoff given in Hz (SP2025 #5(b), FA2024 #5(c), FA2023
#5) and the analog cutoff of a given H_d (FA2019 #9(b)) are in sampling/ad-filter-da; the sampled
spectrum alone at T_N and 2T_N (FA2024 #5(b)) is in sampling/sampled-spectrum."""

import random
from fractions import Fraction as F

from ece310 import fmt, mc, svg
from ece310 import w6_sampling as w
from ece310.w7_analog import chain_svg

FS = [30, 40, 48, 50, 60, 72, 80, 90, 100, 120, 150, 160, 180, 200, 240, 300, 360, 400, 480, 500, 600, 800, 900,
      1000, 1200, 1500, 1600, 2000, 2400, 3000]
BS = [50, 100, 150, 200, 250, 300, 400, 500, 600, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
CS = [F(1, 3), F(1, 3), F(1, 3), F(1, 4), F(1, 2), F(2, 3), F(1, 6), F(1, 5), F(2, 5), F(3, 4), F(3, 5)]
CS_HP = [F(1, 3), F(1, 3), F(1, 4), F(1, 2), F(1, 6), F(1, 5), F(2, 5)]
FRACS = sorted({F(k, d) for d in (2, 3, 4, 5, 6, 8, 10, 12) for k in range(1, 2 * d)})
U = sorted({F(k, d) for d in (2, 3, 4, 5, 6, 8, 10, 12) for k in range(1, d + 1)})
NAMED = [("fa25", 60, F(1, 3), 20, "FA2025 Midterm 2 #5"), ("fa25", 300, F(1, 3), 100, "SP2025 Midterm 2 #5(a)"),
         ("fa19", None, F(1, 3), 500, "FA2019 Midterm 2 #9")]


def _nice(q, d=12):
    return F(q).denominator <= d


def _regimes(lp, c):
    """(name, lo, hi, lo_open, hi_open) for b = 2BT."""
    if lp:
        return [("apart", c, F(1), True, False), ("stop", F(1), 2 - c, True, True), ("pass", 2 - c, F(2), True, True)]
    return [("apart", c, F(1), True, False), ("pass", F(1), F(2), True, True)]


def _in(b, lo, hi, lo_open, hi_open):
    return (lo < b if lo_open else lo <= b) and (b < hi if hi_open else b <= hi)


def _xd_nice(q):
    """X_d values a student can write as a short decimal: denominator 1, 2 or 4."""
    return 4 % F(q).denominator == 0


def _choose_u(sp, c, b):
    """A query point: inside the overlap when there is one (X_d is constant there, the sum of two opposite
    ramps), unless the cutoff already lies in it (then Y_d(w_c) asks for that value): then one copy only."""
    X0, Yc = sp.value(0), sp.value(c)
    if b > 1 and not 2 - b < c:
        pool = [u for u in U if 2 - b < u <= 1 and u != c]
        if 1 in pool and random.random() < 0.3 and _xd_nice(sp.value(1)) and sp.value(1) not in (X0, Yc):
            return F(1)
    elif b > 1:
        pool = [u for u in U if u < 2 - b and u != c]
    else:
        inside = [u for u in U if u < b and u != c]
        outside = [u for u in U if b <= u <= 1 and u != c]
        pool = inside if (inside and (random.random() < 0.75 or not outside)) else outside
    X0, Yc = sp.value(0), sp.value(c)
    pool = [u for u in pool if _xd_nice(sp.value(u)) and sp.value(u) not in (X0, Yc)]
    return random.choice(pool) if pool else None


def _attempt(kind, lp, c, reg, fs=None, B=None):
    name, lo, hi, lo_o, hi_o = reg
    bs = [x for x in FRACS if _in(x, lo, hi, lo_o, hi_o) and x != 2 - c and c < x < 2]
    b = random.choice(bs)
    if kind == "fa19":
        B = F(B or random.choice(BS))
        fs = 2 * B / b
        if fs.denominator != 1 or (fs % 10 and fs > 100):
            return None
    else:
        fs = F(fs or random.choice(FS))
        B = b * fs / 2
        if B.denominator != 1:
            return None
    if kind == "fa25" and (((2 - c) * fs / 2).denominator != 1 or (fs / 2).denominator != 1):
        return None
    if kind == "hpf" and (fs / 2).denominator != 1:
        return None
    return fs, B, b


def _pick():
    while True:
        kind = random.choices(["fa25", "fa19", "hpf"], weights=[45, 30, 25])[0]
        lp = kind != "hpf"
        src = None
        fs0 = B0 = None
        if random.random() < 0.05:
            kind, fs0, c, B0, src = random.choice(NAMED)
            lp = True
        else:
            c = random.choice(CS if lp else CS_HP)
        regs = _regimes(lp, c)
        if src and kind == "fa25":
            b = 2 * F(B0) / fs0
            reg = next(r for r in regs if _in(b, *r[1:]))
        else:
            reg = random.choices(regs, weights=[30, 30, 40] if lp else [50, 50])[0]
        for _ in range(300):
            if src and kind == "fa25":
                got = (F(fs0), F(B0), 2 * F(B0) / fs0)
            else:
                got = _attempt(kind, lp, c, reg, fs0, B0)
            if got is None:
                continue
            fs, B, b = got
            sp = w.Spec("tri", 1, B, fs)
            Yd = sp.value(c)
            if not (_xd_nice(Yd) and _nice(Yd / fs, 12)) or Yd == 0:
                continue
            u1 = _choose_u(sp, c, b)
            if u1 is None:
                continue
            return kind, lp, fs, c, B, b, reg[0], sp, u1, src


# ----------------------------------------------------------------------------- figures
def _xa_svg(kind, B):
    if kind == "fa19":
        lab = svg.fmt_pi(2 * B)
        ticks = [(-1, "−" + lab), (1, lab), (0, "0")]
    else:
        ticks = [(-1, "−2πB"), (1, "2πB"), (0, "0")]
    pts = [(-1.45, 0), (-1, 0), (0, 1), (1, 0), (1.45, 0)]
    return svg.plot([pts], xlim=(-1.5, 1.5), ylim=(0, 1.3), xticks=ticks, yticks=[(1, "1")], xlabel="Ω",
                    ylabel="Xₐ(Ω)", width=520, height=200, max_width=480, title="triangle spectrum of x_a(t)")


def _hd_svg(lp, c):
    cf = float(c)
    if lp:
        pts = [(-1, 0), (-cf, 0), (-cf, 1), (cf, 1), (cf, 0), (1, 0)]
    else:
        pts = [(-1, 1), (-cf, 1), (-cf, 0), (cf, 0), (cf, 1), (1, 1)]
    ticks = [(-1, "−π"), (1, "π"), (0, "0"), (cf, svg.fmt_pi(c)), (-cf, svg.fmt_pi(-c))]
    fig = svg.plot([pts], xlim=(-1.15, 1.15), ylim=(0, 1.35), xticks=ticks, yticks=[(1, "1")],
                   xlabel="ω", ylabel="Hd(ω)", width=520, height=200, max_width=480,
                   title="ideal digital filter, one period")
    return fig.replace(">Hd(ω)<", '>H<tspan font-size="10" dy="4">d</tspan><tspan dy="-4">(ω)</tspan><', 1)


def _sol_svg(sp, lp, c):
    xd = [(float(x), float(y)) for x, y in sp.polyline(F(-1), F(1))]
    br = sorted({x for x, _ in sp.polyline(F(-1), F(1))} | {c, -c})
    if lp:
        inner = [(float(x), float(sp.value(x))) for x in br if -c <= x <= c]
        yd = [(-1, 0), (-float(c), 0)] + inner + [(float(c), 0), (1, 0)]
    else:
        left = [(float(x), float(sp.value(x))) for x in br if x <= -c]
        right = [(float(x), float(sp.value(x))) for x in br if x >= c]
        yd = left + [(-float(c), 0), (float(c), 0)] + right
    top = float(sp.H)
    xt = [(-1, "−π"), (1, "π"), (0, "0"), (float(c), svg.fmt_pi(c)), (-float(c), svg.fmt_pi(-c))]
    e = sp.b if sp.b <= 1 else 2 - sp.b
    if e not in (c, 1) and abs(e - c) > F(1, 5):
        xt += [(float(e), svg.fmt_pi(e)), (-float(e), svg.fmt_pi(-e))]
    yv = sp.value(c)
    yt = [(top, w.plain_dec(sp.H))] + ([(float(yv), w.plain_dec(yv))] if abs(float(yv) - top) > 0.12 * top else [])
    return svg.plot([xd, yd], xlim=(-1.15, 1.15), ylim=(0, 1.3 * top), xticks=xt, yticks=yt, xlabel="ω",
                    ylabel="", dashed=[True, False], colors=[svg.ACCENT, svg.ACCENT2], width=560, height=230,
                    max_width=540, labels=[(-0.82, 1.18 * top, "Xd dashed"), (0.82, 1.18 * top, "Yd solid")],
                    title="X_d (dashed) and Y_d (solid) on one period")


# ----------------------------------------------------------------------------- generate
def _hz(q):
    return w.tex_dec(q)


def generate(data):
    kind, lp, fs, c, B, b, reg, sp, u1, src = _pick()
    p = data["params"]
    ans = data["correct_answers"]
    p.update({"kind": kind, "lp": lp, "fs": str(fs), "c": str(c), "B": str(B), "u1": str(u1), "regime": reg,
              "named": src is not None, "src_html": src or ""})
    p["is_fa25"], p["is_fa19"], p["is_hpf"] = kind == "fa25", kind == "fa19", kind == "hpf"
    p["chain_html"] = chain_svg(True)
    p["xa_svg_html"] = _xa_svg(kind, B)
    p["hd_svg_html"] = _hd_svg(lp, c)
    ct = w.tex_pi(c)
    if lp:
        p["hd_tex"] = (rf"H_d(\omega) = \begin{{cases}} 1, & \lvert\omega\rvert \le {ct} \\ 0, & {ct} \lt \lvert\omega\rvert \le \pi"
                       r"\end{cases}")
    else:
        p["hd_tex"] = (rf"H_d(\omega) = \begin{{cases}} 0, & \lvert\omega\rvert \lt {ct} \\ 1, & {ct} \le \lvert\omega\rvert \le \pi"
                       r"\end{cases}")
    p["kind_word"] = "low-pass" if lp else "high-pass"
    p["T_tex"] = sp.tex_T()
    p["B_tex"] = _hz(B)
    p["Omax_tex"] = w.tex_pi(2 * B)
    p["wc_tex"] = ct
    p["u1_tex"] = w.tex_pi(u1)

    # (a)
    p["H_work_tex"] = rf"\tfrac{{1}}{{T}} = {_hz(fs)}"
    p["b_work_tex"] = rf"\omega_B = \Omega_B T = 2\pi\cdot {_hz(B)}\cdot {sp.tex_T()} = {w.tex_pi(b)}"
    p["central_tex"] = w.central_tex(sp)
    p["X1_line"] = w.value_tex(sp, u1, rf"X_d\!\left({w.tex_pi(u1)}\right)")
    p["overlap"] = b > 1
    p["ov_tex"] = w.tex_pi(2 - b)
    ans["Xd0"] = float(sp.H)
    ans["wB"] = float(b)
    ans["Xd1"] = float(sp.value(u1))

    # (b)
    Yd = sp.value(c)
    Ya = Yd / fs
    p["Yd_line"] = w.value_tex(sp, c, rf"Y_d\!\left({ct}\right) = X_d\!\left({ct}\right)")
    p["Omc_tex"] = rf"\Omega_c = \frac{{\omega_c}}{{T}} = {ct}\cdot {_hz(fs)} = {w.tex_pi(c * fs)}"
    p["Ya_tex"] = rf"Y_a(\Omega_c) = T\,Y_d(\omega_c) = \frac{{{w.tex_dec(Yd)}}}{{{_hz(fs)}}} = {w.tex_dec(Ya)}"
    yes = reg in ("apart", "stop")
    p["lti_choices"] = mc.yes_no(yes)
    p["lti_html"] = {
        "apart": rf"The copies do not overlap ($\omega_B = {w.tex_pi(b)} \le \pi$): no aliasing at all, so the chain acts as "
                 rf"$H_a(\Omega) = H_d(\Omega T)$ on the whole band of $x_a$, an ideal analog {p['kind_word']} filter with cutoff "
                 rf"$\tfrac{{\omega_c}}{{T}}$: <b>Yes</b>.",
        "stop": rf"The copies overlap ($\omega_B = {w.tex_pi(b)} \gt \pi$), but the copy centred at $2\pi$ starts at "
                rf"$2\pi - \omega_B = {w.tex_pi(2 - b)} \gt \omega_c = {ct}$: the overlap lies in the stopband, where $H_d$ deletes it. "
                r"$Y_d$ is the central copy cut at $\omega_c$, and $y_a(t)$ is $x_a(t)$ low-pass filtered at $\tfrac{\omega_c}{T}$: "
                r"<b>Yes</b>, although $X_d$ itself is aliased.",
        "pass": (rf"The copy centred at $2\pi$ starts at $2\pi - \omega_B = {w.tex_pi(2 - b)} \lt \omega_c = {ct}$: aliased content "
                 r"lies inside the passband and reaches $y_a(t)$ (it is why $Y_d(\omega_c)$ includes a second copy): <b>No</b>."
                 if lp else
                 rf"The copies overlap ($\omega_B = {w.tex_pi(b)} \gt \pi$) starting at $2\pi - \omega_B = {w.tex_pi(2 - b)} \lt \pi$, "
                 r"inside the high-pass band $\omega_c \le \lvert\omega\rvert \le \pi$; also the D/A rebuilds only "
                 r"$\lvert\Omega\rvert \lt \tfrac{\pi}{T} \lt 2\pi B$: <b>No</b>."),
    }[reg]
    ans["Ydc"] = float(Yd)
    ans["Omc"] = float(c * fs)
    ans["Yac"] = float(Ya)

    # (c)
    if kind == "fa25":
        By, Bx = (2 - c) * fs / 2, fs / 2
        p["c_html"] = (rf"With bandwidth $B$ the copy centred at $0$ ends at $2\pi BT$ and the copy centred at $2\pi$ starts at "
                       rf"$2\pi - 2\pi BT$. <b>No aliasing in $Y_d$</b>: that copy must stay out of the passband "
                       rf"$\lvert\omega\rvert \le {ct}$,"
                       rf"$$2\pi - 2\pi B\cdot {sp.tex_T()} \ge {ct} \iff B \le \frac{{2\pi - {ct}}}{{2\pi T}}"
                       rf" = {_hz(By)}\ \text{{Hz}}.$$ (The key states it as $B \lt {_hz(By)}$ Hz: at the bound the copy reaches the "
                       rf"cutoff with height $0$.) Enter $\mathbf{{{_hz(By)}}}$. <b>No aliasing in $X_d$</b> needs the copies apart: "
                       rf"$2\pi BT \le \pi$, $B \le \tfrac{{1}}{{2T}} = \mathbf{{{_hz(Bx)}}}$ Hz, the Nyquist bound. The filter buys "
                       rf"${_hz(By - Bx)}$ Hz: aliasing is harmless where $H_d$ deletes it.")
    elif kind == "hpf":
        By = Bx = fs / 2
        p["c_html"] = (rf"The copy centred at $2\pi$ starts at $2\pi - 2\pi BT$. The high-pass band $\lvert\omega\rvert \ge {ct}$ reaches "
                       r"all the way to $\pi$, exactly where the copies meet first, so <b>any</b> overlap lands in the passband: no "
                       r"aliasing in $Y_d$ requires $2\pi - 2\pi BT \ge \pi$, i.e. the Nyquist bound "
                       rf"$B \le \tfrac{{1}}{{2T}} = \mathbf{{{_hz(By)}}}$ Hz, the same as for $X_d$ ($2\pi BT \le \pi$): "
                       rf"$\mathbf{{{_hz(Bx)}}}$ Hz. The low-pass trick $B \le \big(1 - \tfrac{{\omega_c}}{{2\pi}}\big)\tfrac{{1}}{{T}} = "
                       rf"{_hz((2 - c) * fs / 2)}$ Hz does not work for a high-pass.")
    else:
        By, Bx = (2 - c) / (2 * B), 1 / (2 * B)
        p["c_html"] = (rf"<b>Nyquist period:</b> $T_{{\text{{NYQ}}}} = \tfrac{{\pi}}{{\Omega_{{\max}}}} = \tfrac{{1}}{{2B}} = "
                       rf"\mathbf{{{fmt.tex_num(Bx)}}}$ s. <b>Largest LTI period:</b> the copy centred at $2\pi$ starts at "
                       rf"$2\pi - \Omega_{{\max}}T$; it must stay out of the passband $\lvert\omega\rvert \le {ct}$:"
                       rf"$$2\pi - {w.tex_pi(2 * B)}\,T \ge {ct} \iff T \le \frac{{2\pi - {ct}}}{{{w.tex_pi(2 * B)}}} = "
                       rf"\mathbf{{{fmt.tex_num(By)}}}\ \text{{s}},$$ "
                       rf"${_hz(2 - c)}$ times the Nyquist period. At that $T$ the chain is still an ideal analog low-pass "
                       rf"filter with cutoff $\tfrac{{\omega_c}}{{T}} = {w.tex_pi(c / By)}$ rad/s.")
    ans["c_y"] = float(By)
    ans["c_x"] = float(Bx)
    p["sol_svg_html"] = _sol_svg(sp, lp, c)
