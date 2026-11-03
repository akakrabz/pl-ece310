"""The DTFT of a finite-length cosine: two periodic sincs (Lecture 22 Sec. 1.1, eqs. 1-10).

x[n] = A cos(w0 n), 0 <= n <= N-1, has X_d(w) = (A/2) C(w - w0) + (A/2) C(w + w0) with
C(w) = e^{-j w (N-1)/2} sin(N w / 2) / sin(w / 2): a main lobe of height A N / 2 at w0, nulls at
w0 +- 2 pi / N (width 4 pi / N), side lobes between the nulls every 2 pi / N.

Two variant types (Midterm 2 family "DFT spectral analysis and resolution"):
  read     Lecture 22 slide Example 1 (centre pi/3, nulls 5pi/24 and 11pi/24, height 64 -> w0 = pi/3,
           N = 16, A = 8), in the wording of SP2021 MT2 #7 ("The DTFT ... is as shown below. Determine
           the constants"): |X_d(w)| on [0, pi] is plotted (dense polyline, ece310.svg.plot) with
           tick labels at the two main-lobe nulls and the peak height; asked w0 (multiple of pi), N, A.
  predict  the forward direction, FA2025 MT2 #8 (four DFT plots of the same two cosines: 16 samples,
           zero-padded to 32, 32 samples, Hamming -- which change makes the lobes taller/narrower) and
           FA2019 MT2 #1(c) (a finite cosine has lobes, not impulses): A, N, w0 given; asked the
           main-lobe height, the two nulls and which change halves the main-lobe width (MC: N' = 2N
           samples, against zero padding, a Hamming window, N/2 or 4N samples, a larger amplitude).

Variants: w0 = pi k / N with k an integer, so X_d(w0) = A N / 2 exactly (w0 is a null of the mirror
lobe C(w + w0)) and the labelled nulls w0 +- 2pi/N are exact zeros of X_d (2 w0 +- 2pi/N is a multiple
of 2pi/N); the mirror lobe moves the plotted maximum by less than 2 % ("the peak height is AN/2 to
the accuracy of the plot"). w0 stays at least 4pi/N away from 0 and pi, so the mirror lobes at -w0 and
2pi - w0 keep their main lobes off the plotted range."""

import cmath
import math
import random
from fractions import Fraction as Fr

from ece310 import fmt, mc, svg, w8_dft as W

RS = [Fr(a, q) for q in (2, 3, 4, 5, 6, 7, 8, 9, 10, 12) for a in range(1, q) if math.gcd(a, q) == 1]
WIDTH, LEFT, RIGHT = 640, 48, 34
PX_PER_PI = WIDTH - LEFT - RIGHT


def _C(w, N):
    """Periodic sinc C(w) = e^{-jw(N-1)/2} sin(Nw/2)/sin(w/2) (limit N cos(Nw/2)/cos(w/2) at its peaks)."""
    s = math.sin(w / 2)
    ratio = N * math.cos(N * w / 2) / math.cos(w / 2) if abs(s) < 1e-12 else math.sin(N * w / 2) / s
    return cmath.exp(-0.5j * w * (N - 1)) * ratio


def _Xd(w, A, w0, N):
    return abs(A / 2 * (_C(w - w0, N) + _C(w + w0, N)))


def _label_px(lab):
    return 7.5 * len(lab)


def _candidates():
    out = []
    for r in RS:
        for N in range(8, 33):
            k = r * N
            if k.denominator != 1:
                continue
            lo, hi = r - Fr(2, N), r + Fr(2, N)
            if r < Fr(4, N) or r > 1 - Fr(4, N):
                continue
            ticks = [(Fr(0), "0"), (lo, svg.fmt_pi(lo)), (hi, svg.fmt_pi(hi)), (Fr(1), "π")]
            ok = True
            for (x0, a), (x1, b) in zip(ticks, ticks[1:]):
                if float(x1 - x0) * PX_PER_PI < (_label_px(a) + _label_px(b)) / 2 + 10:
                    ok = False
            if ok and _mirror_ok(r, N):
                out.append((r, N))
    return out


def _shape(A, r, N):
    """Peak of |X_d| over the main lobe (dense search) and its location, as multiples of pi."""
    w0 = math.pi * float(r)
    best, at = 0.0, w0
    for i in range(-400, 401):
        w = w0 + i * (2 * math.pi / N) / 400
        v = _Xd(w, A, w0, N)
        if v > best:
            best, at = v, w
    return best, at / math.pi


def _figure(p, A, r, N):
    w0 = math.pi * float(r)
    H = A * N // 2
    lo, hi = r - Fr(2, N), r + Fr(2, N)
    xs = sorted({i / 720 for i in range(721)} | {float(lo), float(hi), float(r)})
    pts = [(x, _Xd(math.pi * x, A, w0, N)) for x in xs]
    ymax = round(1.2 * H, 3)
    p["ymax"] = ymax
    p["fig_html"] = svg.plot([pts], xlim=(0, 1), ylim=(0, ymax),
                             xticks=[(0, "0"), (float(lo), svg.fmt_pi(lo)), (float(hi), svg.fmt_pi(hi)), (1, "π")],
                             yticks=[(H, str(H))], guides=[((0, H), (float(r), H))],
                             xlabel="ω", ylabel="|Xd(ω)|", width=WIDTH, max_width=WIDTH, height=260,
                             title="Magnitude of the DTFT of a finite-length cosine, 0 ≤ ω ≤ π")


def _mc_width(N, A):
    correct = rf"Keep $N' = {2 * N}$ samples of the cosine instead of ${N}$."
    wrong = [
        rf"Zero-pad the ${N}$ samples to length ${2 * N}$.",
        rf"Keep $N' = {4 * N}$ samples of the cosine instead of ${N}$.",
        rf"Keep $N' = {N // 2}$ samples of the cosine instead of ${N}$." if N % 2 == 0 else
        rf"Zero-pad the ${N}$ samples to length ${4 * N}$.",
        rf"Multiply the ${N}$ samples by a Hamming window.",
        rf"Double the amplitude to ${2 * A}$, keeping the ${N}$ samples.",
    ]
    return mc.choices(correct, random.sample(wrong, len(wrong)), n=4)


def _mirror_ok(r, N):
    peak, at = _shape(2, r, N)
    return abs(peak - N) < 0.02 * N and abs(at - float(r)) < 0.02 * 4 / N


CANDS = _candidates()


def generate(data):
    p, c = data["params"], data["correct_answers"]
    while True:
        r = random.choice(sorted({rr for rr, _ in CANDS}))
        N = random.choice([n for rr, n in CANDS if rr == r])
        A = random.randint(2, 12)
        if (A * N) % 2 == 0:
            break
    k = int(r * N)
    H = A * N // 2
    lo, hi = r - Fr(2, N), r + Fr(2, N)
    peak, at = _shape(A, r, N)
    if not (abs(peak - H) < 0.02 * H and abs(at - float(r)) < 0.02 * 4 / N):
        raise RuntimeError("mirror lobe moves the peak by more than 2 %")   # excluded by w0 = pi k / N
    read = random.random() < 0.55
    p["read"], p["predict"] = read, not read
    p["N"], p["A"], p["H"], p["k"] = N, A, H, k
    p["Nm1"], p["N2"] = N - 1, 2 * N
    p["w0_tex"] = W.tex_pi_frac(r)
    p["w0n_tex"] = (W.tex_pi_frac(r) + "n") if r.denominator != 1 else W.tex_pi_frac(r) + " n"
    p["lo_tex"], p["hi_tex"] = W.tex_pi_frac(lo), W.tex_pi_frac(hi)
    p["width_tex"] = W.tex_pi_frac(hi - lo)
    p["half_tex"] = W.tex_pi_frac(Fr(2, N))
    p["lo_plain"], p["hi_plain"] = svg.fmt_pi(lo), svg.fmt_pi(hi)
    p["sum_tex"] = W.tex_pi_frac(lo + hi)
    p["width_over_pi_tex"] = fmt.tex_num(hi - lo)
    p["peak_exact"] = round(peak, 2)
    p["peak_pct"] = f"{100 * abs(peak - H) / H:.2f}"
    p["r"] = [r.numerator, r.denominator]
    p["lohi"] = [[lo.numerator, lo.denominator], [hi.numerator, hi.denominator]]
    if read:
        _figure(p, A, r, N)
        c["w0"] = float(r)
        c["N"] = N
        c["A"] = float(A)
    else:
        p["width_choices"] = _mc_width(N, A)
        c["Hp"] = float(H)
        c["lo"] = float(lo)
        c["hi"] = float(hi)
