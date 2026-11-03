"""Two- and three-tap FIRs, causal or not: |H_d|, principal phase, number of +-pi jumps, group delay
(Lectures 15-16). One variant type per past Midterm 2 problem of this shape (new numbers, same
statement style):

  fa25   non-causal LCCDE  y[n] = x[n+a] +- x[n-b]            FA2025 #3: y[n] = x[n+3] - x[n]
  sp25   LCCDE  y[n] = x[n] + c x[n-d] + x[n-2d]              SP2025 #3: x[n] + 2x[n-3] + x[n-6]
  fa24   h[n] = delta[n+d] + c delta[n] + delta[n-d]          FA2024 #3: delta[n+2] + delta[n] + delta[n-2]
  sp23   h[n] = delta[n] +- delta[n-N]                         SP2023 #5: delta[n] + delta[n-6]
  sp21b  h[n] = g delta[n+a] -+ g delta[n-b], g in {2, 3}     SP2021 #3: 3 delta[n+1] - 3 delta[n-7]
  sp21a  H_d(w) = g e^{-jMw} sin(Kw) given in closed form     SP2021 #2: e^{-j4w} sin 2w

Model: H_d(w) = (j if J) e^{-j m w} A(w) with A real: 2g cos(dw) or 2g sin(dw) for two taps at
n1 < n2 (m = (n1 + n2)/2, d = (n2 - n1)/2; equal signs give cos, opposite signs give j sin),
c + 2cos(dw) for three symmetric taps, g sin(Kw) for the closed form.
Asked:
  (a) |H_d(w0)| exactly (rational or one surd);
  (b) the principal phase (multiple of pi) at two frequencies with A != 0 and answer != +-pi, at least
      one with A < 0 (the +-pi correction, whenever A takes negative values) and one whose angle
      -m w (+ pi/2) (+ pi if A < 0) leaves (-pi, pi] (a 2 pi wrap; 70 %: at a point with A > 0);
      fa24 has m = 0 and a real H_d (FA2024 #3(a)): the real values H_d(w1), H_d(w2) instead, one < 0;
  (c) the number of +-pi jumps in (-pi, pi): sign changes of A only (a zero that A touches without
      changing sign, e.g. 2 + 2cos 3w, is no jump; 2 pi wraps are no jumps);
  (d) the group delay m in samples (negative = advance)."""

import math
import random
from fractions import Fraction as F

import numpy as np
import sympy as sp
from ece310 import fmt
from ece310.w5_freqresp import dec, plain, principal, rat, tex_pi, tex_real

PH_FREQS = [F(k, 12) for k in (2, 3, 4, 6, 8, 9, 10)]          # pi/6 ... 5pi/6
MAG_FREQS = [F(0)] + PH_FREQS + [F(1)]
TYPES = [("fa25", 24), ("sp25", 16), ("fa24", 14), ("sp23", 16), ("sp21b", 16), ("sp21a", 14)]
GRID_N = 4000                       # grid t_i = -1 + (2i+1)/N never hits a zero (odd/4000)

SOURCE = {
    "fa25": "FA2025 Midterm 2 #3 ($y[n] = x[n+3] - x[n]$)",
    "sp25": "SP2025 Midterm 2 #3 ($y[n] = x[n] + 2x[n-3] + x[n-6]$)",
    "fa24": r"FA2024 Midterm 2 #3 ($h[n] = \delta[n+2] + \delta[n] + \delta[n-2]$)",
    "sp23": r"SP2023 Midterm 2 #5 ($h[n] = \delta[n] + \delta[n-6]$)",
    "sp21b": r"SP2021 Midterm 2 #3 ($3\delta[n+1] - 3\delta[n-7]$)",
    "sp21a": r"SP2021 Midterm 2 #2 ($e^{-j4\omega}\sin 2\omega$)",
}


# ----------------------------------------------------------------------------- LaTeX pieces
def kw_tex(k):
    """k*omega as a trig argument: \\omega, 3\\omega, \\tfrac{3\\omega}{2}."""
    k = F(k)
    if k == 1:
        return r"\omega"
    if k.denominator == 1:
        return f"{k.numerator}\\omega"
    return r"\tfrac{" + ("" if k.numerator == 1 else str(k.numerator)) + r"\omega}{" + str(k.denominator) + "}"


def trig_tex(fn, k):
    a = kw_tex(k)
    return rf"\{fn}" + (" " if not a.startswith("\\") else "") + a


def ej_tex(x):
    """e^{-j x omega} (x a Fraction): '' for 0, e^{-j\\omega}, e^{j3\\omega}, e^{-j\\frac{3}{2}\\omega}."""
    x = F(x)
    if x == 0:
        return ""
    s = "-" if x > 0 else ""
    a = abs(x)
    if a == 1:
        body = r"\omega"
    elif a.denominator == 1:
        body = f"{a.numerator}\\omega"
    else:
        body = r"\frac{" + str(a.numerator) + "}{" + str(a.denominator) + r"}\omega"
    return f"e^{{{s}j{body}}}"


def ramp_tex(m):
    """-m*omega in LaTeX: -3\\omega, \\tfrac{3}{2}\\omega, \\omega."""
    m = F(m)
    a = abs(m)
    body = r"\omega" if a == 1 else (f"{a.numerator}\\omega" if a.denominator == 1 else
                                     r"\tfrac{" + str(a.numerator) + "}{" + str(a.denominator) + r"}\omega")
    return ("-" if m > 0 else "") + body


def x_tex(s):
    """x[n+s] in LaTeX."""
    return "x[n]" if s == 0 else (f"x[n+{s}]" if s > 0 else f"x[n-{-s}]")


def d_tex(n):
    """delta[n - n0] in LaTeX."""
    return r"\delta[n]" if n == 0 else (rf"\delta[n-{n}]" if n > 0 else rf"\delta[n+{-n}]")


def simple(x):
    """Exam-nice exact value: rational or a rational multiple of one square root (sqrt 2, sqrt 3)."""
    x = sp.nsimplify(x)
    if x.is_rational:
        return True
    return bool(sp.nsimplify(sp.expand(x * x)).is_rational)


def nice(x):
    """Exact and printable: a sum of rationals and rational multiples of sqrt 2, sqrt 3, sqrt 6."""
    x = sp.nsimplify(x)
    for t in sp.Add.make_args(sp.expand(x)):
        if not (t.is_rational or any((t / sp.sqrt(d)).is_rational for d in (2, 3, 6))):
            return False
    return True


# ----------------------------------------------------------------------------- the system
class System:
    def __init__(self, typ):
        self.typ = typ
        self.taps = None             # {n: coefficient} for the FIR types
        self.closed = None           # (g, M, K) for sp21a

    # A(w) = sum coef * fn(k w)
    def A_exact(self, t):
        f = sp.cos if self.fn == "cos" else sp.sin
        return sp.nsimplify(sum((rat(c) * (f(sp.pi * rat(k) * rat(t)) if k != 0 else 1) for c, k in self.terms),
                                sp.Integer(0)))

    def A_num(self, w):
        f = np.cos if self.fn == "cos" else np.sin
        w = np.asarray(w, float)
        return sum(float(c) * (f(float(k) * w) if k != 0 else 1.0) for c, k in self.terms)


def _two_tap(sysm, n1, g1, n2, g2):
    sysm.m = F(n1 + n2, 2)
    d = F(n2 - n1, 2)
    if g2 == g1:
        sysm.J, sysm.fn, sysm.terms = F(0), "cos", [(F(2 * g1), d)]
    else:
        sysm.J, sysm.fn, sysm.terms = F(1, 2), "sin", [(F(2 * g1), d)]
    sysm.d = d


def _draw_system():
    typ = random.choices([t for t, _ in TYPES], weights=[w for _, w in TYPES])[0]
    s = System(typ)
    if typ == "fa25":
        while True:
            a, b = random.choice([1, 2, 3, 3, 4, 5]), random.choice([0, 0, 1, 2, 3])
            if a != b:
                break
        sg = random.choice([-1, -1, 1])
        if random.random() < 0.7:          # advanced term first, as in the exam
            s.taps = {-a: 1, b: sg}
            s.sys_tex = "y[n] = " + fmt.tex_sum([(1, x_tex(a)), (sg, x_tex(-b))])
        else:
            s.taps = {-a: sg, b: 1}
            s.sys_tex = "y[n] = " + fmt.tex_sum([(1, x_tex(-b)), (sg, x_tex(a))])
        _two_tap(s, -a, s.taps[-a], b, s.taps[b])
    elif typ == "sp25":
        d = random.choice([1, 2, 3, 3])
        c = random.choice([2, 2, 1, -1, -2, 1, 3])
        s.taps = {0: 1, d: c, 2 * d: 1}
        s.m, s.J, s.fn, s.terms, s.d, s.c = F(d), F(0), "cos", [(F(c), F(0)), (F(2), F(d))], F(d), c
        s.sys_tex = "y[n] = " + fmt.tex_sum([(1, "x[n]"), (c, x_tex(-d)), (1, x_tex(-2 * d))])
    elif typ == "fa24":
        d = random.choice([1, 2, 2, 3, 4])
        c = random.choice([1, 1, -1, -2])
        s.taps = {-d: 1, 0: c, d: 1}
        s.m, s.J, s.fn, s.terms, s.d, s.c = F(0), F(0), "cos", [(F(c), F(0)), (F(2), F(d))], F(d), c
        s.sys_tex = "h[n] = " + fmt.tex_sum([(1, d_tex(-d)), (c, d_tex(0)), (1, d_tex(d))])
    elif typ == "sp23":
        N = random.choice([2, 3, 4, 5, 6, 6])
        sg = random.choice([1, 1, -1])
        s.taps = {0: 1, N: sg}
        _two_tap(s, 0, 1, N, sg)
        s.sys_tex = "h[n] = " + fmt.tex_sum([(1, d_tex(0)), (sg, d_tex(N))])
    elif typ == "sp21b":
        while True:
            a, b = random.choice([1, 1, 2, 3]), random.choice([2, 3, 4, 5, 6, 7])
            if a != b:
                break
        g = random.choice([2, 3])
        sg = random.choice([-1, -1, -1, 1])
        s.taps = {-a: g, b: sg * g}
        _two_tap(s, -a, g, b, sg * g)
        s.sys_tex = "h[n] = " + fmt.tex_sum([(g, d_tex(-a)), (sg * g, d_tex(b))])
    else:   # sp21a
        g, M, K = random.choice([1, 1, 1, 2]), random.choice([1, 2, 3, 4, 4]), random.choice([1, 2, 2, 3])
        s.closed = (g, M, K)
        s.m, s.J, s.fn, s.terms, s.d = F(M), F(0), "sin", [(F(g), F(K))], F(K)
        s.sys_tex = r"H_d(\omega) = " + ("" if g == 1 else str(g)) + ej_tex(M) + trig_tex("sin", K)
    return s


# ----------------------------------------------------------------------------- zeros and jumps
def _zeros(sysm):
    """(crossings, touches): zeros of A on (-1, 1) (multiples of pi, Fractions), verified exactly."""
    t = -1 + (2 * np.arange(GRID_N) + 1) / GRID_N
    a = sysm.A_num(np.pi * t)
    cross, touch = [], []
    for i in range(GRID_N - 1):
        if (a[i] > 0) != (a[i + 1] > 0):
            lo, hi = t[i], t[i + 1]
            for _ in range(60):
                mid = (lo + hi) / 2
                if (sysm.A_num(np.pi * mid) > 0) == (a[i] > 0):
                    lo = mid
                else:
                    hi = mid
            cross.append(F(lo).limit_denominator(60))
    for i in range(1, GRID_N - 1):
        if abs(a[i]) < 1e-4 and abs(a[i]) <= abs(a[i - 1]) and abs(a[i]) <= abs(a[i + 1]) \
                and (a[i - 1] > 0) == (a[i + 1] > 0) == (a[i] > 0):
            z = F(float(t[i])).limit_denominator(60)
            if z not in touch:
                touch.append(z)
    for z in cross + touch:
        if sysm.A_exact(z) != 0:
            raise ValueError(f"zero {z} not exact")
    return sorted(cross), sorted(touch)


# ----------------------------------------------------------------------------- frequencies
def _phase_cands(sysm):
    out = []
    for f in PH_FREQS:
        for t in (f, -f):
            A = sysm.A_exact(t)
            if A == 0 or not nice(A):
                continue
            neg = bool(A < 0)
            ramp = -sysm.m * t + sysm.J
            raw = ramp + (1 if neg else 0)
            ph = principal(raw)
            if ph == 1:
                continue
            out.append({"t": t, "A": A, "neg": neg, "ramp": ramp, "raw": raw, "ph": ph, "wraps": raw != ph})
    return out


def _pick_phase(cands):
    """Two frequencies: at least one with A < 0 (if any), at least one whose angle needs a 2 pi
    wrap (raw = ramp + pi[A < 0] outside (-1, 1]); preferably (70 %) a wrap at a point with A > 0,
    where it cannot be avoided by choosing -pi instead of +pi."""
    if len(cands) < 2 or not any(c["wraps"] for c in cands):
        return None
    any_neg = any(c["neg"] for c in cands)
    any_pos = any(not c["neg"] for c in cands)
    pos_wrap = any(c["wraps"] and not c["neg"] for c in cands) and random.random() < 0.7
    both = any_neg and any_pos and random.random() < 0.5
    for _ in range(300):
        n_negf = random.choice([0, 0, 1])
        pos = [c for c in cands if c["t"] > 0]
        negs = [c for c in cands if c["t"] < 0]
        if len(pos) < 2 - n_negf or len(negs) < n_negf:
            continue
        pick = random.sample(pos, 2 - n_negf) + random.sample(negs, n_negf)
        if abs(pick[0]["t"]) == abs(pick[1]["t"]) or pick[0]["ph"] == pick[1]["ph"]:
            continue
        if any_neg and not any(c["neg"] for c in pick):
            continue
        if both and all(c["neg"] for c in pick):
            continue
        if not any(c["wraps"] for c in pick):
            continue
        if pos_wrap and not any(c["wraps"] and not c["neg"] for c in pick):
            continue
        return sorted(pick, key=lambda c: c["t"])
    return None


def _pick_values(sysm, avoid=()):
    """fa24: two frequencies in [0, pi] with nonzero nice real values, at least one negative."""
    cands = [(t, sysm.A_exact(t)) for t in MAG_FREQS if t not in avoid]
    cands = [(t, A) for t, A in cands if A != 0 and simple(A)]
    negs = [c for c in cands if c[1] < 0]
    if not negs or len(cands) < 2:
        return None
    first = random.choice(negs)
    second = random.choice([c for c in cands if c[0] != first[0]])
    return sorted([first, second], key=lambda c: c[0])


def _pick_mag(sysm, avoid):
    cands = [t for t in MAG_FREQS if t not in avoid]
    cands = [(t, sysm.A_exact(t)) for t in cands]
    cands = [(t, A) for t, A in cands if A != 0 and simple(A)]
    if not cands:
        return None
    rational = [c for c in cands if sp.nsimplify(c[1]).is_rational]
    surd = [c for c in cands if not sp.nsimplify(c[1]).is_rational]
    pool = surd if (surd and random.random() < 0.45) else (rational or surd)
    return random.choice(pool)


def _draw():
    while True:
        sysm = _draw_system()
        if sysm.typ == "fa24":
            vals = _pick_values(sysm)
            if not vals:
                continue
            mag = _pick_mag(sysm, avoid={vals[0][0], vals[1][0]})
            if not mag:
                continue
            return sysm, None, vals, mag
        if any(sysm.A_num(math.pi * float(s * f)) < -1e-9 for f in PH_FREQS for s in (1, -1)) and \
                not any(c["neg"] for c in _phase_cands(sysm)):
            continue                       # A < 0 occurs, but at no usable frequency
        pick = _pick_phase(_phase_cands(sysm))
        if not pick:
            continue
        mag = _pick_mag(sysm, avoid={abs(c["t"]) for c in pick})
        if not mag:
            continue
        return sysm, pick, None, mag


# ----------------------------------------------------------------------------- generate
def _A_tex(sysm):
    if len(sysm.terms) == 1:
        c, k = sysm.terms[0]
        return fmt.tex_sum([(c, trig_tex(sysm.fn, k))])
    return fmt.tex_sum([(c, "" if k == 0 else trig_tex(sysm.fn, k)) for c, k in sysm.terms])


def _H_tex(sysm):
    """Final form (j) e^{-jmw} A(w)."""
    lin = ej_tex(sysm.m)
    if len(sysm.terms) == 1:
        c, k = sysm.terms[0]
        ct = "" if c == 1 else ("-" if c == -1 else fmt.tex_num(c))
        jt = (r"j\," if lin else "j") if sysm.J else ""
        return ct + jt + lin + trig_tex(sysm.fn, k)
    A = _A_tex(sysm)
    return (lin + r"\left(" + A + r"\right)") if lin else A


def _def_tex(sysm):
    """sum_n h[n] e^{-j w n} written out."""
    return fmt.tex_sum([(v, ej_tex(n)) for n, v in sorted(sysm.taps.items())])


def _pair_tex(sysm):
    """The bracket after pulling out e^{-j m w}."""
    out = []
    for n, v in sorted(sysm.taps.items()):
        out.append((v, ej_tex(n - sysm.m)))
    return fmt.tex_sum(out)


def _row_phase(sysm, c):
    t, A = c["t"], c["A"]
    a_val = tex_real(A) + ("" if sp.nsimplify(A).is_rational else rf" \approx {dec(A)}")
    lin = -sysm.m * t
    steps = [rf"{tex_pi(lin)}{r' + \tfrac{\pi}{2}' if sysm.J else ''}" + ("" if not sysm.J else rf" = {tex_pi(c['ramp'])}")]
    if c["neg"]:
        steps.append(rf"{tex_pi(c['ramp'])} + \pi = {tex_pi(c['raw'])}")
    if c["raw"] != c["ph"]:
        k = (c["ph"] - c["raw"]) / 2
        steps.append(rf"{tex_pi(c['raw'])} {'+' if k > 0 else '-'} {tex_pi(2 * abs(k))} = {tex_pi(c['ph'])}")
    return {"w_tex": tex_pi(t), "A_line": rf"A\!\left({tex_pi(t)}\right) = {a_val} " + (r"\lt 0" if c["neg"] else r"\gt 0"),
            "steps": r" \;\to\; ".join(steps), "neg": c["neg"]}


def _list_tex(ts):
    return ",\\ ".join(tex_pi(t) for t in ts)


def generate(data):
    sysm, pick, vals, mag = _draw()
    p, cr = data["params"], data["correct_answers"]
    typ = sysm.typ
    p["typ"] = typ
    for t, _ in TYPES:
        p["is_" + t] = typ == t
    p["is_lccde"] = typ in ("fa25", "sp25")
    p["is_h"] = typ in ("fa24", "sp23", "sp21b")
    p["is_closed"] = typ == "sp21a"
    p["is_ramp"] = typ != "fa24"
    p["source_html"] = SOURCE[typ]
    p["sys_tex"] = sysm.sys_tex
    p["noncausal"] = bool(sysm.taps and min(sysm.taps) < 0)
    p["m_zero"] = sysm.m == 0
    if sysm.taps:
        p["spec"] = {"taps": [[n, v] for n, v in sorted(sysm.taps.items())]}
        p["h_tex"] = "h[n] = " + fmt.tex_sum([(v, d_tex(n)) for n, v in sorted(sysm.taps.items())])
        p["def_tex"] = _def_tex(sysm)
        p["pair_tex"] = _pair_tex(sysm)
        p["lin_tex"] = ej_tex(sysm.m)
        p["m_zero"] = sysm.m == 0
        p["euler_tex"] = (r"e^{j\theta} + e^{-j\theta} = 2\cos\theta" if sysm.fn == "cos"
                          else r"e^{j\theta} - e^{-j\theta} = 2j\sin\theta")
    else:
        g, M, K = sysm.closed
        p["spec"] = {"closed": {"g": g, "M": M, "K": K}}
    p["H_tex"] = _H_tex(sysm)
    p["A_tex"] = _A_tex(sysm)
    p["has_j"] = bool(sysm.J)
    p["m_tex"] = fmt.tex_num(sysm.m)
    p["ramp_tex"] = ramp_tex(sysm.m) if sysm.m != 0 else "0"
    p["unit_tex"] = (r"j\," if sysm.J else "") + (ej_tex(sysm.m) or "1")

    # (a) magnitude
    t0, A0 = mag
    p["f0"] = plain(t0)
    p["f0_tex"] = tex_pi(t0)
    p["A0_tex"] = tex_real(A0)
    p["mag0_tex"] = tex_real(sp.Abs(A0))
    p["mag0_irr"] = not sp.nsimplify(A0).is_rational
    p["mag0_dec"] = dec(sp.Abs(A0))
    p["A0_sub_tex"] = fmt.tex_sum([(c, "" if k == 0 else rf"\{sysm.fn}\!\left({tex_pi(k * t0)}\right)")
                                   for c, k in sysm.terms])
    cr["mag"] = float(sp.Abs(A0))

    # (b) phases or values
    if pick:
        rows = [_row_phase(sysm, c) for c in pick]
        for i, (c, row) in enumerate(zip(pick, rows), start=1):
            p[f"f{i}"] = plain(c["t"])
            p[f"f{i}_tex"] = tex_pi(c["t"])
            p[f"ph{i}_plain"] = plain(c["ph"])
            cr[f"ph{i}"] = float(c["ph"])
            row["ph_tex"] = tex_pi(c["ph"])
            row["i"] = i
        p["ph_rows"] = rows
    else:
        vrows = []
        for i, (t, A) in enumerate(vals, start=1):
            p[f"f{i}"] = plain(t)
            p[f"f{i}_tex"] = tex_pi(t)
            cr[f"hv{i}"] = float(A)
            sub = fmt.tex_sum([(c, "" if k == 0 else rf"\cos\!\left({tex_pi(k * t)}\right)") for c, k in sysm.terms])
            vrows.append({"i": i, "w_tex": tex_pi(t), "sub": sub, "val": tex_real(A),
                          "irr": not sp.nsimplify(A).is_rational, "dec": dec(A)})
        p["val_rows"] = vrows

    # (c) jumps
    cross, touch = _zeros(sysm)
    cr["njumps"] = len(cross)
    p["cross"] = [plain(z) for z in cross]
    p["touch"] = [plain(z) for z in touch]
    p["cross_tex"] = _list_tex(cross)
    p["touch_tex"] = _list_tex(touch)
    p["has_cross"] = bool(cross)
    p["has_touch"] = bool(touch)
    p["no_zero"] = not cross and not touch
    if sysm.fn == "cos" and len(sysm.terms) == 2:
        c = sysm.c
        lo, hi = c - 2, c + 2
        p["range_tex"] = rf"{c - 2} \le A(\omega) \le {c + 2}"
        p["sign_word"] = ("never negative" if lo == 0 else "never positive" if hi == 0 else
                          "never zero" if lo > 0 or hi < 0 else "")
    else:
        p["range_tex"], p["sign_word"] = "", ""
    p["is_3tap"] = len(sysm.terms) == 2
    p["njumps"] = len(cross)

    # (d) group delay
    cr["tau"] = float(sysm.m)
    p["tau_tex"] = fmt.tex_num(sysm.m)
    p["advance"] = sysm.m < 0
    p["adv_tex"] = fmt.tex_num(-sysm.m)
