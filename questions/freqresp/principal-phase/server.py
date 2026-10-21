"""Principal phase in (-pi, pi], not (-pi/2, pi/2] (Lecture 16, sections 2-4: the principal angle;
the +-pi of a negative real factor; 2 pi wraps of a linear phase).

(a), (b)  H_d(w_1) = a + jb, H_d(w_2) = c + jd with small integers. (a) always has Re < 0 (Q2 or Q3:
          tan^{-1}(b/a) alone is off by pi); (b) is in another quadrant, or on an axis (on the negative
          real axis the angle is pi: pi and -pi are both accepted).
(c)       H_d(w) = e^{-jMw} A(w), A(w) = alpha + beta cos(Kw) real, M in {1, 2, 3}; the principal
          phase at three frequencies (multiples of pi/6 or pi/4, one possibly negative): at least
          one with A < 0 (+-pi), at least one whose raw angle -Mw (+ pi) needs a 2 pi wrap.
          Never A(w) = 0, never a principal angle of exactly +-pi.
Answers are in radians (floats). grade() accepts |sub - ans| <= tol modulo 2 pi, but only for
submissions inside [-pi - tol, pi + tol]: so -pi is accepted for pi, while ans +- 2 pi is not;
tol = 1e-3 |ans| + 1e-3, the element's own relabs band (rtol = atol = 1e-3), so PrairieLearn's
"incorrect" test value (1.1 tol away) is still wrong."""

import math
import random
from fractions import Fraction as F

import prairielearn as pl
import sympy as sp
from ece310 import fmt
from ece310.w5_freqresp import cos_pi, principal, rat, tex_pi, tex_real

FREQS = [F(k, 12) for k in (2, 3, 4, 6, 8, 9, 10, 12)]       # pi/6, pi/4, ..., pi
NAMES = ["za", "zb", "w1", "w2", "w3"]
RTOL = ATOL = 1e-3


def _rad(t):
    return float(t) * math.pi


def _r4(x):
    return f"{x:.4f}"


# ----------------------------------------------------------------------------- (a), (b)
def _cplx_tex(a, b):
    if b == 0:
        return str(a)
    bt = ("" if abs(b) == 1 else str(abs(b))) + "j"
    if a == 0:
        return ("-" if b < 0 else "") + bt
    return f"{a} {'-' if b < 0 else '+'} {bt}"


def _quad(a, b):
    if b == 0:
        return "pos" if a > 0 else "neg"
    if a == 0:
        return "up" if b > 0 else "down"
    return {(True, True): "Q1", (False, True): "Q2", (False, False): "Q3", (True, False): "Q4"}[(a > 0, b > 0)]


def _draw_z(kind):
    if kind == "neg":
        return -random.randint(1, 5), 0
    if kind in ("up", "down"):
        return 0, random.randint(1, 4) * (1 if kind == "up" else -1)
    while True:
        a, b = random.randint(-4, 4), random.randint(-4, 4)
        if a != 0 and b != 0 and _quad(a, b) == kind:
            return a, b


def _z_angle(a, b):
    """(angle in radians, steps LaTeX list, naive arctan or None)."""
    q = _quad(a, b)
    if q == "neg":
        return math.pi, [rf"\text{{negative real axis}}", r"\angle = \pi\ (\equiv -\pi)"], None
    if q in ("up", "down"):
        s = 1 if q == "up" else -1
        rel = r"\gt" if s > 0 else r"\lt"
        return s * math.pi / 2, [rf"\text{{Re}} = 0,\ \text{{Im}} {rel} 0",
                                 rf"\angle = {'' if s > 0 else '-'}\tfrac{{\pi}}{{2}} \approx {_r4(s * math.pi / 2)}"], None
    r = F(b, a)
    base = math.atan(float(r))
    at = rf"\tan^{{-1}}\!\left({fmt.tex_num(r)}\right) \approx {_r4(base)}"
    if a > 0:
        return base, [rf"\text{{Re}} \gt 0\ (\text{{{q}}})", rf"\angle = {at}"], None
    if b > 0:
        ang = base + math.pi
        return ang, [rf"\text{{Re}} \lt 0,\ \text{{Im}} \gt 0\ (\text{{Q2}})", rf"{at}\ (\text{{Q4: wrong}})",
                     rf"\angle = {_r4(base)} + \pi \approx {_r4(ang)}"], base
    ang = base - math.pi
    return ang, [rf"\text{{Re}} \lt 0,\ \text{{Im}} \lt 0\ (\text{{Q3}})", rf"{at}\ (\text{{Q1: wrong}})",
                 rf"\angle = {_r4(base)} - \pi \approx {_r4(ang)}"], base


# ----------------------------------------------------------------------------- (c)
def _draw_A():
    fam = random.choices(["one2", "cos2", "c_cos", "gen"], weights=[25, 20, 30, 25])[0]
    if fam == "one2":
        return F(random.choice([1, -1])), F(2), 1
    if fam == "cos2":
        return F(0), F(random.choice([2, 1, -2])), 2
    if fam == "c_cos":
        return random.choice([F(1, 2), F(-1, 2)]), F(1), random.choice([1, 1, 2])
    while True:
        al = random.choice([F(0), F(1, 2), F(-1, 2), F(1), F(-1)])
        be = F(random.choice([1, 2, -1, -2]))
        K = random.choice([1, 2])
        if abs(al) < abs(be):
            return al, be, K


def _A_exact(al, be, K, t):
    return sp.nsimplify(rat(al) + rat(be) * cos_pi(K * t))


def _A_tex(al, be, K):
    return fmt.tex_sum([(al, ""), (be, r"\cos\omega" if K == 1 else r"\cos(2\omega)")])


def _lin_tex(M):
    return r"e^{-j\omega}" if M == 1 else f"e^{{-j{M}\\omega}}"


def _pick_freqs(al, be, K, M):
    cands = []
    for f in FREQS:
        for s in (1, -1):
            t = s * f
            A = _A_exact(al, be, K, t)
            if A == 0:
                continue
            neg = bool(A < 0)
            raw = -M * t + (1 if neg else 0)
            ph = principal(raw)
            if ph == 1:
                continue
            cands.append((t, neg, raw, ph))
    for _ in range(200):
        n_neg = random.choice([0, 0, 1])                      # how many negative frequencies
        pos = [c for c in cands if c[0] > 0]
        negs = [c for c in cands if c[0] < 0]
        if len(pos) < 3 - n_neg or len(negs) < n_neg:
            return None
        pick = random.sample(pos, 3 - n_neg) + random.sample(negs, n_neg)
        if len({abs(c[0]) for c in pick}) < 3:
            continue
        if not any(c[1] for c in pick) or all(c[1] for c in pick):
            continue
        if not any(c[2] != c[3] for c in pick):
            continue
        if len({c[3] for c in pick}) < 3:
            continue
        pick.sort(key=lambda c: c[0])
        return pick
    return None


def generate(data):
    p, c = data["params"], data["correct_answers"]
    # (a), (b)
    k1 = random.choice(["Q2", "Q3"])
    k2 = random.choices(["Q1", "Q4", "Q2" if k1 == "Q3" else "Q3", "neg", "up", "down"],
                        weights=[25, 25, 25, 12, 7, 6])[0]
    zs = [_draw_z(k1), _draw_z(k2)]
    if random.random() < 0.3:
        zs.reverse()
    zrows = []
    for name, (a, b), lab in zip(("za", "zb"), zs, ("1", "2")):
        ang, steps, naive = _z_angle(a, b)
        c[name] = ang
        zrows.append({"lab": lab, "z_tex": _cplx_tex(a, b), "steps": " \\;\\to\\; ".join(steps), "ans": _r4(ang)})
        p[f"{name}_ab"] = [a, b]
        p[f"{name}_tex"] = _cplx_tex(a, b)
    p["zrows"] = zrows

    # (c)
    while True:
        al, be, K = _draw_A()
        M = random.choice([1, 2, 3, 2, 3])
        pick = _pick_freqs(al, be, K, M)
        if pick:
            break
    p["A"] = {"alpha": [al.numerator, al.denominator], "beta": int(be), "K": K, "M": M}
    p["A_tex"] = _A_tex(al, be, K)
    p["lin_tex"] = _lin_tex(M)
    p["M"] = M
    p["M_tex"] = "" if M == 1 else str(M)
    wrows = []
    for i, (t, neg, raw, ph) in enumerate(pick, start=1):
        name = f"w{i}"
        c[name] = _rad(ph)
        p[f"f{i}"] = [t.numerator, t.denominator]
        p[f"f{i}_tex"] = tex_pi(t)
        A = _A_exact(al, be, K, t)
        lin = -M * t
        A_val = tex_real(A)
        Aline = (rf"A\!\left({tex_pi(t)}\right) = " + fmt.tex_sum([(al, ""), (be, rf"\cos\!\left({tex_pi(K * t)}\right)")])
                 + rf" = {A_val}" + ("" if A.is_rational else rf" \approx {float(A):.4g}") + (r" \lt 0" if neg else r" \gt 0"))
        steps = [rf"-{M if M > 1 else ''}\omega = {tex_pi(lin)}"]
        if neg:
            steps.append(rf"{tex_pi(lin)} + \pi = {tex_pi(raw)}")
        if raw != ph:
            k = (ph - raw) / 2
            steps.append(rf"{tex_pi(raw)} {'+' if k > 0 else '-'} {tex_pi(2 * abs(k))} = {tex_pi(ph)}")
        steps = [r" \;\to\; ".join(steps) + rf"\ (\approx {_r4(_rad(ph))})"]
        wrows.append({"i": i, "w_tex": tex_pi(t), "A_line": Aline, "steps": r" \;\to\; ".join(steps),
                      "neg": neg, "wrapped": raw != ph})
    p["wrows"] = wrows


def _ok(sub, ans):
    tol = RTOL * abs(ans) + ATOL
    try:
        sub = float(sub)
        d = (sub - ans + math.pi) % (2 * math.pi) - math.pi
        # NaN / inf fail every comparison
        return bool(-math.pi - tol <= sub <= math.pi + tol and abs(d) <= tol)
    except Exception:
        return False


def grade(data):
    changed = False
    for name in NAMES:
        if name in data["format_errors"] or data["submitted_answers"].get(name) is None:
            continue
        ps = data["partial_scores"].get(name)
        if ps is None:
            continue
        ps["score"] = 1.0 if _ok(data["submitted_answers"][name], data["correct_answers"][name]) else 0.0
        changed = True
    if changed:
        pl.set_weighted_score_data(data)
