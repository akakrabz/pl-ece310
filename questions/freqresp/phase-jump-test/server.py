"""Which phase responses have +-pi jumps? (Lecture 16, slide 6 concept check; the sign of the real
amplitude decides, a zero that is only touched is no jump, 2 pi wraps are no jumps.)

Four systems H_k(w) = (unit-modulus factor) * A_k(w) with A_k real, drawn from shapes of the slide and
of past Midterm 2 problems (new numbers):
  quad   e^{-jMw}(a0 + a1 cos w + a2 cos 2w)        L16 slide 6 H2, H4, H5 (2 + cos 2w; 1 - 2cos w - cos 2w;
                                                   3 - 2cos w - cos 2w, touches 0 at w = 0)
  cosK   e^{-jMw} cos(Kw)                           slide H1 (cos 4w)
  sinK   j e^{-jMw} sin(Kw)                         slide H3 (j e^{-j3w} sin 3w)
  even   (b + cos w) e^{+-j|w|}, |w| <= pi          FA2023 MT2 #4(a) ((1/2 + cos w) e^{j|w|})
  fir3a  h[n] = delta[n+d] + c delta[n] + delta[n-d] FA2024 MT2 #3 (delta[n+2] + delta[n] + delta[n-2])
  fir3b  y[n] = x[n] + c x[n-d] + x[n-2d]           SP2025 MT2 #3 (x[n] + 2x[n-3] + x[n-6], touches only)
  fir2   h[n] = delta[n] + delta[n-N]                SP2023 MT2 #5 (delta[n] + delta[n-6])
Each system is in one category: "jump" (A changes sign on (-pi, pi)), "touch" (A has zeros but never
changes sign: 3 - 2cos w - cos 2w, 2 + 2cos 3w, 1 + cos 2w, (1 + cos w) e^{j|w|}), "range" (A never 0:
2 + cos 2w, |a1| + |a2| < |a0| for the quad shape; 3 + 2cos dw; (2 + cos w) e^{j|w|}).
The four always contain 1-3 jump systems, a touch system whenever a non-jump one is present and a
range system whenever two are. Asked:
  (a) pl-checkbox (net-correct): the systems whose phase has +-pi jumps on (-pi, pi);
  (b) for one designated system (any category: 60 % jump, 30 % touch, 10 % range) the number of +-pi
      jumps on (-pi, pi) (integer, 0 allowed);
  (c) the smallest w in (0, pi) where its phase jumps, as a multiple of pi, or an empty box if there is
      none (allow-blank, correct answer ""); for jump systems only variants with a nice location.
Counts and locations come from a dense grid of A (never on a zero) with bisection, verified exactly
with SymPy; quad categories are also derived exactly from the quadratic in c = cos w."""

import math
import random
from fractions import Fraction as F

import numpy as np
import sympy as sp
from ece310 import fmt
from ece310.w5_freqresp import plain, rat, tex_pi

GRID_N = 4000
FIR_KINDS = ("fir3a", "fir3b", "fir2")
KINDS_FOR = {
    "jump": ["quad", "quad", "cosK", "sinK", "even", "fir3a", "fir3b", "fir2"],
    "touch": ["quad", "quad", "even", "fir3a", "fir3b", "fir2"],
    "range": ["quad", "quad", "even", "fir3a", "fir3b"],
}
NICE_T = [F(0), F(1, 3), F(1, 2), F(2, 3), F(1), F(1, 6), F(1, 4), F(3, 4), F(5, 6)]


# ----------------------------------------------------------------------------- small LaTeX helpers
def kw(k):
    return r"\omega" if k == 1 else f"{k}\\omega"


def cos_tex(k, fn="cos"):
    return rf"\{fn}({kw(k)})"


def ej(M):
    return "" if M == 0 else (r"e^{-j\omega}" if M == 1 else f"e^{{-j{M}\\omega}}")


def d_tex(n):
    return r"\delta[n]" if n == 0 else (rf"\delta[n-{n}]" if n > 0 else rf"\delta[n+{-n}]")


def x_tex(s):
    return "x[n]" if s == 0 else (f"x[n+{s}]" if s > 0 else f"x[n-{-s}]")


def acos_t(c):
    """arccos(c)/pi for the nice values of c (exact), else None."""
    table = {sp.Integer(1): F(0), sp.sqrt(3) / 2: F(1, 6), sp.sqrt(2) / 2: F(1, 4), sp.Rational(1, 2): F(1, 3),
             sp.Integer(0): F(1, 2), -sp.Rational(1, 2): F(2, 3), -sp.sqrt(2) / 2: F(3, 4), -sp.sqrt(3) / 2: F(5, 6),
             sp.Integer(-1): F(1)}
    c = sp.nsimplify(c)
    for k, v in table.items():
        if sp.simplify(c - k) == 0:
            return v
    return None


# ----------------------------------------------------------------------------- the systems
class System:
    """kind, spec (JSON), A (float fn of w), A_sym (sympy fn of t = w/pi), texts."""

    def zeros(self):
        """(crossings, touches) on (-1, 1) as Fractions (exactly verified), endpoint zero flag."""
        t = -1 + (2 * np.arange(GRID_N) + 1) / GRID_N
        a = self.A(np.pi * t)
        cross, touch = [], []
        for i in range(GRID_N - 1):
            if (a[i] > 0) != (a[i + 1] > 0):
                lo, hi = t[i], t[i + 1]
                for _ in range(60):
                    mid = (lo + hi) / 2
                    if (self.A(np.pi * mid) > 0) == (a[i] > 0):
                        lo = mid
                    else:
                        hi = mid
                cross.append(F(lo).limit_denominator(60))
        for i in range(1, GRID_N - 1):
            if abs(a[i]) < 1e-4 and abs(a[i]) <= abs(a[i - 1]) and abs(a[i]) <= abs(a[i + 1]) \
                    and (a[i - 1] > 0) == (a[i] > 0) == (a[i + 1] > 0):
                z = F(float(t[i])).limit_denominator(60)
                if z not in touch:
                    touch.append(z)
        exact_c = [z for z in cross if self.A_sym(z) == 0]
        exact_t = [z for z in touch if self.A_sym(z) == 0]
        end = abs(float(self.A(np.pi))) < 1e-12
        return sorted(cross), sorted(exact_t), end, len(exact_c) == len(cross), len(exact_t) == len(touch)

    def analyse(self):
        cross, touch, end, ok_c, ok_t = self.zeros()
        self.cross, self.touch, self.end_zero = cross, touch, end
        self.cross_exact, self.touch_exact = ok_c, ok_t
        self.njumps = len(cross)
        self.cat = "jump" if cross else ("touch" if (touch or end) else "range")
        pos = [z for z in cross if z > 0]
        self.first = pos[0] if pos else None
        # designatable: jumps at exact (nice) locations, or no jump at all
        self.nice = (not cross) or (ok_c and all(z.denominator <= 12 for z in cross))
        return self


def _quad_class(a0, a1, a2):
    """Exact category of a0 + a1 cos w + a2 cos 2w from Q(c) = 2a2 c^2 + a1 c + (a0 - a2), |c| <= 1."""
    al, be, ga = 2 * a2, a1, a0 - a2
    Q = lambda c: al * c * c + be * c + ga
    if al == 0:
        if be == 0:
            return None
        r = F(-ga, be)
        return "jump" if -1 < r < 1 else ("touch" if abs(r) == 1 else "range")
    D = be * be - 4 * al * ga
    if D < 0:
        return "range"
    if D == 0:
        r = F(-be, 2 * al)
        return "touch" if -1 <= r <= 1 else "range"
    roots = [(-be + s * math.sqrt(D)) / (2 * al) for s in (1, -1)]
    cross = [r for r in roots if -1 + 1e-12 < r < 1 - 1e-12]
    if cross:
        return "jump"
    if Q(1) == 0 or Q(-1) == 0:
        return "touch"
    return "range"


QUADS = {"jump": [], "touch": [], "range": []}
for _a0 in range(-4, 5):
    for _a1 in range(-3, 4):
        for _a2 in range(-2, 3):
            if _a1 == 0 and _a2 == 0:
                continue
            _c = _quad_class(_a0, _a1, _a2)
            if _c == "range" and not abs(_a1) + abs(_a2) < abs(_a0):
                continue          # range cases only where the bound |a1| + |a2| < |a0| shows it at once
            if _c == "jump":      # a pair of nice points of opposite sign must show it
                _v = [_a0 + _a1 * math.cos(math.pi * t) + _a2 * math.cos(2 * math.pi * t) for t in NICE_T]
                if not (max(_v) > 1e-9 and min(_v) < -1e-9):
                    continue
            if _c and math.gcd(math.gcd(_a0, _a1), _a2) == 1:
                QUADS[_c].append((_a0, _a1, _a2))


def make(kind, cat, idx):
    s = System()
    s.kind, s.idx = kind, idx
    lab = f"H_{idx}"
    if kind == "quad":
        a0, a1, a2 = random.choice(QUADS[cat])
        M = random.choice([0, 1, 2, 2, 3])
        s.spec = {"kind": kind, "a": [a0, a1, a2], "M": M}
        s.A = lambda w, a=(a0, a1, a2): a[0] + a[1] * np.cos(w) + a[2] * np.cos(2 * w)
        s.A_sym = lambda t, a=(a0, a1, a2): sp.nsimplify(a[0] + a[1] * sp.cos(sp.pi * rat(t)) + a[2] * sp.cos(2 * sp.pi * rat(t)))
        s.A_tex = fmt.tex_sum([(a0, ""), (a1, cos_tex(1)), (a2, cos_tex(2))])
        s.unit_tex = ej(M)
        s.def_tex = rf"{lab}(\omega) = " + (ej(M) + r"\left(" + s.A_tex + r"\right)" if M else s.A_tex)
        s.abc = (a0, a1, a2)
    elif kind == "cosK":
        K, M = random.choice([2, 3, 4]), random.choice([0, 0, 1, 2])
        s.spec = {"kind": kind, "K": K, "M": M}
        s.A = lambda w, K=K: np.cos(K * w)
        s.A_sym = lambda t, K=K: sp.cos(K * sp.pi * rat(t))
        s.A_tex = cos_tex(K)
        s.unit_tex = ej(M)
        s.def_tex = rf"{lab}(\omega) = " + ej(M) + cos_tex(K)
    elif kind == "sinK":
        K, M = random.choice([1, 2, 3, 3]), random.choice([1, 2, 3])
        s.spec = {"kind": kind, "K": K, "M": M}
        s.A = lambda w, K=K: np.sin(K * w)
        s.A_sym = lambda t, K=K: sp.sin(K * sp.pi * rat(t))
        s.A_tex = cos_tex(K, "sin")
        s.unit_tex = r"j\," + ej(M)
        s.def_tex = rf"{lab}(\omega) = j\," + ej(M) + cos_tex(K, "sin")
    elif kind == "even":
        b = random.choice({"jump": [F(1, 2), F(1, 2), F(-1, 2), F(0)], "touch": [F(1), F(-1)],
                           "range": [F(3, 2), F(2), F(-2), F(-3, 2)]}[cat])
        sg = random.choice([1, -1])
        s.spec = {"kind": kind, "b": plain(b), "s": sg}
        s.A = lambda w, b=float(b): b + np.cos(w)
        s.A_sym = lambda t, b=b: sp.nsimplify(rat(b) + sp.cos(sp.pi * rat(t)))
        s.A_tex = fmt.tex_sum([(b, ""), (1, cos_tex(1))])
        s.unit_tex = r"e^{j|\omega|}" if sg == 1 else r"e^{-j|\omega|}"
        s.def_tex = rf"{lab}(\omega) = \left({s.A_tex}\right){s.unit_tex},\quad |\omega| \le \pi"
        s.b = b
    elif kind in ("fir3a", "fir3b"):
        c = random.choice({"jump": [1, 1, -1], "touch": [2, -2, 2], "range": [3, -3, 3]}[cat])
        d = random.choice([1, 2, 2, 3] if kind == "fir3a" else [1, 2, 3, 3])
        taps = {-d: 1, 0: c, d: 1} if kind == "fir3a" else {0: 1, d: c, 2 * d: 1}
        m = 0 if kind == "fir3a" else d
        s.spec = {"kind": kind, "taps": [[n, v] for n, v in sorted(taps.items())]}
        s.A = lambda w, c=c, d=d: c + 2 * np.cos(d * w)
        s.A_sym = lambda t, c=c, d=d: sp.nsimplify(c + 2 * sp.cos(d * sp.pi * rat(t)))
        s.A_tex = fmt.tex_sum([(c, ""), (2, cos_tex(d))])
        s.unit_tex = ej(m)
        if kind == "fir3a":
            s.def_tex = rf"h_{idx}[n] = " + fmt.tex_sum([(1, d_tex(-d)), (c, d_tex(0)), (1, d_tex(d))])
        else:
            s.def_tex = "y[n] = " + fmt.tex_sum([(1, "x[n]"), (c, x_tex(-d)), (1, x_tex(-2 * d))])
        mid = fmt.tex_sum([(1, f"e^{{j{kw(d)}}}"), (c, ""), (1, f"e^{{-j{kw(d)}}}")])
        s.derive_tex = (rf"{lab}(\omega) = " + (ej(m) + r"\left(" + mid + r"\right) = " if m else mid + " = ")
                        + (ej(m) + r"\left(" + s.A_tex + r"\right)" if m else s.A_tex))
        s.c, s.d = c, d
    else:   # fir2
        N = random.choice([2, 3, 4, 5, 6, 6]) if cat == "jump" else 1
        s.spec = {"kind": kind, "taps": [[0, 1], [N, 1]]}
        half = F(N, 2)
        s.A = lambda w, h=float(half): 2 * np.cos(h * w)
        s.A_sym = lambda t, h=half: sp.nsimplify(2 * sp.cos(sp.pi * rat(h) * rat(t)))
        hN = kw(N // 2) if N % 2 == 0 else r"\tfrac{" + str(N) + r"}{2}\omega"
        s.A_tex = rf"2\cos({hN})" if N % 2 == 0 else rf"2\cos\left({hN}\right)"
        lin = f"e^{{-j{hN}}}"
        s.unit_tex = lin
        s.def_tex = rf"h_{idx}[n] = " + fmt.tex_sum([(1, d_tex(0)), (1, d_tex(N))])
        s.derive_tex = (rf"{lab}(\omega) = 1 + e^{{-j{kw(N)}}} = {lin}\left(e^{{j{hN}}} + e^{{-j{hN}}}\right) = {lin}\cdot {s.A_tex}")
        s.N = N
    s.analyse()
    return s


# ----------------------------------------------------------------------------- explanations
def _val_tex(v):
    v = sp.nsimplify(v)
    return sp.latex(v).replace(r"\frac", r"\tfrac")


def _where(s):
    parts = [tex_pi(z) for z in s.touch] + ([r"\pm\pi"] if s.end_zero else [])
    return ",\\ ".join(parts)


def why(s):
    """One-line HTML reason for the category (the slide's speed tests)."""
    cat, k = s.cat, s.kind
    zs = ",\\ ".join(tex_pi(z) for z in s.cross)
    if k == "quad":
        a0, a1, a2 = s.abc
        if cat == "range":
            b = abs(a1) + abs(a2)
            return (rf"$|{fmt.tex_sum([(a1, cos_tex(1)), (a2, cos_tex(2))])}| \le {b} \lt {abs(a0)}$, so $A$ stays between "
                    rf"${a0 - b}$ and ${a0 + b}$: never $0$, no sign change.")
        if cat == "jump":
            pts = [(t, s.A_sym(t)) for t in NICE_T]
            for t1, v1 in pts:
                for t2, v2 in pts:
                    if v1 > 0 and v2 < 0:
                        return (rf"$A({tex_pi(t1)}) = {_val_tex(v1)} \gt 0$ and $A({tex_pi(t2)}) = {_val_tex(v2)} \lt 0$: two points of "
                                r"opposite sign, so $A$ changes sign in between.")
            raise ValueError("no sign pair")
        c = sp.Symbol("c")
        Qc = sp.factor(2 * a2 * c ** 2 + a1 * c + (a0 - a2))
        sign = "\\ge 0" if any(s.A_sym(t) > 0 for t in NICE_T) else "\\le 0"
        where = _where(s)
        return (rf"With $c = \cos\omega \in [-1, 1]$ and $\cos(2\omega) = 2c^2 - 1$: $A = {sp.latex(Qc)} {sign}$ on $[-1, 1]$, "
                rf"zero only at $\omega = {where}$: it touches $0$ there without changing sign.")
    if k == "cosK":
        return rf"$\cos(K\omega)$ with $K = {s.spec['K']}$ changes sign at every zero, $\omega = {zs}$."
    if k == "sinK":
        return rf"$\sin(K\omega)$ with $K = {s.spec['K']}$ changes sign at every zero, $\omega = {zs}$ (including $\omega = 0$)."
    if k == "even":
        b = s.b
        lo, hi = b - 1, b + 1
        rng = rf"${fmt.tex_num(b)} + \cos\omega$ lies between ${fmt.tex_num(lo)}$ and ${fmt.tex_num(hi)}$"
        tail = r" (the factor $e^{\pm j|\omega|}$ is continuous: a kink at $0$, no jump)."
        if cat == "jump":
            return rng + rf", with $0$ strictly inside: it changes sign at $\omega = {zs}$" + tail
        if cat == "touch":
            at = r"\omega = \pm\pi" if b == 1 else r"\omega = 0"
            return rng + rf": it only touches $0$ (at ${at}$, where $\cos\omega$ is extreme), no sign change" + tail
        return rng + ": never $0$" + tail
    if k in ("fir3a", "fir3b"):
        c = s.c
        rng = rf"$A = {s.A_tex}$ lies between ${c - 2}$ and ${c + 2}$"
        if cat == "jump":
            return rng + rf": it changes sign at $\omega = {zs}$."
        if cat == "touch":
            where = _where(s)
            return rng + rf": it is {'never negative' if c > 0 else 'never positive'} and only touches $0$ (at $\omega = {where}$)."
        return rng + ": never $0$."
    if s.cat == "jump":
        return rf"$\cos\left(\tfrac{{{s.N}\omega}}{{2}}\right)$ changes sign at $\omega = {zs}$."
    return r"$\cos\left(\tfrac{\omega}{2}\right) \gt 0$ on $(-\pi, \pi)$; its only zeros are the endpoints $\pm\pi$: no jump inside."


# ----------------------------------------------------------------------------- composition
def _compose():
    while True:
        kj = random.choices([1, 2, 3], weights=[2, 3, 2])[0]
        cats = ["jump"] * kj
        n_non = 4 - kj
        if n_non >= 1:
            cats.append("touch")
        if n_non >= 2:
            cats.append("range")
        if n_non == 3:
            cats.append(random.choice(["touch", "range"]))
        random.shuffle(cats)
        kinds = [random.choice(KINDS_FOR[c]) for c in cats]
        if "quad" not in kinds or not any(k in FIR_KINDS for k in kinds):
            continue
        if kinds.count("quad") > 2 or any(kinds.count(k) > 1 for k in set(kinds) - {"quad"}):
            continue
        systems = [make(k, c, i + 1) for i, (k, c) in enumerate(zip(kinds, cats))]
        if any(s.cat != c for s, c in zip(systems, cats)):
            continue                      # the numeric analysis must agree with the exact category
        if len({s.def_tex.split("=", 1)[1] for s in systems}) < 4:
            continue
        weights = {"jump": 6, "touch": 3, "range": 1}
        cand = [s for s in systems if s.nice and (s.cat != "jump" or s.first is not None or s.kind == "sinK")]
        if not cand:
            continue
        des = random.choices(cand, weights=[weights[s.cat] for s in cand])[0]
        return systems, des


def generate(data):
    systems, des = _compose()
    p, cr = data["params"], data["correct_answers"]
    p["systems"] = []
    for s in systems:
        is_fir = s.kind in FIR_KINDS
        p["systems"].append({
            "idx": s.idx, "spec": s.spec, "def_tex": s.def_tex, "is_fir": is_fir, "is_lccde": s.kind == "fir3b",
            "jumps": s.cat == "jump", "cat": s.cat, "why_html": why(s),
            "form_tex": rf"H_{s.idx}(\omega) = " + ((s.unit_tex + r"\cdot ") if s.unit_tex else "") + "A(\\omega)",
            "A_tex": s.A_tex, "derive_tex": getattr(s, "derive_tex", ""), "unit_tex": s.unit_tex or "1",
            "njumps": s.njumps,
        })
    p["jump_choices"] = [{"text": f"System {s.idx}", "correct": s.cat == "jump"} for s in systems]
    p["k"] = des.idx
    p["des_cat"] = des.cat
    p["des_why_html"] = why(des)
    p["des_n"] = des.njumps
    p["des_cross_tex"] = ",\\ ".join(tex_pi(z) for z in des.cross)
    p["des_has_cross"] = bool(des.cross)
    p["des_touch_tex"] = _where(des)
    p["des_has_touch"] = bool(des.touch or des.end_zero)
    p["has_loc"] = des.first is not None
    p["cross"] = [plain(z) for z in des.cross]
    p["touch"] = [plain(z) for z in des.touch]
    cr["njumps"] = des.njumps
    if des.first is not None:
        cr["wjump"] = float(des.first)
        p["loc_tex"] = tex_pi(des.first)
        p["loc_plain"] = plain(des.first)
    else:
        cr["wjump"] = ""
        p["loc_why"] = ("its only jump is at $\\omega = 0$" if des.cross else "it has no jump at all")
