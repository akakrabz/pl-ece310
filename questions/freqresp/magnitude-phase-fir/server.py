"""Magnitude and phase of a (anti)symmetric FIR filter by pulling out a linear-phase factor
(Lecture 16 section 3; SP2025 Midterm 2 #3, FA2024 #3, SP2023 #5; HW6 #6).

h[n] is a short symmetric (or antisymmetric) sequence from a vetted bank, scaled by s in {1, -1, 2}
and shifted so that its centre alpha is an integer or half-integer between -2 and 3. Then
    symmetric:     H_d(w) = e^{-j alpha w} A(w),     A(w) = sum_i h_i cos((c - i) w)
    antisymmetric: H_d(w) = j e^{-j alpha w} A(w),   A(w) = sum_i h_i sin((c - i) w)
(c = (L - 1)/2, the centre measured from the first tap). Asked: alpha; A(w) (pl-symbolic-input in w,
graded numerically in grade() because SymPy cannot always prove cos/exponential identities);
|H_d| at two frequencies; the principal phase at a frequency where A < 0 (the +-pi jump trap,
never exactly +-pi); the group delay (= alpha, jumps ignored); the smallest zero of |H_d| in (0, pi].
Every bank entry has a sign change of A and a zero in (0, pi] that is a multiple of pi/12."""

import cmath
import math
import random
from fractions import Fraction as F

import prairielearn as pl
import prairielearn.sympy_utils as psu
import sympy as sp
from ece310 import fmt
from ece310.w5_freqresp import dec, plain, principal, rat, tex_pi, tex_r, tex_real

# (taps, kind): A(w) changes sign on [-pi, pi] and |H_d| has a zero in (0, pi] at a multiple of pi/12
BANK = [
    ([1, 0, 1], "sym"), ([1, 1, 1], "sym"), ([1, -1, 1], "sym"),
    ([1, 0, 1, 0, 1], "sym"), ([1, 0, 0, 0, 1], "sym"), ([1, 0, 0, 0, 0, 0, 1], "sym"),
    ([1, 1, 0, 1, 1], "sym"), ([1, 0, -1, 0, 1], "sym"),
    ([1, 0, 0, 1], "sym"), ([1, 1, 1, 1], "sym"), ([1, 2, 2, 1], "sym"), ([1, -1, -1, 1], "sym"),
    ([1, 0, -1], "anti"), ([1, 0, 0, 0, -1], "anti"), ([1, 0, 0, -1], "anti"),
    ([1, 1, -1, -1], "anti"), ([1, -1, 1, -1], "anti"), ([1, 0, 0, 0, 0, 0, -1], "anti"),
    ([1, 2, 0, -2, -1], "anti"),
]
MAG_FREQS = [F(0), F(1, 6), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(5, 6), F(1)]
PHASE_FREQS = [F(1, 6), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(5, 6)]
W_TEST = [0.37, 1.23, 2.11, -0.81, 2.93, -2.47, 0.05, -1.66]


# ----------------------------------------------------------------------------- A(w) pieces
def _terms(h, kind):
    """[(coef, k)] with A(w) = sum coef * cos(k w) (sym) or coef * sin(k w) (anti), k > 0 or k = 0."""
    L = len(h)
    c = F(L - 1, 2)
    out = {}
    for i, v in enumerate(h):
        k = c - i
        if v == 0:
            continue
        if kind == "sym":
            out[abs(k)] = out.get(abs(k), 0) + v            # cos is even
        else:
            if k == 0:
                continue
            out[abs(k)] = out.get(abs(k), 0) + (v if k > 0 else -v)   # sin is odd
    return sorted(((v, k) for k, v in out.items() if v != 0), key=lambda t: t[1])


def _A_exact(terms, kind, t):
    """A(t pi) exactly."""
    fn = sp.cos if kind == "sym" else sp.sin
    return sp.nsimplify(sum((v * fn(sp.pi * rat(k) * rat(t)) for v, k in terms), sp.Integer(0)))


def _A_num(terms, kind, w):
    fn = math.cos if kind == "sym" else math.sin
    return sum(v * fn(float(k) * w) for v, k in terms)


def _kw_tex(k):
    """k*omega as an argument: \\omega, 2\\omega, \\tfrac{\\omega}{2}, \\tfrac{3\\omega}{2}."""
    k = F(k)
    if k == 1:
        return r"\omega"
    if k.denominator == 1:
        return f"{k.numerator}\\omega"
    p = "" if k.numerator == 1 else str(k.numerator)
    return r"\tfrac{" + p + r"\omega}{" + str(k.denominator) + "}"


def _A_tex(terms, kind):
    fn = r"\cos" if kind == "sym" else r"\sin"
    return fmt.tex_sum((v, "" if k == 0 else f"{fn} {_kw_tex(k)}" if F(k).denominator == 1 else f"{fn}{_kw_tex(k)}")
                       for v, k in terms)


def _A_sym(terms, kind):
    fn = "cos" if kind == "sym" else "sin"
    parts = [fmt.sym_num(v) if k == 0 else f"{fmt.sym_num(v)}*{fn}({fmt.sym_num(k)}*w)" for v, k in terms]
    return " + ".join(parts)


def _A_exp_sym(terms, kind):
    """The same A(w) written with exponentials (an equivalent form for the Test button / probes)."""
    parts = []
    for v, k in terms:
        if k == 0:
            parts.append(fmt.sym_num(v))
        elif kind == "sym":
            parts.append(f"({F(v) / 2})*(exp(j*{fmt.sym_num(k)}*w) + exp(-j*{fmt.sym_num(k)}*w))")
        else:
            parts.append(f"({F(v) / 2})*(exp(j*{fmt.sym_num(k)}*w) - exp(-j*{fmt.sym_num(k)}*w))/j")
    return " + ".join(parts)


def _exp_tex(n):
    """e^{-j w n} for the definition line."""
    if n == 0:
        return ""
    k = "" if abs(n) == 1 else str(abs(n))
    return f"e^{{{'-' if n > 0 else ''}j{k}\\omega}}"


def _pair_tex(h, kind):
    """The bracket after pulling out the centre: h_c + h(e^{jkw} + e^{-jkw}) + ..."""
    L = len(h)
    c = F(L - 1, 2)
    out = []
    for i in range(L // 2 + (L % 2)):
        k = c - i
        v = h[i]
        if v == 0:
            continue
        if k == 0:
            out.append((v, ""))
            continue
        kt = _kw_tex(k)
        op = "+" if kind == "sym" else "-"
        out.append((v, rf"\left(e^{{j{kt}}} {op} e^{{-j{kt}}}\right)"))
    return fmt.tex_sum(out)


# ----------------------------------------------------------------------------- generate
def _zeros(terms, kind):
    """Zeros of A on (0, pi] as Fractions (multiples of pi/12), exact."""
    zs = []
    for k in range(1, 13):
        t = F(k, 12)
        if abs(_A_num(terms, kind, math.pi * float(t))) < 1e-9 and _A_exact(terms, kind, t) == 0:
            zs.append(t)
    return zs


def _simple(x):
    """2 if x is rational, 1 if x^2 is rational (one square root), 0 otherwise."""
    x = sp.nsimplify(x)
    if x.is_rational:
        return 2
    return 1 if sp.nsimplify(sp.expand(x * x)).is_rational else 0


def _draw():
    while True:
        taps, kind = random.choice(BANK)
        s = random.choice([1, 1, -1, 2])
        h = [s * v for v in taps]
        L = len(h)
        c = F(L - 1, 2)
        alpha_choices = [a for a in (F(x, 2) for x in range(-4, 7)) if (a - c).denominator == 1]
        alpha = random.choice(alpha_choices)
        start = int(alpha - c)
        terms = _terms(h, kind)
        zeros = _zeros(terms, kind)
        if not zeros:
            continue
        j = F(1, 2) if kind == "anti" else F(0)
        negc, posc = [], []
        for f in PHASE_FREQS + [-f for f in PHASE_FREQS]:
            Av = float(_A_exact(terms, kind, f))
            if abs(Av) < 1e-12:
                continue
            neg = Av < 0
            ph = principal(-alpha * f + j + (1 if neg else 0))
            if ph not in (F(1), F(-1)):
                (negc if neg else posc).append((f, ph, neg))
        if not negc:
            continue
        cands = negc if (not posc or random.random() < 0.7) else posc     # usually the +-pi trap, not always
        pos = [cd for cd in cands if cd[0] > 0]
        f1, ph1, neg1 = random.choice(pos if pos and random.random() < 0.85 else cands)
        mags = [f for f in MAG_FREQS if f != abs(f1) and _A_exact(terms, kind, f) != 0]
        if len(mags) < 2:
            continue
        rational = [f for f in mags if _simple(_A_exact(terms, kind, f)) == 2]
        surd = [f for f in mags if _simple(_A_exact(terms, kind, f)) >= 1]
        if len(rational) >= 2 and random.random() < 0.7:
            pool = rational
        elif len(surd) >= 2 and random.random() < 0.9:
            pool = surd
        else:
            pool = mags
        fa, fb = sorted(random.sample(pool, 2))
        return h, start, kind, alpha, terms, zeros, f1, ph1, neg1, fa, fb


def generate(data):
    h, start, kind, alpha, terms, zeros, f1, ph1, neg1, fa, fb = _draw()
    p = data["params"]
    c = data["correct_answers"]
    p["h"] = h
    p["start"] = start
    p["kind"] = kind
    p["is_anti"] = kind == "anti"
    p["h_tex"] = fmt.tex_seq(h, start)
    p["sym_word"] = "symmetric" if kind == "sym" else "antisymmetric"
    p["n_first"], p["n_last"] = start, start + len(h) - 1
    p["form_tex"] = (r"e^{-j\alpha\omega}A(\omega)" if kind == "sym" else r"j\,e^{-j\alpha\omega}A(\omega)")
    jmark = F(1, 2) if kind == "anti" else F(0)

    # (a) alpha, (b) A
    p["alpha"] = plain(alpha)
    p["alpha_tex"] = fmt.tex_num(alpha)
    p["alpha_zero"] = alpha == 0
    p["lin_tex"] = (r"e^{-j" + (fmt.tex_num(alpha) if alpha != 1 else "") + r"\omega}" if alpha > 0 else
                    (r"e^{j" + (fmt.tex_num(-alpha) if alpha != -1 else "") + r"\omega}" if alpha < 0 else ""))
    if alpha.denominator != 1:
        p["lin_tex"] = (r"e^{-j\frac{" + str(abs(alpha.numerator)) + r"}{2}\omega}" if alpha > 0 else
                        r"e^{j\frac{" + str(abs(alpha.numerator)) + r"}{2}\omega}")
    ns = list(range(start, start + len(h)))
    p["def_tex"] = fmt.tex_sum((v, _exp_tex(n)) for n, v in zip(ns, h))
    p["pair_tex"] = _pair_tex(h, kind)
    p["euler_tex"] = (r"e^{jk\omega}+e^{-jk\omega} = 2\cos k\omega" if kind == "sym"
                      else r"e^{jk\omega}-e^{-jk\omega} = 2j\sin k\omega")
    p["A_tex"] = _A_tex(terms, kind)
    p["A_sym"] = _A_sym(terms, kind)
    p["A_exp_sym"] = _A_exp_sym(terms, kind)
    # No correct_answers["A"]: with it, pl-symbolic-input runs its own sympy `equals` check, which
    # times out on equivalent exponential forms under load ("did not converge" format error).
    # grade() decides numerically; the worked solution shows A(w).
    p["terms"] = [[plain(v), plain(k)] for v, k in terms]
    p["jpre_tex"] = "j" if kind == "anti" else ""

    # (c) magnitudes
    rows = []
    for name, f in (("mag_a", fa), ("mag_b", fb)):
        Av = _A_exact(terms, kind, f)
        m = sp.Abs(Av)
        c[name] = float(m)
        if _simple(m):
            a_tex = "= " + tex_real(Av)
            m_tex = "= " + tex_r(m)
        else:
            a_tex, m_tex = r"\approx " + dec(Av), r"\approx " + dec(m)
        rows.append({"w_tex": tex_pi(f), "A_val_tex": a_tex, "mag_tex": m_tex,
                     "irr": bool(_simple(m) == 1), "dec": dec(m)})
    p["fa"], p["fb"] = plain(fa), plain(fb)
    p["fa_tex"], p["fb_tex"] = tex_pi(fa), tex_pi(fb)
    p["mag_rows"] = rows

    # (d) phase at f1 (usually A < 0 there: the +-pi trap; sometimes A > 0: no pi)
    A1 = _A_exact(terms, kind, f1)
    raw = -alpha * f1 + jmark + (1 if neg1 else 0)
    p["f1"] = plain(f1)
    p["f1_tex"] = tex_pi(f1)
    p["A1_neg"] = neg1
    p["A1_tex"] = ("= " + tex_real(A1)) if _simple(A1) else (r"\approx " + dec(A1))
    p["rule_tex"] = r"-\alpha\omega" + (r" + \tfrac{\pi}{2}" if kind == "anti" else "") + (r" \pm \pi" if neg1 else "")
    p["lin_at_f1_tex"] = tex_pi(-alpha * f1)
    p["raw_ph_tex"] = tex_pi(-alpha * f1 + jmark) + (r" + \pi" if neg1 else "")
    p["raw_ph_sum_tex"] = tex_pi(raw)
    p["ph_tex"] = tex_pi(ph1)
    p["ph_plain"] = plain(ph1)
    p["ph_needs_wrap"] = raw != ph1
    p["show_sum"] = neg1 or raw != ph1
    p["alt_tex"] = tex_pi(principal(raw + 1))     # the other decision about the +-pi
    c["phase"] = float(ph1)

    # (e) group delay, (f) first zero
    c["alpha"] = float(alpha)
    c["tau"] = float(alpha)
    c["wzero"] = float(zeros[0])
    p["zeros"] = [plain(z) for z in zeros]
    p["zeros_tex"] = ",\\ ".join(tex_pi(z) for z in zeros)
    p["wzero_tex"] = tex_pi(zeros[0])
    p["wzero_plain"] = plain(zeros[0])
    p["several_zeros"] = len(zeros) > 1


# ----------------------------------------------------------------------------- grading of A(w)
def _A_true(params, w):
    fn = math.cos if params["kind"] == "sym" else math.sin
    return sum(float(F(v)) * fn(float(F(k)) * w) for v, k in params["terms"])


def _matches(expr, params):
    syms = {s.name: s for s in expr.free_symbols}
    if set(syms) - {"w"}:
        return False
    for w in W_TEST:
        try:
            val = complex(expr.subs(syms["w"], w).evalf(30)) if syms else complex(expr.evalf(30))
        except (TypeError, ValueError, ZeroDivisionError, OverflowError):
            return False
        tru = _A_true(params, w)
        if not abs(val - tru) <= 1e-8 * max(1.0, abs(tru)):
            return False
    return True


def grade(data):
    sub = data["submitted_answers"].get("A")
    if sub is None or "A" in data["format_errors"]:
        return
    try:
        expr = (psu.json_to_sympy(sub, allow_complex=True) if isinstance(sub, dict)
                else psu.convert_string_to_sympy(sub, ["w"], allow_complex=True))
        ok = _matches(expr, data["params"])
    except Exception:
        ok = False
    data["partial_scores"]["A"] = {"score": 1.0 if ok else 0.0, "weight": 1}
    pl.set_weighted_score_data(data)


def test(data):
    tt = data["test_type"]
    if tt == "correct":
        data["raw_submitted_answers"]["A"] = data["params"]["A_sym"]
        data["partial_scores"]["A"] = {"score": 1.0, "weight": 1}
    elif tt == "incorrect":
        data["raw_submitted_answers"]["A"] = "-(" + data["params"]["A_sym"] + ")"
        data["partial_scores"]["A"] = {"score": 0.0, "weight": 1}
    if tt in ("correct", "incorrect"):
        pl.set_weighted_score_data(data)
