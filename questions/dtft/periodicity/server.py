"""Periodicity of discrete-time sinusoids and equivalent frequencies (Lecture 13 §1; slide 6 concept
check; Singer and Munson Fig. 2.4; the cos(17pi n/3) = cos(5pi n/3) = cos(pi n/3) example of the notes).

(a) x[n] built from 2-3 complex exponentials / sinusoids. Periodic? (Yes/No) and the fundamental
    period N0 (integer, left EMPTY when x is not periodic; the box is shown in every variant).
    Families
      sum   2-3 terms with frequencies omega = r*pi (r rational): N0 = lcm of the reduced periods;
            some terms are written 2*pi*(k/N)*n with k/N not in lowest terms (trap: period N).
      prod  a product (e^{ja n} e^{jb n}, cos(a n) e^{jb n}, cos(a n) cos(b n), cos^2): expand into
            exponentials first; N0 = lcm over the expanded frequencies (trap: lcm of the factors).
      irr   one rational term plus one term with omega/2pi irrational (cos(n/3), e^{j2n}, ...): not
            periodic.
    N0 is computed exactly: x[n] = sum_i c_i e^{j 2 pi f_i n} with distinct f_i (mod 1) and c_i != 0
    is N-periodic iff every f_i N is an integer (linear independence of distinct exponentials).
(b) checkbox: which listed sequences equal cos(omega0 n) for all n (cos(a n) = cos(b n) for all n
    iff a = +-b mod 2 pi; phase-shifted options cos(. + 2pi), cos(. + pi), sin(. + pi/2)).
(c) omega1 > pi: the frequency in [0, pi] with the same cosine samples (a multiple of pi).
"""

import math
import random
from fractions import Fraction as F

import prairielearn as pl
from ece310 import fmt, mc
from ece310.w4_dtft import fold, mod2, plain, same_cos, tex_pi, tex_pi_n

RAT = [F(1, 6), F(1, 4), F(1, 3), F(2, 5), F(1, 2), F(3, 5), F(2, 3), F(3, 4), F(4, 5), F(5, 6), F(3, 8),
       F(5, 8), F(1, 5), F(1, 8), F(7, 6), F(5, 4), F(3, 2), F(7, 8)]
IRR = [(1, 2), (1, 3), (2, 3), (3, 4), (2, 1), (3, 2), (1, 1), (2, 5)]       # omega = p/q rad


def period(r):
    """Fundamental period of e^{j r pi n} (r rational): denominator of r/2 reduced mod 1."""
    return (F(r) / 2 - math.floor(F(r) / 2)).denominator


def lcm_list(xs):
    out = 1
    for v in xs:
        out = out * v // math.gcd(out, v)
    return out


# ----------------------------------------------------------------------------- display of terms
def arg_tex(t):
    """LaTeX of the argument omega*n of a term."""
    if t["irr"]:
        p, q = t["irr"]
        num = "n" if p == 1 else f"{p}n"
        return num if q == 1 else r"\tfrac{" + num + "}{" + str(q) + "}"
    if t.get("dft"):
        k, N = t["dft"]
        return rf"2\pi\cdot\tfrac{{{k}}}{{{N}}}\,n"
    return tex_pi_n(t["r"])


def body_tex(t):
    a = arg_tex(t)
    if t["kind"] == "exp":
        return "e^{j" + a + "}"
    return rf"\{t['kind']}\left({a}\right)"


def term_freq_tex(t):
    """omega of the term in LaTeX (for the solution)."""
    if t["irr"]:
        p, q = t["irr"]
        return str(p) if q == 1 else rf"\tfrac{{{p}}}{{{q}}}"
    return tex_pi(t["r"])


def ratio_tex(t):
    """omega/(2 pi) in lowest terms (or the irrational value)."""
    if t["irr"]:
        v = F(t["irr"][0], 2 * t["irr"][1])
        num, den = v.numerator, v.denominator
        return rf"\tfrac{{{num}}}{{\pi}}" if den == 1 else rf"\tfrac{{{num}}}{{{den}\pi}}"
    v = F(t["r"]) / 2
    return fmt.tex_num(v)


# ----------------------------------------------------------------------------- part (a)
def _term(kind, r=None, irr=None, amp=1, dft=None):
    return {"kind": kind, "r": r, "irr": irr, "amp": amp, "dft": dft}


def _dft_form(r):
    """(k, N) with k/N = r/2 but not in lowest terms (N <= 32), or None."""
    v = F(r) / 2
    gs = [g for g in (2, 3, 4) if v.denominator * g <= 32]
    if not gs:
        return None
    g = random.choice(gs)
    return (v.numerator * g, v.denominator * g)


def _gen_sum():
    while True:
        m = random.choice([2, 2, 3])
        rs = random.sample(RAT, m)
        if any(same_cos(a, b) for i, a in enumerate(rs) for b in rs[i + 1:]):
            continue
        Ns = [period(r) for r in rs]
        N0 = lcm_list(Ns)
        if N0 > 120 or N0 == max(Ns):
            continue
        terms = []
        for r in rs:
            kind = random.choice(["cos", "cos", "sin", "exp"])
            amp = random.choice([1, 1, 2, 3, -1, -2, F(1, 2)])
            dft = _dft_form(r) if random.random() < 0.3 else None
            terms.append(_term(kind, r=r, amp=amp, dft=dft))
        prod = 1
        for v in Ns:
            prod *= v
        unreduced = lcm_list([t["dft"][1] if t["dft"] else period(t["r"]) for t in terms])
        if prod == N0 and unreduced == N0:
            continue                                     # no trap
        return terms, N0, {"product": prod, "unreduced": unreduced, "max": max(Ns)}


def _expanded_freqs(kind, a, b):
    """Frequencies (multiples of pi, mod 2) of the product's exponentials, with nonzero coefficients."""
    if kind == "expexp":
        fs = [a + b]
    elif kind == "cosexp":
        fs = [b + a, b - a]
    else:                                                  # cos cos: 1/2 cos((a-b)n) + 1/2 cos((a+b)n)
        fs = [a + b, -(a + b), a - b, b - a]
    return sorted({mod2(f) for f in fs})


def _gen_prod():
    while True:
        kind = random.choice(["expexp", "cosexp", "coscos", "cossq"])
        a, b = random.sample(RAT, 2)
        if kind == "cossq":
            b = a
        if kind in ("cosexp", "coscos", "cossq") and mod2(a) in (0, 1):
            continue
        if kind == "expexp" and mod2(a + b) == 0:
            continue                                       # e^{jan} e^{-jan} = 1: a constant
        fs = _expanded_freqs("coscos" if kind == "cossq" else kind, a, b)
        N0 = lcm_list([period(f) for f in fs])
        factors = lcm_list([period(a), period(b)])
        if N0 > 60 or factors == N0 and random.random() < 0.75:
            continue
        amp = random.choice([1, 1, 2, 3])
        prod = period(a) * period(b)
        return kind, a, b, amp, fs, N0, {"factors": factors, "product": prod, "max": max(period(a), period(b))}


def _gen_irr():
    r = random.choice(RAT)
    p, q = random.choice(IRR)
    k1, k2 = random.choice(["cos", "sin", "exp"]), random.choice(["cos", "sin", "exp"])
    terms = [_term(k1, r=r, amp=random.choice([1, 2, 3, -1])), _term(k2, irr=(p, q), amp=random.choice([1, 2, -1]))]
    random.shuffle(terms)
    return terms


def _sum_tex(terms):
    return fmt.tex_sum((t["amp"], body_tex(t)) for t in terms)


def _prod_tex(kind, a, b, amp):
    pre = "" if amp == 1 else str(amp)
    A = {"r": a, "irr": None, "dft": None}
    B = {"r": b, "irr": None, "dft": None}
    if kind == "expexp":
        return pre + "e^{j" + arg_tex(A) + "}\\,e^{j" + arg_tex(B) + "}"
    if kind == "cosexp":
        return pre + rf"\cos\left({arg_tex(A)}\right)e^{{j{arg_tex(B)}}}"
    if kind == "cossq":
        return pre + rf"\cos^2\left({arg_tex(A)}\right)"
    return pre + rf"\cos\left({arg_tex(A)}\right)\cos\left({arg_tex(B)}\right)"


def exp_tex(f):
    """e^{j f pi n} with the sign in front of j."""
    return "e^{" + ("-" if f < 0 else "") + "j" + tex_pi_n(abs(f)) + "}"


def _prod_expansion_tex(kind, a, b, amp):
    c = F(amp)
    if kind == "expexp":
        return fmt.tex_sum([(c, exp_tex(a + b))])
    if kind == "cosexp":
        parts = []
        for f in (b + a, b - a):
            parts.append((c / 2, exp_tex(f) if f != 0 else ""))
        return fmt.tex_sum(parts)
    if kind == "cossq":
        return fmt.tex_sum([(c / 2, ""), (c / 2, rf"\cos\left({tex_pi_n(2 * a)}\right)")])
    parts = []
    for f in (a - b, a + b):
        parts.append((c / 2, rf"\cos\left({tex_pi_n(abs(f))}\right)" if f != 0 else ""))
    return fmt.tex_sum(parts)


# ----------------------------------------------------------------------------- part (b)
R0 = [F(1, 5), F(2, 5), F(3, 5), F(4, 5), F(1, 3), F(2, 3), F(1, 4), F(3, 4), F(1, 6), F(5, 6), F(3, 8), F(5, 8)]


def _opt_truth(func, a, ph, r0):
    """Is func(a pi n + ph pi) equal to cos(r0 pi n) for every n?"""
    if func == "cos":
        return mod2(ph) == 0 and same_cos(a, r0)
    return mod2(ph - F(1, 2)) == 0 and same_cos(a, r0)          # sin(theta + pi/2) = cos(theta)


def _opt_tex(func, a, ph):
    arg = tex_pi_n(a)
    if ph != 0:
        arg += (" + " if ph > 0 else " - ") + tex_pi(abs(ph))
    return rf"$\{func}\left({arg}\right)$"


def _gen_options(r0):
    true_fr = [r0 + 2, 2 - r0, -r0, r0 - 2, r0 + 4, 4 - r0, -r0 - 2, r0 - 4]
    false_fr = [r0 + 1, 1 - r0, 2 * r0, r0 + 3, 3 - r0, r0 + F(1, 2), 1 + 2 * r0, r0 / 2 + 1]
    false_fr = [f for f in false_fr if not same_cos(f, r0) and f != 0]
    while True:
        nt = random.choice([2, 3, 3, 4])
        opts = [("cos", a, F(0)) for a in random.sample(true_fr, nt)]
        if random.random() < 0.4:
            opts.append(random.choice([("cos", r0, F(2)), ("sin", r0, F(1, 2)), ("sin", 2 - r0, F(1, 2))]))
        opts += [("cos", a, F(0)) for a in random.sample(false_fr, 6 - len(opts))]
        if random.random() < 0.35:                          # swap in a phase trap
            opts[-1] = random.choice([("cos", r0, F(1)), ("sin", r0, F(-1, 2)), ("sin", r0 + 1, F(1, 2))])
        texts = [_opt_tex(*o) for o in opts]
        if len(set(texts)) == len(texts):
            break
    out = []
    for o, t in zip(opts, texts):
        out.append({"text": t, "correct": _opt_truth(o[0], o[1], o[2], r0), "func": o[0], "a": plain(o[1]),
                    "ph": plain(o[2])})
    random.shuffle(out)
    return out


def _why_option(o, r0):
    a, ph = F(o["a"]), F(o["ph"])
    red = mod2(a)
    s = rf"{tex_pi(a)}" + (rf" \equiv {tex_pi(red)} \pmod{{2\pi}}" if red != a else "")
    if o["func"] == "sin":
        lead = (r"$\sin(\theta + \tfrac{\pi}{2}) = \cos\theta$; " if ph == F(1, 2)
                else r"$\sin(\theta - \tfrac{\pi}{2}) = -\cos\theta$: the negative of $x_1[n]$. ")
        if ph != F(1, 2):
            return lead
        return lead + (f"frequency ${s}$, " + ("which is $\\pm\\omega_0$: same samples." if same_cos(a, r0)
                                               else "which is not $\\pm\\omega_0$."))
    if ph == 1:
        return r"$\cos(\theta + \pi) = -\cos\theta$: the negative of $x_1[n]$."
    if ph == 2:
        return r"Adding $2\pi$ to the phase changes nothing."
    if same_cos(a, r0):
        return f"${s}$, which is $\\pm\\omega_0$ modulo $2\\pi$: the same samples."
    if same_cos(a, r0 + 1):
        return f"${s}$, which is $\\pm\\omega_0 + \\pi$: this is $(-1)^n\\cos(\\omega_0 n)$."
    return f"${s}$, which is neither $\\omega_0$ nor $-\\omega_0$ modulo $2\\pi$."


# ----------------------------------------------------------------------------- part (c)
R1 = [F(7, 6), F(4, 3), F(5, 4), F(7, 5), F(5, 3), F(7, 4), F(9, 5), F(11, 6), F(13, 6), F(7, 3), F(9, 4),
      F(12, 5), F(8, 3), F(11, 4), F(14, 5), F(17, 6), F(19, 6), F(10, 3), F(13, 4), F(17, 5), F(11, 3),
      F(15, 4), F(19, 5), F(23, 6), F(25, 6), F(13, 3), F(17, 4), F(22, 5), F(14, 3), F(19, 4), F(29, 6)]


def generate(data):
    p = data["params"]
    fam = random.choices(["sum", "prod", "irr"], weights=[45, 25, 30])[0]
    p["family"] = fam
    rows = []
    if fam == "sum":
        terms, N0, traps = _gen_sum()
        p["x_tex"] = _sum_tex(terms)
        for t in terms:
            rows.append({"term_tex": body_tex(t), "w_tex": term_freq_tex(t), "ratio_tex": ratio_tex(t),
                         "N": str(period(t["r"])),
                         "note_html": (rf"(written as $2\pi\cdot\tfrac{{{t['dft'][0]}}}{{{t['dft'][1]}}}$: reduce the fraction first)"
                                  if t["dft"] else "")})
        p["is_periodic"] = True
        p["terms"] = [{"kind": t["kind"], "amp": plain(t["amp"]), "r": plain(t["r"]), "irr": None} for t in terms]
        p["N_list_tex"] = ", ".join(r["N"] for r in rows)
        p["traps"] = traps
    elif fam == "prod":
        kind, a, b, amp, fs, N0, traps = _gen_prod()
        p["x_tex"] = _prod_tex(kind, a, b, amp)
        p["expansion_tex"] = _prod_expansion_tex(kind, a, b, amp)
        for f in fs:
            if f == 0:
                rows.append({"term_tex": r"\text{constant}", "w_tex": "0", "ratio_tex": "0", "N": "1", "note_html": ""})
            else:
                rows.append({"term_tex": "e^{j" + tex_pi_n(f) + "}", "w_tex": tex_pi(f),
                             "ratio_tex": fmt.tex_num(f / 2), "N": str(period(f)), "note_html": ""})
        p["is_periodic"] = True
        p["prod"] = {"kind": kind, "a": plain(a), "b": plain(b), "amp": amp}
        p["N_list_tex"] = ", ".join(r["N"] for r in rows)
        p["factor_note_html"] = (
            rf"(The factor alone has period ${period(a)}$; squaring creates the frequencies $0$ and $2\omega$, so the "
            rf"fundamental period of the square can be shorter.)" if kind == "cossq" else
            rf"(The factors alone have periods ${period(a)}$ and ${period(b)}$; their least common multiple is <i>a</i> "
            rf"period, but the product's frequencies combine, and the fundamental period can be shorter.)")
        p["traps"] = traps
    else:
        terms = _gen_irr()
        p["x_tex"] = _sum_tex(terms)
        for t in terms:
            rows.append({"term_tex": body_tex(t), "w_tex": term_freq_tex(t), "ratio_tex": ratio_tex(t),
                         "N": (str(period(t["r"])) if not t["irr"] else r"\text{none}"),
                         "note_html": "(irrational)" if t["irr"] else ""})
        irr_t = next(t for t in terms if t["irr"])
        rat_t = next(t for t in terms if not t["irr"])
        p["irr_w_tex"] = term_freq_tex(irr_t)
        p["irr_ratio_tex"] = ratio_tex(irr_t)
        p["is_periodic"] = False
        p["terms"] = [{"kind": t["kind"], "amp": plain(t["amp"]), "r": (plain(t["r"]) if t["r"] is not None else None),
                       "irr": list(t["irr"]) if t["irr"] else None} for t in terms]
        p["traps"] = {"rational_period": period(rat_t["r"]),
                      "two_pi_over_w": round(2 * math.pi * irr_t["irr"][1] / irr_t["irr"][0])}
        N0 = None
    p["rows"] = rows
    p["is_prod"] = fam == "prod"
    p["is_sum"] = fam == "sum"
    p["is_irr"] = fam == "irr"
    p["periodic_choices"] = mc.yes_no(N0 is not None)
    p["N0_tex"] = str(N0) if N0 is not None else ""
    data["correct_answers"]["N0"] = N0 if N0 is not None else ""

    # (b)
    r0 = random.choice(R0)
    opts = _gen_options(r0)
    p["r0"] = plain(r0)
    p["w0_tex"] = tex_pi(r0)
    p["cos0_tex"] = rf"\cos\left({tex_pi_n(r0)}\right)"
    p["same_choices"] = [{"text": o["text"], "correct": o["correct"]} for o in opts]
    p["same_why"] = [{"text": o["text"], "verdict": "same" if o["correct"] else "different",
                      "why_html": _why_option(o, r0)} for o in opts]
    p["same_opts"] = [{"func": o["func"], "a": o["a"], "ph": o["ph"], "correct": o["correct"]} for o in opts]

    # (c)
    r1 = random.choice(R1)
    ans = fold(r1)
    red = mod2(r1)
    p["r1"] = plain(r1)
    p["w1_tex"] = tex_pi(r1)
    p["w1_n_tex"] = tex_pi_n(r1)
    k = (r1 - red) / 2
    p["sub_tex"] = (rf"{tex_pi(r1)} - {tex_pi(2 * k)} = {tex_pi(red)}" if k else tex_pi(red))
    p["needs_flip"] = red > 1
    p["flip_tex"] = rf"2\pi - {tex_pi(red)} = {tex_pi(2 - red)}"
    p["w2_tex"] = tex_pi(ans)
    p["w2_plain"] = plain(ans)
    data["correct_answers"]["w2"] = float(ans)


def grade(data):
    """Non-periodic x: the N0 box must stay empty, and the empty box only counts with the answer No."""
    if data["params"].get("is_periodic"):
        return
    ps = data["partial_scores"]
    if "N0" not in ps:
        return
    sub = data["submitted_answers"].get("N0")
    blank = sub is None or (isinstance(sub, str) and sub.strip() == "")
    per_score = (ps.get("periodic") or {}).get("score") or 0
    ps["N0"]["score"] = per_score if blank else 0
    if blank and not per_score:
        ps["N0"]["feedback"] = "An empty box only counts together with the answer No."
    pl.set_weighted_score_data(data)
