"""Geometric sums and roots of unity, "two tools the z-transform will need" (Lecture 2 section 4 on
the study site; HW1 #4 roots of unity; HW3 #1 z-transforms as geometric sums). Every variant asks:

  (a) one numeric sum (exact rational answer), from one of five families
        fin    sum_{n=n1}^{n2} a^n                (n1 != 0, 3-6 terms)
        inf    sum_{n=n1}^{inf} a^n, |a| < 1      (n1 != 0)
        left   sum_{n=-inf}^{n2} a^n, |a| > 1
        two    sum_n a^{|n|}, or sum_n [a^n u[n] + b^n u[-n]] / [... u[-n-1]] (HW3 #1(c), (d))
        roots  sum_{n=0}^{N-1} e^{j 2 pi k n / N}  (N or 0; the exponent is shown reduced)
  (b) a closed form in n (pl-symbolic-input, n an integer), base different from (a)
        G1 sum_{k=0}^{n} a^k (n >= 0), G2 sum_{k=n}^{inf} a^k, G3 sum_{k=0}^{n} a^k b^{n-k} (n >= 0),
        G4 sum_{k=-inf}^{n} a^k, G5 sum_{k=m}^{n} a^k (n >= m)
  (c) for which a a parameterized series converges (six fixed options in terms of |a| and c).
"""

import random
from fractions import Fraction as Fr

from ece310 import answers, fmt
from ece310.w1_signals import plain_frac

SMALL = [Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(-1, 3), Fr(2, 3), Fr(-2, 3), Fr(1, 4), Fr(-1, 4), Fr(3, 4), Fr(-3, 4)]
BIG = [Fr(2), Fr(-2), Fr(3), Fr(-3), Fr(3, 2), Fr(-3, 2), Fr(4)]


def T(x, paren=False):
    return fmt.tex_num(x, paren=paren)


def P(a, e="n"):
    return fmt.tex_pow(a, e)


def D(num, den):
    return fmt.tex_frac(num, den)


def _nice(S, maxden=64, maxnum=600):
    return S != 0 and S.denominator <= maxden and abs(S.numerator) <= maxnum


def _close_frac(x, y):
    return abs(float(x) - float(y)) > 0.02 * max(1.0, abs(float(y)))


# ----------------------------------------------------------------------------- part (a)
def _a_fin():
    while True:
        a = random.choice(SMALL + BIG)
        n1 = random.choice([-2, -1, 1, 2])
        N = random.randint(3, 6)
        n2 = n1 + N - 1
        S = sum(a ** n for n in range(n1, n2 + 1))
        if not _nice(S) or abs(a ** n1) > 16 or abs(a ** (n2 + 1)) > 81:
            continue
        steps = (rf"\sum_{{n={n1}}}^{{{n2}}} {P(a)} = \dfrac{{a^{{{n1}}} - a^{{{n2 + 1}}}}}{{1 - a}}"
                 rf" = {D(T(a ** n1) + ' - ' + T(a ** (n2 + 1), paren=True), T(1 - a))} = {T(S)}")
        words = (rf"A finite geometric sum: first term minus first omitted term, over one minus the ratio "
                 rf"(or pull out $a^{{{n1}}}$: $a^{{{n1}}}\big(1 + a + \dots + a^{{{N - 1}}}\big) = "
                 rf"a^{{{n1}}}\,\tfrac{{1 - a^{{{N}}}}}{{1 - a}}$). Here $a = {T(a)}$ and there are ${N}$ terms.")
        return a, S, rf"\sum_{{n={n1}}}^{{{n2}}} {P(a)}", steps, words, {"fam": "fin", "a": str(a), "lo": n1, "hi": n2}


def _a_inf():
    while True:
        a = random.choice(SMALL)
        n1 = random.choice([-2, -1, 1, 2, 3])
        S = a ** n1 / (1 - a)
        if _nice(S):
            break
    steps = (rf"\sum_{{n={n1}}}^{{\infty}} {P(a)} = \dfrac{{a^{{{n1}}}}}{{1 - a}} = {D(T(a ** n1), T(1 - a))} = {T(S)}")
    words = (rf"$\lvert a\rvert = {T(abs(a))} \lt 1$, so the series converges. Pull out the first term: "
             rf"$\sum_{{n={n1}}}^{{\infty}} a^n = a^{{{n1}}}\sum_{{m=0}}^{{\infty}} a^m$ with $m = n - ({n1})$.")
    return a, S, rf"\sum_{{n={n1}}}^{{\infty}} {P(a)}", steps, words, {"fam": "inf", "a": str(a), "lo": n1, "hi": None}


def _a_left():
    while True:
        a = random.choice(BIG)
        n2 = random.choice([-2, -1, 0, 1, 2])
        S = a ** (n2 + 1) / (a - 1)
        if _nice(S):
            break
    steps = (rf"\sum_{{n=-\infty}}^{{{n2}}} {P(a)} = \sum_{{m={-n2}}}^{{\infty}} {P(1 / a, 'm')}"
             rf" = \dfrac{{(1/a)^{{{-n2}}}}}{{1 - 1/a}} = {D(T(a ** n2), T(1 - 1 / a))} = {T(S)}")
    words = (rf"Substitute $m = -n$: the sum runs over $m \ge {-n2}$ with ratio $1/a = {T(1 / a)}$, and "
             rf"$\lvert 1/a\rvert \lt 1$ because $\lvert a\rvert = {T(abs(a))} \gt 1$.")
    return a, S, rf"\sum_{{n=-\infty}}^{{{n2}}} {P(a)}", steps, words, {"fam": "left", "a": str(a), "lo": None, "hi": n2}


def _a_two():
    if random.random() < 0.5:
        a = random.choice(SMALL)
        S = (1 + a) / (1 - a)
        tex = rf"\sum_{{n=-\infty}}^{{\infty}} {P(a, r'\lvert n\rvert')}"
        steps = (rf"{tex} = \underbrace{{\sum_{{n=0}}^{{\infty}} a^n}}_{{n \ge 0}} + \underbrace{{\sum_{{m=1}}^{{\infty}} a^m}}_{{n = -m \le -1}}"
                 rf" = \dfrac{{1}}{{1 - a}} + \dfrac{{a}}{{1 - a}} = \dfrac{{1 + a}}{{1 - a}} = {D(T(1 + a), T(1 - a))} = {T(S)}")
        words = (rf"Split at $n = 0$ and count the $n = 0$ term once; both halves have ratio $a = {T(a)}$, "
                 rf"$\lvert a\rvert \lt 1$ (HW3 #1(d): $(\tfrac14)^{{\lvert n\rvert}}$).")
        return a, S, tex, steps, words, {"fam": "two_abs", "a": str(a), "lo": None, "hi": None}
    while True:
        a, b = random.choice(SMALL), random.choice(BIG)
        incl = random.random() < 0.6           # b^n u[-n] includes n = 0; b^n u[-n-1] does not
        left = (b / (b - 1)) if incl else (1 / (b - 1))
        S = 1 / (1 - a) + left
        if _nice(S):
            break
    uarg = "-n" if incl else "-n - 1"
    tex = rf"\sum_{{n=-\infty}}^{{\infty}}\left({P(a)}\,u[n] + {P(b)}\,u[{uarg}]\right)"
    top = 0 if incl else -1
    steps = (rf"{tex} = \sum_{{n=0}}^{{\infty}} {P(a)} + \sum_{{n=-\infty}}^{{{top}}} {P(b)}"
             rf" = \dfrac{{1}}{{1 - a}} + " + (r"\dfrac{1}{1 - 1/b}" if incl else r"\dfrac{1/b}{1 - 1/b}")
             + rf" = {T(1 / (1 - a))} + {T(left, paren=True)} = {T(S)}")
    words = (rf"Split the sum: $u[n]$ keeps $n \ge 0$ and $u[{uarg}]$ keeps $n \le {top}$, so the $n = 0$ sample "
             + (r"is counted in both parts, which is right here: $x[0] = 1 + 1 = 2$" if incl
                else r"belongs to the first part only")
             + rf". The right part needs $\lvert a\rvert = {T(abs(a))} \lt 1$, the left part $\lvert b\rvert = {T(abs(b))} \gt 1$ "
             rf"(HW3 #1(c): $3^n u[-n] + 2^{{-n}}u[n]$).")
    return a, S, tex, steps, words, {"fam": "two_pair", "a": str(a), "b": str(b), "incl": incl, "lo": None, "hi": None}


def _ordinal(N):
    """3rd, 4th, 12th."""
    suf = "th" if 10 <= N % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(N % 10, "th")
    return f"{N}{suf}"


def _a_roots():
    N = random.choice([3, 4, 5, 6, 8, 10, 12])
    if random.random() < 0.35:
        k = N * random.choice([1, -1, 2])
    else:
        k = random.choice([kk for kk in range(-N + 1, 2 * N) if kk % N != 0])
    r = Fr(2 * k, N)
    p_, q_ = abs(r.numerator), r.denominator
    ex = ("-" if r < 0 else "") + "j" + ("" if p_ == 1 else str(p_)) + r"\pi n" + ("" if q_ == 1 else f"/{q_}")
    tex = rf"\sum_{{n=0}}^{{{N - 1}}} e^{{{ex}}}"
    one = r.denominator == 1 and r.numerator % 2 == 0
    S = Fr(N) if one else Fr(0)
    wtex = "e^{" + ex.replace(" n", "") + "}"
    if one:
        steps = rf"{tex} = \sum_{{n=0}}^{{{N - 1}}} 1 = {N}"
        words = (rf"The ratio is $w = {wtex} = 1$ (a multiple of $2\pi$), so every one of the ${N}$ terms equals $1$; "
                 r"the formula $\frac{1 - w^N}{1 - w}$ does not apply when $w = 1$.")
    else:
        steps = (rf"{tex} = \dfrac{{1 - w^{{{N}}}}}{{1 - w}},\quad w = {wtex} \ne 1,\quad w^{{{N}}} = e^{{j2\pi\cdot({k})}} = 1"
                 rf"\ \Rightarrow\ \text{{sum}} = 0")
        words = (rf"Every term $w^n$ is one of the {_ordinal(N)} roots of unity: the ratio $w = {wtex}$ is not $1$, "
                 rf"but $w^{{{N}}} = 1$, so the finite geometric sum $\frac{{1 - w^{{{N}}}}}{{1 - w}}$ is $0$.")
    return None, S, tex, steps, words, {"fam": "roots", "N": N, "k": k}


# ----------------------------------------------------------------------------- part (b)
BASES_B = [Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(-1, 3), Fr(2, 3), Fr(1, 4), Fr(3, 4), Fr(2), Fr(-2), Fr(3)]


def _b_closed(avoid):
    fam = random.choice(["G1", "G2", "G3", "G4", "G5"])
    pool = [x for x in BASES_B if x not in avoid]
    if fam == "G1":
        a = random.choice(pool)
        c0, c1 = 1 / (1 - a), -a / (1 - a)
        sum_tex, valid = rf"\sum_{{k=0}}^{{n}} {P(a, 'k')}", r"n \ge 0"
        form_tex = rf"\dfrac{{1 - a^{{n+1}}}}{{1 - a}} = {D('1 - ' + P(a, 'n+1'), T(1 - a))}"
        spec = {"fam": fam, "a": str(a)}
        terms = [(c0, ""), (c1, a)]
    elif fam == "G2":
        a = random.choice([x for x in pool if abs(x) < 1])
        sum_tex, valid = rf"\sum_{{k=n}}^{{\infty}} {P(a, 'k')}", r"\text{all } n"
        form_tex = rf"\dfrac{{a^{{n}}}}{{1 - a}} = {D(P(a), T(1 - a))}"
        spec = {"fam": fam, "a": str(a)}
        terms = [(1 / (1 - a), a)]
    elif fam == "G3":
        a, b = random.sample(pool, 2)
        sum_tex, valid = rf"\sum_{{k=0}}^{{n}} {P(a, 'k')}\,{P(b, 'n-k')}", r"n \ge 0"
        form_tex = (rf"b^n\sum_{{k=0}}^{{n}}\left(\tfrac{{a}}{{b}}\right)^k = \dfrac{{a^{{n+1}} - b^{{n+1}}}}{{a - b}}"
                    rf" = {D(P(a, 'n+1') + ' - ' + P(b, 'n+1'), T(a - b))}")
        spec = {"fam": fam, "a": str(a), "b": str(b)}
        terms = [(a / (a - b), a), (-b / (a - b), b)]
    elif fam == "G4":
        a = random.choice([x for x in pool if abs(x) > 1])
        sum_tex, valid = rf"\sum_{{k=-\infty}}^{{n}} {P(a, 'k')}", r"\text{all } n"
        form_tex = rf"\dfrac{{a^{{n}}}}{{1 - 1/a}} = \dfrac{{a^{{n+1}}}}{{a - 1}} = {D(P(a, 'n+1'), T(a - 1))}"
        spec = {"fam": fam, "a": str(a)}
        terms = [(a / (a - 1), a)]
    else:
        a = random.choice(pool)
        m = random.choice([1, 2, -1])
        sum_tex, valid = rf"\sum_{{k={m}}}^{{n}} {P(a, 'k')}", rf"n \ge {m}"
        form_tex = (rf"\dfrac{{a^{{{m}}} - a^{{n+1}}}}{{1 - a}} = {D(T(a ** m) + ' - ' + P(a, 'n+1'), T(1 - a))}")
        spec = {"fam": fam, "a": str(a), "m": m}
        terms = [(a ** m / (1 - a), ""), (-a / (1 - a), a)]
    sym = fmt.sym_sum((c, fmt.sym_pow(base) if base != "" else "") for c, base in terms)
    tex = fmt.tex_sum((c, P(base) if base != "" else "") for c, base in terms)
    return fam, sum_tex, valid, form_tex, tex, sym, spec


# ----------------------------------------------------------------------------- part (c)
OPTS = ["lt_inv", "lt_c", "gt_inv", "gt_c", "all", "none"]


def _opt_text(key, c):
    z = r"\lvert a\rvert"
    return {"lt_inv": rf"${z} \lt \tfrac{{1}}{{{c}}}$", "lt_c": rf"${z} \lt {c}$",
            "gt_inv": rf"${z} \gt \tfrac{{1}}{{{c}}}$", "gt_c": rf"${z} \gt {c}$",
            "all": r"for every $a \ne 0$", "none": r"for no value of $a$"}[key]


# (base, exponent, lo, hi, correct option, reason)
SERIES = [
    ("a/c", "n", 0, None, "lt_c"),
    ("a/c", "n", None, 0, "gt_c"),
    ("a/c", "abs", None, None, "lt_c"),
    ("a/c", "n", None, None, "none"),
    ("ca", "-n", 1, None, "gt_inv"),
    ("ca", "n", 0, None, "lt_inv"),
    ("ca", "n", None, -1, "gt_inv"),
    ("a/c", "n", -3, 4, "all"),
    ("ca", "-abs", None, None, "gt_inv"),
    ("a/c", "-n", 0, None, "gt_c"),
    ("ca", "abs", None, None, "lt_inv"),
]


def _c_series():
    base, ex, lo, hi, ans = random.choice(SERIES)
    c = random.choice([2, 3, 4])
    btex = rf"\left(\tfrac{{a}}{{{c}}}\right)" if base == "a/c" else rf"\left({c}a\right)"
    etex = {"n": "n", "-n": "-n", "abs": r"\lvert n\rvert", "-abs": r"-\lvert n\rvert"}[ex]
    lo_t = r"-\infty" if lo is None else str(lo)
    hi_t = r"\infty" if hi is None else str(hi)
    tex = rf"\sum_{{n={lo_t}}}^{{{hi_t}}} {btex}^{{{etex}}}"
    r = rf"\tfrac{{\lvert a\rvert}}{{{c}}}" if base == "a/c" else rf"{c}\lvert a\rvert"
    if ans == "all":
        why = r"A sum of finitely many finite terms: it is finite for every $a \ne 0$ (no convergence question at all)."
    elif ans == "none":
        why = (rf"The terms with $n \to \infty$ need ${r} \lt 1$ and the terms with $n \to -\infty$ need ${r} \gt 1$; "
               r"both cannot hold, so the sum converges for no $a$ (like $x[n] = 1$ for all $n$, which has no z-transform).")
    else:
        need_small = (ex == "n" and hi is None) or ex == "abs"
        if base == "a/c":
            rr = r if need_small else rf"\tfrac{{{c}}}{{\lvert a\rvert}}"
        else:
            rr = r if need_small else rf"\tfrac{{1}}{{{c}\lvert a\rvert}}"
        res = _opt_text(ans, c).strip("$")
        why = (r"Along the infinite end" + (r"s (both sides)" if ex in ("abs", "-abs") else "")
               + rf" the terms have magnitude $\left({rr}\right)^{{m}}$ with $m = 0, 1, 2, \dots$ (up to a constant factor), "
               + rf"so the geometric series converges iff ${rr} \lt 1$, i.e. ${res}$.")
    choices = [{"text": _opt_text(k, c), "correct": k == ans} for k in OPTS]
    spec = {"base": base, "exp": ex, "lo": lo, "hi": hi, "c": c}
    return tex, why, choices, spec


# ----------------------------------------------------------------------------- generate
def generate(data):
    p = data["params"]
    fam = random.choice(["fin", "inf", "left", "two", "roots"])
    a, S, tex, steps, words, spec = {"fin": _a_fin, "inf": _a_inf, "left": _a_left, "two": _a_two,
                                     "roots": _a_roots}[fam]()
    p["sumA_tex"], p["sumA_steps_tex"], p["sumA_words_html"], p["sumA_spec"] = tex, steps, words, spec
    p["sumA_frac"] = plain_frac(S)
    avoid = set()
    if a is not None:
        avoid = {a, 1 / a}
    if spec.get("b"):
        avoid |= {Fr(spec["b"]), 1 / Fr(spec["b"])}
    bfam, sum_tex, valid, form_tex, ctex, sym, bspec = _b_closed(avoid)
    p["sumB_tex"], p["sumB_valid_tex"], p["sumB_form_tex"], p["sumB_closed_tex"] = sum_tex, valid, form_tex, ctex
    p["sumB_spec"] = bspec
    p["sumB_sym"] = sym
    ctex_c, why, choices, cspec = _c_series()
    p["conv_tex"], p["conv_why_html"], p["conv_choices"], p["conv_spec"] = ctex_c, why, choices, cspec
    data["correct_answers"]["S"] = float(S)
    data["correct_answers"]["closed"] = answers.symbolic(sym)
