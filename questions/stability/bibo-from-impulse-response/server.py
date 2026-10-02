"""BIBO stability straight from h[n] (Lecture 4: stable <=> sum |h[n]| < infinity; HW4 #2, #5;
FA2025 #1a, FA2024 #1b, FA2019 #1c "bounded is not enough").

Families: A a^n u[n-k], A a^|n|, A n a^n u[n], A/n u[n-1], cos / r^n cos(pi n/2) / r^n sin(pi n/2),
finite sequences and pulses, two-sided A a^n u[n] + B b^n u[-n-1], delta + A a^n u[n-1] (HW4 #5 h2),
C u[n] + D a^n u[n] (HW4 #5 h1), A b^n u[-n]. The student answers Yes/No and, only if the system is
stable, enters sum |h[n]| exactly (the box is always shown and must be left blank for an unstable
system, so the presence of the box does not give the answer away)."""

import random
from fractions import Fraction as F

import prairielearn as pl
from ece310 import fmt, mc, zt

STAB_A = [F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(2, 3), F(-2, 3), F(3, 4), F(-3, 4), F(1, 4), F(-1, 4)]
UNSTAB_A = [F(2), F(-2), F(3, 2), F(-3, 2), F(1), F(-1), F(5, 4), F(3)]
BIG_B = [F(2), F(-2), F(3), F(-3), F(3, 2), F(-3, 2), F(4), F(-4), F(5, 2), F(4, 3), F(5, 4)]
COEF = [F(1), F(2), F(3), F(4), F(-1), F(-2), F(-3), F(1, 2), F(3, 2)]


def _u(k):
    return f"u[{zt.u_arg_right(k)}]"


def _pow(a, e="n"):
    return fmt.tex_pow(a, e)


def _term(A, body):
    return fmt.tex_sum([(A, body)])


def _num(x):
    return fmt.tex_num(x)


def _abs(x):
    return abs(x)


def fam_rexp(stable):
    A, k = random.choice(COEF), random.choice([0, 0, 1, 2, 3, -1, -2])
    a = random.choice(STAB_A if stable else UNSTAB_A)
    body = (_pow(a) + r"\," if a != 1 else "") + _u(k)
    h_tex = _term(A, body)
    spec = [{"t": "rexp", "A": str(A), "a": str(a), "k": k}]
    if not stable:
        why = (fr"$\lvert h[n]\rvert = {_num(abs(A))}\cdot{_pow(abs(a))}$ for $n \ge {k}$ does not decay "
               r"(the base has magnitude $\ge 1$), so the terms do not even go to $0$ and "
               r"$\sum_n \lvert h[n]\rvert = \infty$.")
        return h_tex, None, None, [], spec, why
    m = abs(a)
    S = abs(A) * m**k / (1 - m)
    signed = A * a**k / (1 - a)
    wrong = [abs(A) / (1 - m)] if k != 0 else []
    why = (fr"$h[n]$ is nonzero for $n \ge {k}$ with $\lvert h[n]\rvert = {_num(abs(A))}\left({_num(m)}\right)^{{n}}$. "
           fr"A geometric series with ratio ${_num(m)} \lt 1$ that starts at $n = {k}$:"
           fr"$$\sum_n \lvert h[n]\rvert = {_num(abs(A))}\cdot\frac{{\left({_num(m)}\right)^{{{k}}}}}{{1 - {_num(m)}}} = {_num(S)}.$$")
    return h_tex, S, signed, wrong, spec, why


def fam_abs(stable):
    A = random.choice(COEF)
    a = random.choice(STAB_A if stable else [F(2), F(3, 2), F(-2), F(-1), F(4, 3)])
    h_tex = _term(A, _pow(a, r"\lvert n\rvert"))
    spec = [{"t": "abs", "A": str(A), "a": str(a)}]
    if not stable:
        why = (r"$\lvert h[n]\rvert = " + _num(abs(A)) + r"\cdot" + _pow(abs(a), r"\lvert n\rvert") + r"$ does not decay as "
               r"$n \to \pm\infty$, so $\sum_n \lvert h[n]\rvert = \infty$.")
        return h_tex, None, None, [], spec, why
    m = abs(a)
    S = abs(A) * (1 + m) / (1 - m)
    signed = A * (1 + a) / (1 - a)
    wrong = [2 * abs(A) / (1 - m)]
    why = (r"Split the two-sided sum at $n = 0$ (count $n = 0$ once):"
           fr"$$\sum_n \lvert h[n]\rvert = {_num(abs(A))}\Big(1 + 2\sum_{{n=1}}^{{\infty}} \left({_num(m)}\right)^{{n}}\Big)"
           fr" = {_num(abs(A))}\Big(1 + \frac{{2\cdot{_num(m)}}}{{1 - {_num(m)}}}\Big) = {_num(S)}.$$")
    return h_tex, S, signed, wrong, spec, why


def fam_nexp(stable):
    A = random.choice([F(1), F(2), F(3), F(-1), F(-2), F(4), F(1, 2)])
    a = random.choice([F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(2, 3), F(-2, 3), F(1, 4), F(-1, 4), F(3, 4)] if stable
                      else [F(1), F(-1), F(2), F(-2), F(3, 2), F(5, 4), F(-3, 2)])
    body = r"n\," + ((_pow(a) + r"\,") if a != 1 else "") + "u[n]"
    h_tex = _term(A, body)
    spec = [{"t": "nexp", "A": str(A), "a": str(a)}]
    if not stable:
        why = (fr"$\lvert h[n]\rvert = {_num(abs(A))}\,n\cdot{_pow(abs(a))}$ grows without bound, so the sum diverges: "
               r"$\sum_n \lvert h[n]\rvert = \infty$.")
        return h_tex, None, None, [], spec, why
    m = abs(a)
    S = abs(A) * m / (1 - m) ** 2
    signed = A * a / (1 - a) ** 2
    wrong = [abs(A) / (1 - m) ** 2]
    why = (r"Use $\sum_{n=0}^{\infty} n\,r^{n} = \dfrac{r}{(1-r)^2}$ for $\lvert r\rvert \lt 1$ "
           r"(the derivative of the geometric series) with $r = " + _num(m) + r"$:"
           fr"$$\sum_n \lvert h[n]\rvert = {_num(abs(A))}\cdot\frac{{{_num(m)}}}{{\left(1 - {_num(m)}\right)^{{2}}}} = {_num(S)}.$$")
    return h_tex, S, signed, wrong, spec, why


def fam_harm(stable):
    A = random.choice([1, 2, 3, 4, 5, -1, -2])
    c = random.choice([0, 1, 2, 3])
    den = "n" if c == 0 else f"n+{c}"
    start = 1 if c == 0 else 0
    h_tex = (r"-" if A < 0 else "") + fr"\dfrac{{{abs(A)}}}{{{den}}}\,{_u(start)}"
    spec = [{"t": "harm", "A": str(A), "c": c, "s": start}]
    why = (fr"$h[n]$ is bounded (it decays to $0$), but only like $1/n$: $\sum_n \lvert h[n]\rvert = {abs(A)}"
           fr"\sum_{{m={start + c}}}^{{\infty}} \frac{{1}}{{m}} = \infty$ (harmonic series). Decaying to zero is not enough.")
    return h_tex, None, None, [], spec, why


_W = {(1, 2): r"\tfrac{\pi}{2}", (1, 3): r"\tfrac{\pi}{3}", (1, 4): r"\tfrac{\pi}{4}", (2, 3): r"\tfrac{2\pi}{3}", (1, 6): r"\tfrac{\pi}{6}"}


def fam_cos(stable):
    A = random.choice([F(1), F(2), F(3), F(4), F(-1), F(-2), F(1, 2)])
    fn = random.choice(["cos", "sin"])
    if not stable:
        w = random.choice(list(_W))
        h_tex = _term(A, fr"\{fn}\!\left({_W[w]}n\right)u[n]")
        spec = [{"t": "trig", "A": str(A), "r": "1", "w": list(w), "fn": fn}]
        why = (fr"$h[n]$ is bounded by ${_num(abs(A))}$ but periodic: it never decays, so infinitely many terms have "
               r"$\lvert h[n]\rvert$ bounded away from $0$ and $\sum_n \lvert h[n]\rvert = \infty$ (bounded is not enough).")
        return h_tex, None, None, [], spec, why
    r = random.choice([F(1, 2), F(1, 3), F(2, 3), F(3, 4), F(1, 4)])
    h_tex = _term(A, _pow(r) + fr"\{fn}\!\left(\tfrac{{\pi}}{{2}}n\right)u[n]")
    spec = [{"t": "trig", "A": str(A), "r": str(r), "w": [1, 2], "fn": fn}]
    if fn == "cos":
        S = abs(A) / (1 - r * r)
        signed = A / (1 + r * r)
        why = (r"$\cos(\tfrac{\pi}{2}n)$ is $1, 0, -1, 0, 1, \dots$: only even $n = 2m$ contribute, with "
               fr"$\lvert h[2m]\rvert = {_num(abs(A))}\left({_num(r * r)}\right)^{{m}}$."
               fr"$$\sum_n \lvert h[n]\rvert = {_num(abs(A))}\sum_{{m=0}}^{{\infty}}\left({_num(r * r)}\right)^{{m}}"
               fr" = \frac{{{_num(abs(A))}}}{{1 - {_num(r * r)}}} = {_num(S)}.$$")
    else:
        S = abs(A) * r / (1 - r * r)
        signed = A * r / (1 + r * r)
        why = (r"$\sin(\tfrac{\pi}{2}n)$ is $0, 1, 0, -1, 0, \dots$: only odd $n = 2m+1$ contribute, with "
               fr"$\lvert h[2m+1]\rvert = {_num(abs(A))}\cdot{_num(r)}\left({_num(r * r)}\right)^{{m}}$."
               fr"$$\sum_n \lvert h[n]\rvert = {_num(abs(A) * r)}\sum_{{m=0}}^{{\infty}}\left({_num(r * r)}\right)^{{m}}"
               fr" = \frac{{{_num(abs(A) * r)}}}{{1 - {_num(r * r)}}} = {_num(S)}.$$")
    wrong = [abs(A) / (1 - r)]
    return h_tex, S, signed, wrong, spec, why


def fam_finite(stable):
    if random.random() < 0.3:
        M, N = random.randint(0, 3), random.randint(2, 5)
        A = random.choice([1, 2, 3, -1, -2])
        lead = "" if A == 1 else ("-" if A == -1 else str(A))
        h_tex = lead + fr"\big(u[{zt.u_arg_right(-M)}] - u[n-{N}]\big)"
        vals, start = [F(A)] * (M + N), -M
        S = abs(A) * (M + N)
        why = (fr"$h[n] = {A}$ for ${-M} \le n \le {N - 1}$ ({M + N} samples) and $0$ otherwise. "
               fr"A finite-length $h$ is always absolutely summable: $\sum_n \lvert h[n]\rvert = {M + N}\cdot{abs(A)} = {_num(S)}$.")
    else:
        while True:
            L = random.randint(3, 6)
            vals = [F(random.randint(-3, 3)) for _ in range(L)]
            vals[0] = F(random.choice([-3, -2, -1, 1, 2, 3]))
            vals[-1] = F(random.choice([-3, -2, -1, 1, 2, 3]))
            if any(v < 0 for v in vals):
                break
        start = random.randint(-2, 2)
        h_tex = fmt.tex_seq(vals, start)
        S = sum(abs(v) for v in vals)
        why = (r"A finite-length $h$ is always absolutely summable (FIR $\Rightarrow$ stable): "
               fr"$\sum_n \lvert h[n]\rvert = {' + '.join(_num(abs(v)) for v in vals if v != 0)} = {_num(S)}$.")
    spec = [{"t": "fin", "values": [str(v) for v in vals], "start": start}]
    signed = sum(vals)
    return h_tex, S, signed, [], spec, why


def fam_two(stable):
    A, B = random.choice([F(1), F(2), F(3), F(-1), F(-2)]), random.choice([F(1), F(2), F(-1), F(3), F(-2)])
    if stable:
        a, b = random.choice(STAB_A), random.choice(BIG_B)
    elif random.random() < 0.6:
        a, b = random.choice(STAB_A), random.choice([F(1, 2), F(-1, 2), F(1, 3), F(2, 3), F(3, 4)])
    else:
        a, b = random.choice([F(2), F(-2), F(3, 2), F(1)]), random.choice(BIG_B)
    t1 = _term(A, (_pow(a) + r"\," if a != 1 else "") + "u[n]")
    t2 = _term(B, _pow(b) + r"\,u[-n-1]")
    h_tex = t1 + (" " if t2.startswith("-") else " + ") + t2
    spec = [{"t": "rexp", "A": str(A), "a": str(a), "k": 0}, {"t": "lexp", "B": str(B), "b": str(b), "m": -1}]
    if not stable:
        if abs(b) < 1:
            why = (fr"The left-sided part ${t2}$ blows up as $n \to -\infty$: for $n = -m$, "
                   fr"$\lvert h[-m]\rvert = {_num(abs(B))}\cdot{_pow(1 / abs(b), 'm')}$ with ${_num(1 / abs(b))} \gt 1$. "
                   r"A left-sided exponential needs $\lvert b\rvert \gt 1$. So $\sum_n \lvert h[n]\rvert = \infty$.")
        else:
            why = (fr"The right-sided part ${t1}$ does not decay as $n \to \infty$ (base of magnitude $\ge 1$), "
                   r"so $\sum_n \lvert h[n]\rvert = \infty$.")
        return h_tex, None, None, [], spec, why
    ma, mb = abs(a), abs(b)
    S = abs(A) / (1 - ma) + abs(B) / (mb - 1)
    signed = A / (1 - a) + B / (b - 1)
    wrong = [abs(A) / (1 - ma) + abs(B) * mb / (mb - 1)]
    why = (r"The two parts never overlap, so add their absolute sums. For $n \le -1$ put $n = -m$, $m \ge 1$:"
           fr"$$\sum_n \lvert h[n]\rvert = {_num(abs(A))}\sum_{{n=0}}^{{\infty}}\left({_num(ma)}\right)^{{n}}"
           fr" + {_num(abs(B))}\sum_{{m=1}}^{{\infty}}\left({_num(1 / mb)}\right)^{{m}}"
           fr" = \frac{{{_num(abs(A))}}}{{1 - {_num(ma)}}} + {_num(abs(B))}\cdot\frac{{{_num(1 / mb)}}}{{1 - {_num(1 / mb)}}} = {_num(S)}.$$")
    return h_tex, S, signed, wrong, spec, why


def fam_hw4(stable):
    if stable:
        A, a = random.choice([F(-3), F(-2), F(-1), F(1), F(2), F(3)]), random.choice(STAB_A)
        t2 = _term(A, _pow(a) + r"\,u[n-1]")
        h_tex = r"\delta[n]" + (" " if t2.startswith("-") else " + ") + t2
        spec = [{"t": "delta", "A": "1", "k": 0}, {"t": "rexp", "A": str(A), "a": str(a), "k": 1}]
        m = abs(a)
        S = 1 + abs(A) * m / (1 - m)
        signed = 1 + A * a / (1 - a)
        wrong = [abs(A) * m / (1 - m)]
        why = (fr"The impulse contributes $\lvert h[0]\rvert = 1$; the exponential starts at $n = 1$ (HW4 #5, $h_2$):"
               fr"$$\sum_n \lvert h[n]\rvert = 1 + {_num(abs(A))}\sum_{{n=1}}^{{\infty}}\left({_num(m)}\right)^{{n}}"
               fr" = 1 + {_num(abs(A))}\cdot\frac{{{_num(m)}}}{{1 - {_num(m)}}} = {_num(S)}.$$")
        return h_tex, S, signed, wrong, spec, why
    C, D, a = random.choice([F(1), F(2), F(3)]), random.choice([F(-1), F(-2), F(-3), F(1)]), random.choice(STAB_A)
    t2 = _term(D, _pow(a) + r"\,u[n]")
    h_tex = _term(C, "u[n]") + (" " if t2.startswith("-") else " + ") + t2
    spec = [{"t": "rexp", "A": str(C), "a": "1", "k": 0}, {"t": "rexp", "A": str(D), "a": str(a), "k": 0}]
    why = (fr"As $n \to \infty$, $h[n] \to {_num(C)}$ (the exponential dies out, the step does not), so infinitely many "
           r"terms have $\lvert h[n]\rvert$ near $" + _num(C) + r"$ and $\sum_n \lvert h[n]\rvert = \infty$ (HW4 #5, $h_1$).")
    return h_tex, None, None, [], spec, why


def fam_left(stable):
    A = random.choice(COEF)
    b = random.choice(BIG_B if stable else [F(1, 2), F(-1, 2), F(1, 3), F(2, 3), F(3, 4), F(1), F(-1)])
    h_tex = _term(A, (_pow(b) + r"\," if b != 1 else "") + "u[-n]")
    spec = [{"t": "lexp", "B": str(A), "b": str(b), "m": 0}]
    if not stable:
        why = (r"$h[n]$ lives on $n \le 0$. Put $n = -m$: $\lvert h[-m]\rvert = " + _num(abs(A)) + r"\cdot"
               + _pow(1 / abs(b), "m") + r"$, which does not decay as $m \to \infty$ because $\lvert b\rvert \le 1$. "
               r"So $\sum_n \lvert h[n]\rvert = \infty$.")
        return h_tex, None, None, [], spec, why
    mb = abs(b)
    S = abs(A) * mb / (mb - 1)
    signed = A * b / (b - 1)
    wrong = [abs(A) / (mb - 1)]
    why = (r"$h[n]$ lives on $n \le 0$. Put $n = -m$, $m \ge 0$ (the $m = 0$ term counts):"
           fr"$$\sum_n \lvert h[n]\rvert = {_num(abs(A))}\sum_{{m=0}}^{{\infty}}\left({_num(1 / mb)}\right)^{{m}}"
           fr" = \frac{{{_num(abs(A))}}}{{1 - {_num(1 / mb)}}} = {_num(S)}.$$")
    return h_tex, S, signed, wrong, spec, why


FAMILIES = [(fam_rexp, 3, 0.55), (fam_abs, 2, 0.6), (fam_nexp, 2, 0.6), (fam_harm, 0.6, 0.0), (fam_cos, 2, 0.5),
            (fam_finite, 2, 1.0), (fam_two, 3, 0.55), (fam_hw4, 2, 0.6), (fam_left, 2, 0.55)]


def generate(data):
    fam, _, pstab = random.choices(FAMILIES, weights=[w for _, w, _ in FAMILIES])[0]
    stable = random.random() < pstab
    h_tex, S, signed, wrong, spec, why = fam(stable)

    p = data["params"]
    p["h_tex"] = h_tex
    p["is_stable"] = stable
    p["family"] = fam.__name__[4:]
    p["spec"] = spec
    p["why_html"] = why
    p["verdict"] = "BIBO stable" if stable else "not BIBO stable"
    p["S_tex"] = _num(S) if stable else ""
    p["S_str"] = str(S) if stable else ""
    p["signed_str"] = str(signed) if stable else ""
    p["wrong_strs"] = [str(w) for w in wrong]
    p["stable_choices"] = mc.yes_no(stable)
    data["correct_answers"]["sumabs"] = float(S) if stable else ""


def grade(data):
    """Unstable h: part (b) has no content of its own (the box must stay empty), so it only counts when the
    box is empty AND part (a) is right; an empty box after answering "stable" earns nothing."""
    if data["params"].get("is_stable"):
        return
    ps = data["partial_scores"]
    if "sumabs" not in ps:
        return
    sub = data["submitted_answers"].get("sumabs")
    blank = sub is None or (isinstance(sub, str) and sub.strip() == "")
    stable_score = (ps.get("stable") or {}).get("score") or 0
    ps["sumabs"]["score"] = stable_score if blank else 0
    if blank and not stable_score:
        ps["sumabs"]["feedback"] = "An empty box only counts together with the answer No in part (a)."
    pl.set_weighted_score_data(data)
