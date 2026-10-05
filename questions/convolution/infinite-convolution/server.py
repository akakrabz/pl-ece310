"""Infinite-length convolution y = x * h (exam family "infinite-length convolution", part (b) of
the convolution problem on 5 of 7 past Midterm 1 exams; HW2 #5(b)-(d); Lecture 4 section 2.2.3).

Three templates (chosen at random):
  "ab"     x = A a^n u[n-k1], h = B b^n u[n-k2] with |a| != |b|   (SP2025 #4b, SP2023 #3b)
           -> y[n] = C_a a^n + C_b b^n for n >= n0 = k1 + k2
  "aa"     the same with a = b                                        (the (n+1) a^n trap)
           -> y[n] = AB (n - n0 + 1) a^n for n >= n0
  "finite" x = {c_0, c_1, c_2} starting at nx, h = B b^n u[n-k]        (HW2 #5(b), FA2023 #3b)
           -> y[n] = sum_i c_i h[n - nx - i]; for n >= n0 + 2 this is C b^n
Student enters n0 (integer), the closed form (pl-symbolic-input in n, without the step factor)
and one sample: y[n0 + 2] ("ab", "aa") or y[n0 + 1] ("finite": a transient sample, the third copy of h
has not started yet, so part (c) is never part (b) evaluated at a point).

The symbolic answer is stored with n declared an integer (convlib.sym_answer_json) so that, e.g.,
(-3)^(-n) is accepted for (-1/3)^n."""

import random
from fractions import Fraction as F

from ece310 import convlib as cl
from ece310 import fmt

DECAY = [F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(1, 4), F(-1, 4), F(2, 3), F(3, 4)]
GROW = [F(2), F(-2), F(3)]


def _base():
    return random.choice(DECAY) if random.random() < 0.8 else random.choice(GROW)


def _nice(*qs, num=24, den=48):
    return all(abs(q.numerator) <= num and q.denominator <= den for q in qs)


def _coef_ab(A, a, k1, B, b, k2):
    """Closed form of A a^n u[n-k1] * B b^n u[n-k2] (a != b) for n >= k1 + k2: (C_a, C_b)."""
    Ca = A * B * a ** (1 - k2) * b ** k2 / (a - b)
    Cb = A * B * a ** k1 * b ** (1 - k1) / (b - a)
    return Ca, Cb


def _gen_ab(p):
    while True:
        a, b = _base(), _base()
        if abs(a) == abs(b):
            continue
        A, B = random.choice([1, 1, 1, 2, -1]), random.choice([1, 1, 1, 2, -1, 3])
        k1, k2 = random.choice([-1, 0, 0, 1, 2]), random.choice([-1, 0, 1, 1, 2, 3])
        if k1 == 0 and k2 == 0:
            continue
        Ca, Cb = _coef_ab(F(A), a, k1, F(B), b, k2)
        n0 = k1 + k2
        y2 = Ca * a ** (n0 + 2) + Cb * b ** (n0 + 2)
        if _nice(Ca, Cb) and _nice(y2, num=100, den=144):
            break
    r = a / b
    AB = F(A * B)
    p.update(template="ab", a=str(a), b=str(b), A=A, B=B, k1=k1, k2=k2, n_first=n0, c_off=2)
    p["x_tex"] = cl.tex_exp_step(A, a, k1)
    p["h_tex"] = cl.tex_exp_step(B, b, k2)
    p["sum_tex"] = (r"\sum_{k=-\infty}^{\infty} "
                    + cl.tex_coef(A, fmt.tex_pow(a, "k")) + r"\,u[" + cl.tex_arg(k1, "k") + r"]\;"
                    + cl.tex_coef(B, fmt.tex_pow(b, "n-k")) + r"\,u[" + _tex_nk(k2) + "]")
    p["range_tex"] = rf"{k1} \le k \le {cl.tex_arg(k2)}"
    p["geo_tex"] = (cl.tex_coef(AB, fmt.tex_pow(b)) + r"\sum_{k=" + str(k1) + "}^{" + cl.tex_arg(k2) + "}"
                    + fmt.tex_pow(r, "k")
                    + " = " + cl.tex_coef(AB, fmt.tex_pow(b)) + r"\cdot\dfrac{" + fmt.tex_pow(r, str(k1))
                    + " - " + fmt.tex_pow(r, cl.tex_arg(k2 - 1)) + "}{1 - " + fmt.tex_num(r, paren=True) + "}")
    p["r_tex"] = fmt.tex_num(r)
    p["y_tex"] = cl.tex_signed_terms([(Ca, fmt.tex_pow(a)), (Cb, fmt.tex_pow(b))])
    p["y0_tex"] = (rf"x[{k1}]\,h[{k2}] = " + fmt.tex_num(A * a ** k1 * B * b ** k2))
    p["yc_tex"] = (cl.tex_signed_terms([(Ca, fmt.tex_pow(a, str(n0 + 2))), (Cb, fmt.tex_pow(b, str(n0 + 2)))])
                   + " = " + fmt.tex_num(y2))
    sym = fmt.sym_sum([(Ca, fmt.sym_pow(a)), (Cb, fmt.sym_pow(b))])
    return n0, sym, y2


def _tex_nk(k2):
    """argument n - k - k2 of the step in h[n - k]."""
    if k2 == 0:
        return "n-k"
    return f"n-k-{k2}" if k2 > 0 else f"n-k+{-k2}"


def _gen_aa(p):
    while True:
        a = _base()
        A, B = random.choice([1, 1, 2, -1]), random.choice([1, 1, 2, -1, 3])
        k1, k2 = random.choice([-1, 0, 1, 2]), random.choice([-1, 0, 1, 2, 3])
        n0 = k1 + k2
        if n0 == 0:
            continue
        y2 = A * B * 3 * a ** (n0 + 2)
        if _nice(y2, num=100, den=144):
            break
    AB = F(A * B)
    p.update(template="aa", a=str(a), b=str(a), A=A, B=B, k1=k1, k2=k2, n_first=n0, c_off=2)
    p["x_tex"] = cl.tex_exp_step(A, a, k1)
    p["h_tex"] = cl.tex_exp_step(B, a, k2)
    p["sum_tex"] = (r"\sum_{k=-\infty}^{\infty} "
                    + cl.tex_coef(A, fmt.tex_pow(a, "k")) + r"\,u[" + cl.tex_arg(k1, "k") + r"]\;"
                    + cl.tex_coef(B, fmt.tex_pow(a, "n-k")) + r"\,u[" + _tex_nk(k2) + "]")
    p["range_tex"] = rf"{k1} \le k \le {cl.tex_arg(k2)}"
    lin = r"\left(" + cl.tex_arg(n0 - 1) + r"\right)"
    p["geo_tex"] = (cl.tex_coef(AB, fmt.tex_pow(a)) + r"\sum_{k=" + str(k1) + "}^{" + cl.tex_arg(k2) + "} 1"
                    + " = " + cl.tex_coef(AB, fmt.tex_pow(a)) + r"\,\big[(" + cl.tex_arg(k2) + ") - (" + str(k1)
                    + r") + 1\big]")
    p["y_tex"] = cl.tex_coef(AB, lin + fmt.tex_pow(a))
    p["y0_tex"] = rf"x[{k1}]\,h[{k2}] = " + fmt.tex_num(A * a ** k1 * B * a ** k2)
    p["yc_tex"] = (cl.tex_coef(AB, r"\cdot 3\cdot " + fmt.tex_pow(a, str(n0 + 2))) + " = " + fmt.tex_num(y2)
                   if AB not in (1, -1) else
                   ("-" if AB < 0 else "") + r"3\cdot " + fmt.tex_pow(a, str(n0 + 2)) + " = " + fmt.tex_num(y2))
    sym = f"{fmt.sym_num(AB)}*(n - ({n0 - 1}))*{fmt.sym_pow(a)}"
    return n0, sym, y2


def _gen_finite(p):
    while True:
        L = 3          # with 3 samples y[n0 + 1] is a transient sample (the last copy has not started)
        c = [random.choice([-2, -1, 1, 2]) if i in (0, L - 1) else random.randint(-2, 2) for i in range(L)]
        nx = random.choice([-1, 0, 0, 1])
        b = random.choice([F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(2), F(-2), F(3), F(1, 4)])
        B = random.choice([1, 1, 2, -1])
        k = random.choice([-1, 0, 0, 1, 2])
        C = B * sum(ci * b ** (-(nx + i)) for i, ci in enumerate(c))
        if C == 0 or not _nice(C, num=40, den=40):
            continue
        n0 = nx + k
        h = lambda m: B * b ** m if m >= k else F(0)          # noqa: E731
        yv = lambda m: sum(ci * h(m - nx - i) for i, ci in enumerate(c))   # noqa: E731
        y1 = yv(n0 + 1)
        if y1 != 0 and _nice(y1, num=100, den=144):
            break
    p.update(template="finite", c=c, nx=nx, b=str(b), B=B, k=k, L=L, n_first=n0, c_off=1, n_full_off=L - 1)
    p["x_tex"] = fmt.tex_seq(c, nx)
    p["h_tex"] = cl.tex_exp_step(B, b, k)
    p["xdelta_tex"] = cl.tex_signed_terms([(ci, cl.tex_delta(nx + i)) for i, ci in enumerate(c)])
    copies = []
    for i, ci in enumerate(c):
        s = nx + i
        copies.append((ci * B, fmt.tex_pow(b, cl.tex_arg(s)) + r"\,u[" + cl.tex_arg(s + k) + "]"))
    p["copies_tex"] = cl.tex_signed_terms(copies)
    p["hcopies_tex"] = cl.tex_signed_terms([(ci, f"h[{cl.tex_arg(nx + i)}]") for i, ci in enumerate(c)])
    p["full_tex"] = rf"n \ge {n0 + L - 1}"
    p["factor_tex"] = (cl.tex_coef(B, fmt.tex_pow(b)) + r"\left("
                       + cl.tex_signed_terms([(ci, fmt.tex_pow(b, str(-(nx + i))) if nx + i != 0 else "")
                                              for i, ci in enumerate(c)])
                       + r"\right)")
    p["y_tex"] = cl.tex_coef(C, fmt.tex_pow(b))
    p["y0_tex"] = rf"x[{nx}]\,h[{k}] = " + fmt.tex_num(c[0] * B * b ** k)
    terms = [(ci, f"h[{n0 + 1 - nx - i}]") for i, ci in enumerate(c) if n0 + 1 - nx - i >= k]
    p["yc_tex"] = cl.tex_signed_terms(terms) + " = " + fmt.tex_num(y1)
    p["n1_lt_full"] = L == 3
    sym = f"{fmt.sym_num(C)}*{fmt.sym_pow(b)}"
    return n0, sym, y1


_DEFAULTS = {  # every param the template mentions exists in every variant (sections hide the unused ones)
    "sum_tex": "", "range_tex": "", "geo_tex": "", "r_tex": "", "k1": 0, "k2": 0,
    "xdelta_tex": "", "hcopies_tex": "", "copies_tex": "", "full_tex": "", "factor_tex": "",
    "nx": 0, "k": 0, "n_full_off": 0, "n1_lt_full": False, "L": 0,
}


def generate(data):
    u = random.random()
    p = data["params"]
    p.update(_DEFAULTS)
    if u < 0.5:
        n0, sym, yc = _gen_ab(p)
    elif u < 0.75:
        n0, sym, yc = _gen_aa(p)
    else:
        n0, sym, yc = _gen_finite(p)
    p["is_ab"] = p["template"] == "ab"
    p["is_aa"] = p["template"] == "aa"
    p["is_finite"] = p["template"] == "finite"
    p["is_exp"] = not p["is_finite"]
    p["y_sym"] = sym
    p["yc_frac"] = fmt.plain_num(yc)
    data["correct_answers"]["n0"] = n0
    data["correct_answers"]["yform"] = cl.sym_answer_json(sym)
    data["correct_answers"]["yc"] = float(yc)
