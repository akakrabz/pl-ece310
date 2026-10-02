"""Step response <-> impulse response of an LTI system (s = h * u = running sum of h; h[n] = s[n] - s[n-1]).

Template "fir" (FA2019 #3(b), Lecture 4): finite h[n] given; enter s[n] on a window that starts one
               sample before h and ends just before the last sample of h (row vector), and the
               steady-state value s[n] for all n >= last index of h (= sum of h).
Template "exp" (FA2023 #4 run with an infinite step response; HW2 #3): s[n] given in closed form,
               either c(1 - a^{n+1})u[n] or (c + d a^n)u[n]; enter h[0] (number) and the closed form
               of h[n] for n >= 1 (pl-symbolic-input, n declared an integer)."""

import random
from fractions import Fraction as F

import numpy as np
import prairielearn as pl
from ece310 import convlib as cl
from ece310 import fmt

BASES = [F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(1, 4), F(-1, 4), F(2, 3), F(-2, 3), F(3, 4), F(-3, 4), F(1, 5)]


def _gen_fir(p, ca):
    L = random.choice([3, 4, 4])
    h = [random.randint(-3, 3) for _ in range(L)]
    h[0], h[-1] = random.choice([-3, -2, -1, 1, 2, 3]), random.choice([-3, -2, -1, 1, 2, 3])
    nh = random.choice([-2, -1, 0, 1])
    s = list(np.cumsum(h).astype(int))
    w0 = nh - 1
    win = [0] + [int(v) for v in s[:-1]]
    steady = int(s[-1])
    p.update(template="fir", h=h, nh=nh, L=L, w0=w0, w_last=nh + L - 2, n_last=nh + L - 1)
    p["h_tex"] = fmt.tex_seq(h, nh)
    p["w_tex"] = r",\ ".join(str(n) for n in range(w0, nh + L - 1))
    p["steps_tex"] = cl.tex_signed_terms([(hv, cl.tex_u(nh + i)) for i, hv in enumerate(h)])
    ns = list(range(nh - 1, nh + L + 1))
    hv = lambda n: h[n - nh] if 0 <= n - nh < L else 0                   # noqa: E731
    sv = lambda n: sum(hv(k) for k in range(nh, n + 1))                  # noqa: E731
    head = "".join(f"<th>$n={n}$</th>" for n in ns)
    rh = "".join(f"<td>${hv(n)}$</td>" for n in ns)
    rs = "".join(f"<td><b>${sv(n)}$</b></td>" for n in ns)
    p["table_html"] = ('<table class="table table-sm table-bordered text-center" style="width:auto">'
                       f"<tr><th></th>{head}<th></th></tr><tr><th>$h[n]$</th>{rh}<td>$\\cdots$</td></tr>"
                       f"<tr><th>$s[n]$</th>{rs}<td>$\\cdots$</td></tr></table>")
    p["sum_tex"] = " + ".join(fmt.tex_num(v, paren=(i > 0)) for i, v in enumerate(h)) + f" = {steady}"
    ca["win"] = pl.to_json(np.array([[float(v) for v in win]]))
    ca["steady"] = float(steady)


def _gen_exp(p, ca):
    while True:
        a = random.choice(BASES)
        if random.random() < 0.5:                       # s[n] = c (1 - a^{n+1}) u[n]
            c = random.choice([F(1), F(2), F(3), F(4), F(6), F(-1), F(-2), F(-3), F(1, 2), F(3, 2), F(-1, 2)])
            d = -c * a
            form = "geo"
        else:                                           # s[n] = (c + d a^n) u[n]
            c = F(random.choice([1, 2, 3, 4, -1, -2, -3]))
            d = F(random.choice([-4, -3, -2, -1, 1, 2, 3, 4]))
            form = "lin"
        h0 = c + d
        K = d * (a - 1) / a                              # h[n] = K a^n for n >= 1
        if h0.denominator <= 12 and K.denominator <= 12 and abs(K.numerator) <= 24 and abs(h0.numerator) <= 24:
            break
    p.update(template="exp", a=str(a), c=str(c), d=str(d), form=form)
    if form == "geo":
        p["s_tex"] = (cl.tex_coef(c, r"\left(1 - " + fmt.tex_pow(a, "n+1") + r"\right)") + r"u[n]")
    else:
        p["s_tex"] = r"\left(" + cl.tex_signed_terms([(c, ""), (d, fmt.tex_pow(a))]) + r"\right)u[n]"
    p["s_inner_tex"] = cl.tex_signed_terms([(c, ""), (d, fmt.tex_pow(a))])
    p["s_inner1_tex"] = cl.tex_signed_terms([(c, ""), (d, fmt.tex_pow(a, "n-1"))])
    p["h0_tex"] = fmt.tex_num(c + d)
    p["s0_tex"] = cl.tex_signed_terms([(c, ""), (d, fmt.tex_pow(a, "0"))])
    p["diff_tex"] = (cl.tex_signed_terms([(d, fmt.tex_pow(a))]) + " - " + r"\left(" +
                     cl.tex_signed_terms([(d, fmt.tex_pow(a, "n-1"))]) + r"\right)")
    p["fact_tex"] = cl.tex_coef(d, fmt.tex_pow(a, "n-1")) + r"\left(" + fmt.tex_num(a) + r" - 1\right)"
    p["hn_tex"] = cl.tex_coef(K, fmt.tex_pow(a))
    p["delta_coef"] = fmt.plain_num(h0 - K)
    p["has_delta"] = h0 != K
    p["h_full_tex"] = (cl.tex_signed_terms([(h0 - K, r"\delta[n]"), (K, fmt.tex_pow(a) + r"\,u[n]")])
                       if h0 != K else cl.tex_coef(K, fmt.tex_pow(a) + r"\,u[n]"))
    p["h0_frac"] = fmt.plain_num(h0)
    sym = f"{fmt.sym_num(K)}*{fmt.sym_pow(a)}"
    p["hn_sym"] = sym
    ca["h0"] = float(h0)
    ca["hn"] = cl.sym_answer_json(sym)


_DEFAULTS = {"h": [], "nh": 0, "L": 0, "w0": 0, "w_last": 0, "n_last": 0, "h_tex": "", "w_tex": "", "steps_tex": "",
             "table_html": "", "sum_tex": "", "a": "0", "c": "0", "d": "0", "form": "", "s_tex": "", "s_inner_tex": "",
             "s_inner1_tex": "", "h0_tex": "", "s0_tex": "", "diff_tex": "", "fact_tex": "", "hn_tex": "",
             "delta_coef": "0", "has_delta": False, "h_full_tex": "", "h0_frac": "0", "hn_sym": ""}


def generate(data):
    p, ca = data["params"], data["correct_answers"]
    p.update(_DEFAULTS)
    if random.random() < 0.5:
        _gen_fir(p, ca)
    else:
        _gen_exp(p, ca)
    p["is_fir"] = p["template"] == "fir"
    p["is_exp"] = p["template"] == "exp"
