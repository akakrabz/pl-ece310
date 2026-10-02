"""Causal H(z) with complex-conjugate poles r e^{+-j theta}: h[n] = r^n (A cos(theta n) + B sin(theta n)) u[n].

Model: Lecture 8 Exercise 2 (conjugate poles e^{+-j pi/3}, h = 3 cos(pi n/3) u[n]) and SP2021
Midterm 1 #5(a) (3z^-1/(1 + z^-2), |z| > 1  ->  3 sin(pi n/2) u[n]).

theta in {pi/2, pi/3, 2pi/3}; r in {1/3, 1/2, 2/3, 3/4, 1, 4/3, 3/2, 2}. For theta = pi/3, 2pi/3 the sine
coefficient is B = k sqrt(3) with k rational, so every coefficient of H and every sample h[n] is
rational while the closed form is exact. Asked: h[0], h[1], h[2] (numbers) and the closed form for
n >= 0 (pl-symbolic-input in an integer n, with cos, sin, pi, sqrt)."""

import random
from fractions import Fraction as F

import prairielearn as pl
from ece310 import fmt, zinv

R_CHOICES = [F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(1), F(4, 3), F(3, 2), F(2)]
# theta: (cos theta, sin theta written as s*sqrt(3) (sqrt3=True) or s, cos(theta n), sin(theta n) / unit)
TH = {
    "pi/2": {"cos": F(0), "sin": F(1), "sqrt3": False, "c": [F(1), F(0), F(-1)], "s": [F(0), F(1), F(0)],
             "tex": r"\tfrac{\pi}{2}", "sym": "pi*n/2", "deg": "90"},
    "pi/3": {"cos": F(1, 2), "sin": F(1, 2), "sqrt3": True, "c": [F(1), F(1, 2), F(-1, 2)], "s": [F(0), F(1, 2), F(1, 2)],
             "tex": r"\tfrac{\pi}{3}", "sym": "pi*n/3", "deg": "60"},
    "2*pi/3": {"cos": F(-1, 2), "sin": F(1, 2), "sqrt3": True, "c": [F(1), F(-1, 2), F(-1, 2)], "s": [F(0), F(1, 2), F(-1, 2)],
               "tex": r"\tfrac{2\pi}{3}", "sym": "2*pi*n/3", "deg": "120"},
}
A_CHOICES = [1, 2, 3, 4, -1, -2, -3, -4]
B_PLAIN = [1, 2, 3, -1, -2, -3, F(1, 2), F(-1, 2), F(3, 2), F(-3, 2)]
K_SQRT3 = [1, 2, -1, -2, F(1, 3), F(-1, 3), F(2, 3), F(-2, 3), F(1, 2), F(-1, 2)]


def _pick():
    while True:
        th = random.choice(list(TH))
        t = TH[th]
        r = random.choice(R_CHOICES)
        kind = random.random()
        A = F(random.choice(A_CHOICES)) if kind > 0.12 else F(0)
        if t["sqrt3"]:
            k = F(random.choice(K_SQRT3)) if kind < 0.88 else F(0)   # B = k sqrt(3)
        else:
            k = F(random.choice(B_PLAIN)) if kind < 0.88 else F(0)   # B = k
        if A == 0 and k == 0:
            continue
        # B sin(theta) as a rational number: sqrt3 family sin(theta) = sqrt(3)/2, B = k sqrt(3) -> 3k/2
        Bs = 3 * k / 2 if t["sqrt3"] else k
        b0 = A
        b1 = r * (-A * t["cos"] + Bs)
        a1 = -2 * r * t["cos"]
        a2 = r * r
        if b1 == 0 and A == 0:
            continue
        h = []
        for n in range(3):
            sn = 3 * k * t["s"][n] if t["sqrt3"] else k * t["s"][n]   # B sin(theta n), rational
            h.append(r**n * (A * t["c"][n] + sn))
        if not all(zinv.nice(x, maxden=18, maxnum=40, bad_dens=()) for x in (b1, a1, a2) + tuple(h)):
            continue
        return th, r, A, k, b0, b1, a1, a2, h


def _B_tex(k, sqrt3):
    """(sign, magnitude LaTeX) of B = k sqrt(3) or B = k."""
    if not sqrt3:
        return (-1 if k < 0 else 1), fmt.tex_num(abs(k))
    m = abs(k)
    if m == 1:
        mag = r"\sqrt{3}"
    elif m.denominator == 1:
        mag = f"{m.numerator}\\sqrt{{3}}"
    else:
        mag = (f"\\tfrac{{{m.numerator}\\sqrt{{3}}}}{{{m.denominator}}}" if m.numerator != 1
               else f"\\tfrac{{\\sqrt{{3}}}}{{{m.denominator}}}")
    return (-1 if k < 0 else 1), mag


def generate(data):
    th, r, A, k, b0, b1, a1, a2, h = _pick()
    t = TH[th]
    sq = t["sqrt3"]
    cos_t = r"\cos\!\left(" + t["tex"] + r"n\right)"
    sin_t = r"\sin\!\left(" + t["tex"] + r"n\right)"
    inner = []
    if A != 0:
        inner.append(("-" if A < 0 else "") + ("" if abs(A) == 1 else fmt.tex_num(abs(A))) + cos_t)
    if k != 0:
        sg, mag = _B_tex(k, sq)
        mag = "" if mag == "1" else mag
        body = mag + ("\\," if mag else "") + sin_t
        inner.append(("-" if sg < 0 else "") + body if not inner else (" - " if sg < 0 else " + ") + body)
    inner_tex = "".join(inner)
    rn = "" if r == 1 else fmt.tex_pow(r, "n")
    closed = (rn + r"\left(" + inner_tex + r"\right)") if (rn and len(inner) > 1) else (rn + r"\," + inner_tex if rn else inner_tex)

    # SymPy string (exact)
    parts = []
    if A != 0:
        parts.append(f"{fmt.sym_num(A)}*cos({t['sym']})")
    if k != 0:
        parts.append(f"{fmt.sym_num(k)}*sqrt(3)*sin({t['sym']})" if sq else f"{fmt.sym_num(k)}*sin({t['sym']})")
    sym = f"{fmt.sym_pow(r, 'n')}*(" + " + ".join(parts) + ")"

    p = data["params"]
    p["b"] = [str(b0), str(b1)]
    p["a"] = ["1", str(a1), str(a2)]
    p["r"], p["theta"], p["A"], p["k"], p["sqrt3"] = str(r), th, str(A), str(k), sq
    p["H_tex"] = fmt.tex_frac(fmt.tex_poly_zinv([b0, b1]), fmt.tex_poly_zinv([1, a1, a2]))
    p["roc_tex"] = r"\lvert z\rvert \gt " + fmt.tex_num(r)
    p["den_tex"] = fmt.tex_poly_zinv([1, a1, a2])
    p["r_tex"], p["r2_tex"] = fmt.tex_num(r), fmt.tex_num(r * r)
    p["cos_th_tex"] = fmt.tex_num(t["cos"])
    p["th_tex"] = t["tex"]
    p["sin_th_tex"] = r"\tfrac{\sqrt{3}}{2}" if sq else "1"
    p["twor_cos_tex"] = fmt.tex_num(2 * r * t["cos"])
    p["rcos_tex"] = fmt.tex_num(r * t["cos"])
    p["rsin_tex"] = (fmt.tex_num(r / 2) + r"\sqrt{3}" if r / 2 != 1 else r"\sqrt{3}") if sq else fmt.tex_num(r)
    p["A_tex"] = fmt.tex_num(A)
    sgB, magB = _B_tex(k, sq)
    p["B_tex"] = ("0" if k == 0 else ("-" if sgB < 0 else "") + magB)
    p["Brsin_tex"] = fmt.tex_num(3 * k / 2 * r if sq else k * r)
    p["b1_tex"] = fmt.tex_num(b1)
    p["closed_tex"] = closed
    p["h0_tex"], p["h1_tex"], p["h2_tex"] = (fmt.tex_num(x) for x in h)
    p["b_tex_num"] = fmt.tex_poly_zinv([b0, b1])
    p["rec1_tex"] = fr"h[1] = {fmt.tex_num(b1)} {'+' if -a1 >= 0 else '-'} {fmt.tex_num(abs(a1), paren=False)}\cdot{fmt.tex_num(h[0], paren=True)} = {fmt.tex_num(h[1])}"
    p["rec2_tex"] = (fr"h[2] = {fmt.tex_num(-a1, paren=False)}\cdot{fmt.tex_num(h[1], paren=True)} - "
                     fr"{fmt.tex_num(a2)}\cdot{fmt.tex_num(h[0], paren=True)} = {fmt.tex_num(h[2])}")
    p["lfilt_tex"] = (fr"h[n] = {fmt.tex_sum([(-a1, 'h[n-1]'), (-a2, 'h[n-2]'), (b0, r'\delta[n]'), (b1, r'\delta[n-1]')])}")

    c = data["correct_answers"]
    c["h0"], c["h1"], c["h2"] = float(h[0]), float(h[1]), float(h[2])
    c["hn"] = zinv.sym_json(sym)


def _samples(params, count=12):
    """Exact h[0..count-1] of the causal recursion y[n] + a1 y[n-1] + a2 y[n-2] = b0 x[n] + b1 x[n-1]."""
    b = [F(x) for x in params["b"]]
    a = [F(x) for x in params["a"]]
    h = []
    for n in range(count):
        v = b[n] if n < len(b) else F(0)
        for k in range(1, len(a)):
            if n - k >= 0:
                v -= a[k] * h[n - k]
        h.append(v)
    return h


def grade(data):
    """h[n] is a sequence, so the closed form is graded by its values at n = 0..11: SymPy's equality
    test cannot see that e.g. cos(4*pi*n/3) = cos(2*pi*n/3) for every integer n. The element's own
    parse still supplies the format errors. NaN/inf-safe: a sample matches only if abs(err) <= tol."""
    name = "hn"
    if name in data["format_errors"]:
        return
    sub = data["submitted_answers"].get(name)
    if not isinstance(sub, dict):
        return
    ok = False
    try:
        from prairielearn import sympy_utils as psu

        expr = psu.json_to_sympy(sub, allow_complex=False)
        syms = list(expr.free_symbols)
        if all(str(v) == "n" for v in syms):
            ok = True
            for k, want in enumerate(_samples(data["params"])):
                got = complex(expr.subs({v: k for v in syms}).evalf(30))
                err = abs(got - float(want))
                if not (err <= 1e-9 * max(1.0, abs(float(want)))):
                    ok = False
                    break
    except Exception:
        ok = False
    ps = data["partial_scores"].setdefault(name, {"weight": 1})
    ps["score"] = 1.0 if ok else 0.0
    pl.set_weighted_score_data(data)
