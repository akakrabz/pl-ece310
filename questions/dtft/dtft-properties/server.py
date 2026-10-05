"""DTFT of a modified signal from the Lecture 14 property table (Table 2), plus one number the
property gives quickly. Lecture 14; HW5 Problems 1(a), 1(c), 1(d), 2(d).

Base pairs (Lecture 14 Table 1 / definition):
  geo  x[n] = b a^n u[n]             <->  X_d(w) = b / (1 - a e^{-jw})
  fin  x[n] = {c0, c1[, c2]} from n=0 <->  X_d(w) = c0 + c1 e^{-jw} + c2 e^{-j2w}
Transformations:
  shift  y = x[n-k]            Y = e^{-jkw} X_d(w)                  number: Y_d(pi)
  mod    y = e^{j w0 n} x[n]   Y = X_d(w - w0),  w0 = +-pi/2         number: Y_d(w0) = X_d(0)
  alt    y = (-1)^n x[n]       Y = X_d(w - pi)                      number: Y_d(0) or Y_d(pi)
  cos    y = x[n] cos(pi n/2)  Y = 1/2 X_d(w - pi/2) + 1/2 X_d(w + pi/2)   number: Y_d(0) or Y_d(pi/2)
  rev    y = x[-n]             Y = X_d(-w)                          number: Parseval integral / pi
  conv   y = x * x             Y = X_d(w)^2                         number: Y_d(0) or Y_d(pi)
The differentiation property n x[n] is deliberately not used: the Lecture 14 table prints it as
-j dX_d/dw, but the correct pair is n x[n] <-> j dX_d/dw.

Every option is a sum of terms coef * prefactor(w) * X'(m w + c)^p, where X' is the base transform
(possibly with modified coefficients); its LaTeX and its numerical value are built from the same
term list, and options that coincide numerically with the answer or with each other are dropped."""

import cmath
import math
import random
from fractions import Fraction as F

from ece310 import fmt

PI_TEX = {F(1, 2): r"\tfrac{\pi}{2}", F(-1, 2): r"-\tfrac{\pi}{2}"}
W_PROBE = [0.31, 1.07, 1.9, 2.6, -0.77, -2.2, 0.0, 3.0]


# ----------------------------------------------------------------------------- bases
def _base_coeffs(base):
    return base["c"] if base["kind"] == "fin" else None


def _alt(base):
    """Base of (-1)^n x[n]: a -> -a, c_n -> (-1)^n c_n."""
    if base["kind"] == "geo":
        return {"kind": "geo", "a": -base["a"], "b": base["b"]}
    return {"kind": "fin", "c": [(-1) ** n * v for n, v in enumerate(base["c"])]}


def _square_samples(base):
    """Base of x[n]^2 (product instead of convolution)."""
    if base["kind"] == "geo":
        return {"kind": "geo", "a": base["a"] ** 2, "b": base["b"] ** 2}
    return {"kind": "fin", "c": [v * v for v in base["c"]]}


def _recip(base):
    return {"kind": "geo", "a": 1 / base["a"], "b": base["b"]}


# ----------------------------------------------------------------------------- arguments theta = m w + c
def _theta_tex(m, c):
    """m w + c with c a multiple of pi (Fraction, tagged 'pi') or plain radians (int, tagged 'rad')."""
    cval, ckind = c
    out = fmt.tex_sum([(m, r"\omega")])
    if cval == 0:
        return out
    ctex = (r"\tfrac{\pi}{2}" if abs(cval) == F(1, 2) else r"\pi") if ckind == "pi" else str(abs(cval))
    if m == -1:                                       # write c - w rather than -w + c
        return (("-" if cval < 0 else "") + ctex) + r" - \omega"
    return out + (" - " if cval < 0 else " + ") + ctex


def _exp_tex(n, m, c):
    """e^{-j n theta} for theta = m w + c ('' for n = 0)."""
    cval, _ = c
    if n == 0:
        return ""
    if cval == 0:
        p = n * m
        if p == 0:
            return ""
        k = "" if abs(p) == 1 else str(abs(p))
        return f"e^{{{'-' if p > 0 else ''}j{k}\\omega}}"
    k = "" if n == 1 else str(n)
    return f"e^{{-j{k}\\left({_theta_tex(m, c)}\\right)}}"


def _theta_val(m, c, w):
    cval, ckind = c
    return m * w + (float(cval) * math.pi if ckind == "pi" else float(cval))


def _X_parts(base, m, c):
    """(numerator-free LaTeX of the base transform at theta, value function)."""
    if base["kind"] == "geo":
        a, b = base["a"], base["b"]
        den = fmt.tex_sum([(1, ""), (-a, _exp_tex(1, m, c))])
        return ("geo", den), (lambda w: float(b) / (1 - float(a) * cmath.exp(-1j * _theta_val(m, c, w))))
    cs = base["c"]
    poly = fmt.tex_sum((v, _exp_tex(n, m, c)) for n, v in enumerate(cs))
    return ("fin", poly), (lambda w: sum(float(v) * cmath.exp(-1j * n * _theta_val(m, c, w)) for n, v in enumerate(cs)))


# ----------------------------------------------------------------------------- terms and options
def _pref_tex(pref):
    if pref is None:
        return ""
    kind, val = pref
    if kind == "ejw":                                  # e^{-j val w}
        return _exp_tex(1, val, (0, "pi"))
    if kind == "const_pi":                             # e^{j val pi}, val = +-1/2
        return "e^{" + ("" if val > 0 else "-") + "j" + PI_TEX[abs(val)] + "}"
    raise ValueError(kind)


def _pref_val(pref, w):
    if pref is None:
        return 1
    kind, val = pref
    if kind == "ejw":
        return cmath.exp(-1j * val * w)
    return cmath.exp(1j * float(val) * math.pi)


def _term(coef, pref, base, m, c, power=1, wrap=False):
    (kind, body), f = _X_parts(base, m, c)
    num_body = _pref_tex(pref)
    if kind == "geo":
        cb = coef * base["b"] ** power
        num = fmt.tex_sum([(cb, num_body)]) if num_body else fmt.tex_num(abs(cb))
        den = body if power == 1 else rf"\left({body}\right)^{{2}}"
        neg = cb < 0 and not num_body                 # the sign goes in front of the fraction
        tex = fmt.tex_frac(num, den)
    else:
        inner = body if (power == 1 and not num_body and coef == 1 and not wrap) else rf"\left({body}\right)"
        if power == 2:
            inner += "^{2}"
        lead = fmt.tex_sum([(coef, num_body)]) if num_body else ("" if coef == 1 else ("-" if coef == -1 else fmt.tex_num(coef)))
        sep = r"\," if lead and lead not in ("-",) and not lead.endswith("}") else ""
        tex = lead + sep + inner
        neg = False
    val = (lambda w, coef=coef, pref=pref, f=f, power=power: float(coef) * _pref_val(pref, w) * f(w) ** power)
    return tex, val, neg


def _option(terms_spec, halves=False):
    """terms_spec: list of (sign, coef, pref, base, m, c, power). halves: write each term as 1/2 . term."""
    texs, vals = [], []
    for i, (sign, coef, pref, base, m, c, power) in enumerate(terms_spec):
        t, v, neg = _term(coef, pref, base, m, c, power, wrap=(halves or sign < 0 or len(terms_spec) > 1))
        shown = -sign if neg else sign                 # v already carries the term's own sign
        if halves:
            t = r"\tfrac{1}{2}" + ("" if t.startswith(r"\left(") else r"\cdot ") + t
            v = (lambda w, v=v: 0.5 * v(w))
        if i == 0:
            texs.append(("-" if shown < 0 else "") + t)
        else:
            texs.append((" - " if shown < 0 else " + ") + t)
        vals.append((lambda w, v=v, sign=sign: sign * v(w)))
    return "".join(texs), (lambda w: sum(v(w) for v in vals))


# ----------------------------------------------------------------------------- exact values at 0, pi/2, pi
def _X_exact(base, w0):
    """X_d(w0) for w0 in {'0', 'pi/2', '-pi/2', 'pi'} as an exact (re, im) pair."""
    e = {"0": (F(1), F(0)), "pi": (F(-1), F(0)), "pi/2": (F(0), F(-1)), "-pi/2": (F(0), F(1))}[w0]   # e^{-j w0}
    if base["kind"] == "geo":
        a, b = base["a"], base["b"]
        dr, di = 1 - a * e[0], -a * e[1]
        d2 = dr * dr + di * di
        return b * dr / d2, -b * di / d2
    re = im = F(0)
    pr, pi_ = F(1), F(0)                      # e^{-j w0 n}
    for v in base["c"]:
        re, im = re + v * pr, im + v * pi_
        pr, pi_ = pr * e[0] - pi_ * e[1], pr * e[1] + pi_ * e[0]
    return re, im


def _energy(base):
    if base["kind"] == "geo":
        return base["b"] ** 2 / (1 - base["a"] ** 2)
    return F(sum(v * v for v in base["c"]))


# ----------------------------------------------------------------------------- variants
def _draw_base():
    if random.random() < 0.5:
        a = random.choice([F(1, 2), F(-1, 2), F(1, 3), F(-1, 3), F(1, 4), F(-1, 4), F(2, 3), F(-2, 3), F(3, 4), F(-3, 4)])
        b = random.choice([F(1), F(1), F(2), F(3), F(-1), F(-2)])
        return {"kind": "geo", "a": a, "b": b}
    while True:
        L = random.choice([2, 3, 3])
        c = [random.choice([-3, -2, -1, 1, 2, 3])] + [random.randint(-3, 3) for _ in range(L - 2)] + [random.choice([-3, -2, -1, 1, 2, 3])]
        x0, xpi = sum(c), sum((-1) ** n * v for n, v in enumerate(c))
        if x0 != 0 and xpi != 0 and x0 != xpi and abs(x0) != abs(xpi):
            return {"kind": "fin", "c": [F(v) for v in c]}


def _build(base, kind):
    """(y_tex, correct spec, distractor specs, numeric (label_tex, suffix, value, slips), solution pieces)."""
    W = (0, "pi")
    geo = base["kind"] == "geo"
    if kind == "shift":
        k = random.choice([-3, -1, 1, 3])
        y_tex = f"x[n{'-' if k > 0 else '+'}{abs(k)}]"
        corr = [(1, 1, ("ejw", k), base, 1, W, 1)]
        dis = [[(1, 1, ("ejw", -k), base, 1, W, 1)],                       # sign of the exponent
               [(1, 1, None, base, 1, (-k, "rad"), 1)]]                    # X_d(w - k): time/frequency mix-up
        if geo:
            dis.append([(1, base["a"] ** k, ("ejw", k), base, 1, W, 1)])   # transform of a^n u[n-k]
        else:
            dis.append([(1, 1, ("ejw", k), base, -1, W, 1)])               # e^{-jkw} X_d(-w)
        xr, _ = _X_exact(base, "pi")
        val = (-1) ** (k % 2) * xr
        INFO.update(k=k, at="pi")
        prop = (r"\text{time shifting: } x[n-k] \;\longleftrightarrow\; e^{-jk\omega}X_d(\omega),\ k = " + str(k))
        num = (r"Y_d(\pi)", "", val, rf"Y_d(\pi) = e^{{-j({k})\pi}}X_d(\pi) = {(-1) ** (k % 2)}\cdot {fmt.tex_num(xr, paren=True)} = {fmt.tex_num(val)}")
        return y_tex, corr, dis, num, prop, "Lecture 14, Table 2 (time shifting)"
    if kind == "mod":
        w0 = random.choice([F(1, 2), F(-1, 2)])
        y_tex = rf"e^{{{'' if w0 > 0 else '-'}j\frac{{\pi}}{{2}}n}}\,x[n]"
        corr = [(1, 1, None, base, 1, (-w0, "pi"), 1)]
        dis = [[(1, 1, None, base, 1, (w0, "pi"), 1)],                     # X_d(w + w0)
               [(1, 1, ("const_pi", w0), base, 1, W, 1)],                  # e^{j w0} X_d(w)
               [(1, 1, None, base, -1, (w0, "pi"), 1)]]                    # X_d(w0 - w)
        x0, _ = _X_exact(base, "0")
        w0t = PI_TEX[w0]
        INFO.update(w0=fmt.plain_num(w0), at="w0")
        num = (rf"Y_d\!\left({w0t}\right)", "", x0, rf"Y_d\!\left({w0t}\right) = X_d\!\left({w0t} - \left({w0t}\right)\right) = X_d(0) = {fmt.tex_num(x0)}")
        prop = r"\text{frequency shifting: } e^{j\omega_0 n}x[n] \;\longleftrightarrow\; X_d(\omega - \omega_0),\ \omega_0 = " + w0t
        return y_tex, corr, dis, num, prop, "Lecture 14, Table 2 (frequency shifting); the same step as HW5 Problem 1(d)"
    if kind == "alt":
        ab = _alt(base)
        y_tex = r"(-1)^{n}\,x[n]"
        corr = [(1, 1, None, ab, 1, W, 1)]
        dis = [[(-1, 1, None, base, 1, W, 1)],                             # -X_d(w)
               [(1, 1, None, base, -1, W, 1)],                             # X_d(-w)
               [(1, 1, None, ab, -1, W, 1)]]                               # X_d(pi - w)
        at = random.choice(["0", "pi"])
        INFO.update(at=at)
        xr, _ = _X_exact(base, "pi" if at == "0" else "0")
        lab = r"Y_d(0)" if at == "0" else r"Y_d(\pi)"
        work = (rf"Y_d(0) = X_d(0 - \pi) = X_d(\pi) = {fmt.tex_num(xr)}" if at == "0"
                else rf"Y_d(\pi) = X_d(\pi - \pi) = X_d(0) = {fmt.tex_num(xr)}")
        num = (lab, "", xr, work)
        prop = (r"(-1)^n = e^{j\pi n}, \text{ so frequency shifting by } \omega_0 = \pi: \; (-1)^n x[n] \;\longleftrightarrow\; X_d(\omega - \pi)")
        return y_tex, corr, dis, num, prop, "Lecture 14, Table 2 (frequency shifting with $\\omega_0 = \\pi$)"
    if kind == "cos":
        h = (F(1, 2), "pi")
        y_tex = r"x[n]\cos\!\left(\tfrac{\pi}{2}n\right)"
        corr = [(1, 1, None, base, 1, (-F(1, 2), "pi"), 1), (1, 1, None, base, 1, h, 1)]
        dis = [("full", [(1, 1, None, base, 1, (-F(1, 2), "pi"), 1), (1, 1, None, base, 1, h, 1)]),   # no 1/2
               ("half", [(1, 1, None, base, 1, (-F(1, 2), "pi"), 1), (-1, 1, None, base, 1, h, 1)]),  # minus sign
               ("full", [(1, 1, None, base, 1, (-F(1, 2), "pi"), 1)])]                                # one term only
        at = random.choice(["0", "pi/2"])
        INFO.update(at=at)
        if at == "pi/2":
            x0, _ = _X_exact(base, "0")
            xp, _ = _X_exact(base, "pi")
            val = (x0 + xp) / 2
            work = (rf"Y_d\!\left(\tfrac{{\pi}}{{2}}\right) = \tfrac{{1}}{{2}}X_d(0) + \tfrac{{1}}{{2}}X_d(\pi) = \tfrac{{1}}{{2}}\left({fmt.tex_num(x0)}\right)"
                    rf" + \tfrac{{1}}{{2}}\left({fmt.tex_num(xp)}\right) = {fmt.tex_num(val)}")
            lab = r"Y_d\!\left(\tfrac{\pi}{2}\right)"
        else:
            xr, xi = _X_exact(base, "pi/2")
            val = xr                                   # (X(-pi/2) + X(pi/2)) / 2 = Re X(pi/2) for real x
            work = (rf"Y_d(0) = \tfrac{{1}}{{2}}X_d\!\left(-\tfrac{{\pi}}{{2}}\right) + \tfrac{{1}}{{2}}X_d\!\left(\tfrac{{\pi}}{{2}}\right)"
                    rf" = \operatorname{{Re}}\,X_d\!\left(\tfrac{{\pi}}{{2}}\right) = {fmt.tex_num(val)}"
                    rf"\quad\left(X_d\!\left(\tfrac{{\pi}}{{2}}\right) = {_cplx_tex(xr, xi)}\right)")
            lab = r"Y_d(0)"
        num = (lab, "", val, work)
        prop = (r"\text{modulation: } x[n]\cos(\omega_0 n) \;\longleftrightarrow\; \tfrac{1}{2}X_d(\omega - \omega_0) + \tfrac{1}{2}X_d(\omega + \omega_0),"
                r"\ \omega_0 = \tfrac{\pi}{2}")
        return y_tex, ("half", corr), dis, num, prop, "Lecture 14, Table 2 (modulation); Euler's formula as in HW5 Problem 1(c)"
    if kind == "rev":
        y_tex = "x[-n]"
        corr = [(1, 1, None, base, -1, W, 1)]
        dis = [[(1, 1, None, base, 1, W, 1)],                              # unchanged
               [(1, 1, None, _alt(base), 1, W, 1)]]                        # (-1)^n confusion
        dis.append([(1, 1, None, _recip(base), 1, W, 1)] if geo else [(-1, 1, None, base, 1, W, 1)])
        en = _energy(base)
        val = 2 * en
        INFO.update(at="energy")
        work = (rf"\int_{{-\pi}}^{{\pi}} |Y_d(\omega)|^2\,d\omega = \int_{{-\pi}}^{{\pi}} |X_d(-\omega)|^2\,d\omega = 2\pi\sum_n |x[n]|^2"
                rf" = 2\pi\cdot {fmt.tex_num(en)} = {fmt.tex_num(val)}\,\pi")
        num = (r"\displaystyle\int_{-\pi}^{\pi} \left|Y_d(\omega)\right|^2 d\omega", r"\cdot\,\pi", val, work)
        prop = r"\text{time reversal: } x[-n] \;\longleftrightarrow\; X_d(-\omega)"
        return y_tex, corr, dis, num, prop, "Lecture 14, Table 2 (time reversal, Parseval's relation); Parseval as in HW5 Problem 2(d)"
    # conv
    y_tex = r"x[n] * x[n]"
    corr = [(1, 1, None, base, 1, W, 2)]
    dis = [[(1, 2, None, base, 1, W, 1)],                                  # 2 X_d(w)
           [(1, 1, None, base, 2, W, 1)],                                  # X_d(2w)
           [(1, 1, None, _square_samples(base), 1, W, 1)]]                 # DTFT of x[n]^2 (multiplied, not convolved)
    at = random.choice(["0", "pi"])
    INFO.update(at=at)
    xr, _ = _X_exact(base, at)
    val = xr * xr
    lab = r"Y_d(0)" if at == "0" else r"Y_d(\pi)"
    wt = "0" if at == "0" else r"\pi"
    work = rf"{lab} = X_d({wt})^2 = {fmt.tex_num(xr, paren=True)}^2 = {fmt.tex_num(val)}"
    num = (lab, "", val, work)
    prop = r"\text{convolution: } x[n] * h[n] \;\longleftrightarrow\; X_d(\omega)H_d(\omega),\ h = x"
    return y_tex, corr, dis, num, prop, "Lecture 14, Table 2 (convolution)"


def _cplx_tex(re, im):
    if im == 0:
        return fmt.tex_num(re)
    jt = ("" if abs(im) == 1 else fmt.tex_num(abs(im))) + "j"
    if re == 0:
        return ("-" if im < 0 else "") + jt
    return f"{fmt.tex_num(re)} {'-' if im < 0 else '+'} {jt}"


def _make(spec):
    if isinstance(spec, tuple) and spec[0] in ("half", "full"):
        return _option(spec[1], halves=(spec[0] == "half"))
    return _option(spec)


def _differs(f, g):
    return max(abs(f(w) - g(w)) for w in W_PROBE) > 1e-6


INFO = {}


def generate(data):
    while True:
        INFO.clear()
        base = _draw_base()
        kind = random.choice(["shift", "mod", "alt", "cos", "rev", "conv"])
        if kind == "cos" and base["kind"] == "fin" and len(base["c"]) > 2:
            continue                                                   # keep the option texts short
        y_tex, corr, dis, num, prop, src = _build(base, kind)
        ctex, cval = _make(corr)
        opts = []
        for d in dis:
            t, v = _make(d)
            if _differs(v, cval) and all(_differs(v, ov) for _, ov in opts) and t != ctex:
                opts.append((t, v))
        if len(opts) == 3:
            break

    p = data["params"]
    p["kind"] = kind
    p["base_kind"] = base["kind"]
    if base["kind"] == "geo":
        a, b = base["a"], base["b"]
        p["a"], p["b"] = fmt.plain_num(a), fmt.plain_num(b)
        lead = "" if b == 1 else ("-" if b == -1 else fmt.tex_num(b))
        p["x_tex"] = lead + fmt.tex_pow(a) + "u[n]"
        p["X_tex"] = ("-" if b < 0 else "") + fmt.tex_frac(fmt.tex_num(abs(b)), fmt.tex_sum([(1, ""), (-a, r"e^{-j\omega}")]))
    else:
        p["c"] = [int(v) for v in base["c"]]
        p["x_tex"] = fmt.tex_seq(base["c"], 0)
        p["X_tex"] = fmt.tex_sum((v, _exp_tex(n, 1, (0, "pi"))) for n, v in enumerate(base["c"]))
    p["y_tex"] = y_tex
    p["Y_tex"] = ctex
    p["yd_choices"] = [{"text": f"${ctex}$", "correct": True}] + [{"text": f"${t}$", "correct": False} for t, _ in opts]
    random.shuffle(p["yd_choices"])
    lab, suffix, val, work = num
    p["num_label"] = lab
    p["num_suffix"] = suffix
    p["has_suffix"] = bool(suffix)
    p["num_work_tex"] = work
    p["prop_tex"] = prop
    p["source_html"] = src
    p["is_fin"] = base["kind"] == "fin"
    p["k"] = INFO.get("k", 0)
    p["w0"] = INFO.get("w0", "0")
    p["at"] = INFO["at"]
    data["correct_answers"]["num"] = float(val)
