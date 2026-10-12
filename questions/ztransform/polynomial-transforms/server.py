"""Finite sequences <-> polynomials in z and z^{-1} (Lecture 6 section 2 definition; Lecture 7 notes
section 1 ROC rule 5 and slide 6 for finite-length ROCs; Lecture 7 convolution property; HW3 #1(a), #1(e),
#2(e), #3(a); FA2025 Midterm 1 #5(b), FA2019 #5).

One variant = an input x[n] (3-4 integer samples, possibly starting at a negative index) and an FIR
system given by its transfer function H(z), a Laurent polynomial such as 2z^2 - 1 + 3z^{-3}
(2-3 terms, usually with missing powers in between). Asked:
  (a) X(z)            pl-symbolic-input in z; SymPy proves Laurent-polynomial identities, so any
                      equivalent form (positive powers over z^k, 1/z, ...) is accepted;
  (b) ROC of X(z)     4 fixed options: all z / all z except 0 / except infinity / except both;
  (c) h[n] from H(z)  first index (integer) + samples (row vector, zeros included);
  (d) y[k]            one output sample of y = x * h, i.e. the coefficient of z^{-k} in X(z)H(z),
                      chosen with >= 2 contributing products, y[k] != 0 and k != 0;
  (e) Y(z0)           z0 = 1 (= sum_n y[n] = X(1)H(1)) or z0 = -1 (alternating sums).

Supports: x is never even about n = 0 (so the z^{+n} sign mix-up is always wrong); the support of h
is never symmetric about n = 0 (so the mirrored first index is always wrong). The ROC family of x
(both ends / only z = 0 / only z = infinity excluded) is drawn first, 40/30/30."""

import random

import numpy as np
import prairielearn as pl
from ece310 import fmt, seq

NZ = [-3, -2, -1, 1, 2, 3]

ROC_OPTS = [
    ("all", r"all $z$"),
    ("no0", r"all $z$ except $z = 0$, i.e. $\lvert z\rvert \gt 0$"),
    ("noinf", r"all $z$ except $z = \infty$, i.e. $\lvert z\rvert \lt \infty$"),
    ("noboth", r"all $z$ except $z = 0$ and $z = \infty$, i.e. $0 \lt \lvert z\rvert \lt \infty$"),
]


def _val(vals, start, n):
    i = n - start
    return vals[i] if 0 <= i < len(vals) else 0


def _draw_x():
    fam = random.choices(["noboth", "no0", "noinf"], weights=[40, 30, 30])[0]
    L = random.choice([3, 3, 4])
    while True:
        vals = [random.randint(-3, 3) for _ in range(L)]
        vals[0], vals[-1] = random.choice(NZ), random.choice(NZ)
        if fam == "noboth":
            start = -1 if L == 3 else random.choice([-2, -1])
        elif fam == "no0":
            start = random.choice([0, 0, 1])
        else:
            start = random.choice([0, 0, -1]) - L + 1
        d = {start + i: v for i, v in enumerate(vals)}
        if all(d.get(n, 0) == d.get(-n, 0) for n in d):     # even about n = 0
            continue
        return fam, vals, start


def _draw_h():
    while True:
        T = random.choice([2, 3, 3])
        ns = sorted(random.sample(range(-3, 4), T))
        span = ns[-1] - ns[0]
        if not 2 <= span <= 5 or ns[0] == -ns[-1]:
            continue
        if ns[0] >= 0 and random.random() < 0.75:           # mostly at least one positive power of z
            continue
        vals = [0] * (span + 1)
        for n in ns:
            vals[n - ns[0]] = random.choice(NZ)
        return vals, ns[0]


def _mistakes_yk(x, nx, h, nh, y, ny, k):
    """Classic wrong values for y[k]: coefficient of z^{+k} (= y[-k]), correlation (h not flipped),
    pointwise product x[k] h[k]."""
    right = _val(y, ny, k)
    cands = [_val(y, ny, -k),
             sum(_val(x, nx, m) * _val(h, nh, m - k) for m in range(nx, nx + len(x))),
             _val(x, nx, k) * _val(h, nh, k)]
    out = []
    for c in cands:
        if c != right and c not in out:
            out.append(c)
    return out


def _sgn(n):
    """(-1)^n as an int, for any integer n."""
    return 1 if n % 2 == 0 else -1


def _zval(vals, start, z0):
    """X(z0) for z0 = +-1 (exact integer)."""
    return sum(vals) if z0 == 1 else sum(v * _sgn(start + i) for i, v in enumerate(vals))


def _mistakes_Y(x, nx, h, nh, z0):
    X0, H0 = _zval(x, nx, z0), _zval(h, nh, z0)
    right = X0 * H0
    if z0 == 1:
        cands = [X0 + H0, X0, H0]                       # added instead of multiplied; forgot a factor
    else:
        cands = [_zval(x, nx, 1) * _zval(h, nh, 1), X0 + H0, X0]   # forgot the (-1)^n; added
    out = []
    for c in cands:
        if c != right and c not in out:
            out.append(c)
    return out


def _draw():
    while True:
        fam, x, nx = _draw_x()
        h, nh = _draw_h()
        yF, ny = seq.conv(x, nx, h, nh)
        y = [int(v) for v in yF]
        ks = []
        for k in range(ny, ny + len(y)):
            if k == 0 or _val(y, ny, k) == 0:
                continue
            prods = [m for m in range(nx, nx + len(x)) if _val(x, nx, m) * _val(h, nh, k - m) != 0]
            if len(prods) >= 2 and len(_mistakes_yk(x, nx, h, nh, y, ny, k)) >= 2:
                ks.append(k)
        z0s = [z0 for z0 in (1, 1, -1) if _zval(x, nx, z0) * _zval(h, nh, z0) != 0
               and len(_mistakes_Y(x, nx, h, nh, z0)) >= 2]
        if ks and z0s:
            return fam, x, nx, h, nh, y, ny, random.choice(ks), random.choice(z0s)


def _laurent_tex(vals, start):
    """sum_n vals[n - start] z^{-n}: 2z^{2} - 1 + 3z^{-3} (ordered by increasing n)."""
    return fmt.tex_sum((v, fmt.zpow_tex(start + i)) for i, v in enumerate(vals))


def _laurent_sym(vals, start):
    terms = [f"{fmt.sym_num(v)}*{fmt.zpow_sym(start + i)}" for i, v in enumerate(vals) if v != 0]
    return " + ".join(terms)


def generate(data):
    fam, x, nx, h, nh, y, ny, k, z0 = _draw()
    ex, eh = nx + len(x) - 1, nh + len(h) - 1
    p = data["params"]
    p["x_vals"], p["x_start"], p["h_vals"], p["h_start"] = x, nx, h, nh
    p["family"], p["k"], p["z0"] = fam, k, z0
    p["x_tex"] = fmt.tex_seq(x, nx)
    p["H_tex"] = _laurent_tex(h, nh)
    p["X_tex"] = _laurent_tex(x, nx)
    p["X_def_tex"] = " + ".join(rf"x[{nx + i}]\,{fmt.zpow_tex(nx + i)}" if nx + i != 0 else "x[0]"
                                for i in range(len(x)) if x[i] != 0)
    p["k_tex"] = str(k)
    p["zk_tex"] = fmt.zpow_tex(k)
    p["z0_tex"] = str(z0)

    # (b) ROC
    if fam == "noboth":
        why = (rf"$x[n]$ has samples at $n \lt 0$ (here down to $n = {nx}$), which give positive powers of $z$ "
               rf"and blow up as $z \to \infty$, and samples at $n \gt 0$ (up to $n = {ex}$), which give negative "
               r"powers of $z$ and blow up at $z = 0$. Both points are excluded: $0 \lt \lvert z\rvert \lt \infty$.")
    elif fam == "no0":
        why = (rf"Every sample sits at $n \ge 0$ and some at $n \gt 0$ (up to $n = {ex}$), so $X(z)$ contains "
               r"negative powers of $z$ only (plus perhaps a constant): they blow up at $z = 0$, while "
               r"$X(\infty) = x[0]$ is finite. ROC: $\lvert z\rvert \gt 0$, all $z$ except $z = 0$.")
    else:
        why = (rf"Every sample sits at $n \le 0$ and some at $n \lt 0$ (down to $n = {nx}$), so $X(z)$ contains "
               r"positive powers of $z$ only (plus perhaps a constant): $X(0) = x[0]$ is finite, but $X(z) \to \infty$ "
               r"as $z \to \infty$. ROC: $\lvert z\rvert \lt \infty$, all $z$ except $z = \infty$.")
    p["roc_why_html"] = why
    p["roc_choices"] = [{"text": t, "correct": key == fam} for key, t in ROC_OPTS]

    # (c) h[n] from H(z)
    rows = []
    for i, v in enumerate(h):
        if v == 0:
            continue
        n = nh + i
        term = fmt.tex_sum([(v, fmt.zpow_tex(n))])
        rows.append(f"<tr><td>${term}$</td><td>$z^{{-n}}$ with $n = {n}$</td><td>$h[{n}] = {fmt.tex_num(v)}$</td></tr>")
    p["h_table_html"] = ('<table class="table table-sm table-bordered text-center" style="width:auto">'
                         "<tr><th>term of $H(z)$</th><th>power</th><th>sample</th></tr>" + "".join(rows) + "</table>")
    p["h_tex"] = fmt.tex_seq(h, nh)
    p["h_last"] = eh
    p["h_len"] = len(h)
    p["h_plain"] = "[" + ", ".join(str(v) for v in h) + "]"
    p["has_gaps"] = any(v == 0 for v in h)

    # (d) y[k]
    ms = [m for m in range(nx, ex + 1) if _val(x, nx, m) * _val(h, nh, k - m) != 0]
    p["yk_sum_tex"] = " + ".join(rf"x[{m}]\,h[{k - m}]" for m in ms)
    p["yk_vals_tex"] = " + ".join(rf"{fmt.tex_num(_val(x, nx, m), paren=True)}\cdot{fmt.tex_num(_val(h, nh, k - m), paren=True)}"
                                  for m in ms)
    p["yk_pairs_tex"] = r",\ ".join(rf"({m},\,{k - m})" for m in ms)
    p["yk_tex"] = str(_val(y, ny, k))
    p["y_tex"] = fmt.tex_seq(y, ny)
    p["Y_tex"] = _laurent_tex(y, ny)

    # (e) Y(z0)
    X0, H0 = _zval(x, nx, z0), _zval(h, nh, z0)
    if z0 == 1:
        p["Xz0_tex"] = " + ".join(fmt.tex_num(v, paren=True) for v in x if v != 0) + f" = {X0}"
        p["Hz0_tex"] = " + ".join(fmt.tex_num(v, paren=True) for v in h if v != 0) + f" = {H0}"
        p["z0_rule_tex"] = r"z^{-n} = 1 \text{ for every } n\text{, so } X(1) = \sum_n x[n]"
    else:
        def alt(vals, start):
            return " + ".join(rf"{fmt.tex_num(v, paren=True)}\cdot(-1)^{{{start + i}}}" for i, v in enumerate(vals) if v != 0)
        p["Xz0_tex"] = alt(x, nx) + f" = {X0}"
        p["Hz0_tex"] = alt(h, nh) + f" = {H0}"
        p["z0_rule_tex"] = r"(-1)^{-n} = (-1)^{n}\text{, so } X(-1) = \sum_n (-1)^n x[n]"
    p["X0"], p["H0"] = X0, H0
    p["Yz0_tex"] = str(X0 * H0)
    p["sum_check_tex"] = str(sum(y)) if z0 == 1 else str(sum(v * _sgn(ny + i) for i, v in enumerate(y)))

    c = data["correct_answers"]
    c["X"] = _laurent_sym(x, nx)
    c["nh"] = nh
    c["h"] = pl.to_json(np.array([[float(v) for v in h]]))
    c["yk"] = int(_val(y, ny, k))
    c["Yz"] = int(X0 * H0)
