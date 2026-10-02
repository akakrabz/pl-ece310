"""z-transform and ROC of a signal made of one or two standard pieces (exam family "z-transform
with ROC", 6 of 7 past Midterm 1 exams; HW3 #1; Lectures 6-7).

Pieces: right-sided exponential with a shifted step A a^n u[n-k] (needs the shift factor a^k),
left-sided exponential B b^n u[-n+m], n a^n u[n], a short finite sequence, and (in a minority of
variants) a damped cosine. Pole magnitudes are distinct and the ROC is never empty, so the ROC is
unambiguous. The student types X(z) (pl-symbolic-input) and picks the ROC (multiple choice: the
correct ROC, 0/infinity membership toggled, side flipped, wrong radius)."""

import random
from fractions import Fraction as F

from ece310 import fmt, mc, zt

MAGS = [F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(4, 3), F(3, 2), F(2), F(3)]
COEFS = [F(1), F(1), F(2), F(3), F(-1), F(-2), F(1, 2)]
COS_W = ["pi/2", "pi/3", "2*pi/3"]
COS_VAL = {"pi/2": F(0), "pi/3": F(1, 2), "2*pi/3": F(-1, 2)}
COS_TEX = {"pi/2": r"\tfrac{\pi}{2}", "pi/3": r"\tfrac{\pi}{3}", "2*pi/3": r"\tfrac{2\pi}{3}"}
LZ = r"\lvert z\rvert"


def _base(exclude=(), allow_one=False):
    pool = [m for m in MAGS + ([F(1)] if allow_one else []) if m not in exclude]
    m = random.choice(pool)
    return -m if (m != 1 and random.random() < 0.3) else m


def _finite(start_choices):
    length = random.choice([2, 3])
    vals = [random.randint(-3, 3) for _ in range(length)]
    vals[0] = random.choice([-3, -2, -1, 1, 2, 3])
    vals[-1] = random.choice([-3, -2, -1, 1, 2, 3])
    return vals, random.choice(start_choices)


# ----------------------------------------------------------------------------- explanation per term
def _right_step(A, a, k):
    t = zt.term_right_exp(A, a, k)
    c = A * a**k
    arg = zt.u_arg_right(k)
    if k != 0 and a != 1:
        rew = t.x_tex + " = " + fmt.tex_sum([(c, fmt.tex_pow(a, arg) + r"\,u[" + arg + "]")])
        powk = fmt.tex_num(a, paren=True) if k == 1 else fmt.tex_pow(a, str(k))
        const = powk if A == 1 else fmt.tex_num(A, paren=True) + r"\cdot " + powk
        const = const if const in (fmt.tex_num(c), fmt.tex_num(c, paren=True)) else const + " = " + fmt.tex_num(c)
        how = (rf"Make the exponent match the step: $a^{{n}}u[n-k] = a^{{k}}\cdot a^{{n-k}}u[n-k]$ with "
               rf"$a = {fmt.tex_num(a)}$, $k = {k}$, so the constant is ${const}$. "
               rf"Then $a^{{n}}u[n] \leftrightarrow \frac{{1}}{{1-az^{{-1}}}}$, and the shift by $k = {k}$ multiplies by ${fmt.zpow_tex(k)}$.")
    elif k != 0:
        rew = t.x_tex
        how = rf"$u[n] \leftrightarrow \frac{{1}}{{1-z^{{-1}}}}$, and the shift by ${k}$ multiplies by ${fmt.zpow_tex(k)}$."
    else:
        rew = t.x_tex
        how = r"Table pair $a^{n}u[n] \leftrightarrow \frac{1}{1-az^{-1}}$."
    if k < 0:
        rnote = (rf"right-sided, pole at $z = {fmt.tex_num(a)}$: ${LZ} \gt {fmt.tex_num(abs(a))}$; it starts at "
                 rf"$n = {k} \lt 0$ (positive powers of $z$), so $z = \infty$ is excluded")
    else:
        rnote = (rf"right-sided, pole at $z = {fmt.tex_num(a)}$: ${LZ} \gt {fmt.tex_num(abs(a))}$; it starts at "
                 rf"$n = {k} \ge 0$, so $z = \infty$ is included")
    return t, {"rew_tex": rew, "how": how, "X_tex": t.X_tex, "roc_tex": t.roc.tex(), "rnote": rnote}


def _left_step(B, b, m):
    t = zt.term_left_exp(B, b, m)
    j = m + 1
    c = -B * b**j
    if j != 0:
        inner = fmt.tex_pow(b, f"n-{j}" if j > 0 else f"n+{-j}") + r"\,u[-(" + (f"n-{j}" if j > 0 else f"n+{-j}") + ")-1]"
        rew = t.x_tex + " = " + fmt.tex_sum([(B * b**j, inner)])
        how = (rf"Rewrite with the standard left-sided pair $-b^{{n}}u[-n-1] \leftrightarrow \frac{{1}}{{1-bz^{{-1}}}}$, "
               rf"$b = {fmt.tex_num(b)}$: the step $u[{zt.u_arg_left(m)}]$ is $u[-n-1]$ shifted by ${j}$, which costs the constant "
               rf"${'b' if j == 1 else 'b^{' + str(j) + '}'}$ and the factor ${fmt.zpow_tex(j)}$; the pair's minus sign gives the constant ${fmt.tex_num(c)}$.")
    else:
        rew = t.x_tex
        how = (rf"Standard left-sided pair $-b^{{n}}u[-n-1] \leftrightarrow \frac{{1}}{{1-bz^{{-1}}}}$ with $b = {fmt.tex_num(b)}$ "
               r"(note the minus sign).")
    if m >= 1:
        rnote = (rf"left-sided, pole at $z = {fmt.tex_num(b)}$: ${LZ} \lt {fmt.tex_num(abs(b))}$; it ends at "
                 rf"$n = {m} \gt 0$ (negative powers of $z$), so $z = 0$ is excluded")
    else:
        rnote = (rf"left-sided, pole at $z = {fmt.tex_num(b)}$: ${LZ} \lt {fmt.tex_num(abs(b))}$; it ends at "
                 rf"$n = {m} \le 0$, so $z = 0$ is included")
    return t, {"rew_tex": rew, "how": how, "X_tex": t.X_tex, "roc_tex": t.roc.tex(), "rnote": rnote}


def _nexp_step(A, a):
    t = zt.term_n_exp(A, a)
    how = (r"Differentiation property: $n\,x[n] \leftrightarrow -z\frac{dX}{dz}$ applied to "
           rf"$a^{{n}}u[n] \leftrightarrow \frac{{1}}{{1-az^{{-1}}}}$ gives $\frac{{az^{{-1}}}}{{(1-az^{{-1}})^{{2}}}}$ (double pole at $z = {fmt.tex_num(a)}$).")
    rnote = rf"right-sided from $n = 0$, double pole at $z = {fmt.tex_num(a)}$: ${LZ} \gt {fmt.tex_num(abs(a))}$"
    return t, {"rew_tex": t.x_tex, "how": how, "X_tex": t.X_tex, "roc_tex": t.roc.tex(), "rnote": rnote}


def _finite_step(vals, start):
    t = zt.term_finite(vals, start)
    first, last = start, start + len(vals) - 1
    how = r"A finite sequence: each $c\,\delta[n-n_0]$ contributes $c\,z^{-n_0}$."
    parts = []
    parts.append(rf"samples at $n \gt 0$ exclude $z = 0$" if last > 0 else r"no samples at $n \gt 0$, so $z = 0$ is included")
    parts.append(rf"samples at $n \lt 0$ exclude $z = \infty$" if first < 0 else r"no samples at $n \lt 0$, so $z = \infty$ is included")
    rnote = "finite length: " + "; ".join(parts)
    return t, {"rew_tex": t.x_tex, "how": how, "X_tex": t.X_tex, "roc_tex": t.roc.tex(), "rnote": rnote}


def _cos_step(A, r, w):
    t = zt.term_cos(A, r, w)
    c = COS_VAL[w]
    how = (r"Table pair $r^{n}\cos(\omega_0 n)u[n] \leftrightarrow \frac{1-r\cos\omega_0\,z^{-1}}{1-2r\cos\omega_0\,z^{-1}+r^{2}z^{-2}}$ "
           rf"with $r = {fmt.tex_num(r)}$, $\cos\omega_0 = \cos{COS_TEX[w]} = {fmt.tex_num(c)}$.")
    rnote = rf"right-sided from $n = 0$, poles at $z = {fmt.tex_num(r)}e^{{\pm j{COS_TEX[w]}}}$ (magnitude ${fmt.tex_num(r)}$): ${LZ} \gt {fmt.tex_num(r)}$"
    return t, {"rew_tex": t.x_tex, "how": how, "X_tex": t.X_tex, "roc_tex": t.roc.tex(), "rnote": rnote}


# ----------------------------------------------------------------------------- ROC distractors
def roc_distractors(roc, radii):
    """Wrong ROCs ordered by priority: 0/infinity toggled, side flipped, then wrong radii."""
    R = zt.ROC
    first, rest = [], []
    k = roc.kind()
    if k == "right":
        r = roc.inner
        first.append(R(inner=r, outer=None, has0=False, hasinf=not roc.hasinf))
        first.append(R(inner=F(0), outer=r, has0=True, hasinf=False))
        for w in [1 / r, F(1)] + list(radii):
            if w != r and w != 0:
                rest.append(R(inner=w, outer=None, has0=False, hasinf=roc.hasinf))
        rest.append(R(inner=F(0), outer=r, has0=False, hasinf=False))
        rest.append(R(inner=F(0), outer=None, has0=False, hasinf=roc.hasinf))
    elif k == "left":
        r = roc.outer
        first.append(R(inner=F(0), outer=r, has0=not roc.has0, hasinf=False))
        first.append(R(inner=r, outer=None, has0=False, hasinf=True))
        for w in [1 / r, F(1)] + list(radii):
            if w != r and w != 0:
                rest.append(R(inner=F(0), outer=w, has0=roc.has0, hasinf=False))
        rest.append(R(inner=r, outer=None, has0=False, hasinf=False))
        rest.append(R(inner=F(0), outer=None, has0=roc.has0, hasinf=False))
    elif k == "two-sided":
        r1, r2 = roc.inner, roc.outer
        first.append(R(inner=r1, outer=None, has0=False, hasinf=True))
        first.append(R(inner=F(0), outer=r2, has0=True, hasinf=False))
        cands = sorted({1 / r1, 1 / r2, r1, r2, F(1)} | set(radii))
        for lo in cands:
            for hi in cands:
                if 0 < lo < hi and (lo, hi) != (r1, r2):
                    rest.append(R(inner=lo, outer=hi, has0=False, hasinf=False))
        rest.append(R(inner=r1, outer=None, has0=False, hasinf=False))
        rest.append(R(inner=F(0), outer=r2, has0=False, hasinf=False))
    else:  # finite shapes
        for h0 in (True, False):
            for hi in (True, False):
                if (h0, hi) != (roc.has0, roc.hasinf):
                    first.append(R(inner=F(0), outer=None, has0=h0, hasinf=hi))
        for w in list(radii) + [F(1)]:
            rest.append(R(inner=w, outer=None, has0=False, hasinf=True))
    random.shuffle(rest)
    return first + rest


def roc_choices(roc, radii, n=5):
    opts = roc_distractors(roc, radii)
    by_text = {}
    for o in [roc] + opts:
        by_text.setdefault("$" + o.tex() + "$", o)
    texts = ["$" + o.tex() + "$" for o in opts]
    ch = mc.choices("$" + roc.tex() + "$", texts, n=n)
    for c in ch:
        c["roc"] = by_text[c["text"]].to_json()
    return ch


def ends_agree(terms, roc):
    """Linearity only promises 'ROC at least the intersection': if samples of different pieces cancel
    at an end of the summed signal, the true ROC can gain z = 0 or z = infinity. True when the summed
    signal's first/last nonzero index agrees with the intersection's 0/infinity membership."""
    x = {n: sum((zt.x_value(t.spec, n) for t in terms), F(0)) for n in range(-20, 21)}
    if roc.inner == 0 and roc.has0 != all(x[n] == 0 for n in range(1, 21)):
        return False
    if roc.outer is None and roc.hasinf != all(x[n] == 0 for n in range(-20, 0)):
        return False
    return True


# ----------------------------------------------------------------------------- generate
def _pick_terms():
    fam = random.choice(["R", "R", "L", "L", "RL", "RL", "RL", "RF", "LF", "NL", "NR", "C"])
    steps = []
    if fam == "R":
        a = _base(allow_one=random.random() < 0.15)
        steps.append(_right_step(random.choice(COEFS), a, random.choice([-2, -1, 1, 2, 3])))
    elif fam == "L":
        steps.append(_left_step(random.choice(COEFS), _base(), random.choice([-1, 0, 1, 2])))
    elif fam in ("RL", "NL", "C"):
        while True:
            a, b = _base(), _base()
            if abs(a) < abs(b):
                break
        if fam == "RL":
            steps.append(_right_step(random.choice(COEFS), a, random.choice([-1, 0, 1, 2])))
            steps.append(_left_step(random.choice(COEFS), b, random.choice([-1, 0, 1])))
        elif fam == "NL":
            steps.append(_nexp_step(random.choice([F(1), F(2), F(-1), F(3)]), a))
            steps.append(_left_step(random.choice(COEFS), b, random.choice([-1, 0])))
        else:
            r = random.choice([F(1, 2), F(2, 3), F(3, 4), F(1)])
            steps.append(_cos_step(random.choice([F(1), F(2), F(-1)]), r, random.choice(COS_W)))
            if random.random() < 0.6:
                bb = random.choice([m for m in MAGS if m > r])
                steps.append(_left_step(random.choice(COEFS), bb if random.random() < 0.7 else -bb, random.choice([-1, 0])))
            else:
                vals, s = _finite([-2, -1])
                steps.append(_finite_step(vals, s))
    elif fam == "RF":
        vals, s = _finite([-2, -1, 0, 1])
        steps.append(_right_step(random.choice(COEFS), _base(allow_one=random.random() < 0.15), random.choice([0, 1, 2])))
        steps.append(_finite_step(vals, s))
    elif fam == "LF":
        vals, s = _finite([-1, 0, 1, 2])
        steps.append(_left_step(random.choice(COEFS), _base(), random.choice([-1, 0])))
        steps.append(_finite_step(vals, s))
    else:  # NR
        a = _base()
        b = _base(exclude=(abs(a),))
        steps.append(_nexp_step(random.choice([F(1), F(2), F(-1), F(3)]), a))
        steps.append(_right_step(random.choice(COEFS), b, random.choice([-1, 0, 1, 2])))
    random.shuffle(steps)
    return fam, steps


def generate(data):
    while True:
        fam, steps = _pick_terms()
        terms = [t for t, _ in steps]
        x_tex, X_tex, X_sym, roc = zt.sum_terms(terms)
        radii = sorted({zt.mag(p) if not isinstance(p, complex) else F(abs(p)).limit_denominator(100) for t in terms for p in t.poles})
        if roc.empty or not ends_agree(terms, roc):
            continue
        choices = roc_choices(roc, radii)
        if len(choices) >= 4:
            break

    p = data["params"]
    p["fam"] = fam
    p["x_tex"] = x_tex
    p["X_tex"] = X_tex
    p["roc_tex"] = roc.tex()
    p["specs"] = [t.spec for t in terms]
    p["two"] = len(terms) == 2
    p["steps"] = [dict(s, i=i + 1) for i, (_, s) in enumerate(steps)]
    p["roc_choices"] = choices
    data["correct_answers"]["X"] = X_sym
