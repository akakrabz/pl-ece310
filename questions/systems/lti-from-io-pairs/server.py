"""What do three input-output pairs of an unknown system prove? (Lecture 3 definitions of
linearity and time-invariance, Lecture 4 "an LTI system is fully described by h[n]"; HW2 #4;
the Midterm 1 True/False statements "the output of any system is x * h" (FA2023 #1(d), FA2024 #1(a)).)

Pair 1 is always a shifted impulse x1 = delta[n - m] (m != 0), so an LTI system would need
h[n] = y1[n + m]. Pairs 2 and 3 have short multi-sample inputs; x3 is related to x2 by a shift, a
scaling, a sum with a (shifted) copy of x1, a scaled shift, or not at all. Four families:

  lti    every pair consistent with h (any relation type)                -> "consistent with an LTI system"
  notL   x3 = c x2 or x3 = x2 + d x1, outputs break that linear relation -> "not linear" (TI undecided)
  notTI  x3 = x2[n - k], y3 != y2[n - k]                                 -> "not time-invariant" (L undecided)
  joint  x3 unrelated / c x2[n - k] / x2 + d x1[n - k]; no linear relation and no pure shift among the
         inputs, but some y_i != x_i * h                                  -> "cannot be LTI, neither property alone fails"

Every draw is classified exactly (rank over Q of the inputs vs inputs+outputs for linearity, pure
shifts for time-invariance, x_i * h for LTI) and redrawn unless the class equals the family, so
exactly one of the four options is right. Part (b), meaningful in every variant: the output an LTI
system with pair 1's h would give for a new input v (first index + row vector)."""

import random
from fractions import Fraction

import numpy as np
import prairielearn as pl
from ece310 import fmt

OPTIONS = [
    ("lti", "The pairs are consistent with an LTI system."),
    ("notL", "The system is not linear (whether it is time-invariant cannot be decided from these pairs)."),
    ("notTI", "The system is not time-invariant (whether it is linear cannot be decided from these pairs)."),
    ("joint", "The system cannot be LTI, but these pairs show neither that it is nonlinear nor that it is time-varying."),
]


# ----------------------------------------------------------------------------- finite sequences (vals, start)
def trim(vals, start):
    vals = list(vals)
    while vals and vals[0] == 0:
        vals.pop(0)
        start += 1
    while vals and vals[-1] == 0:
        vals.pop()
    return (vals, start) if vals else ([], 0)


def shift(s, k):
    return (list(s[0]), s[1] + k) if s[0] else s


def scale(s, c):
    return trim([c * v for v in s[0]], s[1])


def add(a, b):
    if not a[0]:
        return b
    if not b[0]:
        return a
    lo = min(a[1], b[1])
    hi = max(a[1] + len(a[0]), b[1] + len(b[0]))
    out = [0] * (hi - lo)
    for vals, st in (a, b):
        for i, v in enumerate(vals):
            out[st - lo + i] += v
    return trim(out, lo)


def conv(a, b):
    if not a[0] or not b[0]:
        return ([], 0)
    out = [0] * (len(a[0]) + len(b[0]) - 1)
    for i, u in enumerate(a[0]):
        for j, v in enumerate(b[0]):
            out[i + j] += u * v
    return trim(out, a[1] + b[1])


def same(a, b):
    return trim(*a) == trim(*b)


def shift_amount(a, b):
    """k with b = a[n - k], or None."""
    a, b = trim(*a), trim(*b)
    if a[0] and a[0] == b[0]:
        return b[1] - a[1]
    return None


def rank(rows):
    """Exact rank over Q of integer vectors."""
    m = [[Fraction(v) for v in r] for r in rows]
    rk, col, ncol = 0, 0, len(m[0]) if m else 0
    while rk < len(m) and col < ncol:
        piv = next((i for i in range(rk, len(m)) if m[i][col] != 0), None)
        if piv is None:
            col += 1
            continue
        m[rk], m[piv] = m[piv], m[rk]
        for i in range(len(m)):
            if i != rk and m[i][col] != 0:
                f = m[i][col] / m[rk][col]
                m[i] = [a - f * b for a, b in zip(m[i], m[rk])]
        rk += 1
        col += 1
    return rk


def vectors(seqs):
    lo = min(s[1] for s in seqs if s[0])
    hi = max(s[1] + len(s[0]) for s in seqs if s[0])
    out = []
    for vals, st in seqs:
        row = [0] * (hi - lo)
        for i, v in enumerate(vals):
            row[st - lo + i] = v
        out.append(row)
    return out


def classify(xs, ys, h):
    vx = vectors(xs)
    vxy = [a + b for a, b in zip(vx, vectors(ys))]
    l_dis = rank(vxy) > rank(vx)
    ti_dis = False
    for i in range(len(xs)):
        for j in range(len(xs)):
            if i != j:
                k = shift_amount(xs[i], xs[j])
                if k is not None and not same(shift(ys[i], k), ys[j]):
                    ti_dis = True
    lti_dis = any(not same(conv(x, h), y) for x, y in zip(xs, ys))
    if not lti_dis:
        return "lti"
    if l_dis and ti_dis:
        return "both"
    if l_dis:
        return "notL"
    if ti_dis:
        return "notTI"
    return "joint"


# ----------------------------------------------------------------------------- drawing
def rseq(lengths=(2, 3), lo=-2, hi=2, starts=(-1, 0, 1, 2)):
    L = random.choice(lengths)
    vals = [random.randint(lo, hi) for _ in range(L)]
    nz = [v for v in range(lo, hi + 1) if v != 0]
    vals[0], vals[-1] = random.choice(nz), random.choice(nz)
    return (vals, random.choice(starts))


def perturb(s):
    vals, st = list(s[0]), s[1]
    i = random.randrange(len(vals))
    for _ in range(10):
        d = random.choice([-2, -1, 1, 2])
        if vals[i] + d != 0:
            vals[i] += d
            break
    return trim(vals, st)


def kshift(s):
    return random.choice([-2, -1, 1, 2, 3])


def draw(family):
    m = random.choice([-1, 1, 2])
    h = rseq(lengths=(2, 3), lo=-2, hi=3, starts=(-1, 0, 1))
    x1 = ([1], m)
    y1 = shift(h, m)
    x2 = rseq()
    if family in ("lti", "joint"):      # same relation mix as the other families, so the relation type alone
        group = random.choice(["shift", "lin", "other"])          # does not reveal the answer
        rel = {"shift": "shift", "lin": random.choice(["scale", "sum"]),
               "other": random.choice(["none", "scaledshift", "sumshift"])}[group]
    else:
        rel = random.choice({"notL": ["scale", "scale", "sum"], "notTI": ["shift"]}[family])
    k = random.choice([-2, -1, 1, 2, 3])
    c = random.choice([-1, 2, -2, 3])
    d = random.choice([1, -1, 2])
    if rel == "shift":
        x3 = shift(x2, k)
    elif rel == "scale":
        x3 = scale(x2, c)
    elif rel == "sum":
        x3 = add(x2, scale(x1, d))
    elif rel == "scaledshift":
        x3 = scale(shift(x2, k), c)
    elif rel == "sumshift":
        x3 = add(x2, scale(shift(x1, k), d))
    else:
        x3 = rseq()
    y2, y3 = conv(x2, h), conv(x3, h)          # LTI outputs
    how = ""
    if family == "notL":
        if random.random() < 0.25:
            y2 = perturb(y2)
        if rel == "scale":
            how = random.choice(["square", "perturb", "same"])
            y3 = {"square": scale(y2, c * c), "perturb": perturb(scale(y2, c)), "same": y2}[how]
        else:
            how = random.choice(["perturb", "noadd"])
            y3 = perturb(add(y2, scale(y1, d))) if how == "perturb" else y2
    elif family == "notTI":
        if random.random() < 0.25:
            y2 = perturb(y2)
        how = random.choice(["nomove", "wrongmove", "perturb", "flipmove"])
        if how == "nomove":
            y3 = y2
        elif how == "wrongmove":
            y3 = shift(y2, random.choice([j for j in (-2, -1, 1, 2, 3) if j != k]))
        elif how == "flipmove":
            y3 = shift(y2, -k)
        else:
            y3 = perturb(shift(y2, k))
    elif family == "joint" and rel in ("shift", "scale", "sum"):
        how = "y2off"                                 # the relation holds between the outputs, but y2 != x2 * h
        y2 = perturb(y2)
        y3 = {"shift": shift(y2, k), "scale": scale(y2, c), "sum": add(y2, scale(y1, d))}[rel]
    elif family == "joint":
        if rel == "none":
            how = random.choice(["perturb3", "perturb2", "corr"])
            if how == "perturb3":
                y3 = perturb(y3)
            elif how == "perturb2":
                y2 = perturb(y2)
            else:                                     # x3 correlated with h instead of convolved (h not flipped)
                y3 = conv(x3, (h[0][::-1], h[1]))
        elif rel == "scaledshift":
            how = random.choice(["noshift", "noscale", "perturb"])
            y3 = {"noshift": scale(y2, c), "noscale": shift(y2, k), "perturb": perturb(y3)}[how]
        else:
            how = "noshift"
            y3 = add(y2, scale(y1, d))                # added y1 without shifting it
    return {"m": m, "h": h, "x": [x1, x2, x3], "y": [y1, y2, y3], "rel": rel, "k": k, "c": c, "d": d, "how": how}


def acceptable(D):
    xs, ys = D["x"], D["y"]
    for s in xs + ys:
        if not s[0] or len(s[0]) > 6 or max(abs(v) for v in s[0]) > 12:
            return False
        if not (-3 <= s[1] and s[1] + len(s[0]) - 1 <= 6):
            return False
    if same(xs[1], xs[2]) or len(xs[2][0]) < 2:
        return False
    if D["rel"] in ("none", "scaledshift", "sumshift", "shift") and rank(vectors(xs)) != 3:
        return False                                  # the solution says "linearly independent"
    if D["rel"] == "none":                            # x3 must not be a scaled shift of x2
        for kk in range(-8, 9):
            for cc in (-3, -2, -1, 1, 2, 3, Fraction(1, 2), Fraction(-1, 2)):
                if same(scale(shift(xs[1], kk), cc), xs[2]):
                    return False
    return True


def tex_s(s):
    return fmt.tex_seq(s[0], s[1])


def generate(data):
    family = random.choice(["lti", "notL", "notTI", "joint"])
    while True:
        D = draw(family)
        if acceptable(D) and classify(D["x"], D["y"], D["h"]) == family:
            break
    while True:                                           # part (b): a new input
        v = rseq(lengths=(2, 3), lo=-2, hi=2, starts=(-1, 0, 1))
        if all(not any(same(scale(shift(xx, kk), cc), v) for kk in range(-9, 10) for cc in (-2, -1, 1, 2))
               for xx in D["x"][1:]):
            break
    h, m = D["h"], D["m"]
    w = conv(v, h)

    p = data["params"]
    xs, ys = D["x"], D["y"]
    p["family"], p["rel"], p["how"] = family, D["rel"], D["how"]
    p["m"] = m
    p["pairs"] = [{"x": list(x[0]), "nx": x[1], "y": list(y[0]), "ny": y[1]} for x, y in zip(xs, ys)]
    p["v"], p["nv"] = v[0], v[1]
    rows = []
    for i, (x, y) in enumerate(zip(xs, ys), start=1):
        rows.append(f"<tr><td style=\"padding-right:0.8em\">Pair {i}:</td><td>$x_{i}[n] = {tex_s(x)}$</td>"
                    f"<td style=\"padding:0 0.8em\">$\\longmapsto$</td><td>$y_{i}[n] = {tex_s(y)}$</td></tr>")
    p["pairs_html"] = '<table class="mx-auto my-2">' + "".join(rows) + "</table>"
    p["x1_tex"] = rf"\delta[n {'-' if m > 0 else '+'} {abs(m)}]"
    p["v_tex"] = tex_s(v)
    p["verdict_choices"] = [{"text": t, "correct": key == family} for key, t in OPTIONS]

    # ---------------- worked solution
    k, c, d, rel = D["k"], D["c"], D["d"], D["rel"]
    cf = lambda v: "-" if v == -1 else (fmt.tex_num(v) + r"\,")          # noqa: E731  coefficient: -x, 2\,x
    a = lambda j: f"n {'-' if j > 0 else '+'} {abs(j)}" if j else "n"     # noqa: E731
    p["h_tex"] = tex_s(h)
    p["m_abs"], p["m_back"] = abs(m), a(-m)
    lines = []
    if rel == "shift":
        lines.append(rf"$x_3[n] = x_2[{a(k)}]$: a pure shift. Time-invariance forces $y_3[n] = y_2[{a(k)}] = {tex_s(shift(ys[1], k))}$; "
                     + ("the given $y_3$ is exactly that." if same(shift(ys[1], k), ys[2]) else
                        r"the given $y_3$ is <b>not</b>, so the system is <b>not time-invariant</b>."))
        lines.append(r"No input is a linear combination of the others (an impulse, a sequence and its shift are linearly "
                     r"independent), so linearity puts no constraint on the outputs: a linear time-varying system can produce "
                     r"any outputs for them.")
    elif rel == "scale":
        lines.append(rf"$x_3[n] = {cf(c)}x_2[n]$. Linearity (homogeneity) forces $y_3[n] = {cf(c)}y_2[n] = "
                     rf"{tex_s(scale(ys[1], c))}$; "
                     + ("the given $y_3$ is exactly that." if same(scale(ys[1], c), ys[2]) else
                        r"the given $y_3$ is <b>not</b>, so the system is <b>not linear</b>."))
        lines.append(r"No input is a shifted copy of another, so time-invariance puts no constraint on the outputs: "
                     r"a time-invariant (nonlinear) system can map each of these inputs, and its shifts, to anything.")
    elif rel == "sum":
        lines.append(rf"$x_3[n] = x_2[n] {'+' if d > 0 else '-'} {'' if abs(d) == 1 else abs(d)}\,x_1[n]$. Linearity forces "
                     rf"$y_3[n] = y_2[n] {'+' if d > 0 else '-'} {'' if abs(d) == 1 else abs(d)}\,y_1[n] = {tex_s(add(ys[1], scale(ys[0], d)))}$; "
                     + ("the given $y_3$ is exactly that." if same(add(ys[1], scale(ys[0], d)), ys[2]) else
                        r"the given $y_3$ is <b>not</b>, so the system is <b>not linear</b>."))
        lines.append(r"No input is a shifted copy of another (an impulse cannot be shifted into a multi-sample sequence, "
                     r"and $x_3$ and $x_2$ have different sums), so time-invariance puts no constraint on the outputs.")
    elif rel == "scaledshift":
        lines.append(rf"$x_3[n] = {cf(c)}x_2[{a(k)}]$: neither a pure shift (that would need the factor $1$) nor a pure "
                     r"scaling (that would need no shift) of an earlier input, and the three inputs are linearly independent. "
                     r"So linearity alone and time-invariance alone each allow any outputs.")
        lines.append(rf"Together they do not: an LTI system must give $y_3[n] = {cf(c)}y_2[{a(k)}] = "
                     rf"{tex_s(scale(shift(ys[1], k), c))}$; "
                     + ("the given $y_3$ is exactly that." if same(scale(shift(ys[1], k), c), ys[2]) else
                        r"the given $y_3$ is <b>not</b>."))
    elif rel == "sumshift":
        lines.append(rf"$x_3[n] = x_2[n] {'+' if d > 0 else '-'} {'' if abs(d) == 1 else abs(d)}\,x_1[{a(k)}]$: it uses a "
                     r"<em>shifted</em> impulse, which is not among the inputs, and the three inputs are linearly independent and "
                     r"not shifts of one another. So linearity alone and time-invariance alone each allow any outputs.")
        lines.append(rf"Together they do not: an LTI system must give $y_3[n] = y_2[n] {'+' if d > 0 else '-'} "
                     rf"{'' if abs(d) == 1 else abs(d)}\,y_1[{a(k)}] = {tex_s(add(ys[1], scale(shift(ys[0], k), d)))}$; "
                     + ("the given $y_3$ is exactly that." if same(add(ys[1], scale(shift(ys[0], k), d)), ys[2]) else
                        r"the given $y_3$ is <b>not</b> (it adds $y_1$ without shifting it)."))
    else:
        lines.append(r"No input is a shift of another and the three inputs are linearly independent, so neither linearity "
                     r"nor time-invariance alone constrains the outputs. Only both together do: an LTI system has "
                     r"$y = x * h$ for one $h$, and pair 1 fixes $h$.")
    checks = []
    for i in (1, 2):
        pred = conv(xs[i], h)
        ok = same(pred, ys[i])
        checks.append(rf"$x_{i + 1} * h = {tex_s(pred)}$ " + (rf"$= y_{i + 1}$ &#10003;" if ok else rf"$\ne y_{i + 1}$"))
    p["why_lines"] = [{"line_html": s} for s in lines]
    p["lti_checks"] = [{"line_html": s} for s in checks]
    concl = {
        "lti": r"Every pair agrees with $y = x * h$, so the data are <b>consistent with an LTI system</b> (they cannot prove "
               r"that the system is LTI: other inputs might misbehave).",
        "notL": r"Conclusion: <b>not linear</b> (hence not LTI); time-invariance is not decided by these pairs.",
        "notTI": r"Conclusion: <b>not time-invariant</b> (hence not LTI); linearity is not decided by these pairs.",
        "joint": r"Conclusion: <b>no LTI system</b> produces all three pairs, yet no relation among the inputs exposes "
                 r"nonlinearity alone or time-variance alone: the system is nonlinear, or time-varying, or both &mdash; "
                 r"these data cannot say which.",
    }[family]
    p["concl_html"] = concl
    p["w_tex"] = tex_s(w)
    p["vstart"], p["hstart"], p["wstart"] = v[1], h[1], w[1]
    p["y1start"] = ys[0][1]

    cr = data["correct_answers"]
    cr["wfirst"] = w[1]
    cr["wrow"] = pl.to_json(np.array([[float(t) for t in w[0]]]))
