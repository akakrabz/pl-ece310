"""X_d(w0) of a short real sequence at two or three standard frequencies, then |X_d| and the principal
phase at one of them (Quiz 2, Oct 2026: "X_d(w0) at a given w0 for a finite sequence"; HW5 #1;
Lecture 13 DTFT definition).

    X_d(w0) = sum_n x[n] e^{-j w0 n},  w0 in {0, pi, pi/2, -pi/2, pi/3, 2pi/3}
    e^{-j pi n} = (-1)^n,  e^{-j pi n/2} = (-j)^n,  e^{-j pi n/3}: period 6, values (1 +- j sqrt3)/2 ...

Exact arithmetic: Re X_d is rational and Im X_d = q sqrt(d) with d = 3 at multiples of pi/3, d = 1
otherwise. Re/Im are numeric inputs (rtol 1e-3); the phase is entered in radians and graded in
grade() on the difference wrapped to (-pi, pi] (so -pi and pi are both accepted when the answer is pi),
with the submission required to lie in [-pi, pi] (an angle off by 2 pi is wrong).

Design: length 3-5, integers in [-3, 3], nonzero end samples, start in [-3, 2] and not 0 in most
variants (the true index matters); the sequence is not even about n = 0 (the sign error e^{+jwn}
changes the answer) and at least one asked frequency has Im X_d != 0. The magnitude/phase frequency
has |X_d| >= 1 and X_d not a positive real number; in about 60 % of the variants Re X_d < 0 (the
quadrant trap of arctan)."""

import math
import random
from fractions import Fraction as Fr

import prairielearn as pl
from ece310 import fmt
from ece310.w1_signals import Rad, cos_pi, sin_pi, tex_pi

FREQS = [Fr(0), Fr(1), Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(2, 3)]
NZ = [-3, -2, -1, 1, 2, 3]


def _unit(t, n):
    """e^{-j t pi n} = (re, im) with re a Fraction and im a Rad."""
    c, s = cos_pi(t * n), sin_pi(t * n)
    return c.q, -s


def _dval(t):
    return 3 if t.denominator == 3 else 1


def _value(x, start, t):
    """X_d(t pi) exactly: (Re Fraction, Im coefficient Fraction, d)."""
    d = _dval(t)
    re, im = Fr(0), Fr(0)
    for i, v in enumerate(x):
        cr, si = _unit(t, start + i)
        assert si.is_zero() or si.d == d
        re += v * cr
        im += v * si.q
    return re, im, d


def _cplx_tex(re, imq, d):
    """a + j b in LaTeX with b = imq sqrt(d)."""
    im = Rad(imq, d)
    if im.is_zero():
        return fmt.tex_num(re)
    mag = Rad(abs(imq), d)
    jt = "j" if mag == Rad(1) else r"j\," + mag.tex()
    if re == 0:
        return ("-" if imq < 0 else "") + jt
    return f"{fmt.tex_num(re)} {'-' if imq < 0 else '+'} {jt}"


def _dec(v):
    v = float(v)
    if abs(v - round(v)) < 1e-12:
        return str(int(round(v)))
    return f"{v:.4f}".rstrip("0")


def _mag_rad(re, imq, d):
    """|X| as a Rad: sqrt(m) for the rational m = re^2 + 3 imq^2 (or + imq^2)."""
    m = re * re + imq * imq * d
    return Rad(Fr(1, m.denominator), m.numerator * m.denominator), m


def _phase_tex(re, imq, d, ph):
    """How the principal phase follows from Re and Im (one line)."""
    if imq == 0:
        return r"X_d \lt 0 \text{ is a negative real number}:\ \angle X_d = \pi \approx 3.1416"
    if re == 0:
        return r"\angle X_d = " + (r"\tfrac{\pi}{2}" if imq > 0 else r"-\tfrac{\pi}{2}") + rf" \approx {_dec(ph)}"
    ratio = Rad(imq / re, d)
    base = rf"\arctan\!\left({ratio.tex()}\right)"
    q = Fr(ph / math.pi).limit_denominator(12)
    val = (rf"= {tex_pi(q)} \approx {_dec(ph)}" if abs(float(q) * math.pi - ph) < 1e-12 else rf"\approx {_dec(ph)}")
    if re > 0:
        return rf"\mathrm{{Re}} \gt 0:\ \angle X_d = {base} {val}"
    quad = r"\mathrm{Re} \lt 0,\ \mathrm{Im} \gt 0" if imq > 0 else r"\mathrm{Re} \lt 0,\ \mathrm{Im} \lt 0"
    pm = "+" if imq > 0 else "-"
    return rf"{quad}:\ \angle X_d = {base} {pm} \pi {val}"


def _draw():
    L = random.randint(3, 5)
    start = random.choice([-3, -2, -2, -1, -1, 0, 1, 1, 2])
    x = [random.randint(-3, 3) for _ in range(L)]
    x[0], x[-1] = random.choice(NZ), random.choice(NZ)
    return x, start


def generate(data):
    while True:
        x, start = _draw()
        d = {start + i: v for i, v in enumerate(x)}
        if all(d.get(n, 0) == d.get(-n, 0) for n in d):        # even about 0: e^{+jwn} gives the same
            continue
        m = random.choice([2, 3])
        ts = random.sample(FREQS, m)
        if not any(t in (Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(2, 3)) for t in ts):
            continue
        ts = [t for t in FREQS if t in ts]                     # canonical order
        vals = [_value(x, start, t) for t in ts]
        if sum(1 for re, im, _ in vals if re == 0 and im == 0) > 1:
            continue
        if not any(im != 0 for _, im, _ in vals):
            continue
        cand = [k for k, (re, im, dd) in enumerate(vals)
                if re * re + im * im * dd >= 1 and not (im == 0 and re > 0)]
        if not cand:
            continue
        cplx = [k for k in cand if vals[k][1] != 0]          # phase pi (negative real) only now and then
        if cplx and random.random() < 0.88:
            cand = cplx
        neg = [k for k in cand if vals[k][0] < 0]
        k = random.choice(neg) if neg and random.random() < 0.6 else random.choice(cand)
        break

    p = data["params"]
    c = data["correct_answers"]
    p["x"], p["start"] = x, start
    p["x_tex"] = fmt.tex_seq(x, start)
    p["ts"] = [[t.numerator, t.denominator] for t in ts]
    p["has3"] = m == 3
    ns = list(range(start, start + len(x)))
    rows = []
    for i, (t, (re, imq, dd)) in enumerate(zip(ts, vals), start=1):
        p[f"w{i}_tex"] = tex_pi(t)
        p[f"X{i}_tex"] = _cplx_tex(re, imq, dd)
        approx = ""
        if dd == 3 and imq != 0:
            jpart = rf"j\,{_dec(abs(float(Rad(imq, dd))))}"
            approx = (rf" \approx {_dec(re)} {'-' if imq < 0 else '+'} {jpart}" if re != 0
                      else rf" \approx {'-' if imq < 0 else ''}{jpart}")
        p[f"X{i}_approx_tex"] = approx
        cells = "".join(f"<td>${_cplx_tex(*_unit_pair(t, n))}$</td>" for n in ns)
        rows.append(f"<tr><th>${_ewn(t)}$</th>{cells}</tr>")
        c[f"re{i}"] = float(re)
        c[f"im{i}"] = float(Rad(imq, dd))
    head = "<tr><th>$n$</th>" + "".join(f"<th>${n}$</th>" for n in ns) + "</tr>"
    xrow = "<tr><th>$x[n]$</th>" + "".join(f"<td>${fmt.tex_num(v)}$</td>" for v in x) + "</tr>"
    p["table_html"] = ('<table class="table table-sm table-bordered text-center" style="width:auto">'
                       + head + xrow + "".join(rows) + "</table>")

    re, imq, dd = vals[k]
    magr, msq = _mag_rad(re, imq, dd)
    ph = math.atan2(float(Rad(imq, dd)), float(re))
    p["k"] = k
    p["wm_tex"] = tex_pi(ts[k])
    p["Xm_tex"] = _cplx_tex(re, imq, dd)
    p["mp_letter"] = "d" if m == 3 else "c"
    p["msq_tex"] = fmt.tex_num(msq)
    p["mag_tex"] = magr.tex()
    p["mag_dec"] = _dec(float(magr))
    p["mag_exact_is_dec"] = magr.is_rational()
    chain = [rf"\sqrt{{{fmt.tex_num(msq)}}}"]
    if magr.tex() != chain[0]:
        chain.append(magr.tex())
    if not magr.is_rational():
        chain.append(r"\approx " + _dec(float(magr)))
    p["mag_chain_tex"] = " = ".join(chain).replace(r"= \approx", r"\approx")
    p["phase_line_tex"] = _phase_tex(re, imq, dd, ph)
    p["ph_is_pi"] = imq == 0 and re < 0
    c["mag"] = float(magr)
    c["ph"] = ph


def _ewn(t):
    """e^{-j t pi n} in LaTeX with the sign folded in: e^{-j\\pi n}, e^{j\\tfrac{\\pi}{2}n}, e^{0}."""
    if t == 0:
        return "e^{-j0n}"
    s = tex_pi(abs(t))
    return "e^{" + ("-" if t > 0 else "") + "j" + s + (" n" if s.endswith(r"\pi") else "n") + "}"


def _unit_pair(t, n):
    cr, si = _unit(t, n)
    return cr, si.q, (si.d if not si.is_zero() else 1)


# ----------------------------------------------------------------------------- phase grading
def _phase_ok(sub, ans):
    tol = 1e-3 * abs(ans) + 1e-3                     # the element's rtol/atol (its Test button probes this)
    d = (sub - ans + math.pi) % (2 * math.pi) - math.pi
    return abs(d) <= tol and abs(sub) <= math.pi + tol


def grade(data):
    sub = data["submitted_answers"].get("ph")
    if sub is None or "ph" in data["format_errors"]:
        return
    try:
        ok = _phase_ok(float(sub), float(data["correct_answers"]["ph"]))
    except (TypeError, ValueError, OverflowError):
        ok = False
    ps = data["partial_scores"].setdefault("ph", {"weight": 1})
    ps["score"] = 1.0 if ok else 0.0
    pl.set_weighted_score_data(data)
