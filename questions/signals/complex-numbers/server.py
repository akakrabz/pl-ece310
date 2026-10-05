"""Complex-number drills with exact answers (Lecture 2, HW1 #4). One of three templates:

  polar : |w| and the principal angle of a product / quotient / power of nice numbers
          (1 +- j, +-sqrt3 +- j, 1 +- j sqrt3, j, 2e^{j pi/3}, ...);
  roots : the N roots of z^N = c (N in {2, 3, 4}): common magnitude and the principal angle of one
          specified root (HW1 #4: z^4 - 1 = 0, z^4 + 1 = 0);
  sum   : real and imaginary parts of A1 e^{j th1} + A2 e^{j th2}.

Angles are entered as multiples of pi, principal value in (-pi, pi]. Magnitudes such as 2 sqrt2 are
irrational, so pl-number-input uses a relative tolerance and the prompt says a decimal is fine.
Everything is computed exactly (Fractions; magnitudes as q * sqrt(2)^e, sums in Q(sqrt d))."""

import math
import random
from fractions import Fraction as Fr

from ece310 import fmt

# ----------------------------------------------------------------------------- exact helpers


def principal(t):
    """Reduce an angle (multiple of pi, Fraction) to (-1, 1]."""
    t = Fr(t)
    while t <= -1:
        t += 2
    while t > 1:
        t -= 2
    return t


def tex_pi(t):
    """Angle t*pi in LaTeX: \\tfrac{3\\pi}{4}, -\\tfrac{\\pi}{6}, \\pi, 0."""
    t = Fr(t)
    if t == 0:
        return "0"
    s = "-" if t < 0 else ""
    p, q = abs(t.numerator), t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return s + (num if q == 1 else r"\tfrac{" + num + "}{" + str(q) + "}")


def tex_exp(t):
    """e^{j t pi} exponent text: e^{j\\pi/3}, e^{-j2\\pi/3}."""
    t = Fr(t)
    s = "-" if t < 0 else ""
    p, q = abs(t.numerator), t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return f"e^{{{s}j{num}" + ("" if q == 1 else f"/{q}") + "}"


class Mag:
    """q * sqrt(2)^e with q a positive Fraction and e in {0, 1}."""

    def __init__(self, q, e=0):
        q, e = Fr(q), int(e)
        q *= 2 ** (e // 2)
        self.q, self.e = q, e % 2

    def __mul__(self, o):
        return Mag(self.q * o.q, self.e + o.e)

    def __truediv__(self, o):
        # 1/sqrt2 = sqrt2/2
        return Mag(self.q / o.q / (2 if o.e else 1), self.e + o.e)

    def __pow__(self, k):
        out = Mag(1)
        for _ in range(k):
            out = out * self
        return out

    def value(self):
        return float(self.q) * (math.sqrt(2) if self.e else 1.0)

    def tex(self):
        if not self.e:
            return fmt.tex_num(self.q)
        p, q = self.q.numerator, self.q.denominator
        top = r"\sqrt{2}" if p == 1 else f"{p}\\sqrt{{2}}"
        return top if q == 1 else r"\tfrac{" + top + "}{" + str(q) + "}"


class QS:
    """a + b sqrt(d), a, b Fractions, d in {2, 3}."""

    def __init__(self, a, b, d):
        self.a, self.b, self.d = Fr(a), Fr(b), d

    def __add__(self, o):
        return QS(self.a + o.a, self.b + o.b, self.d)

    def scale(self, c):
        return QS(self.a * c, self.b * c, self.d)

    def value(self):
        return float(self.a) + float(self.b) * math.sqrt(self.d)

    def is_rational(self):
        return self.b == 0

    def tex(self):
        parts = []
        if self.a != 0:
            parts.append(fmt.tex_num(self.a))
        if self.b != 0:
            p, q = abs(self.b.numerator), self.b.denominator
            top = rf"\sqrt{{{self.d}}}" if p == 1 else f"{p}\\sqrt{{{self.d}}}"
            t = top if q == 1 else r"\tfrac{" + top + "}{" + str(q) + "}"
            if parts:
                parts.append(("- " if self.b < 0 else "+ ") + t)
            else:
                parts.append(("-" if self.b < 0 else "") + t)
        return " ".join(parts) if parts else "0"


# cos / sin of t*pi for t a multiple of 1/6 or 1/4, exactly
def cos_sin(t, d):
    t = principal(t)
    table6 = {Fr(0): (QS(1, 0, 3), QS(0, 0, 3)), Fr(1, 6): (QS(0, Fr(1, 2), 3), QS(Fr(1, 2), 0, 3)),
              Fr(1, 3): (QS(Fr(1, 2), 0, 3), QS(0, Fr(1, 2), 3)), Fr(1, 2): (QS(0, 0, 3), QS(1, 0, 3)),
              Fr(2, 3): (QS(Fr(-1, 2), 0, 3), QS(0, Fr(1, 2), 3)), Fr(5, 6): (QS(0, Fr(-1, 2), 3), QS(Fr(1, 2), 0, 3)),
              Fr(1): (QS(-1, 0, 3), QS(0, 0, 3))}
    table4 = {Fr(0): (QS(1, 0, 2), QS(0, 0, 2)), Fr(1, 4): (QS(0, Fr(1, 2), 2), QS(0, Fr(1, 2), 2)),
              Fr(1, 2): (QS(0, 0, 2), QS(1, 0, 2)), Fr(3, 4): (QS(0, Fr(-1, 2), 2), QS(0, Fr(1, 2), 2)),
              Fr(1): (QS(-1, 0, 2), QS(0, 0, 2))}
    tab = table6 if d == 3 else table4
    c, s = tab[abs(t)]
    return c, (s if t >= 0 else s.scale(-1))


# ----------------------------------------------------------------------------- template: polar
# (rectangular LaTeX, magnitude, angle / pi)
BASE = [
    ("1 + j", Mag(1, 1), Fr(1, 4)), ("1 - j", Mag(1, 1), Fr(-1, 4)),
    ("-1 + j", Mag(1, 1), Fr(3, 4)), ("-1 - j", Mag(1, 1), Fr(-3, 4)),
    (r"\sqrt{3} + j", Mag(2), Fr(1, 6)), (r"\sqrt{3} - j", Mag(2), Fr(-1, 6)),
    (r"-\sqrt{3} + j", Mag(2), Fr(5, 6)), (r"-\sqrt{3} - j", Mag(2), Fr(-5, 6)),
    (r"1 + j\sqrt{3}", Mag(2), Fr(1, 3)), (r"1 - j\sqrt{3}", Mag(2), Fr(-1, 3)),
    (r"-1 + j\sqrt{3}", Mag(2), Fr(2, 3)), (r"-1 - j\sqrt{3}", Mag(2), Fr(-2, 3)),
    ("2j", Mag(2), Fr(1, 2)), ("-j", Mag(1), Fr(-1, 2)), ("-2", Mag(2), Fr(1)),
    (r"2e^{j\pi/3}", Mag(2), Fr(1, 3)), (r"3e^{-j\pi/4}", Mag(3), Fr(-1, 4)),
    (r"e^{j3\pi/4}", Mag(1), Fr(3, 4)), (r"2e^{-j5\pi/6}", Mag(2), Fr(-5, 6)),
]


def _gen_polar(p):
    while True:
        form = random.choice(["prod", "quot", "pow", "powprod", "powquot", "prodquot"])
        z = random.sample(BASE, 3)
        k = random.choice([2, 3]) if form in ("powprod", "powquot") else random.choice([3, 4, 5])

        def par(b, power=False):
            return b[0] if ("e^" in b[0] and not power) else f"\\left({b[0]}\\right)"

        def mul(a_tex, b_tex):
            """a times b, with a cdot when b starts with a digit (never "2e^{j pi/3}3e^{-j pi/4}")."""
            return a_tex + (r" \cdot " if b_tex[0].isdigit() else "") + b_tex

        if form == "prod":
            w_tex = mul(par(z[0]), par(z[1]))
            used, mag, ang = [z[0], z[1]], z[0][1] * z[1][1], z[0][2] + z[1][2]
            mag_steps = rf"{z[0][1].tex()}\cdot {z[1][1].tex()}"
            ang_steps = rf"{tex_pi(z[0][2])} + \left({tex_pi(z[1][2])}\right)"
        elif form == "quot":
            w_tex = fmt.tex_frac(z[0][0], z[1][0])
            used, mag, ang = [z[0], z[1]], z[0][1] / z[1][1], z[0][2] - z[1][2]
            mag_steps = fmt.tex_frac(z[0][1].tex(), z[1][1].tex())
            ang_steps = rf"{tex_pi(z[0][2])} - \left({tex_pi(z[1][2])}\right)"
        elif form == "pow":
            w_tex = f"{par(z[0], True)}^{{{k}}}"
            used, mag, ang = [z[0]], z[0][1] ** k, k * z[0][2]
            mag_steps = rf"\left({z[0][1].tex()}\right)^{{{k}}}"
            ang_steps = rf"{k}\cdot\left({tex_pi(z[0][2])}\right)"
        elif form == "powprod":
            w_tex = mul(f"{par(z[0], True)}^{{{k}}}", par(z[1]))
            used, mag, ang = [z[0], z[1]], (z[0][1] ** k) * z[1][1], k * z[0][2] + z[1][2]
            mag_steps = rf"\left({z[0][1].tex()}\right)^{{{k}}}\cdot {z[1][1].tex()}"
            ang_steps = rf"{k}\cdot\left({tex_pi(z[0][2])}\right) + \left({tex_pi(z[1][2])}\right)"
        elif form == "powquot":
            w_tex = fmt.tex_frac(f"{par(z[0], True)}^{{{k}}}", z[1][0])
            used, mag, ang = [z[0], z[1]], (z[0][1] ** k) / z[1][1], k * z[0][2] - z[1][2]
            mag_steps = fmt.tex_frac(rf"\left({z[0][1].tex()}\right)^{{{k}}}", z[1][1].tex())
            ang_steps = rf"{k}\cdot\left({tex_pi(z[0][2])}\right) - \left({tex_pi(z[1][2])}\right)"
        else:
            w_tex = fmt.tex_frac(mul(par(z[0]), par(z[1])), z[2][0])
            used, mag, ang = z, z[0][1] * z[1][1] / z[2][1], z[0][2] + z[1][2] - z[2][2]
            mag_steps = fmt.tex_frac(rf"{z[0][1].tex()}\cdot {z[1][1].tex()}", z[2][1].tex())
            ang_steps = rf"{tex_pi(z[0][2])} + \left({tex_pi(z[1][2])}\right) - \left({tex_pi(z[2][2])}\right)"
        pa = principal(ang)
        mv = mag.value()
        if pa == 0 or abs(mv - 1) < 1e-12 or mv > 200 or mv < 1 / 50:
            continue
        break
    p["w_tex"] = w_tex
    p["factors_tex"] = r",\qquad ".join(b[0] if "e^" in b[0] else rf"{b[0]} = {b[1].tex()}\,{tex_exp(b[2])}"
                                       for b in used)
    p["mag_steps_tex"] = mag_steps
    p["mag_tex"] = mag.tex()
    p["ang_steps_tex"] = ang_steps
    p["raw_ang_tex"] = tex_pi(ang)
    p["ang_tex"] = tex_pi(pa)
    p["needs_reduction"] = pa != ang
    p["mag_dec"] = f"{mv:.6g}"
    p["mag_irr"] = mag.e == 1
    p["ang_frac"] = fmt.plain_num(pa)
    p["kind"] = "polar"
    p["exact_mag"] = [mag.q.numerator, mag.q.denominator, mag.e]
    p["exact_ang"] = [pa.numerator, pa.denominator]
    p["raw_ang"] = [ang.numerator, ang.denominator]
    return mv, float(pa)


# ----------------------------------------------------------------------------- template: roots
THETAS = [Fr(0), Fr(1), Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(-1, 3), Fr(2, 3), Fr(-2, 3), Fr(1, 4), Fr(-1, 4),
          Fr(3, 4), Fr(-3, 4), Fr(1, 6), Fr(-1, 6), Fr(5, 6), Fr(-5, 6)]

SELECTORS = {
    "minpos": "the root with the smallest <b>positive</b> principal angle",
    "max": "the root with the largest principal angle",
    "negnear": "the root whose principal angle is <b>negative</b> and closest to $0$",
}


def _c_tex(R, th):
    """c = R e^{j th pi}: rectangular when it is an integer combination, else polar."""
    if th == 0:
        return str(R)
    if th == 1:
        return f"-{R}"
    if th in (Fr(1, 2), Fr(-1, 2)):
        return ("" if th > 0 else "-") + (f"{R}j" if R != 1 else "j")
    if R % 2 == 0 and random.random() < 0.5 and th.denominator in (3, 6):
        c, s = cos_sin(th, 3)
        re, im = c.scale(R), s.scale(R)
        im_t = im.tex()
        sign, mag_t = ("-", im_t[1:]) if im_t.startswith("-") else ("+", im_t)
        return f"{re.tex()} {sign} j" + ("" if mag_t == "1" else mag_t)
    return (f"{R}" if R != 1 else "") + tex_exp(th)


def _gen_roots(p):
    while True:
        N = random.choice([2, 3, 4])
        if random.random() < 0.2:
            R, rho_tex, rho = 2, (r"\sqrt{2}" if N == 2 else rf"\sqrt[{N}]{{2}}"), 2 ** (1 / N)
        else:
            rr = random.choice([1, 1, 2, 2, 3])
            if rr == 3 and N == 4:
                rr = 2
            R, rho_tex, rho = rr ** N, str(rr), float(rr)
        th = random.choice(THETAS)
        angs = sorted(principal((th + 2 * k) / N) for k in range(N))
        sel = random.choice(list(SELECTORS))
        if sel == "minpos":
            cand = [a for a in angs if a > 0]
            target = min(cand)
        elif sel == "max":
            target = max(angs)
        else:
            cand = [a for a in angs if a < 0]
            if not cand:
                continue
            target = max(cand)
        if target == 0:
            continue
        break
    p["N"] = N
    p["Nm1"] = N - 1
    p["k_list_tex"] = ", ".join(str(k) for k in range(N))
    p["rho_irr"] = R == 2
    p["R_pre"] = "" if R == 1 else f"{R}\\,"
    p["rho_pre"] = "" if rho_tex == "1" else rho_tex + "\\;"
    p["c_tex"] = _c_tex(R, th)
    p["R"] = R
    p["c_polar_tex"] = (f"{R}" if R != 1 else "") + tex_exp(th) if th != 0 else str(R)
    p["c_chain_tex"] = (p["c_tex"] if p["c_tex"] == p["c_polar_tex"]
                        else p["c_tex"] + " = " + p["c_polar_tex"])
    p["th_tex"] = tex_pi(th)
    p["rho_tex"] = rho_tex
    p["sel_html"] = SELECTORS[sel]
    p["sel"] = sel
    p["roots_tex"] = r",\quad ".join(tex_pi(a) for a in angs)
    p["gen_angle_tex"] = rf"\left({tex_pi(th)} + 2\pi k\right)/{N}"
    p["ang_tex"] = tex_pi(target)
    p["mag_dec"] = f"{rho:.6g}"
    p["ang_frac"] = fmt.plain_num(target)
    p["kind"] = "roots"
    p["theta"] = [th.numerator, th.denominator]
    p["root_angles"] = [[a.numerator, a.denominator] for a in angs]
    return rho, float(target)


def cplx_tex(re, im):
    """a + jb in LaTeX for a single exact term (im is rational or a multiple of sqrt d)."""
    rt, it = re.tex(), im.tex()
    if it == "0":
        return rt
    neg = it.startswith("-")
    mag = it[1:] if neg else it
    jt = "j" if mag == "1" else rf"j\,{mag}"
    if rt == "0":
        return ("-" if neg else "") + jt
    return f"{rt} {'-' if neg else '+'} {jt}"


# ----------------------------------------------------------------------------- template: sum
def _gen_sum(p):
    while True:
        d = random.choice([3, 3, 2])
        grid = [Fr(k, 6) for k in range(-5, 7)] if d == 3 else [Fr(k, 4) for k in range(-3, 5)]
        grid = [t for t in grid if t not in (0,)]
        A1, A2 = random.choice([1, 2, 3, 4]), random.choice([1, 2, 3, 4])
        t1, t2 = random.sample(grid, 2)
        c1, s1 = cos_sin(t1, d)
        c2, s2 = cos_sin(t2, d)
        re = c1.scale(A1) + c2.scale(A2)
        im = s1.scale(A1) + s2.scale(A2)
        rv, iv = re.value(), im.value()
        if abs(rv) < 1e-9 or abs(iv) < 1e-9 or abs(rv - iv) < 1e-9 or abs(rv + iv) < 1e-9:
            continue
        if t2 > 0 and t1 > 0:          # make at least one exponent negative most of the time
            if random.random() < 0.7:
                continue
        break

    def term(A, t):
        return (f"{A}" if A != 1 else "") + tex_exp(t)

    p["w_tex"] = term(A1, t1) + " + " + term(A2, t2)
    p["terms"] = [[A1, t1.numerator, t1.denominator], [A2, t2.numerator, t2.denominator]]

    def euler(A, t, c, s):
        re, im = c.scale(A), s.scale(A)
        cs = rf"\cos\left({tex_pi(t)}\right) + j\sin\left({tex_pi(t)}\right)"
        return rf"{term(A, t)} = {A if A != 1 else ''}\left({cs}\right) = {cplx_tex(re, im)}"

    p["euler1_tex"] = euler(A1, t1, c1, s1)
    p["euler2_tex"] = euler(A2, t2, c2, s2)
    p["re_tex"] = re.tex()
    p["im_tex"] = im.tex()
    p["re_dec"] = f"{rv:.6g}"
    p["im_dec"] = f"{iv:.6g}"
    p["re_irr"] = not re.is_rational()
    p["im_irr"] = not im.is_rational()
    p["kind"] = "sum"
    return rv, iv


TEMPLATE_KEYS = ["N", "Nm1", "k_list_tex", "mag_irr", "re_irr", "im_irr", "rho_irr", "R_pre", "rho_pre", "R", "ang_frac", "ang_steps_tex", "ang_tex", "c_chain_tex", "c_tex", "euler1_tex",
                 "euler2_tex", "factors_tex", "gen_angle_tex", "im_dec", "im_tex", "mag_dec", "mag_steps_tex",
                 "mag_tex", "needs_reduction", "raw_ang_tex", "re_dec", "re_tex", "rho_tex", "roots_tex",
                 "sel_html", "th_tex", "w_tex"]


def generate(data):
    kind = random.choice(["polar", "roots", "sum"])
    p = data["params"]
    for key in TEMPLATE_KEYS:          # every key the template mentions exists (inactive sections too)
        p[key] = ""
    p["is_polar"] = kind == "polar"
    p["is_roots"] = kind == "roots"
    p["is_sum"] = kind == "sum"
    if kind == "sum":
        rv, iv = _gen_sum(p)
        data["correct_answers"]["re"] = rv
        data["correct_answers"]["im"] = iv
    else:
        mv, av = _gen_polar(p) if kind == "polar" else _gen_roots(p)
        data["correct_answers"]["mag"] = mv
        data["correct_answers"]["ang"] = av
