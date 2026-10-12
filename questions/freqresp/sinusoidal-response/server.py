"""LTI response to a constant, a complex exponential and a real sinusoid, one frequency at a time
(Lecture 15; SP2025 Midterm 2 #4, FA2024 #4, SP2023 #6, FA2023 #4(b), FA2019 #7; HW6 #6-#7).

The system is given as a closed-form H_d(w) (exam style), a short real h[n], an FIR H(z) or a
first-order causal LCCDE:
  absw     |w| e^{+-j pi sin w}        real h (FA2024 #4)          H_d(0) = 0
  cosdelay (1 + cos w) e^{-jkw}        real h, k in {1, 2, 3}      zero at pi
  jw       j w e^{+-j pi sin w}        real h (SP2025 #4)          jumps at +-pi: never used at pi
  wcos     w e^{+-j pi cos w}          complex h (SP2023 #6)       jumps at +-pi: never used at pi
  halfcos  (b + cos w) e^{+-j|w|}      complex h (FA2023 #4), b in {1/2, 1}
  fir      short (anti)symmetric real h[n] (or H(z) = sum h[k] z^-k), H_d = e^{-j alpha w} (j) A(w)
  lccde    y[n] - a y[n-1] = g (x[n] + q x[n-1]), a = +-1/2: the denominator 1 - a e^{-jw} has a
           pi-rational angle (+-pi/6) only at w = pi/3 (a = 1/2) or 2pi/3 (a = -1/2), so the sinusoid
           sits there and the exponential is (-1)^n.
Input: c + a1 E + S, E = e^{j w1 n} (w1 possibly written outside [-pi, pi], or as (-1)^n, j^n,
(-j)^n), S = A cos(w2 n + theta) or A sin(w2 n + theta). Output written
y = D + B E + C+ e^{j w2 n} + C- e^{-j w2 n}; asked D, |B|, angle B, |C+|, angle C+, |C-|, angle C-
(angles as principal multiples of pi; exactly pi accepts 1 and -1 in grade()) and "is h real?".
Every magnitude is exact (rational, rational * pi, or a + b sqrt d) and every angle pi-rational;
H_d is never zero at w1 or +-w2 (so every angle is defined)."""

import cmath
import math
import random
from fractions import Fraction as F

import prairielearn as pl
import sympy as sp
from ece310 import fmt, mc
from ece310.w5_freqresp import Polar, cos_pi, dec, plain, principal, rat, sin_pi, tex_ej, tex_pi, tex_pi_frac, tex_r, tex_real

ANGLES = ("B_ang", "Cp_ang", "Cm_ang")
SIXTHS = [F(k, 6) for k in range(-5, 7)]
TWELFTHS = [F(k, 12) for k in range(-11, 13)]


def simple(x) -> bool:
    """Exact and 'exam nice': rational, rational * pi, or a sum of rationals and rational * sqrt(2|3|6)."""
    x = sp.nsimplify(x)
    if x.is_rational or (x / sp.pi).is_rational:
        return True
    for t in sp.Add.make_args(sp.expand(x)):
        if not (t.is_rational or any((t / sp.sqrt(d)).is_rational for d in (2, 3, 6))):
            return False
    return True


def _sgn_ex(s, fn):
    """e^{+-j pi fn w} exponent text."""
    return ("" if s == 1 else "-") + r"j\pi" + fn


# ----------------------------------------------------------------------------- system families
class System:
    kind = ""
    real = True
    cont_pi = True
    allowed = SIXTHS

    def polar(self, f):            # -> Polar or None (not exact / not nice)
        raise NotImplementedError

    def sub_tex(self, f):          # H_d at w = f pi, with the value substituted
        raise NotImplementedError

    def mid_tex(self, f):          # an optional intermediate step (polar forms of the parts)
        return ""

    def formula_num(self, w):      # closed forms: the formula evaluated at w as written (no reduction)
        return None


class AbsW(System):
    kind = "absw"

    def __init__(self):
        self.s = random.choice([1, -1])
        ex = _sgn_ex(self.s, r"\sin\omega")
        self.H_tex = rf"|\omega|\,e^{{{ex}}}"
        self.domain = r",\qquad |\omega| \le \pi"
        self.spec = {"kind": self.kind, "s": self.s}
        nx = _sgn_ex(-self.s, r"\sin\omega")
        self.real_html = (rf"$H_d(-\omega) = |{{-\omega}}|\,e^{{{_sgn_ex(self.s, r'\sin(-\omega)')}}} = |\omega|\,e^{{{nx}}}"
                          r" = H_d^*(\omega)$ (even magnitude, odd phase), so $h[n]$ is real.")

    def polar(self, f):
        sn = sp.nsimplify(sin_pi(f))
        if not sn.is_rational:
            return None
        return Polar(abs(rat(f)) * sp.pi, self.s * F(int(sn.p), int(sn.q)))

    def sub_tex(self, f):
        return rf"\left|{tex_pi(f)}\right| e^{{{_sgn_ex(self.s, r'\sin(' + tex_pi(f) + ')')}}}"

    def formula_num(self, w):
        return abs(w) * cmath.exp(1j * self.s * math.pi * math.sin(w))


class JW(System):
    kind = "jw"
    cont_pi = False

    def __init__(self):
        self.s = random.choice([1, -1])
        self.H_tex = rf"j\omega\,e^{{{_sgn_ex(self.s, r'\sin\omega')}}}"
        self.domain = r",\qquad |\omega| \le \pi"
        self.spec = {"kind": self.kind, "s": self.s}
        self.real_html = (rf"$H_d^*(\omega) = -j\omega\,e^{{{_sgn_ex(-self.s, r'\sin\omega')}}}$ and "
                          rf"$H_d(-\omega) = -j\omega\,e^{{{_sgn_ex(-self.s, r'\sin\omega')}}}$ agree, so $h[n]$ is real.")

    def polar(self, f):
        sn = sp.nsimplify(sin_pi(f))
        if not sn.is_rational:
            return None
        if f == 0:
            return Polar(0)
        return Polar(abs(rat(f)) * sp.pi, (F(1, 2) if f > 0 else F(-1, 2)) + self.s * F(int(sn.p), int(sn.q)))

    def sub_tex(self, f):
        return rf"j\left({tex_pi(f)}\right) e^{{{_sgn_ex(self.s, r'\sin(' + tex_pi(f) + ')')}}}"

    def formula_num(self, w):
        return 1j * w * cmath.exp(1j * self.s * math.pi * math.sin(w))


class WCos(System):
    kind = "wcos"
    real = False
    cont_pi = False

    def __init__(self):
        self.s = random.choice([1, -1])
        self.H_tex = rf"\omega\,e^{{{_sgn_ex(self.s, r'\cos\omega')}}}"
        self.domain = r",\qquad |\omega| \le \pi"
        self.spec = {"kind": self.kind, "s": self.s}
        self.real_html = (rf"$H_d(-\omega) = -\omega\,e^{{{_sgn_ex(self.s, r'\cos\omega')}}}$, but "
                          rf"$H_d^*(\omega) = \omega\,e^{{{_sgn_ex(-self.s, r'\cos\omega')}}}$: they differ, so $h[n]$ is "
                          "<b>not</b> real, and a sinusoid must be split into its two exponentials.")

    def polar(self, f):
        cs = sp.nsimplify(cos_pi(f))
        if not cs.is_rational:
            return None
        if f == 0:
            return Polar(0)
        return Polar(abs(rat(f)) * sp.pi, (1 if f < 0 else 0) + self.s * F(int(cs.p), int(cs.q)))

    def sub_tex(self, f):
        return rf"\left({tex_pi(f)}\right) e^{{{_sgn_ex(self.s, r'\cos(' + tex_pi(f) + ')')}}}"

    def formula_num(self, w):
        return w * cmath.exp(1j * self.s * math.pi * math.cos(w))


class HalfCos(System):
    kind = "halfcos"
    real = False
    allowed = sorted(set(SIXTHS) | {F(k, 4) for k in range(-3, 5)})

    def __init__(self):
        self.b = random.choice([F(1, 2), F(1, 2), F(1)])
        self.s = random.choice([1, -1])
        self.cont_pi = self.b != 1
        R = fmt.tex_sum([(self.b, ""), (1, r"\cos\omega")])
        ex = r"e^{j|\omega|}" if self.s == 1 else r"e^{-j|\omega|}"
        self.H_tex = rf"\left({R}\right){ex}"
        self.domain = r",\qquad |\omega| \le \pi"
        self.spec = {"kind": self.kind, "s": self.s, "b": plain(self.b)}
        self.real_html = (rf"The phase factor ${ex}$ is even in $\omega$, so $H_d(-\omega) = H_d(\omega)$, which is not "
                          r"$H_d^*(\omega)$: $h[n]$ is <b>not</b> real (a real $h$ has an odd phase).")

    def polar(self, f):
        R = sp.nsimplify(rat(self.b) + cos_pi(f))
        if not simple(R):
            return None
        P = Polar.from_real(R)
        if P.zero:
            return P
        return Polar(P.r, P.t + self.s * abs(f))

    def sub_tex(self, f):
        ex = tex_ej(self.s * abs(f))
        return rf"\left({fmt.tex_num(self.b)} + \cos\left({tex_pi(f)}\right)\right){ex}"

    def formula_num(self, w):
        return (float(self.b) + math.cos(w)) * cmath.exp(1j * self.s * abs(w))


class CosDelay(System):
    kind = "cosdelay"
    allowed = sorted(set(SIXTHS) | {F(k, 4) for k in range(-3, 5)})

    def __init__(self):
        self.k = random.choice([1, 2, 3])
        kt = "" if self.k == 1 else str(self.k)
        self.H_tex = rf"\left(1 + \cos\omega\right)e^{{-j{kt}\omega}}"
        self.domain = ""
        self.spec = {"kind": self.kind, "k": self.k}
        self.real_html = (rf"$H_d(-\omega) = \left(1 + \cos\omega\right)e^{{j{kt}\omega}} = H_d^*(\omega)$ (the factor "
                          r"$1 + \cos\omega$ is real and even, the phase odd), so $h[n]$ is real.")

    def polar(self, f):
        R = sp.nsimplify(1 + cos_pi(f))
        if not simple(R):
            return None
        if R == 0:
            return Polar(0)
        return Polar(R, -self.k * f)

    def sub_tex(self, f):
        return rf"\left(1 + \cos\left({tex_pi(f)}\right)\right){tex_ej(-self.k * f)}"

    def formula_num(self, w):
        return (1 + math.cos(w)) * cmath.exp(-1j * self.k * w)


FIR_BANK = [([1, 0, 1], "sym"), ([1, 1, 1], "sym"), ([1, -1, 1], "sym"), ([1, 2, 1], "sym"), ([2, 1, 2], "sym"),
            ([1, 0, 0, 1], "sym"), ([1, 0, 0, 0, 1], "sym"), ([1, 1], "sym"), ([1, 0, -1], "anti"),
            ([1, -1], "anti"), ([1, 0, 0, -1], "anti"), ([1, 0, 0, 0, -1], "anti")]


class FIR(System):
    kind = "fir"
    allowed = sorted(set(SIXTHS) | {F(k, 4) for k in range(-3, 5)})

    def __init__(self):
        taps, self.sym = random.choice(FIR_BANK)
        s = random.choice([1, 1, -1, 2])
        self.h = [s * v for v in taps]
        self.as_Hz = random.random() < 0.35
        self.start = 0 if self.as_Hz else random.choice([-1, 0, 0, 1, 2])
        L = len(self.h)
        self.c = F(L - 1, 2)
        self.alpha = self.start + self.c
        self.spec = {"kind": self.kind, "h": self.h, "start": self.start}
        if self.as_Hz:
            self.sys_tex = rf"H(z) = {fmt.tex_poly_zinv(self.h)}"
        else:
            self.sys_tex = rf"h[n] = {fmt.tex_seq(self.h, self.start)}"
        self.real_html = "The impulse response has real samples, so $h[n]$ is real."
        fn = r"\cos" if self.sym == "sym" else r"\sin"
        terms = {}
        for i, v in enumerate(self.h):
            k = self.c - i
            if v == 0:
                continue
            if self.sym == "sym":
                terms[abs(k)] = terms.get(abs(k), 0) + v
            elif k != 0:
                terms[abs(k)] = terms.get(abs(k), 0) + (v if k > 0 else -v)
        self.terms = sorted(((v, k) for k, v in terms.items() if v), key=lambda t: t[1])

        def kw(k):
            if k == 1:
                return r"\omega"
            if k.denominator == 1:
                return f"{k.numerator}\\omega"
            return r"\tfrac{" + ("" if k.numerator == 1 else str(k.numerator)) + r"\omega}{" + str(k.denominator) + "}"
        self.A_tex = fmt.tex_sum((v, "" if k == 0 else f"{fn}{kw(k)}" if k.denominator != 1 else f"{fn} {kw(k)}")
                                 for v, k in self.terms)
        lin = ("e^{" + ("-" if self.alpha > 0 else "") + "j" + (fmt.tex_num(abs(self.alpha)) if abs(self.alpha) != 1 else "")
               + r"\omega}") if self.alpha else ""
        self.H_tex = (r"j\," if self.sym == "anti" else "") + lin + rf"\left({self.A_tex}\right)"
        self.domain = ""

    def A(self, f):
        fn = sp.cos if self.sym == "sym" else sp.sin
        return sp.nsimplify(sum((v * fn(sp.pi * rat(k) * rat(f)) for v, k in self.terms), sp.Integer(0)))

    def polar(self, f):
        a = self.A(f)
        if not simple(a):
            return None
        P = Polar.from_real(a)
        if P.zero:
            return P
        return Polar(P.r, P.t - self.alpha * f + (F(1, 2) if self.sym == "anti" else 0))

    def sub_tex(self, f):
        lin = tex_ej(-self.alpha * f)
        a = self.A(f)
        if not lin and self.sym != "anti":
            return tex_r_signed(a)
        return (r"j\," if self.sym == "anti" else "") + (lin + r"\cdot" if lin else "") + r"\left(" + tex_r_signed(a) + r"\right)"


class LCCDE(System):
    kind = "lccde"

    def __init__(self):
        self.a = random.choice([F(1, 2), F(-1, 2)])
        self.q = random.choice([0, -1])
        self.g = random.choice([F(1), F(2), F(1, 2), F(3, 2)])
        self.f0 = F(1, 3) if self.a > 0 else F(2, 3)
        self.allowed = [F(0), F(1), self.f0, -self.f0]
        lhs = fmt.tex_sum([(1, "y[n]"), (-self.a, "y[n-1]")])
        rhs = fmt.tex_sum([(self.g, "x[n]"), (self.g * self.q, "x[n-1]")])
        self.sys_tex = f"{lhs} = {rhs}"
        num = fmt.tex_sum([(1, ""), (self.q, r"e^{-j\omega}")])
        den = fmt.tex_sum([(1, ""), (-self.a, r"e^{-j\omega}")])
        gt = "" if self.g == 1 else fmt.tex_num(self.g) + r"\,"
        self.H_tex = gt + fmt.tex_frac(num, den)
        self.domain = ""
        self.spec = {"kind": self.kind, "a": plain(self.a), "q": self.q, "g": plain(self.g)}
        self.real_html = "The difference equation has real coefficients (and the recursion runs on real numbers), so $h[n]$ is real."
        self.pole_tex = fmt.tex_num(self.a)

    def polar(self, f):
        a, q, g = self.a, self.q, self.g
        if f == 0:
            return Polar.from_real(rat(g * (1 + q) / (1 - a)))
        if f == 1:
            return Polar.from_real(rat(g * (1 - q) / (1 + a)))
        if abs(f) != self.f0:
            return None
        num, den = self.parts(f)
        return Polar(rat(g) * num.r / den.r, num.t - den.t)

    def parts(self, f):
        """Polar forms of 1 + q e^{-jw} and 1 - a e^{-jw} at w = +-f0 pi."""
        if self.q == 0:
            num = Polar(1, 0)
        else:   # 1 - e^{-jw} = 2j sin(w/2) e^{-jw/2}
            P = Polar.from_real(sp.nsimplify(2 * sin_pi(f / 2)))
            num = Polar(P.r, P.t + F(1, 2) - f / 2)
        dang = F(1, 6) if (self.a > 0) == (f > 0) else F(-1, 6)    # angle of 1 - a e^{-jw} at +-f0
        return num, Polar(sp.sqrt(3) / 2, dang)

    def mid_tex(self, f):
        if f in (0, 1):
            return ""
        num, den = self.parts(f)
        gt = "" if self.g == 1 else fmt.tex_num(self.g) + r"\cdot"
        return gt + fmt.tex_frac(num.tex(), den.tex())

    def sub_tex(self, f):
        e = tex_ej(-f) if f != 0 else "1"
        num = fmt.tex_sum([(1, ""), (self.q, e)]) if f != 0 else fmt.tex_num(1 + self.q)
        den = fmt.tex_sum([(1, ""), (-self.a, e)]) if f != 0 else fmt.tex_num(1 - self.a)
        gt = "" if self.g == 1 else fmt.tex_num(self.g) + r"\cdot"
        return gt + fmt.tex_frac(num, den)


def tex_r_signed(x):
    return tex_real(sp.nsimplify(x))


FAMILIES = [(AbsW, 14), (JW, 12), (WCos, 13), (HalfCos, 14), (CosDelay, 9), (FIR, 22), (LCCDE, 16)]

# input pieces ------------------------------------------------------------------------------------
COEFS = [  # (r, t, LaTeX prefix for a coefficient multiplying an exponential)
    (1, F(0), ""), (2, F(0), "2"), (3, F(0), "3"), (1, F(1), "-"), (2, F(1), "-2"),
    (1, F(1, 2), "j"), (2, F(1, 2), "2j"), (1, F(-1, 2), "-j"),
    (2, F(1, 4), r"2e^{j\frac{\pi}{4}}"), (3, F(-1, 3), r"3e^{-j\frac{\pi}{3}}"), (4, F(1, 6), r"4e^{j\frac{\pi}{6}}"),
]
THETAS = [F(0), F(1, 4), F(-1, 4), F(1, 3), F(-1, 3), F(1, 6), F(-1, 6), F(1, 2)]


def _cos_arg_tex(f2, th):
    base = tex_pi(f2) + "n"
    return base + ("" if th == 0 else (" + " + tex_pi(th) if th > 0 else " - " + tex_pi(-th)))


def _exp_written_tex(fw):
    """e^{j fw pi n} with the frequency as written (fw may lie outside (-1, 1])."""
    return tex_ej(fw, "n")


def _draw():
    while True:
        cls = random.choices([c for c, _ in FAMILIES], weights=[w for _, w in FAMILIES])[0]
        sysm = cls()
        pol = {}
        for f in sysm.allowed:
            P = sysm.polar(f)
            if P is not None:
                pol[f] = P
        H0 = sysm.polar(F(0))
        if H0 is None or H0.t not in (0, 1):
            continue
        # sinusoid frequency in (0, 1): H nonzero at +-f2
        f2s = [f for f in pol if 0 < f < 1 and -f in pol and not pol[f].zero and not pol[-f].zero]
        if not f2s:
            continue
        f2 = random.choice(f2s)
        f1s = [f for f in pol if f not in (0, f2, -f2) and not pol[f].zero and (f != 1 or sysm.cont_pi)]
        if not f1s:
            continue
        f1 = random.choice(f1s)
        return sysm, pol, f1, f2


def generate(data):
    sysm, pol, f1, f2 = _draw()
    p = data["params"]
    cr = data["correct_answers"]
    closed = sysm.kind in ("absw", "jw", "wcos", "halfcos", "cosdelay")
    p["spec"] = sysm.spec
    p["closed"] = closed
    p["is_fir"] = sysm.kind == "fir"
    p["is_lccde"] = sysm.kind == "lccde"
    p["H_tex"] = sysm.H_tex
    p["domain_tex"] = sysm.domain
    p["sys_tex"] = getattr(sysm, "sys_tex", "")
    p["fir_Hz"] = sysm.kind == "fir" and sysm.as_Hz
    p["fir_h"] = sysm.kind == "fir" and not sysm.as_Hz

    # ---- the input
    c0 = random.choice([1, 2, 3, 4, 5, -1, -2, -3])
    r1, t1, a1_tex = random.choice(COEFS)
    # how the exponential is written
    if f1 == 1:
        form = random.choice(["m1", "m1", "exp"])
        fw = random.choice([F(1), F(-1), F(3)]) if form == "exp" else F(1)
    elif abs(f1) == F(1, 2) and random.random() < 0.45:
        form, fw = ("jn" if f1 > 0 else "mjn"), f1
    else:
        form = "exp"
        out_prob = 0.45 if closed else 0.25
        fw = f1 + random.choice([2, -2]) if random.random() < out_prob else f1
    if form == "exp" and abs(fw) > 1 and principal(fw) != f1:
        fw = f1
    E_tex = {"m1": "(-1)^n", "jn": "j^n", "mjn": "(-j)^n"}.get(form) or _exp_written_tex(fw)
    p["E_written_outside"] = form == "exp" and abs(fw) > 1
    A = random.choice([2, 4, 6, 1, 3])
    th = random.choice(THETAS)
    trig = random.choice(["cos", "sin"])
    S_tex = (str(A) if A != 1 else "") + (r"\cos" if trig == "cos" else r"\sin") + r"\left(" + _cos_arg_tex(f2, th) + r"\right)"
    ex_term = a1_tex + E_tex if a1_tex in ("", "-") else a1_tex + r"\," + E_tex
    pieces = [str(c0), ex_term, S_tex]
    if random.random() < 0.5:
        pieces = [pieces[0], pieces[2], pieces[1]]
    x_tex = pieces[0]
    for pc in pieces[1:]:
        x_tex += (" - " + pc[1:]) if pc.startswith("-") else (" + " + pc)
    p["x_tex"] = x_tex
    p["E_tex"] = E_tex
    p["w2_tex"] = tex_pi(f2)
    p["y_form_tex"] = (rf"D + B\,{E_tex} + C_{{+}}\,{tex_ej(f2, 'n')} + C_{{-}}\,{tex_ej(-f2, 'n')}")
    p["inp"] = {"c": c0, "a1": [r1, plain(t1)], "fw": plain(fw), "f1": plain(f1), "A": A, "theta": plain(th),
                "trig": trig, "f2": plain(f2)}

    # ---- answers
    H0, H1, Hp, Hm = pol[F(0)], pol[f1], pol[f2], pol[-f2]
    D = sp.nsimplify(c0 * (H0.r if H0.t == 0 else -H0.r))
    cr["D"] = float(D)
    B = Polar(r1, t1) * H1
    thp = th - (F(1, 2) if trig == "sin" else 0)
    Cp = Polar(sp.Rational(A, 2), thp) * Hp
    Cm = Polar(sp.Rational(A, 2), -thp) * Hm
    for name, P in (("B", B), ("Cp", Cp), ("Cm", Cm)):
        cr[f"{name}_mag"] = float(P.r)
        cr[f"{name}_ang"] = float(P.t)
        p[f"{name}_tex"] = P.tex()
        p[f"{name}_mag_tex"] = tex_r(P.r)
        p[f"{name}_mag_dec"] = dec(P.r)
        p[f"{name}_irr"] = not P.r.is_rational
        p[f"{name}_ang_tex"] = tex_pi(P.t)
        p[f"{name}_ang_plain"] = plain(P.t)
        p[f"{name}_pi"] = abs(P.t) == 1
    p["real_choices"] = mc.yes_no(sysm.real)
    p["real_word"] = "real" if sysm.real else "not real"
    p["real_html"] = sysm.real_html
    p["is_real"] = sysm.real
    # exact values for the checker / probes
    p["Hvals"] = {k: [str(P.r), plain(P.t)] for k, P in (("0", H0), ("f1", H1), ("f2", Hp), ("mf2", Hm))}

    # ---- worked solution pieces
    def ev(f, P):
        mid = sysm.mid_tex(f)
        sub = sysm.sub_tex(f)
        final = P.tex()
        steps = [x for x in (sub, mid) if x and x != final]
        return rf"H_d\!\left({tex_pi(f)}\right) = " + " = ".join(steps + [final])
    p["ev_rows"] = [{"line": ev(F(0), H0)}, {"line": ev(f1, H1)}, {"line": ev(f2, Hp)}, {"line": ev(-f2, Hm)}]
    p["c0"] = c0
    p["H0_tex"] = H0.tex()
    p["D_tex"] = tex_r_signed(D)
    p["f1_tex"] = tex_pi(f1)
    if form == "m1":
        p["E_note_html"] = r"$(-1)^n = e^{j\pi n}$: frequency $\pi$."
    elif form == "jn":
        p["E_note_html"] = r"$j^n = e^{j\frac{\pi}{2}n}$: frequency $\tfrac{\pi}{2}$."
    elif form == "mjn":
        p["E_note_html"] = r"$(-j)^n = e^{-j\frac{\pi}{2}n}$: frequency $-\tfrac{\pi}{2}$."
    elif fw == -1:
        p["E_note_html"] = (r"$e^{-j\pi n} = e^{j\pi n} = (-1)^n$: frequency $\pi$ (or $-\pi$; for this system "
                            r"$H_d(-\pi) = H_d(\pi)$, so either endpoint gives the same value).")
    elif fw != f1:
        note = (rf"${_exp_written_tex(fw)} = {_exp_written_tex(f1)}$ because $e^{{j2\pi n}} = 1$: the frequency is "
                rf"${tex_pi(f1)}$")
        raw = sysm.formula_num(math.pi * float(fw))
        if sysm.kind == "cosdelay":
            note += (r". This $H_d$ is a $2\pi$-periodic formula, so plugging in $\omega = " + tex_pi(fw)
                     + r"$ directly gives the same value.")
        elif raw is not None:          # the other closed forms are stated for |w| <= pi only
            same = abs(raw - H1.value()) < 1e-9 * max(1.0, abs(raw))
            note += r", inside $[-\pi, \pi]$, where the formula for $H_d$ is stated" + (
                r" (plugging in $\omega = " + tex_pi(fw) + r"$ happens to give the same value here, but the formula is "
                r"only guaranteed on $[-\pi, \pi]$)." if same else
                r" (plugging $\omega = " + tex_pi(fw) + r"$ into the formula would give a different, wrong value).")
        else:
            note += r"; $H_d(\omega)$ is $2\pi$-periodic, so either form gives the same value."
        p["E_note_html"] = note
    else:
        p["E_note_html"] = rf"${E_tex}$: frequency ${tex_pi(f1)}$."
    p["a1_polar_tex"] = Polar(r1, t1).tex()
    p["half_tex"] = tex_r(sp.Rational(A, 2))
    th_tex = tex_pi_frac(th)
    half = fmt.tex_num(F(A, 2)) + r"\," if A != 2 else ""
    At = str(A) if A != 1 else ""
    if trig == "cos":
        p["split_tex"] = rf"{At}\cos\theta = {half}e^{{j\theta}} + {half}e^{{-j\theta}}"
    else:
        p["split_tex"] = (rf"{At}\sin\theta = \frac{{{A}}}{{2j}}\left(e^{{j\theta}} - e^{{-j\theta}}\right)"
                          rf" = {half}e^{{j(\theta - \frac{{\pi}}{{2}})}} + {half}e^{{-j(\theta - \frac{{\pi}}{{2}})}}")
    p["theta_tex"] = _cos_arg_tex(f2, th)
    p["thp_tex"] = tex_pi(thp)
    p["Cp_in_tex"] = Polar(sp.Rational(A, 2), thp).tex()
    p["Cm_in_tex"] = Polar(sp.Rational(A, 2), -thp).tex()
    p["Hp_tex"], p["Hm_tex"], p["H1_tex"] = Hp.tex(), Hm.tex(), H1.tex()
    p["y1_tex"] = (rf"{tex_r_signed(D)} + " if D != 0 else "") + rf"\left({B.tex()}\right){E_tex}"
    p["y2_tex"] = rf"+ \left({Cp.tex()}\right){tex_ej(f2, 'n')}"
    p["y3_tex"] = rf"+ \left({Cm.tex()}\right){tex_ej(-f2, 'n')}"
    p["fam"] = sysm.kind
    p["pole_tex"] = getattr(sysm, "pole_tex", "")


def grade(data):
    """An angle of exactly pi has two spellings: accept -1 as well as 1."""
    changed = False
    for name in ANGLES:
        corr = data["correct_answers"].get(name)
        ps = data["partial_scores"].get(name)
        if corr is None or ps is None or abs(abs(corr) - 1) > 1e-12 or name in data["format_errors"]:
            continue
        try:
            v = float(data["submitted_answers"].get(name))
        except (TypeError, ValueError):
            continue
        if abs(v + corr) <= 1e-6 + 1e-3 and (ps.get("score") or 0) < 1:
            ps["score"] = 1.0
            ps["feedback"] = "Correct: an angle of -pi and an angle of pi are the same angle."
            changed = True
    if changed:
        pl.set_weighted_score_data(data)
