"""Sinusoids and complex exponential sequences (Lecture 2: Euler's identities, sinusoids
A sin(w0 n + theta), exponentials B a^n with complex a). Every variant asks all three parts:

  (a) Euler: A cos(w0 n + theta) (or A sin) = c e^{j w0 n} + c* e^{-j w0 n}: |c| and the principal
      angle of c (c = (A/2) e^{j theta}, resp. (A/2j) e^{j theta}; negative A adds pi);
  (b) phasor addition: A1 cos(w0 n + th1) + A2 cos/sin(w0 n + th2) = A cos(w0 n + theta), A >= 0.
      Vetted families with exact A and pi-rational theta (phasors psi1, psi2 = psi1 + delta):
        eq      A1 = A2, delta in {+-pi/3, +-pi/2, +-2pi/3}  -> A = 2 A1 |cos(delta/2)|
        sqrt3   amplitudes k and k sqrt3 at right angles     -> A = 2k
        cancel  A2 = 2 A1, delta = +-2pi/3                    -> A = sqrt3 A1, theta = psi1 +- pi/2
        opp     delta = pi, A1 != A2                          -> A = |A1 - A2|
  (c) x[n] = B a^n with a = r e^{j phi} given in rectangular or polar form: |x[n]| as a function of n
      (pl-symbolic-input, n an integer) and Re, Im of x[k] (rational by construction).
Angles are entered as multiples of pi, principal values in (-pi, pi]; variants whose angle would be
+-pi (two spellings) or 0 are not generated."""

import cmath
import math
import random
from fractions import Fraction as Fr

from ece310 import answers, fmt
from ece310.w1_signals import Rad, cos_pi, plain_frac, principal, sin_pi, tex_ej, tex_phase_term, tex_pi

W0 = [Fr(1, 6), Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4), Fr(5, 6)]
G6 = [Fr(k, 6) for k in range(-5, 7)]
G4 = [Fr(k, 4) for k in range(-3, 5)]
G64 = sorted(set(G6 + G4))


def _sinusoid(coef_tex, fn, w, th):
    return rf"{coef_tex}\{fn}\!\left({tex_pi(w)}n{tex_phase_term(th)}\right)"


def _differs(x, y):
    return abs(x - y) > 1e-2 * max(1.0, abs(y))


def _polar(c, t):
    """c e^{j t pi} in LaTeX (c a Rad or a number): 'e^{..}' for c = 1, '-e^{..}' for c = -1, just c for t = 0."""
    c = c if isinstance(c, Rad) else Rad(c)
    e = tex_ej(t)
    if not e:
        return c.tex()
    coef = c.coef_tex()
    return coef + (r"\," if coef not in ("", "-") else "") + e


def _coef_before(x):
    """Coefficient x in front of a factor that always follows: '', '-', or '3' plus a thin space."""
    t = Rad(x).coef_tex()
    return t + (r"\," if t not in ("", "-") else "")


def _distinct(cands, right):
    out = []
    for m in cands:
        if _differs(m, right) and all(abs(m - o) > 1e-9 for o in out):
            out.append(m)
    return out


# ----------------------------------------------------------------------------- (a) Euler
def _part_a(p, w):
    while True:
        fn = random.choice(["cos", "sin"])
        A = random.choice([2, 3, 4, 5, 6, 8, -2, -4, -6])
        th = random.choice([t for t in G64 if t != 1])      # no "+ pi" inside the given sinusoid
        raw = th - (Fr(1, 2) if fn == "sin" else 0) + (1 if A < 0 else 0)
        ang = principal(raw)
        if ang in (0, 1):
            continue
        wrong = [float(-ang), float(principal(th))]
        if fn == "sin":
            wrong.append(float(principal(th + Fr(1, 2) + (1 if A < 0 else 0))))
        if A < 0:
            wrong.append(float(principal(th - (Fr(1, 2) if fn == "sin" else 0))))
        if len(_distinct(wrong, float(ang))) >= 2:
            break
    mag = Fr(abs(A), 2)
    x_tex = _sinusoid(str(A) if A not in (1, -1) else ("-" if A < 0 else ""), fn, w, th)
    p["xa_tex"] = x_tex
    p["a_fn"], p["a_A"], p["a_th"] = fn, A, [th.numerator, th.denominator]
    e_th = (r"\," + tex_ej(th)) if th else ""
    e_mth = (r"\," + tex_ej(-th)) if th else ""
    half = Fr(A, 2)
    arg = rf"\omega_0 n{tex_phase_term(th)}"
    ep = rf"e^{{j({arg})}}" if th else r"e^{j\omega_0 n}"
    em = rf"e^{{-j({arg})}}" if th else r"e^{-j\omega_0 n}"
    if fn == "cos":
        split = (rf"{A}\cdot\dfrac{{{ep} + {em}}}{{2}}"
                 rf" = \underbrace{{{_polar(half, th)}}}_{{c}}\,e^{{j\omega_0 n}}"
                 + (" - " if half < 0 else " + ") + _coef_before(abs(half)) + (tex_ej(-th) + r"\," if th else "")
                 + r"e^{-j\omega_0 n}")
        c_first = rf"c = {_polar(half, th)}"
    else:
        split = (rf"{A}\cdot\dfrac{{{ep} - {em}}}{{2j}}"
                 rf" = \underbrace{{\dfrac{{{A}}}{{2j}}{e_th}}}_{{c}}\,e^{{j\omega_0 n}} - "
                 rf"\dfrac{{{A}}}{{2j}}{e_mth}\,e^{{-j\omega_0 n}}")
        c_first = (rf"c = \dfrac{{{A}}}{{2j}}{e_th} = {_coef_before(half)}e^{{-j\pi/2}}{e_th}"
                   rf"\quad\left(\tfrac{{1}}{{j}} = -j = e^{{-j\pi/2}}\right)")
    p["a_split_tex"] = split
    p["a_c_tex"] = c_first
    p["a_raw_tex"] = tex_pi(raw)
    p["a_neg"] = A < 0
    p["a_absA"] = abs(A)
    p["a_reduced"] = raw != ang
    p["a_mag_tex"] = fmt.tex_num(mag)
    p["a_ang_tex"] = tex_pi(ang)
    p["a_ang_frac"] = plain_frac(ang)
    p["a_raw"] = [raw.numerator, raw.denominator]
    return float(mag), float(ang)


# ----------------------------------------------------------------------------- (b) phasor sum
def _part_b(p, w):
    while True:
        fam = random.choice(["eq", "eq", "sqrt3", "cancel", "opp"])
        swap = False
        if fam == "eq":
            d = random.choice([Fr(1, 3), Fr(-1, 3), Fr(1, 2), Fr(-1, 2), Fr(2, 3), Fr(-2, 3)])
            psi1 = random.choice(G6 if d.denominator == 3 else G4)
            A0 = random.choice([1, 2, 3, 4])
            A1, A2 = Rad(A0), Rad(A0)
            half = d / 2
            cs = cos_pi(half)
            Aex = (Rad(2 * A0) * cs).abs()
            th = psi1 + half + (0 if cs.q > 0 else 1)
            off = None
        elif fam == "sqrt3":
            d = random.choice([Fr(1, 2), Fr(-1, 2)])
            psi1 = random.choice(G6)
            k = random.choice([1, 2, 3])
            swap = random.random() < 0.4
            A1, A2 = (Rad(k, 3), Rad(k)) if swap else (Rad(k), Rad(k, 3))
            Aex = Rad(2 * k)
            off = d.numerator * (Fr(1, 6) if swap else Fr(1, 3))     # angle of A1 +- jA2 (unreduced)
            th = psi1 + off
        elif fam == "cancel":
            d = random.choice([Fr(2, 3), Fr(-2, 3)])
            psi1 = random.choice(G6)
            k = random.choice([1, 2, 3])
            A1, A2 = Rad(k), Rad(2 * k)
            Aex = Rad(k, 3)
            th = psi1 + (Fr(1, 2) if d > 0 else Fr(-1, 2))     # 1 + 2e^{+-j2pi/3} = +-j sqrt3
            off = None
        else:
            d = Fr(1)
            psi1 = random.choice(G64)
            a1, a2 = random.sample([1, 2, 3, 4, 5], 2)
            A1, A2 = Rad(a1), Rad(a2)
            Aex = Rad(abs(a1 - a2))
            th = psi1 + (0 if a1 > a2 else 1)
            off = None
        psi2 = principal(psi1 + d)
        raw = th
        th = principal(th)
        fn2 = random.choice(["cos", "sin"])
        th1, th2 = principal(psi1), principal(psi2 + (Fr(1, 2) if fn2 == "sin" else 0))
        P = float(A1) * cmath.exp(1j * math.pi * float(psi1)) + float(A2) * cmath.exp(1j * math.pi * float(psi2))
        assert abs(abs(P) - float(Aex)) < 1e-9 and abs(cmath.exp(1j * math.pi * float(th)) * float(Aex) - P) < 1e-9, \
            (fam, psi1, d, P, float(Aex), th)
        if th in (0, 1) or th1 == 0 and th2 == 0 or 1 in (th1, th2):
            continue
        wrong_th = [float(-th), float(th1)]
        if fn2 == "sin":   # sine read as a cosine: phasor angle th2 instead of th2 - pi/2
            Pw = float(A1) * cmath.exp(1j * math.pi * float(th1)) + float(A2) * cmath.exp(1j * math.pi * float(th2))
            if abs(Pw) > 1e-9:
                wrong_th.append(cmath.phase(Pw) / math.pi)
        if len(_distinct(wrong_th, float(th))) < 2:
            continue
        break
    t1 = _sinusoid(A1.coef_tex(), "cos", w, th1)
    t2 = _sinusoid(A2.coef_tex(), fn2, w, th2)
    p["xb_tex"] = t1 + " + " + t2
    p["b_fam"] = fam
    p["b_terms"] = [[float(A1), "cos", th1.numerator, th1.denominator], [float(A2), fn2, th2.numerator, th2.denominator]]

    ph = _polar
    conv = ""
    if fn2 == "sin":
        conv = (rf"First write the sine as a cosine: $\sin(\phi) = \cos(\phi - \tfrac{{\pi}}{{2}})$, so the second term is "
                rf"${_sinusoid(A2.coef_tex(), 'cos', w, psi2)}$. ")
    p["b_conv_html"] = conv
    p["b_phasors_tex"] = rf"P = {ph(A1, psi1)} + {ph(A2, psi2)}"
    mid = psi1 + d / 2
    jA2 = "j" if A2 == Rad(1) else "j" + A2.tex()
    if fam == "eq":
        why = (rf"equal amplitudes: factor out the average angle, $P = {ph(A1, mid)}\left({tex_ej(-d / 2)} + {tex_ej(d / 2)}\right)"
               rf" = {(A1.tex() + r'\cdot ') if A1 != Rad(1) else ''}2\cos\!\left({tex_pi(d / 2)}\right){tex_ej(mid)}"
               rf" = {ph(Aex, mid)}$")
    elif fam == "sqrt3":
        why = (rf"the phasors are at right angles: $P = {tex_ej(psi1)}\left({A1.tex()} {'+' if d > 0 else '-'} {jA2}\right)$, "
               rf"and ${A1.tex()} {'+' if d > 0 else '-'} {jA2} = {ph(Aex, off)}$, so $P = {ph(Aex, raw)}$")
    elif fam == "cancel":
        why = (rf"$P = {ph(A1, psi1)}\left(1 + 2{tex_ej(d)}\right)$ and "
               rf"$1 + 2{tex_ej(d)} = 1 + 2\left(-\tfrac{{1}}{{2}} {'+' if d > 0 else '-'} j\tfrac{{\sqrt{{3}}}}{{2}}\right) = "
               rf"{'' if d > 0 else '-'}j\sqrt{{3}} = \sqrt{{3}}\,{tex_ej(Fr(1, 2) if d > 0 else Fr(-1, 2))}$: the real parts cancel, "
               rf"so $P = {ph(Aex, raw)}$")
    else:
        diff = A1.q - A2.q
        e1 = tex_ej(psi1)
        why = (rf"the phasors point in opposite directions, so the amplitudes subtract: "
               rf"$P = ({A1.tex()} - {A2.tex()}){e1} = {ph(Rad(diff), psi1)}$")
        if diff < 0:
            why += (rf", and with $-{fmt.tex_num(-diff)} = {_polar(Rad(-diff), Fr(1))}$ this is "
                    rf"$P = {_coef_before(-diff)}e^{{j\pi}}{(r'\,' + e1) if e1 else ''} = {ph(Aex, raw)}$")
    p["b_why_html"] = why
    p["b_A_tex"] = Aex.tex()
    p["b_A_dec"] = f"{float(Aex):.6g}"
    p["b_A_irr"] = not Aex.is_rational()
    p["b_th_tex"] = tex_pi(th)
    p["b_raw_tex"] = tex_pi(raw)
    p["b_reduced"] = raw != th
    p["b_th_frac"] = plain_frac(th)
    p["b_sum_amp_tex"] = (rf"{A1.tex()} + {A2.tex()} = {fmt.tex_num(A1.q + A2.q)}" if A1.is_rational() and A2.is_rational()
                          else rf"{A1.tex()} + {A2.tex()} \approx {float(A1) + float(A2):.4g}")
    return float(Aex), float(th)


# ----------------------------------------------------------------------------- (c) B a^n
RECT = [("1 + j", Rad(1, 2), Fr(1, 4)), ("1 - j", Rad(1, 2), Fr(-1, 4)), ("-1 + j", Rad(1, 2), Fr(3, 4)),
        ("-1 - j", Rad(1, 2), Fr(-3, 4)), (r"\tfrac{1 + j}{2}", Rad(Fr(1, 2), 2), Fr(1, 4)),
        (r"\tfrac{1 - j}{2}", Rad(Fr(1, 2), 2), Fr(-1, 4)), (r"\tfrac{-1 + j}{2}", Rad(Fr(1, 2), 2), Fr(3, 4)),
        (r"\sqrt{3} + j", Rad(2), Fr(1, 6)), (r"\sqrt{3} - j", Rad(2), Fr(-1, 6)), (r"-\sqrt{3} + j", Rad(2), Fr(5, 6)),
        (r"1 + j\sqrt{3}", Rad(2), Fr(1, 3)), (r"1 - j\sqrt{3}", Rad(2), Fr(-1, 3)), (r"-1 + j\sqrt{3}", Rad(2), Fr(2, 3)),
        (r"\tfrac{1 + j\sqrt{3}}{4}", Rad(Fr(1, 2)), Fr(1, 3)), (r"\tfrac{-1 + j\sqrt{3}}{4}", Rad(Fr(1, 2)), Fr(2, 3)),
        (r"\tfrac{\sqrt{3} - j}{4}", Rad(Fr(1, 2)), Fr(-1, 6)), ("2j", Rad(2), Fr(1, 2)), (r"-\tfrac{j}{2}", Rad(Fr(1, 2)), Fr(-1, 2)),
        ("-3j", Rad(3), Fr(-1, 2))]
R_POLAR = [Fr(1, 2), Fr(2), Fr(3, 2), Fr(1, 3), Fr(3), Fr(2, 3)]


def _exp_tex(phi):
    """e^{j\\pi/3}, e^{-j3\\pi/4}."""
    s = "-" if phi < 0 else ""
    pn, q = abs(phi.numerator), phi.denominator
    return f"e^{{{s}j" + ("" if pn == 1 else str(pn)) + r"\pi" + ("" if q == 1 else f"/{q}") + "}"


def _r_sym(r):
    if r.d == 1:
        return fmt.sym_num(r.q)
    return f"({r.q.numerator}*sqrt({r.d})/{r.q.denominator})"


def _part_c(p):
    while True:
        if random.random() < 0.55:
            a_tex, r, phi = random.choice(RECT)
            polar = False
        else:
            rr = random.choice(R_POLAR)
            phi = random.choice([t for t in G64 if t not in (0, 1)])
            r = Rad(rr)
            a_tex = (fmt.tex_num(rr) if rr != 1 else "") + _exp_tex(phi)
            polar = True
        B = random.choice([1, 2, 3, -1, -2, 4])
        k = random.choice([2, 3, 4, 5, -1, -2])
        rk = r ** k
        re, im = Rad(B) * rk * cos_pi(k * phi), Rad(B) * rk * sin_pi(k * phi)
        if not (re.is_rational() and im.is_rational()):
            continue
        mk = abs(B) * float(rk)
        if mk > 64 or mk < Fr(1, 64) or any(v.q.denominator > 32 or abs(v.q.numerator) > 128 for v in (re, im)):
            continue
        break
    Bc = {1: "", -1: "-"}.get(B, str(B))
    p["xc_tex"] = rf"{Bc}\left({a_tex}\right)^{{n}}"
    p["c_polar"] = polar
    p["c_rect"] = not polar
    p["c_a_tex"] = a_tex
    p["c_r_tex"] = r.tex()
    p["c_phi_tex"] = tex_pi(phi)
    p["c_absB"] = abs(B)
    p["c_k"] = k
    p["c_B"] = B
    mag_sym = (f"{abs(B)}*" if abs(B) != 1 else "") + f"{_r_sym(r)}**n"
    body = rf"\left({r.tex()}\right)^{{n}}" if (r.d != 1 or r.q.denominator != 1) else rf"{r.tex()}^{{n}}"
    p["c_mag_tex"] = fmt.tex_sum([(abs(B), body)])
    kphi = k * phi
    Brk = rk * Rad(B)
    p["c_xk_polar_tex"] = _polar(Brk, kphi) + ((" = " + _polar(Brk, principal(kphi))) if principal(kphi) != kphi else "")
    p["c_re_tex"], p["c_im_tex"] = re.tex(), im.tex()
    p["c_mag_sym"] = mag_sym
    p["c_rect_steps_tex"] = (rf"\lvert {a_tex}\rvert = {r.tex()},\quad \angle\left({a_tex}\right) = {tex_pi(phi)},\quad "
                             rf"{a_tex} = {r.tex() if r != Rad(1) else ''}{_exp_tex(phi)}") if not polar else ""
    return mag_sym, float(re), float(im)


# ----------------------------------------------------------------------------- generate
def generate(data):
    p = data["params"]
    w = random.choice(W0)
    p["w0_tex"] = tex_pi(w)
    p["w0"] = [w.numerator, w.denominator]
    wb = random.choice([x for x in W0 if x != w])
    cmag, cang = _part_a(p, w)
    A, th = _part_b(p, wb)
    mag_sym, re, im = _part_c(p)
    c = data["correct_answers"]
    c["cmag"], c["cang"] = cmag, cang
    c["amp"], c["phase"] = A, th
    c["xmag"] = answers.symbolic(mag_sym)
    c["xre"], c["xim"] = re, im
