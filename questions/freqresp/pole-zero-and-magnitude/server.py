"""Match pole-zero plots to the shape of |H_d(w)| using the distance rule
|H_d(w)| = |K| prod_k |e^{jw} - z_k| / prod_k |e^{jw} - p_k| (study site, concepts/poles-and-zeros,
"In Unit 3"; Lecture 16 Exercise 1: 1 - z^-2 with zeros at +-1 gives 2|sin w|). Enrichment built on
Lecture 9 (poles and zeros) and Lecture 16 (magnitude response); not an exam family.

Four systems, each of one of four shapes (shapes may repeat, at least two distinct, at least one
low-pass or high-pass), randomized but well-separated configurations:
  lowpass   zero(s) at z = -1, real pole p in {1/2, 3/5, 3/4, 4/5} (a double zero adds a pole at z = 0)
  highpass  zero(s) at z = 1, real pole -p
  bandpass  zeros at z = +-1, poles r e^{+-j theta}, r in {3/4, 4/5}
  bandstop  zeros e^{+-j theta} on the unit circle, poles r e^{+-j theta} just inside (a notch at theta)
For one low-pass or high-pass system the gain K of H(z) = K prod(1 - z_k z^-1) / prod(1 - p_k z^-1)
is stated and |H_d(0)| or |H_d(pi)| is asked (its pole and zero values are labelled on its plot, above the
axis next to their symbols; j and -j are labelled outside the unit circle).
pl-matching with partial credit: each correct match is worth 1/4."""

import cmath
import math
import random
from fractions import Fraction as F

from ece310 import fmt, svg

TYPES = [("lowpass", "low-pass"), ("highpass", "high-pass"), ("bandpass", "band-pass"),
         ("bandstop", "band-stop (notch)")]
THETAS = [F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4)]


def _c(r, t):
    """r e^{j t pi} as a complex number."""
    return float(r) * cmath.exp(1j * math.pi * float(t))


def _make(kind):
    if kind in ("lowpass", "highpass"):
        p = random.choice([F(1, 2), F(3, 5), F(3, 4), F(4, 5)])
        m = random.choice([1, 1, 2])
        zr = F(-1) if kind == "lowpass" else F(1)
        pr = p if kind == "lowpass" else -p
        # (1 -+ z^-1)^2 / (1 - p z^-1) = (z -+ 1)^2 / (z (z - p)): a double zero brings a pole at z = 0
        poles = [complex(pr)] + ([0j] if m == 2 else [])
        return dict(kind=kind, poles=poles, zeros=[complex(zr)] * m, real_poles=[pr], real_zeros=[zr] * m)
    r, t = random.choice([F(3, 4), F(4, 5)]), random.choice(THETAS)
    if kind == "bandpass":
        return dict(kind=kind, poles=[_c(r, t), _c(r, -t)], zeros=[1 + 0j, -1 + 0j], theta=t, r=r)
    return dict(kind=kind, poles=[_c(r, t), _c(r, -t)], zeros=[_c(1, t), _c(1, -t)], theta=t, r=r)


W, LIM = 260, 1.6


def _fig(s, label):
    """svg.pole_zero with +-1 ticks only; the j / -j labels are drawn just outside the unit circle and
    the real pole/zero values above the axis next to their symbol (inside the circle), so that no
    label sits on the dashed circle."""
    fig = svg.pole_zero(s["poles"], s["zeros"], lim=LIM, ticks=[(-1, "−1"), (1, "1")], iticks=[],
                        width=W, max_width=270, title=f"Pole-zero plot of system {label}")
    sc = (W - 52) / (2 * LIM)                      # the frame of svg.pole_zero: 26 px margins

    def X(x):
        return 26 + (float(x) + LIM) * sc

    def Y(y):
        return W - 26 - (float(y) + LIM) * sc
    xz, yz = X(0), Y(0)
    extra = []
    for v, lab, dy in ((1, "j", -7), (-1, "−j", 17)):
        extra.append(f'<line x1="{xz - 3:.1f}" y1="{Y(v):.1f}" x2="{xz + 3:.1f}" y2="{Y(v):.1f}" stroke="currentColor"/>')
        extra.append(f'<text x="{xz + 5:.1f}" y="{Y(v) + dy:.1f}" font-size="12">{lab}</text>')
    values = sorted(set((s.get("real_poles") or []) + (s.get("real_zeros") or [])))
    for v in values:
        if abs(v) == 1:
            continue
        x = X(v)
        extra.append(f'<line x1="{x:.1f}" y1="{yz - 3:.1f}" x2="{x:.1f}" y2="{yz + 3:.1f}" stroke="currentColor"/>')
        anchor, tx = ("end", x - 10) if v > 0 else ("start", x + 10)
        extra.append(f'<text x="{tx:.1f}" y="{yz - 9:.1f}" text-anchor="{anchor}" font-size="12">{svg.fmt_num(v)}</text>')
    return fig.replace("</svg>", "".join(extra) + "</svg>")


def tex_pi(t):
    t = F(t)
    p, q = t.numerator, t.denominator
    num = r"\pi" if p == 1 else f"{p}\\pi"
    return num if q == 1 else r"\tfrac{" + num + "}{" + str(q) + "}"


def _why(s):
    k = s["kind"]
    zw = "a double zero" if len(s["zeros"]) > 1 else "the zero"
    origin = (" (the pole at $z = 0$ that comes with the double zero has distance 1 from every point of the circle "
              "and changes nothing)") if len(s["zeros"]) > 1 else ""
    if k == "lowpass":
        return (f"{zw} at $z = -1$ makes $|H_d(\\pi)| = 0$; the pole at $z = {fmt.tex_num(s['real_poles'][0])}$, on the "
                "positive real axis and nearer $z = 1$ than $z = -1$, makes $|H_d|$ largest at $\\omega = 0$"
                f"{origin}: <b>low-pass</b>.")
    if k == "highpass":
        return (f"{zw} at $z = 1$ makes $|H_d(0)| = 0$; the pole at $z = {fmt.tex_num(s['real_poles'][0])}$, on the "
                "negative real axis and nearer $z = -1$ than $z = 1$, makes $|H_d|$ largest at $\\omega = \\pi$"
                f"{origin}: <b>high-pass</b>.")
    if k == "bandpass":
        return ("the zeros at $z = \\pm 1$ give $|H_d(0)| = |H_d(\\pi)| = 0$ (as in Lecture 16, Exercise 1), and the poles close "
                f"to the unit circle at angles $\\pm{tex_pi(s['theta'])}$ make a peak near $\\omega = {tex_pi(s['theta'])}$: "
                "<b>band-pass</b>.")
    return (f"the zeros <i>on</i> the unit circle at angles $\\pm{tex_pi(s['theta'])}$ give an exact null at "
            f"$\\omega = {tex_pi(s['theta'])}$; the poles just inside at the same angles nearly cancel their effect at every "
            "other frequency: a narrow <b>band-stop (notch)</b> filter.")


def generate(data):
    names = [k for k, _ in TYPES]
    while True:
        kinds = [random.choice(names) for _ in range(4)]
        if len(set(kinds)) >= 2 and ({"lowpass", "highpass"} & set(kinds)):
            break
    systems = [_make(k) for k in kinds]
    p = data["params"]
    p["stmts"] = [{"match": s["kind"], "label": i + 1, "fig_html": _fig(s, i + 1)} for i, s in enumerate(systems)]
    p["type_opts"] = [{"name": k, "text": t} for k, t in TYPES]
    p["systems"] = [{"kind": s["kind"], "poles": [[z.real, z.imag] for z in s["poles"]],
                     "zeros": [[z.real, z.imag] for z in s["zeros"]]} for s in systems]
    p["why"] = [{"label": i + 1, "text_html": _why(s)} for i, s in enumerate(systems)]

    # numeric part: one low-pass / high-pass system with a stated gain
    i0 = random.choice([i for i, s in enumerate(systems) if s["kind"] in ("lowpass", "highpass")])
    s = systems[i0]
    K = random.choice([F(1, 10), F(1, 5), F(1, 4), F(1, 2), F(1), F(2)])
    zeros, poles = s["real_zeros"], s["real_poles"]

    def mag_at(zv):
        num = 1
        for z in zeros:
            num *= (1 - z * zv)
        den = 1
        for q in poles:
            den *= (1 - q * zv)
        return abs(K * num / den)
    m0, mpi = mag_at(F(1)), mag_at(F(-1))
    at_pi = random.random() < (0.8 if m0 == 0 else 0.2)          # usually the end where |H_d| is not 0
    val = mpi if at_pi else m0
    p["num_label"] = i0 + 1
    p["num_origin"] = len(zeros) > 1
    p["K_tex"] = fmt.tex_num(K)
    p["w0_tex"] = r"\pi" if at_pi else "0"
    p["at_pi"] = at_pi
    zinv = "-1" if at_pi else "1"          # z^{-1} at z = e^{j w0}
    num_tex = "".join(r"\left(1 - " + fmt.tex_num(z, paren=True) + r"\cdot\left(" + zinv + r"\right)\right)" for z in zeros)
    den_tex = "".join(r"\left(1 - " + fmt.tex_num(q, paren=True) + r"\cdot\left(" + zinv + r"\right)\right)" for q in poles)
    kt = "" if K == 1 else fmt.tex_num(K) + r"\cdot"
    p["num_work_tex"] = (rf"\left|H_d({p['w0_tex']})\right| = \left|{kt}"
                         + fmt.tex_frac(num_tex, den_tex) + rf"\right| = {fmt.tex_num(val)}")
    p["H_num_tex"] = (("" if K == 1 else fmt.tex_num(K) + r"\,")
                      + fmt.tex_frac("".join(fmt.tex_factor(z) for z in zeros), "".join(fmt.tex_factor(q) for q in poles)))
    p["K"] = str(K)
    p["num_zeros"] = [str(z) for z in zeros]
    p["num_poles"] = [str(q) for q in poles]
    data["correct_answers"]["mag"] = float(val)
