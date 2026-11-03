"""LTI response to a sum of a constant, a sinusoid, a complex exponential and (-1)^n, with the
system given only by plots of |H_d(w)| and angle H_d(w) (Lecture 15: eigenfunctions, the cosine rule
for a real h[n]; Lecture 16: reading magnitude and phase).

The system is always the one of FA2025 Midterm 2 #2: |H_d| piecewise constant with two stacked
levels, L1 on |w| <= w1 and L2 on w1 < |w| <= w2 (0 beyond; or L2 all the way to pi), edges
multiples of pi/4, pi/6 or pi/8, and angle H_d = -M w on the passband (M in {1, 2, 3}, drawn as
a line, M w2 <= 3pi/2). Two figures built with ece310.svg.
The input reproduces the past Midterm 2 input shapes, one variant type each (new numbers):
  fa25  A0 + A1 cos(w1 n) + B j^n + C (-1)^n            FA2025 #2(b): 4 + cos(pi n/6) + 3j^n + (-1)^n
  fa24  A0 + A1 cos(w1 n) + B j^n                       FA2024 #4(b): 3 + cos(pi n/6) + j^n
  fa19  A0 + A1 cos(w1 n) + B e^{j w2 n}                FA2019 #7:    3 + 4cos(pi n/4) + e^{j pi n/2}
  sp25  A0 + B e^{j w2 n} + A1 sin(w1 n + th) + C(-1)^n  SP2025 #4:    5 + 2e^{j pi n/4} + sin(pi n/2 + pi/4) + (-1)^n
  sp23  A0 + B e^{j w2 n} + A1 sin(w1 n + th)           SP2023 #6:    3 + e^{j pi n/3} + sin(pi n/2 + pi/4)
(j^n is sometimes (-j)^n; e^{j w2 n} may have a negative frequency.) No input frequency sits on a
band edge; the sinusoid and the exponential are never both removed.
Asked, in the stated form y = D + E trig(w1 n + phi) + F (exp term) [+ G (-1)^n]:
  (a) D = A0 L1; (b) E >= 0 and phi (principal, multiple of pi; 0 and 0 if removed);
  (c) |F| and angle F (0 and 0 if removed); (d) G = C H_d(pi) (0 unless the outer level reaches pi:
  then C L2 (-1)^M) -- only for the shapes with a (-1)^n term.
No phase answer is exactly +-pi."""

import random
from fractions import Fraction as F

from ece310 import fmt, svg
from ece310.w4_dtft import step_svg
from ece310.w5_freqresp import plain, principal, tex_ej, tex_pi

SHAPES = [("fa25", 30), ("fa24", 15), ("fa19", 18), ("sp25", 22), ("sp23", 15)]
SOURCE = {
    "fa25": r"FA2025 Midterm 2 #2(b) ($x[n] = 4 + \cos\left(\tfrac{\pi}{6}n\right) + 3j^n + (-1)^n$)",
    "fa24": r"FA2024 Midterm 2 #4(b) ($x[n] = 3 + \cos\left(\tfrac{\pi n}{6}\right) + j^n$)",
    "fa19": r"FA2019 Midterm 2 #7 ($x[n] = 3 + 4\cos\left(\tfrac{\pi}{4}n\right) + e^{j\frac{\pi}{2}n}$)",
    "sp25": r"SP2025 Midterm 2 #4 ($x[n] = 5 + 2e^{j\frac{\pi}{4}n} + \sin\left(\tfrac{\pi}{2}n + \tfrac{\pi}{4}\right) + (-1)^n$)",
    "sp23": r"SP2023 Midterm 2 #6 ($x[n] = 3 + e^{j\frac{\pi}{3}n} + \sin\left(\tfrac{\pi}{2}n + \tfrac{\pi}{4}\right)$)",
}
LEVELS = [(2, 1), (2, 1), (3, 1), (3, 2), (4, 2), (4, 1), (1, 2), (1, 3)]
SIN_FREQS = [F(k, 12) for k in (2, 3, 4, 6, 8, 9, 10)]
EXP_FREQS = [F(k, 12) for k in (2, 3, 4, 6, 8, 9)]


# ----------------------------------------------------------------------------- the filter
def _draw_filter(to_pi):
    while True:
        q = random.choice([4, 6, 8])
        ks = range(1, q)
        w1 = F(random.choice(list(ks)), q)
        if to_pi:
            w2 = F(1)
        else:
            w2 = F(random.choice([k for k in ks if F(k, q) > w1]), q) if w1 < F(q - 1, q) else None
            if w2 is None:
                continue
        if w1 < F(1, 6) or w2 - w1 < F(1, 6) or (w2 != 1 and 1 - w2 < F(1, 6)):
            continue
        Ms = [M for M in (1, 1, 2, 3) if M * w2 <= F(3, 2)]
        M = random.choice(Ms)
        L1, L2 = random.choice(LEVELS)
        return {"w1": w1, "w2": w2, "M": M, "L1": L1, "L2": L2, "to_pi": to_pi}


def band(fl, t):
    """1 (inner level), 2 (outer level) or 0 (stopped) for |t| (multiple of pi), never on an edge."""
    a = abs(t)
    if a < fl["w1"]:
        return 1
    if a < fl["w2"] or (fl["to_pi"] and a <= 1):
        return 2
    return 0


def gain(fl, t):
    return {1: fl["L1"], 2: fl["L2"], 0: 0}[band(fl, t)]


def on_edge(fl, t):
    return abs(t) in (fl["w1"], fl["w2"])


# ----------------------------------------------------------------------------- figures
def mag_svg(fl):
    w1, w2, L1, L2 = fl["w1"], fl["w2"], fl["L1"], fl["L2"]
    if fl["to_pi"]:
        pts = [(-1, L2), (-w1, L2), (-w1, L1), (w1, L1), (w1, L2), (1, L2)]
        xt = [F(-1), -w1, w1, F(1)]
    else:
        pts = [(-1, 0), (-w2, 0), (-w2, L2), (-w1, L2), (-w1, L1), (w1, L1), (w1, L2), (w2, L2), (w2, 0), (1, 0)]
        xt = [F(-1), -w2, -w1, w1, w2, F(1)]
    return step_svg([pts], max(L1, L2), xt, sorted({L1, L2}), ylabel="|Hd(ω)|", height=210,
                    title="Magnitude response")


def phase_svg(fl):
    """-M w on |w| <= w2 drawn as a line, 0 outside (as in FA2025 #2); tick labels in a row at the
    bottom so the dashed guides at the band edges do not run through them."""
    w2, M = fl["w2"], fl["M"]
    top = float(M * w2)
    width, height = 560, 230
    f = svg._Frame((-1.12, 1.12), (-1.3 * top, 1.3 * top), width, height, left=48, right=34, top=30, bottom=30)
    out = svg._open(width, height, 560, "Phase response")
    xt = [(float(x), svg.fmt_pi(x)) for x in ([F(-1), -w2, w2, F(1)] if w2 != 1 else [F(-1), F(1)])]
    yt = [(top, svg.fmt_pi(M * w2)), (-top, svg.fmt_pi(-M * w2))]
    out += svg._axes(f, xt, yt, "ω", "∠Hd(ω)", False, xlabel_row=f.Y(-1.3 * top) + 16)
    dash = 'stroke="currentColor" stroke-opacity="0.55" stroke-dasharray="4 4" stroke-width="1"'
    guides = [((0, top), (-float(w2), top)), ((0, -top), (float(w2), -top))]
    if w2 != 1:
        guides += [((-float(w2), 0), (-float(w2), top)), ((float(w2), 0), (float(w2), -top))]
    for (a, b), (c, d) in guides:
        out.append(f'<line x1="{f.X(a):.1f}" y1="{f.Y(b):.1f}" x2="{f.X(c):.1f}" y2="{f.Y(d):.1f}" {dash}/>')
    segs = [[(-float(w2), top), (float(w2), -top)]]
    if w2 != 1:
        segs += [[(-1.0, 0.0), (-float(w2), 0.0)], [(float(w2), 0.0), (1.0, 0.0)]]
    for pts in segs:
        d = " ".join(f"{f.X(x):.1f},{f.Y(y):.1f}" for x, y in pts)
        out.append(f'<polyline points="{d}" fill="none" stroke="{svg.ACCENT}" stroke-width="2.4" stroke-linejoin="round"/>')
    out.append("</svg>")
    return "".join(out)


# ----------------------------------------------------------------------------- the input
def _arg_tex(t, th=F(0), style="n"):
    """w n + th in LaTeX: \\tfrac{\\pi}{6}n, \\tfrac{\\pi}{2}n + \\tfrac{\\pi}{4}."""
    s = tex_pi(t) + ("" if tex_pi(t).endswith("}") else " ") + "n"
    if th > 0:
        s += " + " + tex_pi(th)
    elif th < 0:
        s += " - " + tex_pi(-th)
    return s


def _draw_input(shape, fl):
    for _ in range(400):
        inp = {"shape": shape}
        if shape in ("fa25", "fa24", "fa19"):
            inp["trig"] = "cos"
            inp["A1"] = random.choice([1, 1, 2] if shape != "fa19" else [2, 3, 4, 4])
            inp["th"] = F(0) if random.random() < 0.8 else random.choice([F(1, 4), F(-1, 4), F(1, 3)])
        else:
            inp["trig"] = "sin"
            inp["A1"] = random.choice([1, 1, 2])
            inp["th"] = random.choice([F(1, 4), F(1, 4), F(-1, 4), F(1, 3), F(1, 6)])
        inp["w1"] = random.choice(SIN_FREQS)
        if shape in ("fa25", "fa24"):
            inp["w2"] = F(1, 2) if random.random() < 0.8 else F(-1, 2)
            inp["exp_form"] = "jn" if inp["w2"] > 0 else "mjn"
        else:
            inp["w2"] = random.choice(EXP_FREQS) * (1 if random.random() < 0.8 else -1)
            inp["exp_form"] = "exp"
        inp["A0"] = random.choice({"fa25": [2, 3, 4, 5], "fa24": [2, 3, 4], "fa19": [2, 3, 4],
                                   "sp25": [3, 4, 5], "sp23": [2, 3]}[shape])
        inp["B"] = random.choice({"fa25": [1, 2, 3, 3], "fa24": [1, 1, 2], "fa19": [1, 1, 2],
                                  "sp25": [1, 2, 2, 3], "sp23": [1, 1, 2]}[shape])
        inp["has_alt"] = shape in ("fa25", "sp25")
        inp["C"] = random.choice([1, 1, 2, -1]) if inp["has_alt"] else 0
        # placement
        if abs(inp["w1"]) == abs(inp["w2"]) or on_edge(fl, inp["w1"]) or on_edge(fl, inp["w2"]):
            continue
        b1, b2 = band(fl, inp["w1"]), band(fl, inp["w2"])
        if b1 == 0 and b2 == 0:
            continue
        if b1 == b2 and random.random() < 0.6:
            continue
        phi = principal(inp["th"] - fl["M"] * inp["w1"]) if b1 else F(0)
        psi = principal(-fl["M"] * inp["w2"]) if b2 else F(0)
        if phi == 1 or psi == 1:
            continue
        return inp, phi, psi
    return None


def _x_tex(inp):
    trig = r"\cos" if inp["trig"] == "cos" else r"\sin"
    A1 = "" if inp["A1"] == 1 else str(inp["A1"])
    S = A1 + trig + r"\left(" + _arg_tex(inp["w1"], inp["th"]) + r"\right)"
    E = _exp_tex(inp)
    Bt = "" if inp["B"] == 1 else str(inp["B"])
    BE = Bt + E
    parts = [str(inp["A0"])]
    if inp["shape"] in ("sp25", "sp23"):
        parts += [BE, S]
    else:
        parts += [S, BE]
    if inp["has_alt"]:
        parts.append({1: "", -1: "-", 2: "2"}.get(inp["C"], str(inp["C"])) + "(-1)^n")
    out = parts[0]
    for pc in parts[1:]:
        out += (" - " + pc[1:]) if pc.startswith("-") else (" + " + pc)
    return out


def _exp_tex(inp):
    return {"jn": "j^n", "mjn": "(-j)^n"}.get(inp["exp_form"]) or tex_ej(inp["w2"], "n")


# ----------------------------------------------------------------------------- generate
def generate(data):
    shape = random.choices([s for s, _ in SHAPES], weights=[w for _, w in SHAPES])[0]
    while True:
        to_pi = shape in ("fa25", "sp25") and random.random() < 0.35
        fl = _draw_filter(to_pi)
        got = _draw_input(shape, fl)
        if got:
            break
    inp, phi, psi = got
    p, cr = data["params"], data["correct_answers"]
    p["shape"] = shape
    p["source_html"] = SOURCE[shape]
    p["filt"] = {"w1": plain(fl["w1"]), "w2": plain(fl["w2"]), "M": fl["M"], "L1": fl["L1"], "L2": fl["L2"],
                 "to_pi": fl["to_pi"]}
    p["inp"] = {"A0": inp["A0"], "A1": inp["A1"], "trig": inp["trig"], "w1": plain(inp["w1"]), "th": plain(inp["th"]),
                "B": inp["B"], "w2": plain(inp["w2"]), "exp_form": inp["exp_form"], "C": inp["C"],
                "has_alt": inp["has_alt"]}
    p["mag_html"] = mag_svg(fl)
    p["phase_html"] = phase_svg(fl)
    p["x_tex"] = _x_tex(inp)
    p["has_alt"] = inp["has_alt"]
    E_tex = _exp_tex(inp)
    trig = r"\cos" if inp["trig"] == "cos" else r"\sin"
    p["trig_tex"] = trig
    p["trig_word"] = "cosine" if inp["trig"] == "cos" else "sine"
    p["w1_tex"] = tex_pi(inp["w1"])
    p["w2_tex"] = tex_pi(inp["w2"])
    p["E_tex"] = E_tex
    p["arg_form_tex"] = _arg_tex(inp["w1"]) + r" + \phi"
    p["y_form_tex"] = (rf"D + E\,{trig}\left({p['arg_form_tex']}\right) + F\,{E_tex}"
                       + (r" + G\,(-1)^n" if inp["has_alt"] else ""))
    M, L1, L2 = fl["M"], fl["L1"], fl["L2"]
    p["M"] = M
    p["M_tex"] = "" if M == 1 else str(M)
    p["M_one"] = M == 1
    p["L1"], p["L2"] = L1, L2
    p["w1e_tex"], p["w2e_tex"] = tex_pi(fl["w1"]), tex_pi(fl["w2"])
    p["to_pi"] = fl["to_pi"]

    # (a) DC
    D = inp["A0"] * L1
    cr["D"] = float(D)
    p["A0"] = inp["A0"]
    p["D"] = D

    # (b) sinusoid
    g1 = gain(fl, inp["w1"])
    E = inp["A1"] * g1
    cr["E"] = float(E)
    cr["phi"] = float(phi)
    p["b1_html"] = {1: f"inside the inner step ($|\\omega| \\lt {tex_pi(fl['w1'])}$), gain {L1}",
                    2: f"on the outer step, gain {L2}", 0: "outside the passband, gain 0"}[band(fl, inp["w1"])]
    p["s_pass"] = g1 != 0
    p["E"] = E
    p["g1"] = g1
    p["A1"] = inp["A1"]
    p["th_tex"] = tex_pi(inp["th"])
    p["th_zero"] = inp["th"] == 0
    raw = inp["th"] - M * inp["w1"]
    p["phi_raw_tex"] = tex_pi(raw)
    p["phi_tex"] = tex_pi(phi)
    p["phi_wrapped"] = raw != phi and g1 != 0
    p["Mw1_tex"] = tex_pi(M * inp["w1"])
    p["phi_plain"] = plain(phi)

    # (c) exponential
    g2 = gain(fl, inp["w2"])
    cr["Fmag"] = float(inp["B"] * g2)
    cr["Fang"] = float(psi)
    p["b2_html"] = {1: f"inside the inner step, gain {L1}", 2: f"on the outer step, gain {L2}",
                    0: "outside the passband, gain 0"}[band(fl, inp["w2"])]
    p["e_pass"] = g2 != 0
    p["g2"] = g2
    p["B"] = inp["B"]
    p["Fmag"] = inp["B"] * g2
    p["psi_raw_tex"] = tex_pi(-M * inp["w2"])
    p["psi_tex"] = tex_pi(psi)
    p["psi_wrapped"] = (-M * inp["w2"]) != psi and g2 != 0
    p["H2_tex"] = (str(g2) if g2 != 1 else "") + tex_ej(-M * inp["w2"]) if g2 else "0"
    p["exp_html"] = {"jn": r"$j^n = e^{j\frac{\pi}{2}n}$: frequency $\tfrac{\pi}{2}$.",
                     "mjn": r"$(-j)^n = e^{-j\frac{\pi}{2}n}$: frequency $-\tfrac{\pi}{2}$, where the phase is $+M\cdot\tfrac{\pi}{2}$ (odd phase).",
                     "exp": rf"${E_tex}$: frequency ${tex_pi(inp['w2'])}$."}[inp["exp_form"]]

    # (d) (-1)^n
    if inp["has_alt"]:
        G = inp["C"] * (L2 * (-1) ** M if fl["to_pi"] else 0)
        cr["G"] = float(G)
        p["G"] = G
        p["C"] = inp["C"]
        p["Hpi_tex"] = (rf"{L2 if L2 != 1 else ''}e^{{-j{M if M != 1 else ''}\pi}} = {L2 * (-1) ** M}" if fl["to_pi"] else "0")

    # the output, assembled
    terms = [str(D)]
    if E:
        terms.append((str(E) if E != 1 else "") + trig + r"\left(" + _arg_tex(inp["w1"], phi) + r"\right)")
    if g2:
        Bm = inp["B"] * g2
        terms.append((str(Bm) if Bm != 1 else "") + (tex_ej(psi) + r"\," if psi != 0 else "") + E_tex)
    if inp["has_alt"] and p["G"]:
        Gv = p["G"]
        terms.append({1: "", -1: "-"}.get(Gv, str(Gv)) + "(-1)^n")
    y = terms[0]
    for t in terms[1:]:
        y += (" - " + t[1:]) if t.startswith("-") else (" + " + t)
    p["y_tex"] = y
