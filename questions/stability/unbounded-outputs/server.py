"""Causal, marginally stable H(z) (simple poles on |z| = 1, optional pole inside, optional zero):
tick every bounded input that produces an unbounded output (Lecture 11 slide 18; HW4 #3c,d;
SP2025 #6, FA2025 #8b, SP2021 #5b). Translated from the "poles" drill of the study site.

Rule: Y = H X; after cancelling common factors, y is unbounded iff Y has a repeated pole on the
unit circle (no pole of H lies outside it). Built-in traps: the angle 2/3 rad versus 2pi/3, the
conjugate of a single complex pole, (1 - z^-2) with a zero of H on +-1 that cancels one of its
poles, and a repeated pole inside the circle."""

import cmath
import math
import random
from fractions import Fraction as F

from ece310 import fmt

TOL = 1e-9

# unit-circle pole groups: angles ("pi", f) = pi*f or ("rad", f) = f radians, and the LaTeX factor
GROUPS = {
    "1": ([("pi", F(0))], r"\left(1 - z^{-1}\right)"),
    "-1": ([("pi", F(1))], r"\left(1 + z^{-1}\right)"),
    "j": ([("pi", F(1, 2))], r"\left(1 - j\,z^{-1}\right)"),
    "±j": ([("pi", F(1, 2)), ("pi", F(-1, 2))], r"\left(1 + z^{-2}\right)"),
    "pi4": ([("pi", F(1, 4))], r"\left(1 - e^{j\pi/4}z^{-1}\right)"),
    "±pi4": ([("pi", F(1, 4)), ("pi", F(-1, 4))], r"\left(1 - \sqrt{2}\,z^{-1} + z^{-2}\right)"),
    "pi3": ([("pi", F(1, 3))], r"\left(1 - e^{j\pi/3}z^{-1}\right)"),
    "±pi3": ([("pi", F(1, 3)), ("pi", F(-1, 3))], r"\left(1 - z^{-1} + z^{-2}\right)"),
    "2pi3": ([("pi", F(2, 3))], r"\left(1 - e^{j2\pi/3}z^{-1}\right)"),
    "±2pi3": ([("pi", F(2, 3)), ("pi", F(-2, 3))], r"\left(1 + z^{-1} + z^{-2}\right)"),
    "rad23": ([("rad", F(2, 3))], r"\left(1 - e^{j2/3}z^{-1}\right)"),
    "±1": ([("pi", F(0)), ("pi", F(1))], r"\left(1 - z^{-2}\right)"),
}
TRAP_GROUPS = ["2pi3", "±2pi3", "rad23"]
OTHER_GROUPS = [k for k in GROUPS if k not in TRAP_GROUPS] + ["±1", "±1"]   # extra weight: the cancellation case


def aval(a):
    return math.pi * float(a[1]) if a[0] == "pi" else float(a[1])


def apt(a):
    return cmath.exp(1j * aval(a))


def neg(a):
    return (a[0], -a[1])


def same(z, w):
    return abs(complex(z) - complex(w)) < TOL


def exp_tex(a, with_n):
    """exponent j\\pi n/4, -j2\\pi n/3, j2n/3 (without n: j\\pi/4 ...)."""
    f = a[1]
    s = "-" if f < 0 else ""
    n, d = abs(f).numerator, abs(f).denominator
    core = ("" if (n == 1 and a[0] == "pi") else str(n)) + (r"\pi" if a[0] == "pi" else "")
    if with_n:
        core += (" n" if a[0] == "pi" else "n")
    return f"{s}j{core}" + (f"/{d}" if d != 1 else "")


def arg_tex(a):
    """argument of cos/sin: \\tfrac{\\pi}{4}n, \\tfrac{2}{3}n."""
    f = abs(a[1])
    n, d = f.numerator, f.denominator
    if a[0] == "pi":
        num = (r"\pi" if n == 1 else f"{n}\\pi")
        return (r"\tfrac{" + num + "}{" + str(d) + "}n") if d != 1 else num + " n"
    return r"\tfrac{" + str(n) + "}{" + str(d) + "}n"


def pt_name(z, a=None):
    """LaTeX name of a pole/zero location."""
    if a is not None:
        if a[0] == "pi" and a[1] == 0:
            return "1"
        if a[0] == "pi" and abs(a[1]) == 1:
            return "-1"
        if a[0] == "pi" and a[1] == F(1, 2):
            return "j"
        if a[0] == "pi" and a[1] == F(-1, 2):
            return "-j"
        return "e^{" + exp_tex(a, False) + "}"
    return fmt.tex_num(z)


# ----------------------------------------------------------------------------- inputs
def in_exp(a, k=0):
    f = a[1]
    if a[0] == "pi" and f == 0:
        tex = "u[n]"
    elif a[0] == "pi" and abs(f) == 1:
        tex, a = "(-1)^{n}u[n]", ("pi", F(1))
    elif a[0] == "pi" and f == F(1, 2):
        tex = "j^{n}u[n]"
    elif a[0] == "pi" and f == F(-1, 2):
        tex = "(-j)^{n}u[n]"
    else:
        tex = "e^{" + exp_tex(a, True) + "}u[n]"
    if k:
        tex = tex[: -len("u[n]")] + f"u[n-{k}]"
    return {"tex": tex, "poles": [(apt(a), pt_name(None, a))], "zeros": [],
            "spec": {"kind": "exp", "ang": [a[0], a[1].numerator, a[1].denominator], "k": k},
            "key": f"exp{aval(a) % (2 * math.pi):.9f}k{k}", "tier": 2 if k else 0}


def in_combo(a, r):
    """p^n u[n] - r^n u[n] with p = +-1 on the unit circle (HW4 #5 h1 shape): X = (p - r) z^-1 / ((1 - p z^-1)(1 - r z^-1)),
    no finite non-zero zero, poles p and r."""
    base = "u[n]" if a[1] == 0 else "(-1)^{n}u[n]"
    return {"tex": base + " - " + fmt.tex_pow(r) + "u[n]", "poles": [(apt(a), pt_name(None, a)), (complex(r), fmt.tex_num(r))],
            "zeros": [], "spec": {"kind": "combo", "ang": [a[0], a[1].numerator, a[1].denominator], "r": str(r)},
            "key": f"combo{a[1]}|{r}", "tier": 1}


def in_cospi():
    """cos(pi n) u[n] = (-1)^n u[n]: one pole at -1 (the two 'poles' e^{+-j pi} coincide and the zero cancels one)."""
    return {"tex": r"\cos(\pi n)\,u[n]", "poles": [(complex(-1), "-1")], "zeros": [],
            "spec": {"kind": "cos", "ang": ["pi", 1, 1]}, "key": "cospi", "tier": 1}


def in_trig(a, fn):
    a = (a[0], abs(a[1]))
    th = aval(a)
    zeros = [(complex(math.cos(th)), "")] if fn == "cos" and abs(math.cos(th)) > 1e-12 else []
    return {"tex": rf"\{fn}\!\left({arg_tex(a)}\right)u[n]", "poles": [(apt(a), pt_name(None, a)), (apt(neg(a)), pt_name(None, neg(a)))],
            "zeros": zeros, "spec": {"kind": fn, "ang": [a[0], a[1].numerator, a[1].denominator]}, "key": f"{fn}{th:.9f}"}


def in_geo(r):
    return {"tex": fmt.tex_pow(r) + "u[n]", "poles": [(complex(r), fmt.tex_num(r))], "zeros": [],
            "spec": {"kind": "geo", "a": str(r)}, "key": f"geo{r}"}


def in_nexp(r):
    return {"tex": "n" + fmt.tex_pow(r) + "u[n]", "poles": [(complex(r), fmt.tex_num(r))] * 2, "zeros": [],
            "spec": {"kind": "nexp", "a": str(r)}, "key": f"nexp{r}"}


def in_delta(k):
    return {"tex": rf"\delta[n-{k}]", "poles": [], "zeros": [], "spec": {"kind": "fir", "c": [0] * k + [1]}, "key": f"d{k}"}


def in_fir(p):
    return {"tex": r"\delta[n] - \delta[n-1]" if p == 1 else r"\delta[n] + \delta[n-1]", "poles": [],
            "zeros": [(complex(p), fmt.tex_num(p))], "spec": {"kind": "fir", "c": [1, -p]}, "key": f"fir{p}"}


# ----------------------------------------------------------------------------- analysis
def analyse(Hpoles, Hzeros, x):
    """Hpoles/Hzeros: surviving lists of (complex, name). Returns (unbounded, notes)."""
    poles = [(z, nm, "H") for z, nm in Hpoles] + [(z, nm, "X") for z, nm in x["poles"]]
    zeros = [(z, nm, "H") for z, nm in Hzeros] + [(z, nm, "X") for z, nm in x["zeros"]]
    notes = []
    for z, nm, src in zeros:
        for i, (p, pn, psrc) in enumerate(poles):
            if same(p, z) and psrc != src:
                notes.append(("cancel", src, psrc, pn))
                poles.pop(i)
                break
    uc = [(p, pn, s) for p, pn, s in poles if abs(abs(p) - 1) < TOL]
    doubles = []
    for i, (p, pn, s) in enumerate(uc):
        if any(j != i and same(p, q) for j, (q, _, _) in enumerate(uc)) and not any(same(p, d) for d, _ in doubles):
            doubles.append((p, pn))
    outside = [p for p, _, _ in poles if abs(p) > 1 + TOL]
    return bool(doubles or outside), doubles, notes


def generate(data):
    for _ in range(10000):
        gkey = random.choice(TRAP_GROUPS) if random.random() < 0.45 else random.choice(OTHER_GROUPS)
        angs, den_tex = GROUPS[gkey]
        inside = random.choice([F(1, 2), F(-1, 2), F(1, 3), F(2, 3), F(3, 4), F(-1, 3), F(-3, 4)]) if random.random() < 0.55 else None
        q = None
        if gkey == "±1" and random.random() < 0.75:
            q = random.choice([F(1), F(-1)])
        elif random.random() < 0.3:
            q = random.choice([F(1, 2), F(-1, 2), F(3, 4), F(1), F(-1)])
        Hp = [(apt(a), pt_name(None, a), a) for a in angs] + ([(complex(inside), fmt.tex_num(inside), None)] if inside is not None else [])
        if q is not None and inside is not None and q == inside:
            continue
        cancelled = None
        Hz = []
        if q is not None:
            hit = [i for i, (z, _, _) in enumerate(Hp) if same(z, q)]
            if hit:
                if gkey != "±1":
                    continue
                cancelled = Hp.pop(hit[0])
            else:
                Hz = [(complex(q), fmt.tex_num(q))]
        uc_left = [a for z, _, a in Hp if a is not None]
        Hpoles = [(z, nm) for z, nm, _ in Hp]

        match, miss, cancel, generic = [], [], [], []
        for a in uc_left:
            match.append(in_exp(a))
            match += [in_exp(a, 1), in_exp(a, 2)]
            if a[0] == "pi" and a[1] in (0, 1):
                match.append(in_combo(a, random.choice([F(1, 2), F(-1, 2), F(1, 3), F(2, 3)])))
                if a[1] == 1:
                    match.append(in_cospi())
                continue
            match += [in_trig(a, "cos"), in_trig(a, "sin")]
            if not any(same(apt(neg(a)), z) for z, _ in Hpoles):
                miss.append(in_exp(neg(a)))
        trap = []
        if gkey in ("2pi3", "±2pi3"):
            trap = [in_exp(("rad", F(2, 3))), in_trig(("rad", F(2, 3)), "cos")]
        if gkey == "rad23":
            trap = [in_trig(("pi", F(2, 3)), "cos"), in_exp(("pi", F(2, 3))), in_trig(("pi", F(2, 3)), "sin")]
        for o in [F(1, 4), F(3, 4), F(1, 3), F(2, 3), F(1, 2), F(1, 6)]:
            a = ("pi", o)
            if not any(same(apt(a), z) or same(apt(neg(a)), z) for z, _ in Hpoles):
                miss.append(in_trig(a, "cos") if random.random() < 0.5 else in_exp(a))
        if cancelled is not None:
            cancel.append(in_exp(cancelled[2]))
        for a in uc_left:
            if a[0] == "pi" and a[1] in (0, 1):
                cancel.append(in_fir(1 if a[1] == 0 else -1))
        if inside is not None:
            cancel.append(in_geo(inside) if random.random() < 0.6 else in_nexp(inside))
        generic += [in_delta(random.randint(1, 3)), in_geo(random.choice([F(1, 2), F(-1, 2), F(2, 3), F(3, 4)])),
                    in_nexp(F(1, 2)), in_exp(("pi", F(0))), in_exp(("pi", F(1))), in_exp(("pi", F(1, 2)))]

        for x in match + trap + miss + cancel + generic:
            unb, doubles, notes = analyse(Hpoles, Hz, x)
            x["unbounded"], x["doubles"], x["notes"] = unb, doubles, notes

        chosen, keys = [], set()

        def take(lst, k, want):
            for c in lst:
                if k <= 0:
                    break
                if c["key"] not in keys and c["unbounded"] == want:
                    keys.add(c["key"])
                    chosen.append(c)
                    k -= 1
            return k

        def shuffled(lst):
            lst = list(lst)
            random.shuffle(lst)
            return lst

        total = random.randint(5, 6)
        target = random.choice([1, 2, 3])
        # resonant inputs: plain exponentials / sinusoids first, then combinations, then shifted copies
        unb_pool = sorted((c for c in match + generic if c["unbounded"]),
                          key=lambda c: c.get("tier", 0) * 0.45 + random.random())
        if take(unb_pool, target, True) > 0:
            continue
        nb = total - target
        nb = take(shuffled(trap), min(1, nb), False) + nb - min(1, nb)
        nb = take(shuffled(miss), min(1, nb), False) + nb - min(1, nb)
        if cancel and random.random() < 0.8:
            nb = take(shuffled(cancel), min(1, nb), False) + nb - min(1, nb)
        nb = take(shuffled(generic + miss), nb, False)
        nb = take(shuffled(cancel + trap), nb, False)
        if nb > 0:
            continue
        n_unb = sum(x["unbounded"] for x in chosen)
        if n_unb != target or len(chosen) != total:
            continue
        break
    else:
        raise RuntimeError("no variant found")
    random.shuffle(chosen)

    # ---- H(z) LaTeX
    num_tex = "1" if q is None else fmt.tex_sum([(1, ""), (-q, "z^{-1}")])
    den = den_tex + (fmt.tex_factor(inside) if inside is not None else "")
    H_tex = fmt.tex_frac(num_tex, den)

    uc_names = r",\ ".join(pt_name(None, a) for a in angs)
    pole_desc = f"on the unit circle: ${uc_names}$"
    if inside is not None:
        pole_desc += f"; inside: ${fmt.tex_num(inside)}$"
    cancel_html = ""
    if cancelled is not None:
        cancel_html = (f" <b>Cancellation first:</b> $1 - z^{{-2}} = (1 - z^{{-1}})(1 + z^{{-1}})$, so the zero of $H$ at "
                       f"$z = {fmt.tex_num(q)}$ cancels the pole there. The only unit-circle pole that is left is "
                       f"$z = {pt_name(None, uc_left[0])}$.")
    elif q is not None:
        cancel_html = f" The zero of $H$ at $z = {fmt.tex_num(q)}$ cancels no pole of $H$."

    def reason(x):
        parts = []
        if x["poles"]:
            pn = sorted({nm for _, nm in x["poles"]}, key=str)
            parts.append("$X(z)$ has pole" + ("s" if len(pn) > 1 else "") + " at $" + r",\ ".join(pn) + "$.")
        else:
            parts.append("$X(z)$ is a polynomial in $z^{-1}$ (no poles).")
        for kind, zsrc, psrc, pn in x["notes"]:
            who_z = "The input's" if zsrc == "X" else "The zero of $H$"
            who_p = "the pole of $H$" if psrc == "H" else "the input's pole"
            parts.append(f"{who_z} {'zero ' if zsrc == 'X' else ''}cancels {who_p} at ${pn}$.")
        if x["doubles"]:
            parts.append("It lands on the pole of $H$ at $" + r",\ ".join(nm for _, nm in x["doubles"])
                         + "$: $Y(z)$ has a <b>double pole on the unit circle</b>, so $y[n]$ grows like $n$.")
        elif inside is not None and sum(same(z, inside) for z, _ in x["poles"]) + 1 >= 2 and \
                not any(k == "cancel" and pn == fmt.tex_num(inside) for k, _, _, pn in x["notes"]):
            parts.append(f"$Y(z)$ has a repeated pole at ${fmt.tex_num(inside)}$, but it is <i>inside</i> the unit circle "
                         f"(a term like $n\\left({fmt.tex_num(inside)}\\right)^n u[n]$ decays): bounded.")
        elif any(abs(abs(z) - 1) < TOL for z, _ in x["poles"]):
            parts.append("No input pole coincides with a surviving pole of $H$ on the unit circle: only simple poles there, bounded.")
        else:
            parts.append("Every pole of $Y(z)$ on the unit circle is simple: bounded.")
        if x["spec"]["kind"] in ("exp", "cos", "sin") and x["spec"]["ang"][0] == "rad" and gkey in ("2pi3", "±2pi3"):
            parts.append(r"Trap: the angle here is $\tfrac{2}{3}$ rad $\approx 38^\circ$, not $\tfrac{2\pi}{3} = 120^\circ$.")
        if gkey == "rad23" and x["spec"]["kind"] in ("exp", "cos", "sin") and x["spec"]["ang"][0] == "pi" and \
                abs(abs(F(x["spec"]["ang"][1], x["spec"]["ang"][2])) - F(2, 3)) == 0:
            parts.append(r"Trap: the pole of $H$ is $e^{j2/3}$ (angle $\tfrac{2}{3}$ rad $\approx 38^\circ$), not $e^{j2\pi/3}$.")
        return " ".join(parts)

    rows = "".join(f"<tr><td>${x['tex']}$</td><td>{'<b>unbounded</b>' if x['unbounded'] else 'bounded'}</td>"
                   f"<td class=\"text-start\">{reason(x)}</td></tr>" for x in chosen)
    table_html = ('<table class="table table-sm table-bordered text-center" style="width:auto">'
                  "<tr><th>input $x[n]$</th><th>output</th><th>why</th></tr>" + rows + "</table>")

    p = data["params"]
    p["H_tex"] = H_tex
    p["pole_desc_html"] = pole_desc
    p["cancel_html"] = cancel_html
    p["table_html"] = table_html
    p["group"] = gkey
    p["inside"] = None if inside is None else str(inside)
    p["zero"] = None if q is None else str(q)
    p["inputs"] = [{"tex": x["tex"], "spec": x["spec"], "unbounded": x["unbounded"]} for x in chosen]
    p["x_choices"] = [{"text": f"$x[n] = {x['tex']}$", "correct": x["unbounded"]} for x in chosen]
