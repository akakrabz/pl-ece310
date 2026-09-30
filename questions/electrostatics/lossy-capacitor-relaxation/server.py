"""Lossy capacitor (Lecture 10, HW4 #4 family): C, G = (sigma/eps) C, tau = eps/sigma, Q(t) = Q0 exp(-t/tau).

The geometry is drawn at random from
    parallel plates:   C = eps A / d
    coax of length l:  C = 2 pi eps l / ln(b/a)
    spherical shells:  C = 4 pi eps a b / (b - a)
and the point of the question is that tau = C/G = eps/sigma does not care which one it was.
"""
import math
import random

EPS0 = 8.8541878128e-12


def generate(data):
    geometry = random.choice(["plates", "coax", "sphere"])
    er = random.choice([2, 2.2, 3, 4, 5, 6])
    sigma_mant, sigma_exp = random.choice([1, 2, 5]), random.choice([-8, -7, -6])
    sigma = sigma_mant * 10.0 ** sigma_exp
    eps = er * EPS0
    Q0 = random.choice([1, 2, 5, 10])                       # nC

    if geometry == "plates":
        A_cm2 = random.choice([10, 20, 50, 100])
        d_mm = random.choice([0.5, 1, 2])
        C = eps * (A_cm2 * 1e-4) / (d_mm * 1e-3)
        geometry_text = (f"A parallel-plate capacitor has plates of area $A = {A_cm2}\\ \\mathrm{{cm^2}}$ separated by "
                         f"$d = {d_mm}\\ \\mathrm{{mm}}$.")
        C_expl = (r"$C = \dfrac{\epsilon A}{d} = \dfrac{" + f"{er}" + r"\,\epsilon_0\,(" + f"{A_cm2}" + r"\times10^{-4})}{" + f"{d_mm}" + r"\times10^{-3}}$")
        dims = {"A_cm2": A_cm2, "d_mm": d_mm}
    elif geometry == "coax":
        a_mm = random.choice([0.5, 1, 2])
        b_mm = a_mm * random.choice([2, 3, 4, 5])
        L_m = random.choice([0.5, 1, 2])
        C = 2 * math.pi * eps * L_m / math.log(b_mm / a_mm)
        geometry_text = (f"A coaxial cable of length $\\ell = {L_m}\\ \\mathrm{{m}}$ has an inner conductor of radius $a = {a_mm}\\ \\mathrm{{mm}}$ "
                         f"and an outer conductor of inner radius $b = {b_mm}\\ \\mathrm{{mm}}$.")
        C_expl = (r"$C = \dfrac{2\pi\epsilon\,\ell}{\ln(b/a)} = \dfrac{2\pi\,(" + f"{er}" + r"\,\epsilon_0)(" + f"{L_m}" + r")}{\ln " + f"{b_mm/a_mm:g}" + r"}$")
        dims = {"a_mm": a_mm, "b_mm": b_mm, "L_m": L_m}
    else:
        a_cm = random.choice([1, 2, 5, 10])
        b_cm = a_cm * random.choice([2, 3, 4])
        a, b = a_cm * 1e-2, b_cm * 1e-2
        C = 4 * math.pi * eps * a * b / (b - a)
        geometry_text = (f"Two concentric metallic spherical shells have radii $a = {a_cm}\\ \\mathrm{{cm}}$ and $b = {b_cm}\\ \\mathrm{{cm}}$.")
        C_expl = (r"$C = \dfrac{4\pi\epsilon\,ab}{b-a} = \dfrac{4\pi\,(" + f"{er}" + r"\,\epsilon_0)(" + f"{a_cm}" + r"\times10^{-2})(" + f"{b_cm}" + r"\times10^{-2})}{(" + f"{b_cm - a_cm}" + r")\times10^{-2}}$")
        dims = {"a_cm": a_cm, "b_cm": b_cm}

    G = sigma / eps * C
    tau = eps / sigma                                       # s
    t1 = float(f"{tau * random.choice([0.5, 1.5, 2, 3]) * 1e6:.3g}")   # us, rounded as displayed
    Qt = Q0 * math.exp(-t1 * 1e-6 / tau)                   # nC
    I0 = Q0 * 1e-9 / tau                                    # A

    data["params"].update({
        "geometry": geometry, "geometry_text": geometry_text, "C_expl": C_expl, **dims,
        "er": er, "sigma": sigma, "sigma_tex": f"{sigma_mant}\\times10^{{{sigma_exp}}}", "Q0": Q0, "t1": t1,
        "C_str": f"{C*1e12:.4g}", "G_str": f"{G*1e9:.4g}", "tau_str": f"{tau*1e6:.4g}",
        "Qt_str": f"{Qt:.4g}", "I0_str": f"{I0*1e6:.4g}",
    })
    data["correct_answers"]["C"] = C * 1e12       # pF
    data["correct_answers"]["G"] = G * 1e9        # nS
    data["correct_answers"]["tau"] = tau * 1e6    # us
    data["correct_answers"]["Qt"] = Qt            # nC
    data["correct_answers"]["I0"] = I0 * 1e6      # uA
