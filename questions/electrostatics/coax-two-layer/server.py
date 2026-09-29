"""Coaxial cable with two concentric dielectric layers (FA26 Exam 1, problem 4 family).

Inner conductor radius a carries +lambda [C/m]; dielectric eps_r1 for a<r<c, eps_r2 for c<r<b;
outer conductor at r = b carries the return charge -lambda.
    D  = lambda / (2 pi r)                          (same in both layers: free charge only)
    E  = lambda / (2 pi eps0 eps_ri r)              (layer by layer)
    V(a)-V(b) = lambda/(2 pi eps0) [ ln(c/a)/eps_r1 + ln(b/c)/eps_r2 ]   > 0
    C' = 2 pi eps0 / [ ln(c/a)/eps_r1 + ln(b/c)/eps_r2 ]                  (series layers)
"""
import math
import random

EPS0 = 8.8541878128e-12


def generate(data):
    a_mm = random.choice([0.5, 1.0, 1.5, 2.0])
    c_mm = a_mm * random.choice([2, 2.5, 3])
    b_mm = c_mm * random.choice([1.5, 2, 3])
    er1, er2 = random.sample([1, 2, 2.2, 3, 4, 5, 6], 2)
    lam_nC = random.choice([1, 2, 4, 5, 8, 10])

    a, c, b = a_mm * 1e-3, c_mm * 1e-3, b_mm * 1e-3
    lam = lam_nC * 1e-9
    r1 = 0.5 * (a + c)          # a point inside layer 1
    r2 = 0.5 * (c + b)          # a point inside layer 2

    D_r1 = lam / (2 * math.pi * r1)                       # C/m^2
    E_r2 = lam / (2 * math.pi * EPS0 * er2 * r2)          # V/m
    Vab = lam / (2 * math.pi * EPS0) * (math.log(c / a) / er1 + math.log(b / c) / er2)   # V, = V(a)-V(b)
    Cp = lam / Vab                                        # F/m

    data["params"].update({
        "a_mm": a_mm, "c_mm": c_mm, "b_mm": b_mm, "er1": er1, "er2": er2, "lam_nC": lam_nC,
        "r1_mm": f"{r1*1e3:.3g}", "r2_mm": f"{r2*1e3:.3g}",
        "D_r1_str": f"{D_r1*1e9:.4g}", "E_r2_str": f"{E_r2:.4g}", "Vab_str": f"{Vab:.4g}", "Cp_str": f"{Cp*1e12:.4g}",
        "lnca": f"{math.log(c/a):.4g}", "lnbc": f"{math.log(b/c):.4g}",
    })
    data["correct_answers"]["D_r1"] = D_r1 * 1e9      # nC/m^2
    data["correct_answers"]["E_r2"] = E_r2            # V/m
    data["correct_answers"]["q_outer"] = -float(lam_nC)  # nC/m
    data["correct_answers"]["Vab"] = Vab              # V
    data["correct_answers"]["Cp"] = Cp * 1e12         # pF/m
