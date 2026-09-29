"""Gauss's law for an infinite slab with a power-law charge profile.

Family of FA26 Exam 1, problem 2:  rho(x) = rho0 * |x/a|^k  for |x| < a, else 0.
Closed forms (derived and checked numerically in the notes):
    charge per unit area           = 2 rho0 a / (k+1)
    E_x inside  (0 < x < a)         = rho0 x^(k+1) / ((k+1) eps0 a^k)
    E_x outside (x > a)             = rho0 a / ((k+1) eps0)
    V(b) with V(0)=0, b >= a        = rho0 a^2 / ((k+2) eps0) - rho0 a b / ((k+1) eps0)
"""
import random

EPS0 = 8.8541878128e-12  # F/m


def generate(data):
    k = random.choice([1, 2, 3])
    rho0_nC = random.choice([2, 4, 5, 8, 10, 16, 20])      # nC/m^3
    a_mm = random.choice([1, 2, 4, 5, 8, 10])               # mm
    b_over_a = random.choice([2, 3, 4])
    xin_frac = random.choice([(1, 2), (1, 4), (3, 4)])      # inside point x = a * num/den

    rho0 = rho0_nC * 1e-9
    a = a_mm * 1e-3
    b = b_over_a * a
    xin = a * xin_frac[0] / xin_frac[1]

    sigma = 2 * rho0 * a / (k + 1)                                   # C/m^2
    E_in = rho0 * xin ** (k + 1) / ((k + 1) * EPS0 * a ** k)          # V/m at x = xin
    E_out = rho0 * a / ((k + 1) * EPS0)                              # V/m for x > a
    V_b = rho0 * a ** 2 / ((k + 2) * EPS0) - rho0 * a * b / ((k + 1) * EPS0)   # V

    data["params"].update({
        "k": k, "rho0_nC": rho0_nC, "a_mm": a_mm, "b_over_a": b_over_a,
        "xin_num": xin_frac[0], "xin_den": xin_frac[1],
        "k1": k + 1, "k2": k + 2,
        "sigma_nC": f"{sigma * 1e9:.4g}", "E_in_str": f"{E_in:.4g}", "E_out_str": f"{E_out:.4g}", "V_b_str": f"{V_b:.4g}",
    })
    data["correct_answers"]["sigma"] = sigma * 1e9      # answer requested in nC/m^2
    data["correct_answers"]["E_in"] = E_in
    data["correct_answers"]["E_out"] = E_out
    data["correct_answers"]["V_b"] = V_b
