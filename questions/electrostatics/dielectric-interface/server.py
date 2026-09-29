"""Boundary conditions at a charge-free interface between two perfect dielectrics.

Family of FA26 Exam 1, problem 3. Interface z = 0; medium 1 (eps_r1) above, medium 2 (eps_r2) below.
Given E1 just above, find E2, D on both sides, P on both sides, and the bound surface charge.
    tangential E continuous:  E2x = E1x, E2y = E1y
    normal D continuous:      eps_r1 E1z = eps_r2 E2z   (no free surface charge)
    P = (eps_r - 1) eps0 E ;  rho_sb = P2z - P1z  (normal from 2 into 1 is +z)
E1z is chosen so that E2z is an integer.
"""
import random
from math import gcd


def generate(data):
    er1, er2 = random.sample([2, 3, 4, 5, 6, 8], 2)
    g = gcd(er1, er2)
    j = random.choice([1, 2, 3])
    E1x = random.choice([-6, -4, -3, -2, 2, 3, 4, 6])
    E1y = random.choice([0, 0, -3, 3, 5])
    E1z = -j * er2 // g * random.choice([1, -1])          # integer, and er1*E1z/er2 is an integer

    E2x, E2y = E1x, E1y
    E2z = er1 * E1z // er2

    D1z_over_eps0 = er1 * E1z
    D2x_over_eps0 = er2 * E2x
    P1x_over_eps0 = (er1 - 1) * E1x
    P2z_over_eps0 = (er2 - 1) * E2z
    rho_sb_over_eps0 = (er2 - 1) * E2z - (er1 - 1) * E1z       # = E1z - E2z

    def sgn(v):
        return f"+ {v}" if v >= 0 else f"- {abs(v)}"

    E1_tex = f"{E1x}\\,\\hat{{x}} " + (f"{sgn(E1y)}\\,\\hat{{y}} " if E1y != 0 else "") + f"{sgn(E1z)}\\,\\hat{{z}}"
    data["params"].update({
        "er1": er1, "er2": er2, "E1_tex": E1_tex, "E1x": E1x, "E1y": E1y, "E1z": E1z,
        "E2z": E2z, "D1z": D1z_over_eps0, "D2x": D2x_over_eps0, "P1x": P1x_over_eps0, "P2z": P2z_over_eps0,
        "rho_sb": rho_sb_over_eps0, "chi1": er1 - 1, "chi2": er2 - 1,
    })
    data["correct_answers"]["E2x"] = float(E2x)
    data["correct_answers"]["E2z"] = float(E2z)
    data["correct_answers"]["D2x"] = float(D2x_over_eps0)
    data["correct_answers"]["D1z"] = float(D1z_over_eps0)
    data["correct_answers"]["P1x"] = float(P1x_over_eps0)
    data["correct_answers"]["P2z"] = float(P2z_over_eps0)
    data["correct_answers"]["rho_sb"] = float(rho_sb_over_eps0)
