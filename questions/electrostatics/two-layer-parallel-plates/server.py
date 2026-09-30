"""Two dielectric slabs between parallel plates at V = 0 (z = 0) and V = Vp (z = z0).  HW4 #1/#3 family.

eps_1 for 0 < z < d, eps_2 for d < z < z0; no free charge on the interface, so D_z is the same in both slabs:
    eps_r1 E1z = eps_r2 E2z
    V(d) = -E1z d,   Vp = V(d) - E2z (z0 - d)
    rho_s(top)    = -D2z = -eps_r2 eps0 E2z   (n = -z out of the top conductor)
    rho_s(bottom) = +D1z
    C/A = eps0 / (d/eps_r1 + (z0-d)/eps_r2)   (two slabs in series)
E1z is given (negative: the field points from the top plate at Vp > 0 toward ground) and chosen so E2z is an integer.
"""
import random
from math import gcd

EPS0 = 8.8541878128e-12


def generate(data):
    er1, er2 = random.sample([1, 2, 3, 4, 5, 6, 8], 2)
    d = random.choice([1, 2, 3])            # m, lower slab thickness
    t2 = random.choice([1, 2, 3])           # m, upper slab thickness
    z0 = d + t2
    k = random.choice([1, 2, 3, 5])
    E1z = -k * er2 // gcd(er1, er2)         # integer; er1*E1z/er2 is an integer too

    E2z = er1 * E1z // er2
    Vd = -E1z * d
    Vp = Vd - E2z * t2
    rho_top = -er2 * E2z                    # coefficient of eps0
    rho_bot = er1 * E1z
    C_per_A = EPS0 / (d / er1 + t2 / er2)   # F/m^2

    data["params"].update({
        "er1": er1, "er2": er2, "d": d, "t2": t2, "z0": z0, "E1z": E1z, "E2z": E2z,
        "Vd": Vd, "Vp": Vp, "rho_top": rho_top, "rho_bot": rho_bot,
        "C_per_A_str": f"{C_per_A*1e12:.4g}",
    })
    data["correct_answers"]["E2z"] = float(E2z)
    data["correct_answers"]["Vd"] = float(Vd)
    data["correct_answers"]["Vp"] = float(Vp)
    data["correct_answers"]["rho_top"] = float(rho_top)
    data["correct_answers"]["C_per_A"] = C_per_A * 1e12    # pF/m^2
