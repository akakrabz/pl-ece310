"""Uniform free charge sheet rho_s = s*eps0 at z = d between grounded plates at z = 0 and z = z0.  HW4 #2 family.

eps_r1 below the sheet, eps_r2 above.  With V(0) = V(z0) = 0 and V(d) = V0 the potential is linear in each region:
    E1z = -V0/d  (0 < z < d),     E2z = +V0/(z0-d)  (d < z < z0)
    jump of D at the sheet:   s*eps0 = eps2 E2z - eps1 E1z   =>   V0 = s / (eps_r1/d + eps_r2/(z0-d))
    plates:  rho_s(0) = +D1z = eps_r1 E1z eps0,   rho_s(z0) = -D2z = -eps_r2 E2z eps0;   they sum to -s*eps0.
"""
import random


def generate(data):
    d = random.choice([1, 2, 3])
    t2 = random.choice([1, 2, 3, 4])
    z0 = d + t2
    er1, er2 = random.sample([1, 2, 3, 4, 6], 2)
    s = random.choice([-12, -8, -6, -4, -3, 3, 4, 6, 8, 12])

    V0 = s / (er1 / d + er2 / t2)
    E1z = -V0 / d
    E2z = V0 / t2
    rho_0 = er1 * E1z            # coefficient of eps0
    rho_z0 = -er2 * E2z

    data["params"].update({
        "d": d, "t2": t2, "z0": z0, "er1": er1, "er2": er2, "s": s, "neg_s": -s,
        "V0_str": f"{V0:.4g}", "E1z_str": f"{E1z:.4g}", "E2z_str": f"{E2z:.4g}",
        "rho0_str": f"{rho_0:.4g}", "rhoz0_str": f"{rho_z0:.4g}",
    })
    data["correct_answers"]["V0"] = V0
    data["correct_answers"]["E1z"] = E1z
    data["correct_answers"]["E2z"] = E2z
    data["correct_answers"]["rho_0"] = rho_0
    data["correct_answers"]["rho_z0"] = rho_z0
