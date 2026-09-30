"""Long solenoid (Lecture 15): H = nI, Psi = mu0 H pi a^2 per turn, L = N Psi / I = mu0 N^2 pi a^2 / l,
W = L I^2 / 2, tau = L / R.  The radius is kept small compared with the length so the long-solenoid model applies.
"""
import math
import random

MU0 = 4e-7 * math.pi


def generate(data):
    N = random.choice([100, 200, 500, 1000])
    L_cm = random.choice([5, 10, 20, 25])
    a_mm = random.choice([m for m in [2, 5, 10, 20] if m <= L_cm * 10 / 5])   # a <= l/5
    I = random.choice([0.5, 1, 2, 5])
    R = random.choice([0.5, 1, 2, 5])

    L_m, a = L_cm * 1e-2, a_mm * 1e-3
    A = math.pi * a**2
    H = N / L_m * I
    B = MU0 * H
    Psi = B * A
    Lind = N * Psi / I
    W = 0.5 * Lind * I**2
    tau = Lind / R

    data["params"].update({
        "N": N, "L_cm": L_cm, "L_m": f"{L_m:g}", "a_mm": a_mm, "I": I, "R": R,
        "H_str": f"{H:.4g}", "B_str": f"{B*1e3:.4g}", "A_str": f"{A:.4g}", "Psi_str": f"{Psi*1e6:.4g}",
        "L_str": f"{Lind*1e6:.4g}", "W_str": f"{W*1e6:.4g}", "w_str": f"{0.5*MU0*H**2:.4g}", "vol_str": f"{A*L_m:.4g}",
        "tau_str": f"{tau*1e6:.4g}",
    })
    data["correct_answers"]["H"] = H
    data["correct_answers"]["Psi"] = Psi * 1e6      # uWb
    data["correct_answers"]["Lind"] = Lind * 1e6    # uH
    data["correct_answers"]["W"] = W * 1e6          # uJ
    data["correct_answers"]["tau"] = tau * 1e6      # us
