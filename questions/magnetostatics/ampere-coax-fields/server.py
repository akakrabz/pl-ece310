"""Ampère's law for a coaxial cable: solid core (radius a, +I uniform), hollow return conductor (b < r < c, -I uniform).

    2 pi r H_phi = I_enc(r):
        r < a:      I r^2/a^2                 -> H = I r / (2 pi a^2)
        a < r < b:  I                         -> H = I / (2 pi r)
        b < r < c:  I (c^2 - r^2)/(c^2 - b^2) -> H = I (c^2 - r^2) / (2 pi r (c^2 - b^2))
        r > c:      0
Part (e) flips the return current to +I, so H(r > c) = 2I/(2 pi r).
"""
import math
import random


def generate(data):
    a_mm = random.choice([1, 2])
    b_mm = a_mm * random.choice([2, 3, 4])
    c_mm = b_mm + a_mm * random.choice([1, 2])
    I = random.choice([2, 5, 10, 20])

    a, b, c = a_mm * 1e-3, b_mm * 1e-3, c_mm * 1e-3
    r1 = a / 2
    r2 = (a + b) / 2
    r3 = (b + c) / 2
    r4 = 2 * c

    H1 = I * r1 / (2 * math.pi * a**2)
    H2 = I / (2 * math.pi * r2)
    Ienc3 = I * (c**2 - r3**2) / (c**2 - b**2)
    H3 = Ienc3 / (2 * math.pi * r3)
    H4 = 0.0
    H4same = 2 * I / (2 * math.pi * r4)

    data["params"].update({
        "a_mm": a_mm, "b_mm": b_mm, "c_mm": c_mm, "I": I,
        "r1_mm": f"{r1*1e3:g}", "r2_mm": f"{r2*1e3:g}", "r3_mm": f"{r3*1e3:g}", "r4_mm": f"{r4*1e3:g}",
        "H1_str": f"{H1:.4g}", "H2_str": f"{H2:.4g}", "H3_str": f"{H3:.4g}", "Ienc3_str": f"{Ienc3:.4g}", "H4same_str": f"{H4same:.4g}",
    })
    data["correct_answers"]["H1"] = H1
    data["correct_answers"]["H2"] = H2
    data["correct_answers"]["H3"] = H3
    data["correct_answers"]["H4"] = H4
    data["correct_answers"]["H4same"] = H4same
