"""Two infinite current sheets Js1 z at x = 0 and Js2 z at x = d (Lecture 13, Example 1 + superposition).

One sheet:  H = (1/2) Js x n  with n the unit normal from the sheet to the field point, i.e. H_y = (Js/2) sgn(x - x_sheet).
    x < 0:      H_y = -(Js1 + Js2)/2
    0 < x < d:  H_y =  (Js1 - Js2)/2
    x > d:      H_y =  (Js1 + Js2)/2
Force per unit area on sheet 2 from sheet 1:  f = Js2 z x B1(d) = Js2 z x (mu0 Js1/2) y = -(mu0 Js1 Js2 / 2) x
    -> parallel currents (same sign) attract, antiparallel repel.
"""
import math
import random

MU0 = 4e-7 * math.pi


def generate(data):
    Js1 = random.choice([-8, -6, -4, -3, 3, 4, 6, 8, 10])
    Js2 = random.choice([-8, -6, -4, -3, -2, 2, 3, 4, 6, 8])
    d = random.choice([1, 2, 3])

    Hleft = -(Js1 + Js2) / 2
    Hmid = (Js1 - Js2) / 2
    Hright = (Js1 + Js2) / 2
    fx = -MU0 * Js1 * Js2 / 2                      # N/m^2

    if Js1 * Js2 > 0:
        force_words = "negative, i.e. toward sheet 1 — parallel currents attract."
    else:
        force_words = "positive, i.e. away from sheet 1 — antiparallel currents repel."
    dir_correct = r"$+\hat{y}$" if Js1 > 0 else r"$-\hat{y}$"
    dir_wrong = r"$-\hat{y}$" if Js1 > 0 else r"$+\hat{y}$"

    data["params"].update({
        "Js1": Js1, "Js2": Js2, "d": d,
        "Hleft_str": f"{Hleft:g}", "Hmid_str": f"{Hmid:g}", "Hright_str": f"{Hright:g}", "fx_str": f"{fx*1e6:.4g}",
        "force_words": force_words, "dir_correct": dir_correct, "dir_wrong": dir_wrong,
    })
    data["correct_answers"]["Hleft"] = float(Hleft)
    data["correct_answers"]["Hmid"] = float(Hmid)
    data["correct_answers"]["Hright"] = float(Hright)
    data["correct_answers"]["fx"] = fx * 1e6       # uN/m^2
