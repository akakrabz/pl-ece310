"""Flux of D through the whole z = 0 plane from two point charges (Lecture 3 challenge question).

A point charge sends half its flux through any infinite plane that does not contain it.
With the plane's normal n = +z:  a charge q at z = +h contributes -q/2 (its lines cross downward),
a charge q at z = -h contributes +q/2 (its lines cross upward).  psi = (q_below - q_above)/2.
Choices are rendered as multiples of Q with exact fractions.
"""
import random
from fractions import Fraction


def _fmt(fr: Fraction) -> str:
    """Format a Fraction as a LaTeX multiple of Q."""
    if fr == 0:
        return "$0$"
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    if fr.denominator == 1:
        mag = "" if fr.numerator == 1 else str(fr.numerator)
        return f"${sign}{mag}Q$"
    num = "" if fr.numerator == 1 else str(fr.numerator)
    return f"${sign}\\dfrac{{{num}Q}}{{{fr.denominator}}}$"


def generate(data):
    q_above = random.choice([-3, -2, -1, 1, 2, 3])
    q_below = random.choice([-3, -2, -1, 1, 2, 3])
    normal_up = random.choice([True, True, False])       # which way the plane's normal points

    psi = Fraction(q_below - q_above, 2)
    if not normal_up:
        psi = -psi

    def qtex(q):
        return ("+" if q > 0 else "-") + ("" if abs(q) == 1 else str(abs(q))) + "Q"

    data["params"].update({
        "q_above_tex": qtex(q_above), "q_below_tex": qtex(q_below),
        "normal_tex": "+\\hat{z}" if normal_up else "-\\hat{z}",
        "q_above": q_above, "q_below": q_below,
        "psi_tex": _fmt(psi).strip("$"),
        "contrib_above": _fmt(Fraction(-q_above, 2) * (1 if normal_up else -1)).strip("$"),
        "contrib_below": _fmt(Fraction(q_below, 2) * (1 if normal_up else -1)).strip("$"),
    })

    # distractors built from the standard mistakes, de-duplicated against the correct value
    candidates = [
        Fraction(q_above + q_below, 2),            # "both contribute +q/2"
        -Fraction(q_above + q_below, 2),           # ... with the sign flipped
        Fraction(q_above + q_below),               # closed-surface answer (enclosed charge)
        -psi,                                      # wrong normal
        Fraction(q_above - q_below, 2),            # signs of above/below swapped
        Fraction(0),
        Fraction(-q_above, 2), Fraction(q_below, 2),   # only one charge counted
        Fraction(q_above), Fraction(-q_below),
        Fraction(q_above - q_below), Fraction(q_below - q_above),
    ]
    distractors = []
    for c in candidates:
        if c != psi and c not in distractors:
            distractors.append(c)
        if len(distractors) == 4:
            break
    choices = [{"text": _fmt(psi), "correct": "true"}] + [{"text": _fmt(d), "correct": "false"} for d in distractors]
    random.shuffle(choices)
    data["params"]["choices"] = choices
