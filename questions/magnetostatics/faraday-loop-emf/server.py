"""Faraday's law for a stationary N-turn rectangular coil in B(t) = (B0 + k t) z  (Lecture 14, Examples 1 and 6).

Circulation C counterclockwise seen from +z  <->  dS = +z dS.
    Psi(t) = (B0 + k t) a b          (per turn)
    emf    = -N dPsi/dt = -N a b k   (constant)
    I      = emf / R                 (positive = counterclockwise from above)
Lenz: a counterclockwise current makes +z field inside the loop; it appears when the +z flux is decreasing (k < 0).
"""
import random


def generate(data):
    N = random.choice([1, 5, 10, 20, 50])
    a_cm, b_cm = random.choice([5, 10, 20]), random.choice([5, 10, 20])
    R = random.choice([1, 2, 4, 5, 10])
    B0 = random.choice([0.2, 0.5, 1, 2])
    k = random.choice([-4, -2, -1, -0.5, 0.5, 1, 2, 4])          # T/s, dB_z/dt
    t1 = random.choice([0, 0.1, 0.2, 0.5])

    A = a_cm * b_cm * 1e-4
    Bt1 = B0 + k * t1
    Psi = Bt1 * A
    emf = -N * A * k
    I = emf / R
    ccw = emf > 0

    sense_choices = [
        {"text": "counterclockwise (in the direction of $C$)", "correct": "true" if ccw else "false"},
        {"text": "clockwise (opposite to $C$)", "correct": "false" if ccw else "true"},
        {"text": "it does not circulate: a stationary loop in a uniform field carries no current", "correct": "false"},
        {"text": "it alternates direction, because the flux changes sign at $t = " + f"{abs(B0/k):g}" + "$ s", "correct": "false"},
    ]
    random.shuffle(sense_choices)
    lenz_choices = [
        {"text": r"$+\hat{z}$", "correct": "true" if ccw else "false"},
        {"text": r"$-\hat{z}$", "correct": "false" if ccw else "true"},
        {"text": r"in the plane of the loop, along the wire", "correct": "false"},
        {"text": r"it is zero, because the induced field cancels $\mathbf{B}$ exactly", "correct": "false"},
    ]
    random.shuffle(lenz_choices)

    k_signed = f"+ {k:g}" if k > 0 else f"- {abs(k):g}"
    if ccw:
        c_words = "positive, so it flows counterclockwise seen from above."
        d_words = (r"$\mathcal{E} > 0$ drives the current along $C$: counterclockwise. The $+\hat{z}$ flux is <em>decreasing</em> "
                   r"($dB_z/dt = " + f"{k:g}" + r"$ T/s), and the coil responds by making its own $+\hat{z}$ field.")
        e_words = (r"A counterclockwise current (seen from $+z$) produces $+\hat{z}$ field inside the loop — Lenz's law: it tries to replace the flux that is being lost.")
    else:
        c_words = "negative, so it actually flows clockwise seen from above."
        d_words = (r"$\mathcal{E} < 0$ drives the current against $C$: clockwise. The $+\hat{z}$ flux is <em>increasing</em> "
                   r"($dB_z/dt = " + f"{k:g}" + r"$ T/s), and the coil responds by making its own $-\hat{z}$ field.")
        e_words = (r"A clockwise current (seen from $+z$) produces $-\hat{z}$ field inside the loop — Lenz's law: it opposes the growth of the flux.")

    data["params"].update({
        "N": N, "plural": "s" if N > 1 else "", "a_cm": a_cm, "b_cm": b_cm, "R": R, "B0": B0, "k": k, "t1": t1,
        "k_signed": k_signed, "k_signed_plain": f"{k:g}", "A_str": f"{A:.4g}", "Bt1_str": f"{Bt1:.4g}",
        "Psi_str": f"{Psi*1e3:.4g}", "emf_str": f"{emf:.4g}", "I_str": f"{I*1e3:.4g}",
        "c_words": c_words, "d_words": d_words, "e_words": e_words,
        "sense_choices": sense_choices, "lenz_choices": lenz_choices,
    })
    data["correct_answers"]["Psi"] = Psi * 1e3    # mWb
    data["correct_answers"]["emf"] = emf          # V
    data["correct_answers"]["I"] = I * 1e3        # mA
