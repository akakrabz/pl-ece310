"""Parallel plates filled with a dielectric: charge fixed (battery removed) vs voltage fixed (battery connected).
HW4 #5 family, plus the Lecture 10 lossy-capacitor concept for the conducting-fluid part.

Initially vacuum with E = E0 z, so rho_s(0) = eps0 E0 and D = eps0 E0.
    charge fixed:   D unchanged,  E = E0/eps_r,  P = D - eps0 E = eps0 E0 (1 - 1/eps_r),   W/W0 = 1/eps_r
    voltage fixed:  E unchanged,  D = eps_r eps0 E0,  P = (eps_r - 1) eps0 E0,             W/W0 = eps_r
Conducting fluid, steady state:
    charge fixed:   the plates discharge through the fluid (no external path): E = D = P = 0
    voltage fixed:  E = V/W is pinned by the battery, so E, D, P are unchanged and J = sigma E flows steadily.
"""
import random


def generate(data):
    mode = random.choice(["Q", "V"])
    E0 = random.choice([2, 3, 4, 5, 6, 9, 12])
    er = random.choice([2, 3, 4, 5, 9, 81])
    sigma = random.choice([0.5, 1, 2, 4])

    if mode == "Q":
        scenario = "The battery is then disconnected, leaving the plates isolated."
        D, E, P = E0, E0 / er, E0 * (1 - 1 / er)       # D, P as coefficients of eps0
        Wratio = 1 / er
        b_expl = (r"The plates are isolated, so $\rho_s$ and therefore $D_z = \rho_s$ cannot change: $D_z = " + f"{E0}" +
                  r"\,\epsilon_0$. Then $E_z = D_z/\epsilon = " + f"{E0}/{er} = {E:.4g}" + r"\ \mathrm{V/m}$ and " +
                  r"$P_z = D_z - \epsilon_0E_z = " + f"{P:.4g}" + r"\,\epsilon_0$: the bound charge on the fluid's faces cancels most of the plate charge.")
        c_expl = (r"$D$ is unchanged and $E$ dropped by $\epsilon_r$, so the stored energy fell to $1/\epsilon_r$ of its value "
                  r"(the fluid was pulled in by the fringing field, doing work on it).")
        d_expl = (r"With the plates isolated there is no external path for current, but the fluid now conducts: $\mathbf{J} = \sigma\mathbf{E}$ carries charge from the positive plate to the negative one until $\mathbf{E} = 0$. "
                  r"Then $\mathbf{D} = \epsilon\mathbf{E} = 0$ and $\mathbf{P} = (\epsilon-\epsilon_0)\mathbf{E} = 0$ as well: the capacitor has discharged through its own shunt conductance, $Q(t) = Q_0e^{-t/\tau}$ with $\tau = \epsilon/\sigma$.")
    else:
        scenario = "The battery stays connected, holding the plate voltage fixed."
        D, E, P = er * E0, E0, (er - 1) * E0
        Wratio = er
        b_expl = (r"The battery pins $V$, hence $E_z = V/W$ is unchanged: $E_z = " + f"{E0}" + r"\ \mathrm{V/m}$. Then $D_z = \epsilon E_z = " +
                  f"{er}\\times{E0}\\,\\epsilon_0 = {D}" + r"\,\epsilon_0$ and $P_z = D_z - \epsilon_0E_z = " + f"{P}" +
                  r"\,\epsilon_0$. The battery had to deliver extra charge: $\rho_s$ rose from $" + f"{E0}" + r"\,\epsilon_0$ to $" + f"{D}" + r"\,\epsilon_0$.")
        c_expl = (r"$E$ is unchanged and $D$ grew by $\epsilon_r$, so the stored energy rose by $\epsilon_r$ "
                  r"(paid for by the battery, which moved the extra charge onto the plates).")
        d_expl = (r"The battery keeps $V$, and therefore $\mathbf{E} = (V/W)\hat{z}$, fixed; $\mathbf{D} = \epsilon\mathbf{E}$ and $\mathbf{P} = (\epsilon - \epsilon_0)\mathbf{E}$ are unchanged too. "
                  r"What changes is that a steady current density $\mathbf{J} = \sigma\mathbf{E}$ now flows through the fluid, supplied by the battery: the device behaves as $C$ in parallel with $G = (\sigma/\epsilon)C$.")

    choice_Q = {"text": r"$\mathbf{E} = \mathbf{D} = \mathbf{P} = 0$: the conduction current $\mathbf{J} = \sigma\mathbf{E}$ carries the plate charge across the gap until the plates are neutral, and nothing replaces it.",
                "correct": "true" if mode == "Q" else "false"}
    choice_V = {"text": r"$\mathbf{E}$, $\mathbf{D}$ and $\mathbf{P}$ are unchanged; a steady current density $\mathbf{J} = \sigma\mathbf{E}$ flows through the gap, supplied by the battery.",
                "correct": "true" if mode == "V" else "false"}
    distractors = [
        {"text": r"$\mathbf{E} = 0$ inside the fluid (it is now a conductor), but $\mathbf{D}$ and $\mathbf{P}$ keep their pre-salt values because the plate charge has not moved.", "correct": "false"},
        {"text": r"$\mathbf{E}$ is unchanged, but $\mathbf{D}$ and $\mathbf{P}$ drop to zero because a conducting medium cannot be polarized.", "correct": "false"},
    ]
    choices = [choice_Q, choice_V] + distractors
    random.shuffle(choices)

    data["params"].update({
        "mode": mode, "E0": E0, "er": er, "sigma": sigma, "scenario": scenario,
        "b_expl": b_expl, "c_expl": c_expl, "d_expl": d_expl, "Wratio_str": f"{Wratio:.4g}", "choices": choices,
    })
    data["correct_answers"]["rho0"] = float(E0)
    data["correct_answers"]["D"] = float(D)
    data["correct_answers"]["E"] = float(E)
    data["correct_answers"]["P"] = float(P)
    data["correct_answers"]["Wratio"] = float(Wratio)
