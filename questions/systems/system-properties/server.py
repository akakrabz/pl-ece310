"""Linear / time-invariant / causal / BIBO-stable table for one system y[n] = T{x}[n]
(exam family "system properties", 7 of 7 past Midterm 1 exams; HW1 #5-6, HW2 #1; Lecture 3).

15 % of the variants come from the verified bank (past exam tables, homework, Lecture 3 list,
practice problems), 85 % from the compositional generator in ece310.systems_bank, whose flags are
derived by combination rules and re-verified numerically by the checker."""

import random

from ece310 import mc, systems_bank

NAMES = ["lin", "ti", "caus", "stab"]


def generate(data):
    if random.random() < 0.15:
        s = systems_bank.bank_system()
        src = s["source"]
        origin = {"kind": "bank", "id": s["id"]}
        src_html = f"This system is from {src} (system-property bank)."
    else:
        s = systems_bank.generate_system()
        origin = {"kind": "gen", "struct": s["struct"]}
        src_html = ("A new system built from the patterns of the past exam tables (FA2025 #2, SP2025 #2, FA2024 #2, "
                    "FA2023 #2, SP2023 #2, SP2021 #3, FA2019 #2) and the Lecture 3 practice list.")
    p = data["params"]
    p["sys_tex"] = s["tex"]
    p["origin"] = origin
    p["src_html"] = src_html
    rows = []
    for name, prop, flag, why in zip(NAMES, systems_bank.PROPS, s["flags"], s["reasons"]):
        p[f"{name}_choices"] = mc.yes_no(flag)
        rows.append({"prop": prop, "ans": "Yes" if flag else "No", "why": why})
    p["rows"] = rows
