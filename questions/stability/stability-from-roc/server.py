"""H(z) with two or three real poles and a stated ROC: is the system causal? BIBO stable? Which ROC
would make it stable? (Lecture 11, Table 1; HW4 #3; FA2023 #1f, FA2024 #1f, SP2025 #1e, FA2025 #1c.)

Poles never lie on the unit circle, so exactly one of the possible ROCs contains |z| = 1. Some
three-pole variants contain a pair +-r, written unfactored as (1 - r^2 z^-2): the pair shares one
circle, so there are only three ROCs. H is proper in z^-1, so the exterior ROC includes infinity."""

import random
from fractions import Fraction as F

from ece310 import fmt, zt

MAGS = [F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(4, 3), F(3, 2), F(2), F(5, 2), F(3)]
ZEROS = [F(1), F(-1), F(2), F(-2), F(1, 2), F(-1, 2), F(3), F(-3)]


def _kind(r):
    if r.inner == 0:
        return "left-sided (anti-causal)"
    if r.outer is None:
        return "right-sided, causal"
    return "two-sided"


def generate(data):
    while True:
        npoles = random.choice([2, 2, 3])
        pair = npoles == 3 and random.random() < 0.3
        mags = random.sample(MAGS, npoles - 1 if pair else npoles)
        if (all(m < 1 for m in mags) or all(m > 1 for m in mags)) and random.random() < 0.6:
            continue  # keep most variants mixed (inside and outside the unit circle)
        if pair:
            rpair, singles = mags[0], mags[1:]
        else:
            rpair, singles = None, mags
        singles = [m if random.random() < 0.6 else -m for m in singles]
        poles = list(singles) + ([rpair, -rpair] if pair else [])
        zero = random.choice(ZEROS) if random.random() < 0.6 else None
        if zero is not None and zero in poles:
            continue
        break

    cands = zt.rocs_for_poles(poles)          # inside smallest, annuli, outside largest (with infinity)
    gi = random.randrange(len(cands))
    given = cands[gi]
    causal = given.outer is None
    stable = given.contains_unit_circle()
    si = next(i for i, r in enumerate(cands) if r.contains_unit_circle())

    # ---- LaTeX of H(z)
    b0 = random.choice([1, 1, 2, 3])
    if zero is None:
        num_tex = str(b0)
    else:
        num_tex = fmt.tex_sum([(b0, ""), (-b0 * zero, "z^{-1}")])
    factors = "".join(fmt.tex_factor(p) for p in singles)
    if pair:
        factors += r"\left(" + fmt.tex_sum([(1, ""), (-rpair * rpair, "z^{-2}")]) + r"\right)"
    H_tex = fmt.tex_frac(num_tex, factors)

    pole_items = [fr"{fmt.tex_num(p)}" for p in singles]
    if pair:
        pole_items.append(r"\pm" + fmt.tex_num(rpair))
    pole_list_tex = r",\ ".join(pole_items)
    mags_sorted = sorted({abs(p) for p in poles})
    mags_tex = r" \lt ".join(fmt.tex_num(m) for m in mags_sorted)

    rows = []
    for i, r in enumerate(cands):
        cls = ' class="table-active"' if i == gi else ""
        unit = "yes" if r.contains_unit_circle() else "no"
        caus = "yes" if r.is_causal_system() else "no"
        note = "&larr; given" if i == gi else ""
        rows.append(f"<tr{cls}><td>${r.tex()}$</td><td>{_kind(r)}</td><td>{unit}</td><td>{caus}</td><td>{note}</td></tr>")
    table_html = ('<table class="table table-sm table-bordered text-center" style="width:auto">'
                  "<tr><th>ROC</th><th>$h[n]$ is</th><th>contains $\\lvert z\\rvert = 1$?</th><th>causal?</th><th></th></tr>"
                  + "".join(rows) + "</table>")

    p = data["params"]
    p["H_tex"] = H_tex
    p["roc_tex"] = given.tex()
    p["pole_list_tex"] = pole_list_tex
    p["mags_tex"] = mags_tex
    p["nrocs"] = len(cands)
    p["pair"] = pair
    p["pair_tex"] = (fr"1 - {fmt.tex_num(rpair * rpair)}z^{{-2}} = {fmt.tex_factor(rpair)}{fmt.tex_factor(-rpair)}"
                     if pair else "")
    p["zero_tex"] = fmt.tex_num(zero) if zero is not None else ""
    p["has_zero"] = zero is not None
    p["stable_roc_tex"] = cands[si].tex()
    p["causal_word"] = "is" if causal else "is not"
    p["stable_word"] = "contains" if stable else "does not contain"
    p["stable_not"] = "" if stable else "not "
    p["given_kind"] = _kind(given)
    p["table_html"] = table_html
    # data for the checker (exact strings)
    p["poles"] = [str(q) for q in poles]
    p["given_roc"] = {"inner": str(given.inner), "outer": None if given.outer is None else str(given.outer)}
    p["causal_choices"] = [{"text": "Yes", "correct": causal}, {"text": "No", "correct": not causal}]
    p["stable_choices"] = [{"text": "Yes", "correct": stable}, {"text": "No", "correct": not stable}]
    p["roc_choices"] = [{"text": f"${r.tex()}$", "correct": i == si, "inner": str(r.inner),
                         "outer": None if r.outer is None else str(r.outer)} for i, r in enumerate(cands)]
