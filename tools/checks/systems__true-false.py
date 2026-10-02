"""Independent check of questions/systems/true-false.

Fixed statements: answer table typed from the official keys (study-site True/False bank, Q1-Q43).
Parameterized statements: recomputed numerically from the stored parameters (summability of the
candidate impulse responses, lfilter simulations, brute-force delta sums, numpy.roots)."""

import numpy as np
from scipy.signal import lfilter

OFFICIAL = {   # sid -> answer (exam, bank number)
    "h_any_system": False,          # FA2024 #1(a), Q1
    "h_conv_any": False,            # FA2023 #1(d), Q2
    "fully_described_lti": True,    # SP2023 #1(a), Q4
    "abs_sum_any_system": False,    # SP2025 #1(b), Q14
    "causal_conv_causal": True,     # FA2025 #1(b), Q5
    "conv_sum_causal": False,       # SP2025 #1(a), Q6
    "right_sided_causal": False,    # SP2023 #1(c), Q7
    "causal_implies_ti": False,     # FA2023 #1(a), Q9
    "tv_not_causal": False,         # FA2019 #1(f), Q10
    "bounded_h_stable": False,      # FA2025 #1(a), Q11
    "unstable_h_unbounded": False,  # FA2024 #1(b), Q12
    "two_sided_never_stable": False,  # SP2025 #1(f), Q15
    "unstable_any_input": False,    # FA2019 #1(j), Q19
    "stable_unbounded_in_out": False,  # FA2025 #1(e), Q16
    "parallel_stable": True,        # SP2025 #1(d), Q21
    "parallel_stable_any": True,    # FA2023 #1(c), Q23
    "parallel_unstable": False,     # FA2025 #1(f), Q24
    "cascade_unstable": False,      # SP2021 #1(c), Q25
    "series_order": True,           # FA2023 #1(e), Q28
    "roc_no_zeros": False,          # FA2023 #1(b), Q30
    "poles_of_sum": False,          # SP2021 #1(b), Q31
    "lccde_finite_poles": True,     # FA2024 #1(d), Q38
    "fir_stable_or_unstable": False,  # FA2019 #1(e), Q41
}


def f(pq):
    return pq[0] / pq[1]


def summable(h_of_n, side, N=400):
    """Is sum |h[n]| finite on the given side? Compare partial sums over 200 and 400 terms."""
    n = np.arange(0, N) if side == "right" else -np.arange(1, N + 1)
    with np.errstate(all="ignore"):
        a = np.abs(h_of_n(n))
    s1, s2 = a[: N // 2].sum(), a.sum()
    return bool(np.isfinite(s2) and s2 - s1 < 1e-6 * max(1.0, s1))


def stable_side(p):
    """Side ('right'/'left') of the summable term for a simple pole p, or None if neither is summable."""
    if summable(lambda n: p ** n.astype(float), "right"):
        return "right"
    if summable(lambda n: p ** n.astype(float), "left"):
        return "left"
    return None


def grows(y):
    m1, m2 = np.max(np.abs(y[: len(y) // 2])), np.max(np.abs(y))
    return (not np.isfinite(m2)) or m2 > 1.5 * m1 + 1e-9


def recompute(sid, prm):
    if sid in OFFICIAL:
        return OFFICIAL[sid]
    if sid == "stable_pole_causality":                   # "must be non-causal / anti-causal"
        side = stable_side(f(prm["p"]))
        return side == "left"
    if sid == "two_poles_two_sided":
        s1, s2 = stable_side(f(prm["p1"])), stable_side(f(prm["p2"]))
        return s1 != s2
    if sid == "no_roc_not_stable":
        return stable_side(f(prm["a"])) is None or stable_side(f(prm["b"])) is None
    if sid == "causal_ramp":
        p = f(prm["p"])
        with np.errstate(all="ignore"):
            y = lfilter([0, 1], [1, -p], np.ones(800))
        return grows(y)
    if sid == "delta_pick":
        n = np.arange(-60, 61)
        hits = n[(2.0 ** n) * (n >= 0) - 2 ** prm["m"] == 0]
        # the sum equals x at the hit indices; with a single hit it pins that sample
        return len(hits) == 1 and int(hits[0]) == prm["mc"] and prm["vc"] == prm["v"]
    if sid == "delta_all":
        n = np.arange(-200, 201)
        arg = 4 * np.cos(2 * n * np.pi + np.pi / 2) - 6 * np.sin(n * np.pi)
        everywhere = bool(np.all(np.abs(arg) < 1e-9))
        return everywhere and prm["Sc"] == prm["S"]
    if sid == "step_pole":
        p = f(prm["p"])
        step_factor = 1 / (1 - 1 / p) if p != 1 else np.inf  # value of 1/(1 - z^-1) at z = p
        return bool(np.isfinite(step_factor) and step_factor != 0)
    if sid == "ztrans_exp_dne":
        w = np.pi * prm["w"][0] / prm["w"][1]
        n = np.arange(-50, 51)
        return bool(np.allclose(np.abs(np.exp(1j * w * n)), 1.0))   # |x| = 1 on both sides: ROC empty
    if sid == "h_an_stable":
        a = f(prm["a"])
        return summable(lambda n: a ** n.astype(float), "right")
    if sid == "unbounded_input_unstable":
        a = f(prm["a"])
        with np.errstate(all="ignore"):
            x = a ** np.arange(0, 800, dtype=float)
        return not grows(x)                                # True iff the input itself is bounded
    if sid == "lccde_poles_distinct":
        N = prm["N"]
        r = np.roots([1] + [0] * (N - 1) + [-1])
        return len(r) == N and min(abs(a - b) for i, a in enumerate(r) for b in r[i + 1:]) > 1e-6
    if sid == "lccde_causal_or_anticausal":
        a = f(prm["a"])
        n = np.arange(-12, 13)
        hc = np.where(n >= 0, a ** n.astype(float), 0.0)
        ha = np.where(n <= -1, -(a ** n.astype(float)), 0.0)
        d = (n == 0).astype(float)
        ok_c = np.allclose(hc[1:] - a * hc[:-1], d[1:])
        ok_a = np.allclose(ha[1:] - a * ha[:-1], d[1:])
        return bool(ok_c and ok_a)
    raise KeyError(sid)


def check(params, correct):
    probs = []
    items = params["statements"]
    if len(items) != 6:
        probs.append(f"{len(items)} statements instead of 6")
    if len({it["sid"] for it in items}) != len(items):
        probs.append("duplicate statement")
    for it in items:
        want = recompute(it["sid"], it["param"])
        graded = correct[it["name"]]["html"].strip()
        if graded != ("True" if want else "False"):
            probs.append(f"{it['name']} ({it['sid']}, {it['param']}): graded {graded}, recomputed {want}")
        if it["ans"] != graded:
            probs.append(f"{it['name']}: answer panel says {it['ans']} but grading uses {graded}")
    return probs


def _key(params, name, text):
    for opt in params[name]:
        if opt["html"].strip() == text:
            return opt["key"]
    raise KeyError(name)


def submissions(params, correct):
    names = [it["name"] for it in params["statements"]]
    right = {nm: correct[nm]["html"].strip() for nm in names}
    flip = {"True": "False", "False": "True"}
    cases = [({nm: _key(params, nm, flip[right[nm]])}, {nm: 0}) for nm in names[:3]]    # flipped answers
    cases.append(({nm: _key(params, nm, flip[right[nm]]) for nm in names}, {nm: 0 for nm in names}))
    cases.append(({nm: _key(params, nm, right[nm]) for nm in names}, {nm: 1 for nm in names}))   # explicit keys
    cases.append(({names[0]: correct[names[0]]["key"]}, {names[0]: 1}))
    cases.append(({nm: _key(params, nm, "False") for nm in names},                              # "all False" guess
                  {nm: (1 if right[nm] == "False" else 0) for nm in names}))
    return cases
