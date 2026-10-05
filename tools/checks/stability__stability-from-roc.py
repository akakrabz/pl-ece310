"""Independent check of questions/stability/stability-from-roc.

Recomputes causality / stability of the given ROC and the stable ROC from the raw pole list with
plain Fractions (no ece310.zt), and maps the element's displayed options back to the server's
option list by their text."""

import html
import re
from fractions import Fraction


def _norm(s):
    return re.sub(r"\s+", "", html.unescape(s or ""))


def _key(params, name, text):
    for opt in params[name]:
        if _norm(opt["html"]) == _norm(text):
            return opt["key"]
    raise KeyError(f"{name}: no option {text!r}")


def _expected(params):
    mags = sorted({abs(Fraction(p)) for p in params["poles"]})
    g = params["given_roc"]
    inner = Fraction(g["inner"])
    outer = None if g["outer"] is None else Fraction(g["outer"])
    causal = outer is None
    stable = inner < 1 and (outer is None or outer > 1)
    lo = max([m for m in mags if m < 1], default=Fraction(0))
    hi = min([m for m in mags if m > 1], default=None)
    return mags, inner, outer, causal, stable, lo, hi


def check(params, correct):
    probs = []
    mags, inner, outer, causal, stable, lo, hi = _expected(params)
    if any(m == 1 for m in mags):
        probs.append("a pole lies on the unit circle")
    # the given ROC must be one of the rings bounded by consecutive pole circles
    bounds = [Fraction(0)] + mags + [None]
    rings = list(zip(bounds[:-1], bounds[1:]))
    if (inner, outer) not in rings:
        probs.append(f"given ROC ({inner}, {outer}) is not a ring between consecutive poles {mags}")
    if _norm(correct["causal"]["html"]) != ("Yes" if causal else "No"):
        probs.append(f"causal answer {correct['causal']['html']} but expected causal={causal}")
    if _norm(correct["stable"]["html"]) != ("Yes" if stable else "No"):
        probs.append(f"stable answer {correct['stable']['html']} but expected stable={stable}")
    chosen = [c for c in params["roc_choices"] if _norm(c["text"]) == _norm(correct["roc"]["html"])]
    if len(chosen) != 1:
        probs.append("cannot identify the correct ROC option")
    else:
        c = chosen[0]
        got = (Fraction(c["inner"]), None if c["outer"] is None else Fraction(c["outer"]))
        if got != (lo, hi):
            probs.append(f"stable ROC option {got} != expected ({lo}, {hi})")
    if len(params["roc_choices"]) != len(mags) + 1:
        probs.append(f"{len(params['roc_choices'])} ROC options for {len(mags)} distinct pole magnitudes")
    return probs


def submissions(params, correct):
    mags, inner, outer, causal, stable, lo, hi = _expected(params)
    yes_no = lambda b: "Yes" if b else "No"  # noqa: E731
    cases = [
        # the correct keys given explicitly (the only equivalent form a multiple choice has)
        ({"causal": _key(params, "causal", yes_no(causal))}, {"causal": 1}),
        ({"stable": _key(params, "stable", yes_no(stable))}, {"stable": 1}),
        # classic mistakes
        ({"causal": _key(params, "causal", yes_no(not causal))}, {"causal": 0}),
        ({"stable": _key(params, "stable", yes_no(not stable))}, {"stable": 0}),
    ]
    by_bounds = {(Fraction(c["inner"]), None if c["outer"] is None else Fraction(c["outer"])): c["text"]
                 for c in params["roc_choices"]}
    exterior = by_bounds[(mags[-1], None)]
    if (lo, hi) != (mags[-1], None):
        # "causal ROC = stable ROC" confusion: picking the exterior of the outermost pole
        cases.append(({"roc": _key(params, "roc", exterior)}, {"roc": 0}))
    if not stable:
        # keeping the given (unstable) ROC
        cases.append(({"roc": _key(params, "roc", by_bounds[(inner, outer)])}, {"roc": 0}))
    if len(mags) >= 2 and (lo, hi) != (mags[0], mags[1]):
        # "the stable one is always the middle annulus"
        cases.append(({"roc": _key(params, "roc", by_bounds[(mags[0], mags[1])])}, {"roc": 0}))
    cases.append(({"roc": _key(params, "roc", by_bounds[(lo, hi)])}, {"roc": 1}))
    return cases
