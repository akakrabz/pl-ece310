"""Independent check of questions/dtft/periodicity.

Every displayed formula (x_tex, the checkbox options, cos(omega1 n)) is converted from its LaTeX to a
numpy expression in n and evaluated on integers, so the check covers exactly what the student reads:
(a) the fundamental period is found by brute force (smallest N <= 400 with x[n + N] = x[n] on a long
stretch of n), and non-periodicity by the absence of any such N; (b) every option is compared with
cos(omega0 n) sample by sample; (c) the answer is checked by sampling and by arccos(cos(omega1))."""

import math
import re

import numpy as np

NS = np.arange(-60, 400)


def _match(s, i, o="{", c="}"):
    depth = 0
    for k in range(i, len(s)):
        if s[k] == o:
            depth += 1
        elif s[k] == c:
            depth -= 1
            if depth == 0:
                return k
    raise ValueError(f"unbalanced in {s!r}")


def tex_to_py(s):
    s = s.replace("$", "").replace(r"\left(", "(").replace(r"\right)", ")").replace(r"\cdot", "*")
    s = s.replace(r"\,", " ").replace(r"\dfrac", r"\tfrac")
    while r"\tfrac{" in s:
        i = s.index(r"\tfrac{")
        j = _match(s, i + 6)
        k = _match(s, j + 1)
        s = s[:i] + "((" + s[i + 7:j] + ")/(" + s[j + 2:k] + "))" + s[k + 1:]
    while "e^{" in s:
        i = s.index("e^{")
        j = _match(s, i + 2)
        inner = s[i + 3:j]
        assert inner.startswith("j"), inner
        s = s[:i] + "EXP(1j*(" + inner[1:] + "))" + s[j + 1:]
    s = s.replace(r"\cos^2", "COSQ").replace(r"\cos", "COS").replace(r"\sin", "SIN").replace(r"\pi", "PI")
    s = re.sub(r"\s+", " ", s)
    # implicit multiplication: digit/letter/closing paren followed by an operand
    s = re.sub(r"(\d)\s*(?=(PI|n|COS|SIN|EXP|\())", r"\1*", s)
    s = re.sub(r"(\)|PI|n)\s*(?=(PI|n|COS|SIN|EXP|\(|\d))", r"\1*", s)
    return s


def ev(tex, n):
    py = tex_to_py(tex)
    env = {"COS": np.cos, "SIN": np.sin, "EXP": np.exp, "PI": math.pi, "n": n,
           "COSQ": lambda v: np.cos(v) ** 2}
    return np.asarray(eval(py, {"__builtins__": {}}, env), dtype=complex) * np.ones_like(n, dtype=complex)


def fundamental_period(x, max_n=400, tol=1e-9):
    """x: samples on NS. Smallest N with x[n+N] = x[n] for all available n, or None."""
    for N in range(1, max_n):
        a, b = x[:-N], x[N:]
        if len(a) < 40:
            break
        if np.max(np.abs(a - b)) < tol:
            return N
    return None


def check(params, correct):
    probs = []
    x = ev(params["x_tex"], NS.astype(float))
    N0 = fundamental_period(x, max_n=len(NS) - 60)
    want = correct["N0"]
    if params["is_periodic"]:
        if N0 is None or N0 != want:
            probs.append(f"x = {params['x_tex']}: brute-force period {N0}, answer {want!r}")
        if correct["periodic"]["html"].strip() != "Yes":
            probs.append("periodic x but the MC answer is not Yes")
    else:
        if N0 is not None:
            probs.append(f"x = {params['x_tex']} repeats after {N0} samples but is marked not periodic")
        if want != "":
            probs.append(f"non-periodic x needs a blank N0, got {want!r}")
        if correct["periodic"]["html"].strip() != "No":
            probs.append("non-periodic x but the MC answer is not No")
        best = min(np.max(np.abs(x[:-N] - x[N:])) for N in range(1, 300))
        if best < 1e-3:
            probs.append(f"x nearly periodic (min mismatch {best:.2e}): not a clear 'not periodic' variant")
    if params["is_periodic"] and len(_n0_mistakes(params, correct)) < 1:
        probs.append(f"no computed classic N0 mistake besides the blank box ({_n0_mistakes(params, correct)})")
    # (b)
    n = np.arange(-60, 61).astype(float)
    ref = np.cos(float(eval(params["r0"].replace("/", "/"), {"__builtins__": {}})) * math.pi * n)
    if np.max(np.abs(ev(params["cos0_tex"], n) - ref)) > 1e-9:
        probs.append(f"cos0_tex {params['cos0_tex']} does not match r0 = {params['r0']}")
    keys = {o["key"] for o in correct["same"]}
    n_true = 0
    for o in params["same"]:
        same = np.max(np.abs(ev(o["html"], n) - ref)) < 1e-9
        n_true += same
        if same != (o["key"] in keys):
            probs.append(f"option {o['html']} identical={same} but marked {o['key'] in keys}")
    if not 2 <= n_true <= 5 or len(params["same"]) != 6:
        probs.append(f"{n_true} identical options out of {len(params['same'])}")
    # (c)
    r1 = float(eval(params["r1"], {"__builtins__": {}}))
    w2 = correct["w2"]
    if not 0 <= w2 <= 1:
        probs.append(f"omega2 = {w2} pi not in [0, pi]")
    if np.max(np.abs(np.cos(w2 * math.pi * n) - np.cos(r1 * math.pi * n))) > 1e-9:
        probs.append(f"cos({w2} pi n) differs from cos({r1} pi n)")
    if abs(math.acos(math.cos(r1 * math.pi)) / math.pi - w2) > 1e-9:
        probs.append("omega2 is not arccos(cos(omega1))")
    if np.max(np.abs(ev(params["w1_n_tex"], n) - r1 * math.pi * n)) > 1e-9:
        probs.append(f"w1_n_tex {params['w1_n_tex']} is not {r1} pi n")
    if not r1 > 1:
        probs.append("omega1 must exceed pi")
    if abs(w2 - 2 / 7) < 1e-3:
        probs.append("omega2 equals the format example 2/7")
    if w2 in (0.0, 1.0):
        probs.append("omega2 at 0 or pi (degenerate)")
    return probs


def _n0_mistakes(params, correct):
    t = params["traps"]
    cands = [t.get("product"), t.get("unreduced"), t.get("factors"), t.get("max")]
    out = []
    for c in cands:
        if c is not None and c != correct["N0"] and c not in out:
            out.append(c)
    return out


def _key(params, name, text):
    return next(o["key"] for o in params[name] if o["html"].strip() == text)


def _frac(q):
    from fractions import Fraction
    return Fraction(q).limit_denominator(100)


def submissions(params, correct):
    cases = []
    yes, no = _key(params, "periodic", "Yes"), _key(params, "periodic", "No")
    if params["is_periodic"]:
        N0 = correct["N0"]
        cases += [({"N0": f" {N0} "}, {"N0": 1}), ({"N0": f"+{N0}"}, {"N0": 1}),
                  ({"N0": ""}, {"N0": 0}),                                   # "not periodic" box
                  ({"periodic": no, "N0": ""}, {"periodic": 0, "N0": 0}),
                  ({"N0": "2.5"}, {"N0": "invalid"})]
        for m in _n0_mistakes(params, correct):                              # product / unreduced / factors / max
            cases.append(({"N0": str(m)}, {"N0": 0}))
    else:
        t = params["traps"]
        cases += [({"N0": "  "}, {"N0": 1, "periodic": 1}),
                  ({"N0": str(t["rational_period"])}, {"N0": 0}),             # period of the rational term
                  ({"N0": str(t["two_pi_over_w"])}, {"N0": 0}),               # 2 pi / omega, rounded
                  ({"periodic": yes, "N0": ""}, {"periodic": 0, "N0": 0}),    # blank box only with "No"
                  ({"periodic": yes, "N0": str(t["rational_period"])}, {"periodic": 0, "N0": 0})]
    # (b) checkbox: exact keys, everything ticked, one true + one false
    keys = sorted(o["key"] for o in correct["same"])
    allk = sorted(o["key"] for o in params["same"])
    wrong = [k for k in allk if k not in keys]
    cases.append(({"same": keys}, {"same": 1}))
    cases.append(({"same": allk}, {"same": max(0.0, (len(keys) - len(wrong)) / len(keys))}))
    cases.append(({"same": [keys[0], wrong[0]]}, {"same": 0.0}))
    cases.append(({"same": keys[:1]}, {"same": 1 / len(keys)}))
    # (c)
    w2 = _frac(correct["w2"])
    r1 = _frac(eval(params["r1"], {"__builtins__": {}}))
    cases += [({"w2": f"{float(w2):.4g}"}, {"w2": 1}),
              ({"w2": f"{2 * w2.numerator}/{2 * w2.denominator}"}, {"w2": 1})]
    red = r1 - 2 * (r1.numerator // (2 * r1.denominator))
    for m in {r1 - 2, 2 - r1, red, r1, float(w2) * math.pi}:
        m = float(m)
        if abs(m - float(w2)) > 0.02 * max(1.0, abs(float(w2))):
            cases.append(({"w2": f"{m:.6g}"}, {"w2": 0}))
    return cases
