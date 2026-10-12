"""Independent check of questions/convolution/convolution-facts.

Numbers: np.convolve with the start index nx + nh, sums and samples read off that array.
Multiple choice: every option's LaTeX is translated into a Python expression in n (x[.], y[.],
delta[.], u[.], running sums) and evaluated on random integer test signals; exactly one option
must equal the stem's operation, computed with np.convolve, at every n."""

import re

import numpy as np


def _y(params):
    y = np.convolve(params["x"], params["h"])
    return y, params["a1"] + params["a2"]


def check(params, correct):
    probs = []
    y, ny = _y(params)
    nl = ny + len(y) - 1
    want = {"nfirst": ny, "nlast": nl, "leny": len(y), "sumy": int(y.sum()), "yfirst": int(y[0]),
            "ylast": int(y[-1]), "yn0": int(y[params["n0"] - ny])}
    for k, v in want.items():
        if correct[k] != v:
            probs.append(f"{k}: correct answer {correct[k]} but numpy gives {v}")
    if y[0] == 0 or y[-1] == 0:
        probs.append("y has a zero end sample")
    if not ny < params["n0"] < nl:
        probs.append(f"n0 = {params['n0']} is not interior")
    if sum(params["x"]) * sum(params["h"]) == sum(params["x"]) + sum(params["h"]):
        pass   # probe skipped in submissions()
    probs += _check_mc(params, correct)
    return probs


# ----------------------------------------------------------------------------- multiple choice
LO, HI = -40, 40


def _tex_to_py(t):
    s = t.strip().strip("$")
    s = s.replace(r"\big(", "(").replace(r"\big)", ")")
    s = re.sub(r"\\sum_\{m=-\\infty\}\^\{\\infty\} x\[m\]", "SALL()", s)
    s = re.sub(r"\\sum_\{m=-\\infty\}\^\{([^}]*)\} x\[m\]", r"SLE(\1)", s)
    s = re.sub(r"\\sum_\{m=([^}]*)\}\^\{\\infty\} x\[m\]", r"SGE(\1)", s)
    s = s.replace(r"\delta[", "D[").replace(r"\,", "*")
    s = re.sub(r"([xyDu])\[([^\]]*)\]", lambda m: f"{m.group(1).upper()}({m.group(2)})", s)
    s = re.sub(r"\)\s*([XYDU(])", r")*\1", s)
    return s


def _evaluator(xv, x0, yv, y0):
    def X(k):
        i = k - x0
        return xv[i] if 0 <= i < len(xv) else 0

    def Y(k):
        i = k - y0
        return yv[i] if 0 <= i < len(yv) else 0

    def D(k):
        return 1 if k == 0 else 0

    def U(k):
        return 1 if k >= 0 else 0

    def SLE(k):
        return sum(X(m) for m in range(LO - 20, k + 1))

    def SGE(k):
        return sum(X(m) for m in range(k, HI + 20))

    def SALL():
        return sum(xv)

    return {"X": X, "Y": Y, "D": D, "U": U, "SLE": SLE, "SGE": SGE, "SALL": SALL}


def _true_and_env(params, rng):
    kind, prm = params["mc_kind"], params["mc_param"]
    xv = [int(v) for v in rng.integers(-3, 4, 5)]
    xv[0], xv[-1] = 2, -1
    x0 = int(rng.integers(-3, 3))
    hv = [int(v) for v in rng.integers(-3, 4, 3)]
    hv[0], hv[-1] = 1, 3
    h0 = int(rng.integers(-2, 3))
    yv = list(np.convolve(xv, hv))
    env = _evaluator(xv, x0, yv, x0 + h0)
    if kind == "delta":
        k = prm["k"]
        g, g0 = [1], k                                   # delta[n - k]
    elif kind == "diff":
        s = prm["s"]
        g, g0 = [1, -1], -s                              # delta[n + s] - delta[n + s - 1]
    elif kind == "step":
        k = prm["k"]
        g, g0 = [1] * (2 * (HI - LO) + 100), k           # u[n - k], long enough on [LO, HI]
    else:
        a, b = prm["a"], prm["b"]                        # x[n - a] * h[n - b]
        w = np.convolve(xv, hv)
        w0 = x0 + a + h0 + b
        return (lambda n: int(w[n - w0]) if 0 <= n - w0 < len(w) else 0), env
    v = np.convolve(xv, g)
    v0 = x0 + g0
    return (lambda n: int(v[n - v0]) if 0 <= n - v0 < len(v) else 0), env


def _check_mc(params, correct):
    probs = []
    rng = np.random.default_rng(sum(map(ord, params["mc_stem_html"])))
    matches = []
    for trial in range(2):
        true, env = _true_and_env(params, rng)
        ok = []
        for opt in params["prop"]:
            expr = _tex_to_py(opt["html"])
            try:
                vals = [eval(expr, {"__builtins__": {}}, dict(env, n=n)) for n in range(LO, HI + 1)]
            except Exception as exc:
                probs.append(f"cannot evaluate option {opt['html']!r} as {expr!r}: {exc}")
                return probs
            ok.append(all(int(v) == true(n) for v, n in zip(vals, range(LO, HI + 1))))
        matches.append(ok)
    good = [all(m[i] for m in matches) for i in range(len(params["prop"]))]
    if sum(good) != 1:
        probs.append(f"{sum(good)} options match the stem: {[o['html'] for o in params['prop']]}")
    elif params["prop"][good.index(True)]["key"] != correct["prop"]["key"]:
        probs.append(f"marked option {correct['prop']['html']} but {params['prop'][good.index(True)]['html']} matches")
    if len(params["prop"]) != 4:
        probs.append(f"{len(params['prop'])} options")
    return probs


# ----------------------------------------------------------------------------- submissions
def submissions(params, correct):
    y, ny = _y(params)
    nl = ny + len(y) - 1
    a1, a2 = params["a1"], params["a2"]
    sx, sh = sum(params["x"]), sum(params["h"])
    n0 = params["n0"]
    hrev = list(params["h"])[::-1]                         # forgot to flip h inside its support
    yrev = np.convolve(params["x"], hrev)
    wrong_n0 = int(yrev[n0 - ny])
    cases = [
        ({"nfirst": f"+{ny}" if ny >= 0 else f" {ny} "}, {"nfirst": 1}),
        ({"sumy": f" {sx * sh} "}, {"sumy": 1}),
        ({"leny": str(len(y) + 1)}, {"leny": 0}),                    # N + M (off by one)
        ({"nlast": str(nl + 1)}, {"nlast": 0}),
        ({"nfirst": "1.5"}, {"nfirst": "invalid"}),
    ]
    if a2 != 0:
        cases.append(({"nfirst": str(a1 - a2)}, {"nfirst": 0}))    # subtracted the start indices
    if sx + sh != sx * sh:
        cases.append(({"sumy": str(sx + sh)}, {"sumy": 0}))       # added the sums
    if y[0] != y[-1]:
        cases.append(({"yfirst": str(int(y[-1])), "ylast": str(int(y[0]))}, {"yfirst": 0, "ylast": 0}))
    if wrong_n0 != int(y[n0 - ny]):
        cases.append(({"yn0": str(wrong_n0)}, {"yn0": 0}))
    for d in (-1, 1):                                            # neighbouring sample (index slip)
        if int(y[n0 - ny + d]) != int(y[n0 - ny]):
            cases.append(({"yn0": str(int(y[n0 - ny + d]))}, {"yn0": 0}))
            break
    wrong_keys = [o["key"] for o in params["prop"] if o["key"] != correct["prop"]["key"]]
    cases.append(({"prop": wrong_keys[0]}, {"prop": 0}))
    cases.append(({"prop": correct["prop"]["key"]}, {"prop": 1}))
    return cases
