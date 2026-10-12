"""Independent check of questions/dtft/symmetry-and-existence.

(a) X_d(w0) is recomputed from w0 alone with this file's own table of the given formulas (numpy):
    reduce w0 by 2 pi into [-pi, pi]; for a real x use X_d(-w) = X_d(w)^* when the reduced frequency is
    negative. For the real families the Hermitian extension is also inverted numerically (x[n] for
    |n| <= 3 must come out real).
(b) H_d(-w) == H_d(w)^* is tested on a frequency grid with this file's own numpy H_d.
(c), (d) each signal is generated here (numpy) and classified from its samples: sum |x[n]| converges
    -> summable; bounded but sum |x|^2 converges -> ordinary finite-energy DTFT; bounded with infinite
    energy -> impulses; growing -> no DTFT."""

import cmath
import math
import re

import numpy as np

F_NUM = {
    "jw": (r"j\omega", lambda w: 1j * w),
    "1pw_e1": (r"(1 + \omega)\,e^{-j\omega}", lambda w: (1 + w) * np.exp(-1j * w)),
    "w2_e2": (r"\omega^2 e^{-j2\omega}", lambda w: w * w * np.exp(-2j * w)),
    "pimw_e1": (r"(\pi - \omega)\,e^{j\omega}", lambda w: (math.pi - w) * np.exp(1j * w)),
    "w_e3": (r"\omega\,e^{-j3\omega}", lambda w: w * np.exp(-3j * w)),
    "w_jw2": (r"\omega + j\omega^2", lambda w: w + 1j * w * w),
    "j1pw": (r"j(1 + \omega)", lambda w: 1j * (1 + w)),
    "jw_pimw": (r"j\omega(\pi - \omega)", lambda w: 1j * w * (math.pi - w)),
    "1pw2_e1": (r"(1 + \omega^2)\,e^{-j\omega}", lambda w: (1 + w * w) * np.exp(-1j * w)),
    "w2_ep1": (r"\omega^2 e^{j\omega}", lambda w: w * w * np.exp(1j * w)),
}

H_NUM = {
    "absw_sin": lambda w: np.abs(w) * np.exp(-1j * np.pi * np.sin(w)),
    "w_cos": lambda w: w * np.exp(1j * np.pi * np.cos(w)),
    "w2_cos": lambda w: w ** 2 * np.exp(1j * np.cos(w)),
    "jw_e2": lambda w: 1j * w * np.exp(-2j * w),
    "w_e1": lambda w: w * np.exp(-1j * w),
    "cos_jsin2": lambda w: np.cos(w) + 1j * np.sin(2 * w),
    "cos_jcos2": lambda w: np.cos(w) + 1j * np.cos(2 * w),
    "e_jw2": lambda w: np.exp(-1j * w ** 2),
    "absw_e3": lambda w: np.abs(w) * np.exp(-3j * w),
    "jsin_e1": lambda w: 1j * np.sin(w) * np.exp(-1j * w),
    "sin_e1": lambda w: np.sin(w) * np.exp(-1j * w),
    "w2_jsin": lambda w: w ** 2 * np.exp(-1j * np.sin(w)),
    "pole_real": lambda w: 1 / (1 - 0.5 * np.exp(-1j * w)),
    "pole_j": lambda w: 1 / (1 - 0.5j * np.exp(-1j * w)),
    "e_jabs": lambda w: np.exp(1j * np.abs(w)),
    "pimabs_e1": lambda w: (np.pi - np.abs(w)) * np.exp(1j * w),
}


def frac(s):
    a, b = (s.split("/") + ["1"])[:2]
    return int(a) / int(b)


def signal(kind, args, n):
    n = n.astype(float)
    with np.errstate(all="ignore"):
        if kind in ("geo", "geo_big"):
            a = frac(args[0])
            return np.where(n >= 0, np.sign(a) ** n * np.abs(a) ** n, 0.0)
        if kind == "nageo":
            a = frac(args[0])
            return np.where(n >= 0, n * np.sign(a) ** n * np.abs(a) ** n, 0.0)
        if kind == "absgeo":
            a = frac(args[0])
            return np.sign(a) ** np.abs(n) * np.abs(a) ** np.abs(n)
        if kind in ("leftgrow_ok", "left_bad"):
            b = frac(args[0])
            return np.where(n <= 0, np.sign(b) ** n * np.abs(b) ** n, 0.0)
        if kind == "finite":
            return ((n >= 0) & (n < int(args[0]))).astype(float)
        if kind == "ngrow":
            a = frac(args[0])
            return np.where(n >= 0, n * np.sign(a) ** n * np.abs(a) ** n, 0.0)
        if kind == "step":
            return (n >= 0).astype(float)
        if kind == "alt_step":
            return np.where(n >= 0, np.cos(np.pi * n), 0.0)
        if kind == "cos":
            return np.cos(frac(args[0]) * np.pi * n)
        if kind == "cexp":
            return np.exp(1j * frac(args[0]) * np.pi * n)
        if kind == "cexp_step":
            return np.where(n >= 0, np.exp(1j * frac(args[0]) * np.pi * n), 0.0)
        if kind == "sinc":
            wc = frac(args[0]) * np.pi
            safe = np.where(n == 0, 1.0, n)
            return np.where(n == 0, wc / np.pi, np.sin(wc * safe) / (np.pi * safe))
    raise KeyError(kind)


def _match(s, i):
    depth = 0
    for k in range(i, len(s)):
        depth += {"{": 1, "}": -1}.get(s[k], 0)
        if depth == 0:
            return k
    raise ValueError(s)


def tex_signal(tex, n):
    """Evaluate the displayed LaTeX of a signal literally (e^{j-a n} would mean e^{j} e^{-a n})."""
    s = tex.replace(r"\ \text{for all } n", "").replace(r"\left(", "(").replace(r"\right)", ")")
    s = s.replace(r"\,", " ").replace(r"\dfrac", r"\tfrac").replace("|n|", "ABS(n)")
    while r"\tfrac{" in s:
        i = s.index(r"\tfrac{")
        j = _match(s, i + 6)
        k = _match(s, j + 1)
        s = s[:i] + "((" + s[i + 7:j] + ")/(" + s[j + 2:k] + "))" + s[k + 1:]
    while "e^{" in s:
        i = s.index("e^{")
        j = _match(s, i + 2)
        inner = s[i + 3:j].strip()
        sign = ""
        if inner.startswith("-"):
            sign, inner = "-", inner[1:].lstrip()
        assert inner.startswith("j"), inner
        rest = inner[1:]
        body = ("1j" + rest) if rest.lstrip()[:1] in "+-" else ("1j*(" + rest + ")")
        s = s[:i] + "EXP(" + sign + "(" + body + "))" + s[j + 1:]
    s = re.sub(r"\^\{([^{}]*)\}", r"**(\1)", s)
    s = re.sub(r"u\[([^\]]*)\]", r"U(\1)", s)
    s = s.replace(r"\cos", "COS").replace(r"\sin", "SIN").replace(r"\pi", "PI")
    s = re.sub(r"(\d)\s*(?=(PI|n|COS|SIN|EXP|ABS|U\(|\())", r"\1*", s)
    s = re.sub(r"(\)|PI|\bn\b)\s*(?=(PI|\bn\b|COS|SIN|EXP|ABS|U\(|\(|\d))", r"\1*", s)
    env = {"EXP": np.exp, "COS": np.cos, "SIN": np.sin, "PI": math.pi, "ABS": np.abs, "n": n,
           "U": lambda m: (np.asarray(m) >= 0).astype(float)}
    with np.errstate(all="ignore"):
        return np.asarray(eval(s, {"__builtins__": {}}, env), dtype=complex) * np.ones_like(n, dtype=complex)


def classify(kind, args):
    def part(N):
        n = np.arange(-N, N + 1)
        with np.errstate(all="ignore"):
            x = np.abs(signal(kind, args, n))
        return x
    x1, x2 = part(1000), part(4000)
    if not np.all(np.isfinite(x2)) or np.max(x2) > 10 * max(1.0, np.max(part(100))):
        return 3
    s1, s2 = x1.sum(), x2.sum()
    if s2 - s1 < 1e-6 * max(1.0, s1):
        return 0
    e1, e2 = (x1 ** 2).sum(), (x2 ** 2).sum()
    return 2 if e2 - e1 < 1e-3 * e1 else 1


def expected_value(params):
    fid = params["fid"]
    w0 = params["w0"][0] / params["w0"][1]
    red = w0
    while red > 1:
        red -= 2
    while red < -1:
        red += 2
    f = F_NUM[fid][1]
    if params["family"] == "per":
        return complex(f(red * math.pi))
    if red < 0:
        return complex(np.conj(f(-red * math.pi)))
    return complex(f(red * math.pi))


def check(params, correct):
    probs = []
    fam, fid = params["family"], params["fid"]
    if params["F_tex"] != F_NUM[fid][0]:
        probs.append(f"F_tex {params['F_tex']} does not match id {fid}")
    w0 = params["w0"][0] / params["w0"][1]
    ranges = {"herm": (-1, 0), "per": None, "both": (1, 2)}
    if fam == "per" and not (2 < abs(w0) < 3 or 1 < abs(w0) < 2):
        probs.append(f"periodicity variant at w0 = {w0} pi")
    if ranges.get(fam) and not ranges[fam][0] <= w0 < ranges[fam][1]:
        probs.append(f"{fam}: w0 = {w0} pi outside {ranges[fam]}")
    if params["is_real"] != (fam != "per"):
        probs.append("is_real flag inconsistent with the family")
    v = expected_value(params)
    if abs(v.real - correct["re"]) > 1e-9 or abs(v.imag - correct["im"]) > 1e-9:
        probs.append(f"X_d(w0) = {v}, answers ({correct['re']}, {correct['im']})")
    if fam != "per":                       # the Hermitian extension inverts to a real sequence
        w = np.linspace(-np.pi, np.pi, 200000)          # even count: no sample on a jump at w = 0
        f = F_NUM[fid][1]
        X = np.where(w >= 0, f(np.abs(w)), np.conj(f(np.abs(w))))
        for n in range(-3, 4):
            xn = np.trapezoid(X * np.exp(1j * w * n), w) / (2 * np.pi)
            if abs(xn.imag) > 1e-6:
                probs.append(f"Hermitian extension gives complex x[{n}] = {xn}")
                break
    if min(abs(abs(correct["re"]) - 2.244), abs(abs(correct["im"]) - 2.244)) < 0.01:
        probs.append("an answer equals the format example 2.244")
    n_wrong = sum(1 for o, e in submissions(params, correct) if "re" in o and 0 in (e.get("re"), e.get("im")))
    if n_wrong < 2:
        probs.append(f"only {n_wrong} classic Re/Im mistakes differ from the answer")
    # (b)
    w = np.linspace(-np.pi, np.pi, 2001)
    H = H_NUM[params["hid"]]
    real = bool(np.max(np.abs(H(-w) - np.conj(H(w)))) < 1e-9)
    if correct["hreal"]["html"].strip() != ("Yes" if real else "No"):
        probs.append(f"H_d {params['hid']}: Hermitian={real} but answer {correct['hreal']['html']}")
    # (c), (d)
    if fam != "per":
        f = F_NUM[fid][1]
        for w0 in (0.0, math.pi):
            if abs(complex(f(w0)).imag) > 1e-12:
                probs.append(f"real x but X_d({w0:.3f}) = {complex(f(w0))} is not real")
    for name, ans in (("s1", "ex1"), ("s2", "ex2")):
        s = params[name]
        nn = np.arange(-25, 26).astype(float)
        nn = nn[nn != 0]
        shown, own = tex_signal(s["tex"], nn), signal(s["kind"], s["args"], nn)
        if not np.allclose(shown, own, rtol=1e-9, atol=1e-12):
            probs.append(f"{name}: the displayed {s['tex']} is not the stored {s['kind']} {s['args']}")
        cat = classify(s["kind"], s["args"])
        want = params[name + "_choices"][cat]["text"].strip()
        if correct[ans]["html"].strip() != want:
            probs.append(f"{name} = {s['tex']}: classified {cat}, answer {correct[ans]['html']}")
    return probs


def _tol(v):
    return 1e-4 + 1e-3 * abs(v)


def submissions(params, correct):
    re_, im_ = correct["re"], correct["im"]
    cases = [({"re": f"{re_:.4g}", "im": f"{im_:.4g}"}, {"re": 1, "im": 1}),
             ({"re": f"{re_:.7f}", "im": f"{im_:.7f}"}, {"re": 1, "im": 1}),
             ({"re": "1/0"}, {"re": "invalid"})]
    f = F_NUM[params["fid"]][1]
    w0 = params["w0"][0] / params["w0"][1] * math.pi
    fam = params["family"]
    if fam == "herm":
        probes = [np.conj(f(w0)), f(-w0), f(w0), -f(-w0)]          # no sign flip / no conj / formula as is / odd
    elif fam == "per":
        red = w0 - 2 * math.pi * round(w0 / (2 * math.pi))
        probes = [f(w0), np.conj(f(-red)), f(-red)]        # forgot periodicity / symmetry used for a complex x
    else:
        probes = [f(2 * math.pi - w0), f(w0 - 2 * math.pi), np.conj(f(w0 - 2 * math.pi)), f(w0)]
    seen = 0
    for p in probes:
        p = complex(p)
        exp = {}
        for name, got, want in (("re", p.real, re_), ("im", p.imag, im_)):
            d = abs(got - want)
            if d > 2 * _tol(want):
                exp[name] = 0
            elif d < 0.4 * _tol(want):
                exp[name] = 1
            else:
                exp = None
                break
        if exp and 0 in exp.values():
            cases.append(({"re": f"{p.real:.6g}", "im": f"{p.imag:.6g}"}, exp))
            seen += 1
    # (b)
    right = correct["hreal"]["html"].strip()
    wrong_key = next(o["key"] for o in params["hreal"] if o["html"].strip() != right)
    cases.append(({"hreal": wrong_key}, {"hreal": 0}))
    cases.append(({"hreal": correct["hreal"]["key"]}, {"hreal": 1}))
    # (c), (d)
    for ans in ("ex1", "ex2"):
        for o in params[ans]:
            if o["key"] != correct[ans]["key"]:
                cases.append(({ans: o["key"]}, {ans: 0}))
                break
    return cases
