"""Independent check of questions/dtft/real-valued-test.

(a) Every displayed option (the element's own html) is parsed from its LaTeX into a numpy function
    by this file's small converter and tested on a frequency grid inside (-pi, pi):
    x real <=> X_d(-w) == X_d(w)^*. As a second test the inverse DTFT x[n], |n| <= 4, is computed by
    numerical integration and must come out real exactly for the options marked correct. The
    solution's displayed X_d(-w) and X_d^*(w) are parsed too and compared with the option.
(b) A real two-sample sequence x = alpha delta[n] + beta delta[n-1] with X_d(w1) = a + jb is built
    here; X_d(-w1), X_d(w2) and the phase come from direct sums / np.angle."""

import math
import re

import numpy as np

FUNCS = {"cos", "sin", "exp", "ABS"}


def _match(s, i):
    depth = 0
    for k in range(i, len(s)):
        depth += {"{": 1, "}": -1}.get(s[k], 0)
        if depth == 0:
            return k
    raise ValueError(s)


def tex_to_py(tex):
    s = tex.replace(r"\left(", "(").replace(r"\right)", ")").replace(r"\,", " ").replace(r"\!", "")
    for fr in (r"\dfrac{", r"\tfrac{"):
        while fr in s:
            i = s.index(fr)
            j = _match(s, i + len(fr) - 1)
            k = _match(s, j + 1)
            s = s[:i] + "((" + s[i + len(fr):j] + ")/(" + s[j + 2:k] + "))" + s[k + 1:]
    while "e^{" in s:
        i = s.index("e^{")
        j = _match(s, i + 2)
        s = s[:i] + "exp(" + s[i + 3:j] + ")" + s[j + 1:]
    s = s.replace(r"|\omega|", "ABS(w)").replace(r"\omega^2", "(w**2)").replace(r"\omega", "w")
    s = re.sub(r"\\(cos|sin)\s*(\d*)\s*w", lambda m: f"{m.group(1)}({m.group(2) or 1}*w)", s)
    if "\\" in s or "{" in s or "|" in s:
        raise ValueError(f"unparsed LaTeX left in {s!r}")
    toks = re.findall(r"\d+|cos|sin|exp|ABS|[A-Za-z]|\*\*|[()+\-*/]", s)
    if "".join(toks) != s.replace(" ", ""):
        raise ValueError(f"untokenized characters in {s!r}")
    out = []
    for t in toks:
        if out:
            prev = out[-1]
            left = prev == ")" or prev[0].isdigit() or (prev[0].isalpha() and prev not in FUNCS)
            right = t == "(" or t[0].isdigit() or t[0].isalpha()
            if left and right:
                out.append("*")
        out.append(t)
    py = "".join("1j" if t == "j" else t for t in out)
    return py


def tex_eval(tex, w):
    env = {"cos": np.cos, "sin": np.sin, "exp": np.exp, "ABS": np.abs, "w": w}
    return np.asarray(eval(tex_to_py(tex), {"__builtins__": {}}, env), dtype=complex) * np.ones_like(w, dtype=complex)


def _formula(html):
    m = re.fullmatch(r"\s*\$X_d\(\\omega\) = (.*)\$\s*", html, re.S)
    if not m:
        raise ValueError(html)
    return m.group(1)


def is_real(tex):
    w = np.linspace(-np.pi, np.pi, 4001)[1:-1]
    return bool(np.max(np.abs(tex_eval(tex, -w) - np.conj(tex_eval(tex, w)))) < 1e-9)


def inverse_real(tex):
    """x[n] for |n| <= 4 by the midpoint rule; True if all imaginary parts vanish."""
    M = 20000
    w = -np.pi + (np.arange(M) + 0.5) * 2 * np.pi / M
    X = tex_eval(tex, w)
    xs = [np.mean(X * np.exp(1j * w * n)) for n in range(-4, 5)]
    return max(abs(v.imag) for v in xs) < 1e-6


def check(params, correct):
    probs = []
    opts = params["real"]
    ck = {o["key"] for o in correct["real"]}
    if len(opts) != 5 or not (1 <= len(ck) <= 4):
        probs.append(f"{len(opts)} options, {len(ck)} correct")
    rows = {r["X_tex"]: r for r in params["sol_rows"]}
    for o in opts:
        tex = _formula(o["html"])
        r1, r2 = is_real(tex), inverse_real(tex)
        if r1 != r2:
            probs.append(f"{tex}: grid test {r1} but inverse DTFT test {r2}")
        if r1 != (o["key"] in ck):
            probs.append(f"{tex}: real={r1} but marked {o['key'] in ck}")
        row = rows.get(tex)
        if row is None:
            probs.append(f"no solution row for {tex}")
            continue
        w = np.linspace(-3.0, 3.0, 61)
        if np.max(np.abs(tex_eval(row["Xm_tex"], w) - tex_eval(tex, -w))) > 1e-9:
            probs.append(f"solution X_d(-w) = {row['Xm_tex']} wrong for {tex}")
        if np.max(np.abs(tex_eval(row["Xc_tex"], w) - np.conj(tex_eval(tex, w)))) > 1e-9:
            probs.append(f"solution X_d^*(w) = {row['Xc_tex']} wrong for {tex}")
        if row["verdict"] != ("real" if r1 else "not real"):
            probs.append(f"solution verdict for {tex}")
    # (b)
    a, b = params["a"], params["b"]
    w1 = params["w1"][0] / params["w1"][1] * math.pi
    w2 = params["w2"][0] / params["w2"][1] * math.pi
    beta = -b / math.sin(w1)
    alpha = a - beta * math.cos(w1)
    X = lambda w: alpha + beta * np.exp(-1j * w)  # noqa: E731
    if abs(X(w1) - complex(a, b)) > 1e-12:
        probs.append("construction failed")
    Xm, X2 = X(-w1), X(w2)
    for name, v in (("rem", Xm.real), ("imm", Xm.imag), ("re2", X2.real), ("im2", X2.imag)):
        if abs(v - correct[name]) > 1e-9:
            probs.append(f"{name}: {correct[name]} != {v}")
    if abs(np.angle(Xm) - correct["ph"]) > 1e-12:
        probs.append(f"phase {correct['ph']} != {np.angle(Xm)}")
    if a == 0 or b == 0 or not (0 < w1 < math.pi) or not abs(w2) > math.pi:
        probs.append("(b) design violated")
    sh = re.findall(r"\\approx (-?[0-9.]+)", params["ph_line_tex"])
    if not sh or abs(float(sh[-1]) - correct["ph"]) > 1e-4:
        probs.append(f"solution phase {sh}")
    return probs


def submissions(params, correct):
    opts = params["real"]
    ck = [o["key"] for o in correct["real"]]
    wrong = [o["key"] for o in opts if o["key"] not in ck]
    c = len(ck)
    cases = [
        ({"real": list(reversed(ck))}, {"real": 1}),
        ({"real": ck + wrong[:1]}, {"real": (c - 1) / c}),
        ({"real": wrong}, {"real": 0}),
        ({"real": []}, {"real": "invalid"}),
    ]
    if c >= 2:
        cases.append(({"real": ck[:1]}, {"real": 1 / c}))
    a, b = params["a"], params["b"]
    ph = correct["ph"]
    cases += [
        ({"rem": f"{a}.0", "imm": f"{2 * -b}/2"}, {"rem": 1, "imm": 1}),
        ({"imm": str(b)}, {"imm": 0}),                               # forgot the conjugate
        ({"rem": str(-a)}, {"rem": 0}),                              # treated X_d as odd
        ({"ph": f"{ph:.4g}"}, {"ph": 1}),
        ({"ph": f"{ph:.6f}"}, {"ph": 1}),
        ({"ph": f"{-ph:.6g}"}, {"ph": 0}),                           # phase of X_d(w1)
        ({"ph": f"{ph + (2 * math.pi if ph < 0 else -2 * math.pi):.6g}"}, {"ph": 0}),
        ({"ph": f"{math.degrees(ph):.6g}"}, {"ph": 0}),
        ({"ph": "x"}, {"ph": "invalid"}),
    ]
    if a < 0:
        cases.append(({"ph": f"{math.atan(-b / a):.6g}"}, {"ph": 0}))  # no quadrant correction
    if correct["im2"] == -b:
        cases.append(({"im2": str(b)}, {"im2": 0}))                  # X_d(2pi - w1) read as X_d(w1)
    else:
        cases.append(({"im2": str(-b)}, {"im2": 0}))                 # conjugated without need
    return cases
