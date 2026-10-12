"""Independent check of questions/signals/complex-exponentials.

The displayed LaTeX of x1, x2, x3 is translated to Python and evaluated with cmath. (a) c is the
24-point DFT coefficient of x1 at the bin of w0 (the conjugate term falls in another bin);
(b) the phasor A e^{j theta} is twice the DFT coefficient of x2 at its frequency bin, found by
numpy.fft; (c) x3[k] is evaluated directly and |x3[n]| compared with the stored SymPy answer."""

import cmath
import math
import re
from fractions import Fraction

import numpy as np
from checklib import N, S as Sym


def _match_brace(s, i):
    depth = 0
    for k in range(i, len(s)):
        if s[k] == "{":
            depth += 1
        elif s[k] == "}":
            depth -= 1
            if depth == 0:
                return k
    raise ValueError(f"unbalanced braces in {s!r}")


def tex_to_py(s):
    s = s.replace(r"\left(", "(").replace(r"\right)", ")").replace(r"\!", "").replace(r"\,", "")
    s = s.replace(r"\cos", "C").replace(r"\sin", "S")
    while r"\tfrac{" in s:                                        # \tfrac{A}{B} -> ((A)/(B))
        i = s.index(r"\tfrac{")
        j = _match_brace(s, i + 6)
        k = _match_brace(s, j + 1)
        s = s[:i] + "((" + s[i + 7:j] + ")/(" + s[j + 2:k] + "))" + s[k + 1:]
    s = re.sub(r"e\^\{(-?)j(\d*)\\pi(?:/(\d+))?\}",
               lambda m: f"E({m.group(1)}J*{m.group(2) or 1}*P/{m.group(3) or 1})", s)
    s = re.sub(r"\\sqrt\{(\d+)\}", r"Q(\1)", s)
    s = s.replace(r"\pi", "P").replace("^{n}", "**(n)").replace("j", "J").replace(" ", "")
    assert "\\" not in s and "{" not in s, s
    s = re.sub(r"(?<=[\d)PJn])(?=[(SCEQPJn])", "*", s)
    s = re.sub(r"(?<=[)PJn])(?=\d)", "*", s)
    return s


ENV = {"C": cmath.cos, "S": cmath.sin, "E": cmath.exp, "Q": math.sqrt, "P": math.pi, "J": 1j}


def seq(tex):
    code = tex_to_py(tex)
    return lambda n: complex(eval(code, {"__builtins__": {}}, dict(ENV, n=n)))


def coef(tex, half=True):
    """(bin m, X[m]/24) of the 24-point DFT of a real sinusoid at its positive-frequency bin."""
    x = seq(tex)
    v = np.array([x(n) for n in range(24)])
    if np.max(np.abs(v.imag)) > 1e-9:
        raise ValueError("sequence is not real")
    X = np.fft.fft(v.real) / 24
    m = int(np.argmax(np.abs(X[1:12]))) + 1
    others = [abs(X[i]) for i in range(24) if i not in (m, 24 - m)]
    if max(others) > 1e-9:
        raise ValueError(f"not a single sinusoid: {np.round(np.abs(X), 6)}")
    return m, complex(X[m])


def ang_pi(z):
    a = cmath.phase(z) / math.pi
    return 1.0 if a <= -1 + 1e-12 else a


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(b))


def check(params, correct):
    probs = []
    # (a)
    m, c = coef(params["xa_tex"])
    w0 = eval(tex_to_py(params["w0_tex"]), {"__builtins__": {}}, ENV)
    if not close(2 * math.pi * m / 24, w0):
        probs.append(f"(a) x1 oscillates at bin {m}, not at w0 = {w0}")
    if not close(abs(c), correct["cmag"]) or not close(ang_pi(c), correct["cang"]):
        probs.append(f"(a) c = {c} (|c| {abs(c)}, angle {ang_pi(c)}) vs {correct['cmag']}, {correct['cang']}")
    if abs(correct["cang"]) < 1e-9 or abs(abs(correct["cang"]) - 1) < 1e-9:
        probs.append("(a) angle 0 or +-pi")
    # (b)
    m2, P = coef(params["xb_tex"])
    P *= 2
    if not close(abs(P), correct["amp"]) or not close(ang_pi(P), correct["phase"]):
        probs.append(f"(b) phasor {P} (A {abs(P)}, theta {ang_pi(P)}) vs {correct['amp']}, {correct['phase']}")
    if abs(correct["phase"]) < 1e-9 or abs(abs(correct["phase"]) - 1) < 1e-9:
        probs.append("(b) phase 0 or +-pi")
    th = Fraction(correct["phase"]).limit_denominator(12)
    if abs(float(th) - correct["phase"]) > 1e-12:
        probs.append(f"(b) theta/pi = {correct['phase']} is not a multiple of 1/12")
    if A_sq := round(correct["amp"] ** 2, 9):
        if abs(A_sq - round(A_sq)) > 1e-9:
            probs.append(f"(b) A^2 = {A_sq} is not an integer")
    if abs(correct["amp"] - sum(t[0] for t in params["b_terms"])) < 1e-6:
        probs.append("(b) A = A1 + A2 (the amplitude-sum mistake would be right)")
    # (c)
    x3 = seq(params["xc_tex"])
    k = params["c_k"]
    v = x3(k)
    if not close(v.real, correct["xre"]) or not close(v.imag, correct["xim"]):
        probs.append(f"(c) x3[{k}] = {v} vs {correct['xre']} + j{correct['xim']}")
    for val in (correct["xre"], correct["xim"]):
        if abs(float(Fraction(val).limit_denominator(64)) - val) > 1e-12:
            probs.append(f"(c) {val} is not a simple fraction")
    expr = Sym(correct["xmag"]["_value"])
    for n in range(-3, 6):
        got = complex(expr.subs(N, n).evalf(30))
        if not close(got.real, abs(x3(n)), 1e-9) or abs(got.imag) > 1e-12:
            probs.append(f"(c) |x3[{n}]| = {abs(x3(n))} but the answer gives {got}")
            break
    if abs(abs(x3(1)) - abs(x3(0))) < 1e-9:
        probs.append("(c) |x3[n]| is constant")
    return probs


# ----------------------------------------------------------------------------- probes
def _num(x):
    return f"{float(x):.6g}"


def _frac(x):
    q = Fraction(x).limit_denominator(64)
    return f"{q.numerator}/{q.denominator}" if abs(float(q) - x) < 1e-12 else _num(x)


def _wrong(cands, right):
    out = []
    for m in cands:
        if abs(m - right) > 2e-2 * max(1.0, abs(right)) and all(abs(m - o) > 1e-9 for o in out):
            out.append(m)
    return out


def _r_strings(r):
    """Two equivalent spellings of r^n, plus the exact string of r (r rational or q*sqrt(d))."""
    q = Fraction(r).limit_denominator(64)
    if abs(float(q) - r) < 1e-12:
        return [f"{q.numerator}^n/{q.denominator}^n", f"({q.denominator}/{q.numerator})^(-n)"], f"({q.numerator}/{q.denominator})"
    for d in (2, 3):
        q = Fraction(r / math.sqrt(d)).limit_denominator(64)
        if abs(float(q) * math.sqrt(d) - r) < 1e-12:
            qs = f"({q.numerator}/{q.denominator})"
            return [f"{qs}^n*{d}^(n/2)", f"{qs}^n*sqrt({d})^n"], f"({qs}*sqrt({d}))"
    raise ValueError(f"unexpected magnitude {r}")


def submissions(params, correct):
    cases = []
    # (a)
    A, fn = params["a_A"], params["a_fn"]
    th = params["a_th"][0] / params["a_th"][1]
    cm, ca = correct["cmag"], correct["cang"]
    cases += [({"cmag": _frac(cm), "cang": params["a_ang_frac"]}, {"cmag": 1, "cang": 1}),
              ({"cmag": f"{cm:.4g}", "cang": f"{ca:.4g}"}, {"cmag": 1, "cang": 1})]
    for m in _wrong([abs(A), abs(A) / 4], cm):                       # c = A (no 1/2); 1/2 applied twice
        cases.append(({"cmag": _num(m)}, {"cmag": 0}))
    if A < 0:
        cases.append(({"cmag": _num(A / 2)}, {"cmag": 0}))           # negative "magnitude"

    def pr(t):
        t = (t + 1) % 2 - 1
        return 1.0 if t <= -1 + 1e-12 else t

    wrong_a = [-ca, pr(th)]                                          # conjugate coefficient; phase taken as is
    if fn == "sin":
        wrong_a.append(pr(th + 0.5 + (1 if A < 0 else 0)))           # 1/j read as +j
    raw = params["a_raw"][0] / params["a_raw"][1]
    wrong_a.append(raw)                                              # not reduced to (-pi, pi]
    for m in _wrong(wrong_a, ca):
        cases.append(({"cang": _num(m)}, {"cang": 0}))
    # (b)
    Aa, ph = correct["amp"], correct["phase"]
    terms = params["b_terms"]
    cases += [({"amp": f"{Aa:.4g}", "phase": params["b_th_frac"]}, {"amp": 1, "phase": 1}),
              ({"amp": _num(Aa), "phase": f"{ph:.4g}"}, {"amp": 1, "phase": 1})]
    for m in _wrong([terms[0][0] + terms[1][0], abs(terms[0][0] - terms[1][0])], Aa):
        cases.append(({"amp": _num(m)}, {"amp": 0}))
    th1 = terms[0][2] / terms[0][3]
    wrong_b = [-ph, th1]
    if terms[1][1] == "sin":                                         # sine treated as a cosine
        Pw = terms[0][0] * cmath.exp(1j * math.pi * th1) + terms[1][0] * cmath.exp(1j * math.pi * terms[1][2] / terms[1][3])
        if abs(Pw) > 1e-9:
            wrong_b.append(ang_pi(Pw))
    for m in _wrong(wrong_b, ph):
        cases.append(({"phase": _num(m)}, {"phase": 0}))
    # (c)
    x3 = seq(params["xc_tex"])
    r = abs(x3(1)) / abs(x3(0))
    Bm = abs(x3(0))
    Bs = "" if abs(Bm - 1) < 1e-12 else f"{round(Bm)}*"
    eqs, rstr = _r_strings(r)
    cases += [({"xmag": Bs + e}, {"xmag": 1}) for e in eqs]
    bad = [f"{Bs}{rstr}^(2*n)", f"{Bs}{rstr}^(-n)", f"{Bs}{rstr}"]
    if Bs:
        bad.append(f"{rstr}^n")
    if params["c_B"] < 0:
        bad.append(f"-{Bs}{rstr}^n")
    cases += [({"xmag": b}, {"xmag": 0}) for b in bad]
    cases.append(({"xmag": "0.5^n"}, {"xmag": "invalid"}))
    re_, im_ = correct["xre"], correct["xim"]
    cases += [({"xre": _frac(re_), "xim": _frac(im_)}, {"xre": 1, "xim": 1}),
              ({"xre": f"{re_:.4g}", "xim": f"{im_:.4g}"}, {"xre": 1, "xim": 1})]
    k = params["c_k"]
    z = x3(k)
    a1 = x3(1) / x3(0)                                               # the base a
    cands_re, cands_im = [], []
    conj = z.conjugate()                                             # conjugate (angle sign) slip
    cands_im.append(conj.imag)
    noB = z / x3(0)                                                  # forgot B
    cands_re.append(noB.real)
    cands_im.append(noB.imag)
    rk_only = x3(0) * abs(a1) * cmath.exp(1j * k * cmath.phase(a1))  # r instead of r^k
    cands_re.append(rk_only.real)
    cands_im.append(rk_only.imag)
    swap = complex(z.imag, z.real)                                   # Re and Im swapped
    cands_re.append(swap.real)
    cands_im.append(swap.imag)
    for m in _wrong(cands_re, re_):
        cases.append(({"xre": _num(m)}, {"xre": 0}))
    for m in _wrong(cands_im, im_):
        cases.append(({"xim": _num(m)}, {"xim": 0}))
    return cases
