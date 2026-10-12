"""Independent check of questions/stability/pole-zero-plot.

The poles, zeros and multiplicities are read back from the SVG the student sees (marker pixel
positions -> z-plane coordinates via the unit circle's radius) and compared with the stored ones.
H(z) = K prod(z - z_k)/prod(z - p_k) is rebuilt with numpy; the gain condition, every generated
SymPy form, the ROC options (all annuli between distinct pole magnitudes), the causal / stable
ROCs and the causality verdicts are recomputed from numpy magnitudes."""

import cmath
import math
import re
from fractions import Fraction

import numpy as np
import sympy
from checklib import S, Z
from scipy.signal import lfilter

PTS = [0.7 * cmath.exp(0.3j), 1.21 * cmath.exp(-0.77j), 2.3 * cmath.exp(1.9j), 0.4 * cmath.exp(-2.8j)]


def fr(s):
    return Fraction(s)


def roots(params, key):
    return [complex(float(fr(a)), float(fr(b))) for a, b in params[key]]


def H(params, z):
    K = float(fr(params["K"]))
    return K * np.prod([z - q for q in roots(params, "zeros")]) / np.prod([z - q for q in roots(params, "poles")])


def svg_points(html):
    """(poles, zeros) as lists of complex numbers with multiplicity, read from the SVG."""
    w = float(re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', html).group(1))
    circ = re.search(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="([\d.]+)" fill="none" stroke="currentColor"', html)
    cx0, cy0, rpx = map(float, circ.groups())

    def val(px, py):
        return complex((px - cx0) / rpx, (cy0 - py) / rpx)
    poles = [val(float(a) + 7, float(b) + 7) for a, b in re.findall(r'<path d="M([\d.]+),([\d.]+) L[\d.]+,[\d.]+ M', html)]
    zeros = [val(float(a), float(b)) for a, b in re.findall(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="7.0" fill="none"', html)]
    mult = [(val(float(a) - 9, float(b) + 7), int(m), col)
            for a, b, col, m in re.findall(r'<text x="([\d.]+)" y="([\d.]+)" font-size="12" fill="(#\w+)">(\d)</text>', html)]
    for loc, m, col in mult:
        lst = poles if col == "#d1501f" else zeros
        i = min(range(len(lst)), key=lambda k: abs(lst[k] - loc))
        lst += [lst[i]] * (m - 1)
    return poles, zeros, w


def same_multiset(a, b, tol=0.02):
    b = list(b)
    for x in a:
        i = next((k for k, y in enumerate(b) if abs(x - y) < tol), None)
        if i is None:
            return False
        b.pop(i)
    return not b


def parse_r(t):
    """\\tfrac{1}{2}, 2, \\sqrt{2}, \\tfrac{\\sqrt{2}}{2}, \\tfrac{3}{2} -> float."""
    t = t.strip()
    m = re.fullmatch(r"\\tfrac\{(.+)\}\{(\d+)\}", t)
    if m:
        return parse_r(m.group(1)) / int(m.group(2))
    m = re.fullmatch(r"(\d*)\\sqrt\{2\}", t)
    if m:
        return (int(m.group(1)) if m.group(1) else 1) * math.sqrt(2)
    return float(int(t))


def parse_roc(html):
    s = html.replace("$", "").strip()
    z = r"\lvert z\rvert"
    if s.startswith(z + r" \gt "):
        return parse_r(s[len(z + r" \gt "):]), math.inf
    if s.startswith(z + r" \lt "):
        return 0.0, parse_r(s[len(z + r" \lt "):])
    lo, mid, hi = s.split(r" \lt ")
    assert mid == z, html
    return parse_r(lo), parse_r(hi)


def parse_num_tex(t):
    t = t.strip()
    sign = -1 if t.startswith("-") else 1
    t = t.lstrip("-")
    m = re.fullmatch(r"\\tfrac\{(\d+)\}\{(\d+)\}", t)
    return sign * (Fraction(int(m.group(1)), int(m.group(2))) if m else Fraction(int(t)))


def close(r1, r2):
    return all((x == y) if math.inf in (x, y) else abs(x - y) < 1e-9 for x, y in zip(r1, r2))


def mags(params):
    return sorted({round(abs(p), 9) for p in roots(params, "poles")})


def expected(params):
    ms = mags(params)
    rocs = [(0.0, ms[0])] + list(zip(ms, ms[1:])) + [(ms[-1], math.inf)]
    on = any(abs(m - 1) < 1e-9 for m in ms)
    stable = None if on else next(r for r in rocs if r[0] < 1 < r[1])
    return rocs, stable


def check(params, correct):
    probs = []
    P, Zs = roots(params, "poles"), roots(params, "zeros")
    sp, sz, _ = svg_points(params["fig_html"])
    if not same_multiset(sp, P) or not same_multiset(sz, Zs):
        probs.append(f"figure shows poles {sp} zeros {sz}, params say {P} {Zs}")
    if not len(Zs) <= len(P) or any(abs(p) < 1e-12 for p in P):
        probs.append("H must be proper in z with no pole at the origin")
    if any(abs(p - q) < 1e-9 for p in P for q in Zs):
        probs.append("a zero sits on a pole (cancellation)")
    if not np.allclose(np.poly(P).imag, 0) or not np.allclose(np.poly(Zs).imag, 0):
        probs.append("poles/zeros are not closed under conjugation")
    K = fr(params["K"])
    if K == 1:
        probs.append("K = 1: the missing-gain probe is meaningless")
    if params["norm"] == "b0":
        if len(Zs) != len(P):
            probs.append("b0 normalization needs as many finite zeros as poles")
        b0 = H(params, 1e7)                                        # = numerator constant term in z^-1 form
        if abs(b0 - float(parse_num_tex(params["K_tex"]))) > 1e-5:
            probs.append(f"b0 = {b0} but the text says {params['K_tex']}")
    else:
        z0 = 1.0 if params["norm"] == "dc" else -1.0
        c = float(parse_num_tex(params["c_tex"]))
        if abs(H(params, z0) - c) > 1e-9 * max(1, abs(c)):
            probs.append(f"H({z0}) = {H(params, z0)} but the text says {c}")
    for key in ("H_sym", "H_pos_sym", "H_cplx_sym"):
        e = S(params[key])
        for zz in PTS:
            got, want = complex(e.subs(Z, zz).evalf(30)), H(params, zz)
            if abs(got - want) > 1e-9 * max(1, abs(want)):
                probs.append(f"{key} = {params[key]} gives {got} at {zz}, H = {want}")
                break
    rocs, stable = expected(params)
    for name in ("roc_c", "roc_s"):
        opts = [o["html"] for o in params[name] if "none" not in o["html"]]
        got = sorted(parse_roc(o) for o in opts)
        if len(got) != len(rocs) or not all(close(a, b) for a, b in zip(got, rocs)):
            probs.append(f"{name} options {got} != possible ROCs {rocs}")
    lo, hi = parse_roc(correct["roc_c"]["html"])
    if hi != math.inf or abs(lo - rocs[-1][0]) > 1e-9:
        probs.append(f"causal ROC marked {correct['roc_c']['html']}")
    rmax = rocs[-1][0]
    if (correct["stab_c"]["html"].strip() == "Yes") != (rmax < 1 - 1e-9):
        probs.append(f"causal stable? marked {correct['stab_c']['html']} with r_max = {rmax}")
    marked = correct["roc_s"]["html"]
    if stable is None:
        if "none" not in marked:
            probs.append(f"stable ROC should be none, marked {marked}")
    else:
        lo, hi = parse_roc(marked) if "none" not in marked else (None, None)
        if lo is None or not close((lo, hi), stable):
            probs.append(f"stable ROC {stable}, marked {marked}")
    k = correct["kind_s"]["html"]
    want = ("no BIBO" if stable is None else "right-sided" if stable[1] == math.inf
            else "anti-causal" if stable[0] == 0 else "two-sided")
    if want not in k:
        probs.append(f"kind of the stable system: want {want!r}, marked {k!r}")
    if want == "anti-causal" and abs(H(params, 1e-9)) > 1e-6:
        probs.append("stable ROC |z| < r_min but H(0) != 0: h[0] != 0, so the system is not anti-causal")
    # (c) wording: "marginally stable" only if every unit-circle pole is simple (outermost pole on the circle)
    on = [q for q in P if abs(abs(q) - 1) < 1e-9]
    rmax_on = abs(max(abs(q) for q in P) - 1) < 1e-9
    repeated = any(sum(1 for r in on if abs(r - q) < 1e-9) > 1 for q in on)
    if params["marginal"] != (rmax_on and not repeated) or params["on_repeated"] != (rmax_on and repeated):
        probs.append(f"marginal-stability flags wrong: marginal={params['marginal']}, on_repeated={params['on_repeated']}")
    if rmax_on:                                           # what the sentence claims about the causal h[n]
        b = [0.0] * (len(P) - len(Zs)) + list(float(fr(params["K"])) * np.real(np.poly(Zs) if Zs else np.array([1.0])))
        a = np.real(np.poly(P))
        x = np.zeros(800)
        x[0] = 1.0
        h = np.abs(lfilter(b, a, x))
        grows = h[400:].max() > 1.5 * h[:400].max()
        if grows != repeated:
            probs.append(f"causal h[n] {'grows' if grows else 'stays bounded'} but on_repeated = {repeated}")
    return probs


def _lin(roots_, sign=-1):
    """prod (z + sign*q) over the roots, typed the way a student would."""
    out = []
    for q in roots_:
        a, b = Fraction(q.real).limit_denominator(16), Fraction(q.imag).limit_denominator(16)
        if a == 0 and b == 0:
            out.append("z")
        elif b == 0:
            out.append(f"(z + ({sign * a}))")
        else:
            out.append(f"(z + ({sign})*(({a}) + ({b})*j))")
    return "*".join(out) or "1"


def submissions(params, correct):
    P, Zs = roots(params, "poles"), roots(params, "zeros")
    K = fr(params["K"])
    d = len(P) - len(Zs)
    num = np.real(np.poly(Zs)) if Zs else np.array([1.0])
    den = np.real(np.poly(P))

    def frs(c):
        return [Fraction(v).limit_denominator(64) for v in c]
    nz = frs(num)
    dz = frs(den)
    # z^-1 form, multiplied out, written with ^ and 1/z^k: K z^{-d} N(z^-1) / D(z^-1)
    n_zinv = " + ".join(f"({c})*z^(-{i + d})" if i + d else f"({c})" for i, c in enumerate(nz) if c)
    d_zinv = " + ".join(f"({c})*z^(-{i})" if i else f"({c})" for i, c in enumerate(dz) if c)
    cases = [
        ({"H": params["H_sym"]}, {"H": 1}),
        ({"H": params["H_pos_sym"]}, {"H": 1}),
        ({"H": params["H_cplx_sym"]}, {"H": 1}),                                   # complex factors with j
        ({"H": f"({K})*({n_zinv})/({d_zinv})"}, {"H": 1}),
        ({"H": params["H_swap_sym"]}, {"H": 0}),                                   # poles and zeros swapped
        ({"H": f"({_lin(Zs)})/({_lin(P)})"}, {"H": 0}),                            # gain K missing
        ({"H": "x + 1"}, {"H": "invalid"}),
    ]
    mirrored = f"({K})*({_lin(Zs, sign=1)})/({_lin(P, sign=1)})"                    # z + p instead of z - p
    zz = PTS[0]
    e = S(mirrored)
    if abs(complex(e.subs(Z, zz).evalf(20)) - H(params, zz)) > 1e-6:
        cases.append(({"H": mirrored}, {"H": 0}))
    if d:                                                                          # forgot the z^{-d}
        n0 = " + ".join(f"({c})*z^(-{i})" if i else f"({c})" for i, c in enumerate(nz) if c)
        cases.append(({"H": f"({K})*({n0})/({d_zinv})"}, {"H": 0}))
    for name in ("roc_c", "stab_c", "roc_s", "kind_s"):
        for o in params[name]:
            cases.append(({name: o["key"]}, {name: 1 if o["key"] == correct[name]["key"] else 0}))
    return cases
