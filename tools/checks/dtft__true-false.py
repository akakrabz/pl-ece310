"""Independent check of questions/dtft/true-false.

Fixed statements: answer table typed from the official keys (study-site Midterm 2 True/False bank,
Q1-Q16) and from the Lecture 13-16 notes, with the source of each entry in a comment.
Parameterized statements: recomputed numerically from the stored parameters with numpy/scipy
(partial sums for summability, truncated DTFT sums, numerical inverse DTFTs by quadrature,
convolutions, scipy.signal.group_delay, incremental phase on a fine grid) -- never through
server.py's code path."""

import re

import numpy as np
from scipy.signal import group_delay as sp_group_delay

OFFICIAL = {   # sid -> answer
    "ztrans_always": False,          # FA2023 MT2 #1(a), bank Q1 (key: False, 2^n u[n])
    "finite_periodic": True,         # FA2024 MT2 #1(b), bank Q3 (key: True)
    "periodic_only_infinite": False,  # FA2019 MT2 #1(e), bank Q4 (key: False)
    "uniqueness": False,             # Lecture 13 notes section 4: the DTFT of an absolutely summable sequence is unique
    "real_hermitian": True,          # FA2023 MT2 #1(c), bank Q9 (key: True)
    "real_X0_Xpi": True,             # Lecture 14 notes: X_d(0), X_d(pi) real for real x (Hermitian symmetry)
    "real_Xd_real": False,           # Lecture 14: real x gives Hermitian, not real, X_d (delta[n-1] -> e^{-jw})
    "Xd_real_x_real": False,         # Lecture 14: j delta[n+1] - j delta[n-1] -> -2 sin w (real), x imaginary
    "conv_product": True,            # Lecture 14 notes Table 2 (convolution)
    "product_product": False,        # Lecture 14 notes Table 2 (windowing: periodic convolution / 2pi)
    "zp_can_change": False,          # FA2019 MT2 #1(a), bank Q12 (key: False)
    "zp_equal": True,                # FA2023 MT2 #1(f), bank Q14; SP2023 MT2 #1(c), Q13 (keys: True)
    "cos_eigen_every": False,        # FA2021 MT2 #1(a), bank Q11 (key: False), reworded "every"
    "exp_eigen": True,               # Lecture 15 notes: e^{jw0 n} -> H_d(w0) e^{jw0 n}
    "real_cos_formula": False,       # FA2019 MT2 #1(d), bank Q10 (key: False: h must be real)
}

RNG = np.random.default_rng(310)


def f(pq):
    return pq[0] / pq[1]


def dtft(x, n0, w):
    n = n0 + np.arange(len(x))
    return np.exp(-1j * np.outer(np.atleast_1d(w), n)) @ np.asarray(x, dtype=complex)


def summable(vals_of_k, N=400):
    """Is sum_{k>=0} |v(k)| finite? (partial sums over N/2 and N terms agree)."""
    k = np.arange(N, dtype=float)
    with np.errstate(all="ignore"):
        a = np.abs(vals_of_k(k))
    s1, s2 = a[: N // 2].sum(), a.sum()
    return bool(np.isfinite(s2) and s2 - s1 < 1e-6 * max(1.0, s1))


def inv_dtft(Xfun, n, M=200001):
    """x[n] = (1/2pi) int_{-pi}^{pi} X(w) e^{jwn} dw by the trapezoid rule on a symmetric grid."""
    w = np.linspace(-np.pi, np.pi, M)
    X = Xfun(w)
    return np.array([np.trapezoid(X * np.exp(1j * w * k), w) / (2 * np.pi) for k in n])


# ----------------------------------------------------------------------------- recomputation
XD = {   # independent numpy versions of the DTFT formulas used by real_from_Xd
    "w2_ecos": lambda w: w ** 2 * np.exp(1j * np.cos(w)),
    "w2_e2w": lambda w: w ** 2 * np.exp(-2j * w),
    "jw": lambda w: 1j * w,
    "w": lambda w: w + 0j,
    "absw_ew2": lambda w: np.abs(w) * np.exp(1j * w ** 2),
    "twocos_ew": lambda w: (2 + np.cos(w)) * np.exp(-1j * w),
    "twocos_eabs": lambda w: (2 + np.cos(w)) * np.exp(1j * np.abs(w)),
    "cos2_esin": lambda w: np.cos(w) ** 2 * np.exp(1j * np.sin(w)),
}

SYS_H = {   # frequency responses of sinusoid_formula, valid on [-pi, pi]
    "absw_esin": lambda w: np.abs(w) * np.exp(-1j * np.pi * np.sin(w)),
    "cos2_esin": lambda w: np.cos(w) ** 2 * np.exp(1j * np.sin(w)),
    "jw_esin": lambda w: 1j * w * np.exp(1j * np.pi * np.sin(w)),
    "w_ecos": lambda w: w * np.exp(1j * np.pi * np.cos(w)),
    "half_eabs": lambda w: (0.5 + np.cos(w)) * np.exp(1j * np.abs(w)),
    "w_esin": lambda w: w * np.exp(1j * np.sin(w)),
}

PHASE_H = {   # frequency responses of phase_jump
    "cos4": lambda w: np.cos(4 * w) + 0j,
    "H2": lambda w: np.exp(-2j * w) * (2 + np.cos(2 * w)),
    "H3": lambda w: 1j * np.exp(-3j * w) * np.sin(3 * w),
    "H4": lambda w: np.exp(-2j * w) * (1 - 2 * np.cos(w) - np.cos(2 * w)),
    "H5": lambda w: np.exp(-2j * w) * (3 - 2 * np.cos(w) - np.cos(2 * w)),
    "sp2023": lambda w: 2 * np.exp(-3j * w) * np.cos(3 * w),
    "sp2025": lambda w: np.exp(-3j * w) * (2 + 2 * np.cos(3 * w)),
    "delay": lambda w: np.exp(-3j * w),
}


def recompute(sid, prm, probs):
    if sid in OFFICIAL:
        return OFFICIAL[sid]
    W = lambda key: np.pi * f(prm[key])                      # noqa: E731
    if sid == "roc_substitute":                               # valid iff sum |x[n]| < inf
        a = f(prm["a"])
        if prm["side"] == "right":                            # a^n u[n]
            return summable(lambda k: a ** k)
        return summable(lambda k: a ** (-(k + 1)))            # -a^n u[-n-1]: n = -1, -2, ...
    if sid == "abs_summable":
        a = f(prm["a"])
        if prm["form"] == "right":
            return summable(lambda k: a ** k)
        if prm["form"] == "anti":                             # a^n u[-n]: n = 0, -1, -2, ...
            return summable(lambda k: a ** (-k))
        return summable(lambda k: a ** k) and summable(lambda k: a ** (-(-k)))
    if sid == "bandlimited_zero":                             # a triangle spectrum of half-width wc, then test the interval
        wc = W("wc")
        Mx = 4000
        n = np.arange(-Mx, Mx + 1)
        with np.errstate(all="ignore"):
            x = np.where(n == 0, wc / (2 * np.pi), (1 - np.cos(wc * n)) / (np.pi * wc * n.astype(float) ** 2))
        lo = np.pi * f(prm["lo"])
        hi = 6 * np.pi if prm["hi"] is None else np.pi * f(prm["hi"])
        w = lo + (hi - lo) * np.linspace(0.01, 0.99, 197)
        Xw = np.abs(np.exp(-1j * np.outer(w, n)) @ x)
        if np.abs(dtft(x, -Mx, 0.0)[0] - 1) > 1e-2:
            probs.append("bandlimited_zero: test spectrum is not 1 at w = 0")
        return bool(Xw.max() < 0.02)
    if sid == "finite_cos_impulses":                          # invert the claimed impulses by sifting
        w0, N = W("w0"), prm["N"]
        n = np.arange(-6, N + 6)
        claim = (np.pi * np.exp(1j * w0 * n) + np.pi * np.exp(-1j * w0 * n)) / (2 * np.pi)
        x = np.cos(w0 * n) if prm["kind"] == "all" else np.where((n >= 0) & (n < N), np.cos(w0 * n), 0.0)
        return bool(np.allclose(claim, x, atol=1e-12))
    if sid == "impulse_pair":
        w0, A = W("w0"), prm["A"]
        n = np.arange(-12, 13)
        claim = sum((re + 1j * im) * np.pi * np.exp(1j * np.pi * f(loc) * n) for re, im, loc in prm["imps"]) / (2 * np.pi)
        x = {"exp": np.exp(1j * w0 * n), "cos": np.cos(w0 * n) + 0j, "sin": np.sin(w0 * n) + 0j,
             "const": A + 0 * n}[prm["sig"]]
        return bool(np.allclose(claim, x, atol=1e-12))
    if sid == "inverse_finite":                               # sample X_d on 64 frequencies and invert
        c, n1, M = prm["coefs"], prm["n1"], prm["M"]
        wm = 2 * np.pi * np.arange(64) / 64
        Xm = sum(ck * np.exp(-1j * wm * k) for k, ck in zip(range(n1, n1 + len(c)), c))
        n = np.arange(-10, 11)
        x = np.array([np.mean(Xm * np.exp(1j * wm * k)) for k in n])
        support = n[np.abs(x) > 1e-9]
        return bool(support.min() >= 0 and support.max() <= M)
    if sid == "fractional_delay":                             # numerical inverse DTFT of e^{-jwd} on [-pi, pi]
        d = f(prm["d"])
        n = np.arange(-6, 7)
        x = inv_dtft(lambda w: np.exp(-1j * w * d), n)
        if float(d).is_integer():
            return bool(np.allclose(x, (n == d).astype(float), atol=1e-4))
        if np.sum(np.abs(x) > 1e-3) < 5:
            probs.append(f"fractional_delay: inverse DTFT for d = {d} is not spread out")
        return False                                          # delta[n - d] is not a sequence for non-integer d
    if sid == "hermitian_values":                             # test on random real sequences
        a, b = np.pi * f(prm["a"]), np.pi * f(prm["b"])
        ok = True
        for _ in range(6):
            x = RNG.standard_normal(7)
            Xa, Xb = dtft(x, -3, a)[0], dtft(x, -3, b)[0]
            rhs = np.conj(Xb) if prm["rel"] == "conj" else Xb
            ok &= bool(abs(Xa - rhs) < 1e-9)
        return ok
    if sid == "real_from_Xd":                                 # numerical inverse DTFT, look at Im x[n]
        x = inv_dtft(XD[prm["id"]], np.arange(-5, 6))
        m = np.abs(x.imag).max()
        if 1e-6 < m < 1e-3:
            probs.append(f"real_from_Xd {prm['id']}: unclear imaginary part {m}")
        return bool(m < 1e-6)
    if sid == "prop_shift":
        k, form = prm["k"], prm["form"]
        x = RNG.standard_normal(6) + 1j * RNG.standard_normal(6)
        w = RNG.uniform(-np.pi, np.pi, 9)
        Y = dtft(x, -2 + k, w)                                # y[n] = x[n - k]: same samples, start moved by k
        X = lambda v: dtft(x, -2, v)                          # noqa: E731
        claim = {"ok": np.exp(-1j * w * k) * X(w), "sign": np.exp(1j * w * k) * X(w), "freq": X(w - k)}[form]
        return bool(np.allclose(Y, claim))
    if sid == "prop_modulation":
        w0, kind, good = W("w0"), prm["kind"], prm["good"]
        x = RNG.standard_normal(6) + 1j * RNG.standard_normal(6)
        n = -2 + np.arange(6)
        w = RNG.uniform(-np.pi, np.pi, 9)
        X = lambda v: dtft(x, -2, v)                          # noqa: E731
        if kind == "cos":
            Y = dtft(x * np.cos(w0 * n), -2, w)
            claim = 0.5 * X(w - w0) + 0.5 * X(w + w0) if good else X(w - w0) + X(w + w0)
        elif kind == "exp":
            Y = dtft(x * np.exp(1j * w0 * n), -2, w)
            claim = X(w - w0) if good else X(w + w0)
        else:
            Y = dtft(x * (-1.0) ** n, -2, w)
            claim = X(w - np.pi) if good else -X(w)
        return bool(np.allclose(Y, claim))
    if sid == "parseval":
        x, n0 = np.array(prm["x"], float), prm["n0"]
        w = np.linspace(-np.pi, np.pi, 4097)
        X = dtft(x, n0, w)
        val = np.trapezoid(np.abs(X) ** 2, w) if prm["kind"] == "energy" else np.trapezoid(X, w)
        claim = prm["K"] * (np.pi if prm["kpi"] else 1.0)
        return bool(abs(val - claim) < 1e-6)
    if sid == "zp_scaled":
        c, form = prm["c"], prm["form"]
        x = RNG.standard_normal(5) + 1j * RNG.standard_normal(5)
        y = np.concatenate([x, np.zeros(5 * (c - 1))])
        w = RNG.uniform(-np.pi, np.pi, 9)
        X = lambda v: dtft(x, 0, v)                           # noqa: E731
        Y = lambda v: dtft(y, 0, v)                           # noqa: E731
        lhs, rhs = {"same": (Y(w), X(w)), "exam": (X(w), Y(c * w)), "upsample": (Y(w), X(c * w))}[form]
        return bool(np.allclose(lhs, rhs))
    if sid == "cos_multiple":                                 # convolve and fit a constant
        w0 = W("w0")
        if prm["kind"] == "sym":
            h, nh = np.array([1.0, prm["c"], 1.0]), -1
        elif prm["kind"] == "delay":
            h, nh = np.array([1.0]), prm["k"]
        else:
            h, nh = np.array([1.0, 0.0, 1.0]), 0
        n = np.arange(-60, 61)
        x = np.cos(w0 * n)
        y = np.convolve(x, h)
        ny = n[0] + nh + np.arange(len(y))
        keep = np.abs(ny) <= 40
        xs = np.cos(w0 * ny[keep])
        cfit = (y[keep] @ xs) / (xs @ xs)
        if abs(cfit) < 1e-6 and np.abs(y[keep]).max() < 1e-9:
            probs.append("cos_multiple: zero output (H_d(w0) = 0)")
        return bool(np.abs(y[keep] - cfit * xs).max() < 1e-9)
    if sid == "sinusoid_formula":                             # split the cosine into its two exponentials
        H, w0 = SYS_H[prm["id"]], W("w0")
        n = np.arange(0, 25)
        y = 0.5 * H(w0) * np.exp(1j * w0 * n) + 0.5 * H(-w0) * np.exp(-1j * w0 * n)
        claim = np.abs(H(w0)) * np.cos(w0 * n + np.angle(H(w0)))
        return bool(np.allclose(y, claim, atol=1e-9))
    if sid == "group_delay":
        h = np.array(prm["h"], float)
        w = np.linspace(-np.pi, np.pi, 1201)[1:-1] + 1e-3 * np.sqrt(2)
        H = dtft(h, 0, w)
        w = w[np.abs(H) > 1e-3]
        _, gd = sp_group_delay((h, [1.0]), w=w)
        return bool(np.allclose(gd, f(prm["tau"]), atol=1e-6))
    if sid == "mag_nonneg":
        h = np.array(prm["h"], float)
        w = np.linspace(-np.pi, np.pi, 4001)
        H = dtft(h, 0, w)
        A = sum(c * np.cos(f(k) * w) for c, k in prm["terms"])
        if not np.allclose(H * np.exp(1j * f(prm["M"]) * w), A, atol=1e-9):
            probs.append(f"mag_nonneg: H_d(w) is not e^(-jMw) A(w) for h = {prm['h']}")
        return bool(np.allclose(np.abs(H), A, atol=1e-9))
    if sid == "phase_jump":                                   # incremental phase between neighbouring grid points
        w = np.linspace(-np.pi, np.pi, 200001)[1:-1] + 1e-6 * np.sqrt(3)
        H = PHASE_H[prm["id"]](w)
        d = np.angle(H[1:] / H[:-1])
        return bool(np.abs(d).max() > np.pi / 2)
    if sid == "principal_phase":
        k, w0 = prm["k"], W("w0")
        return bool(abs(np.angle(np.exp(-1j * k * w0)) - (-k * w0)) < 1e-9)
    raise KeyError(sid)


# ----------------------------------------------------------------------------- display slips
# Typesetting slips that KaTeX accepts but students would misread; checked in every generated variant.
SLIPS = [
    (re.compile(r"e\^\{\s*[+-]?\s*j\s*[-+\u2212]"), "sign right after j in an exponent (e^{j-...})"),
    (re.compile(r"(?:\+|(?<![<>])-|\u2212)\s*(?:-|\u2212)(?![-\u2212>])"), "double sign (+ -, - -)"),
    (re.compile(r"(?<![\d.\w}\]])1\s*(?:\\[,;:!]\s*)?(?:e\^\{|\\cos|\\sin|\\delta|\\pi|\\omega)"), "coefficient 1 written out"),
    (re.compile(r"\}1\s*[a-z]\b"), "factor 1 glued to a symbol"),
    (re.compile(r"(?<![\d.\w])[01]\s*T(?![\w])"), "0T or 1T"),
    (re.compile(r"(?<![\d.\w])-\s*0(?![\d.])"), "minus zero"),
    (re.compile(r"\\left\(\s*-?\d+\s*\\right\)\s*\\(?:cos|sin)"), "parenthesised integer before cos/sin"),
    (re.compile(r"e\^\{-?j0\}"), "e^{j0}"),
]
TEXT_SLIPS = [(re.compile(r"(?<![\d.])1\$?\s+samples\b"), "1 samples")]
MATH_RX = re.compile(r"\$\$(.+?)\$\$|\$(.+?)\$", re.S)


def display_slips(sid, html):
    out = []
    for a, b in MATH_RX.findall(html):
        e = a or b
        for rx, what in SLIPS:
            m = rx.search(e)
            if m:
                out.append(f"{sid}: display slip ({what}): ...{e[max(0, m.start() - 20):m.end() + 20]}...")
    for rx, what in TEXT_SLIPS:
        if rx.search(html):
            out.append(f"{sid}: display slip ({what})")
    return out


def check(params, correct):
    probs = []
    items = params["statements"]
    if len(items) != 6:
        probs.append(f"{len(items)} statements instead of 6")
    if len({it["sid"] for it in items}) != len(items):
        probs.append("duplicate statement")
    for it in items:
        probs += display_slips(it["sid"], it["stmt_html"]) + display_slips(it["sid"], it["why"])
        want = recompute(it["sid"], it["param"], probs)
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
    cases = []
    for nm in names[:3]:                    # one statement flipped: that one scores 0, the other five 1 (total 5/6)
        cases.append(({nm: _key(params, nm, flip[right[nm]])}, {m: (0 if m == nm else 1) for m in names}))
    cases.append(({nm: _key(params, nm, flip[right[nm]]) for nm in names}, {nm: 0 for nm in names}))    # all flipped
    cases.append(({nm: _key(params, nm, right[nm]) for nm in names}, {nm: 1 for nm in names}))         # all correct
    cases.append(({names[0]: correct[names[0]]["key"]}, {names[0]: 1}))
    cases.append(({nm: _key(params, nm, "True") for nm in names},                                     # "all True" guess
                  {nm: (1 if right[nm] == "True" else 0) for nm in names}))
    return cases
