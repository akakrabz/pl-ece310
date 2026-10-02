"""Independent check of questions/systems/system-properties.

Every system is re-implemented here as a numpy callable T(x, n) -> y[n] (bank systems by id, generated
systems from their structure with this file's own atom implementations) and its four properties are
measured numerically, without looking at the stored flags:
  linear          superposition with random real, then complex, inputs and scalars;
  time-invariant  T{x[n - n0]}[n] == y[n - n0] for random inputs and several n0;
  causal          perturbing x[k] for k > n0 never changes y[n0] (n0 = -8..8);
  BIBO stable     max |y| over bounded inputs (|x| <= 1, incl. 0, delta, 1, (-1)^n, j^n, u[n], u[-n])
                  stays finite and does not grow when the input window grows from ~25 to 400 samples.
"""

import json

import numpy as np

# ============================================================================ signals


class Sig:
    """Finite signal: vals[i] = x[start + i], zero elsewhere."""

    def __init__(self, vals, start):
        self.v = np.asarray(vals, dtype=complex)
        self.s = int(start)

    def __call__(self, idx):
        idx = np.asarray(idx, dtype=np.int64)
        out = np.zeros(idx.shape, dtype=complex)
        i = idx - self.s
        m = (i >= 0) & (i < len(self.v))
        out[m] = self.v[i[m]]
        return out

    def cum(self, idx):
        """sum_{k <= idx} x[k]."""
        c = np.concatenate([[0], np.cumsum(self.v)])
        i = np.clip(np.asarray(idx, dtype=np.int64) - self.s + 1, 0, len(self.v))
        return c[i]

    def total(self):
        return self.v.sum()


def conv(x, n, h, klo, khi):
    """sum_{k=klo}^{khi} h[k] x[n-k], with the infinite side truncated beyond the input support."""
    n = np.asarray(n, dtype=np.int64)
    span = len(x.v) + int(np.abs(n).max()) + abs(x.s) + 5
    lo = klo if klo is not None else -span
    hi = khi if khi is not None else span
    ks = np.arange(lo, hi + 1)
    with np.errstate(all="ignore"):
        hk = np.nan_to_num(h(ks).astype(complex))
    full = np.convolve(x.v, hk)
    return Sig(full, x.s + lo)(n)


def r(a):
    return np.real(a)


# ============================================================================ bank systems (by id)
def _u(k):
    return (np.asarray(k) >= 0).astype(float)


BANK = {
    "abs_n_gain": lambda x, n: np.abs(n) * x(n),
    "abs_diff": lambda x, n: np.abs(x(n) - x(n - 1)),
    "conv_2n_uneg": lambda x, n: conv(x, n, lambda k: 2.0 ** np.minimum(k, 0) * _u(-k), None, 0),
    "conv_jn": lambda x, n: conv(x, n, lambda k: (1j) ** (k % 4) * _u(k), 0, None),
    "abs_idx_aff": lambda x, n: 2 * x(np.abs(n)) + 10,
    "exp_plus1": lambda x, n: np.exp(x(n) + 1),
    "prod_next": lambda x, n: x(n) * x(n + 1),
    "inv_abs_gain": lambda x, n: x(n) / (np.abs(n) + 1),
    "sin_plus_x0": lambda x, n: np.sin(x(n)) + x(np.zeros_like(n)),
    "log_gain": lambda x, n: np.log(np.abs(n) + 1) * x(n),
    "conv_u_np1": lambda x, n: conv(x, n, lambda k: _u(k + 1), -1, None),
    "plus3": lambda x, n: x(n) + 3,
    "conv_m1n": lambda x, n: conv(x, n, lambda k: np.where(k % 2 == 0, 1.0, -1.0) * _u(k), 0, None),
    "div_x2": lambda x, n: x(n) / x(np.full_like(n, 2)),
    "cos2_gain": lambda x, n: np.cos(np.pi * n / 2) ** 2 * x(n),
    "cos_shift_gain": lambda x, n: x(n) * np.cos(np.pi * (n - 2) / 3),
    "x3_prod": lambda x, n: x(np.full_like(n, 3)) * x(n),
    "cgain": lambda x, n: (0.8 + 0.8j) ** n * x(n),
    "idx_absn_plus_n": lambda x, n: x(np.abs(n) + n),
    "clip": lambda x, n: np.clip(r(x(n)), -2, 2) + 0j if np.allclose(np.imag(x(n)), 0) else np.full(n.shape, np.nan),
    "window": lambda x, n: x(n) * ((n >= 0) & (n <= 3)),
    "half_abs_gain": lambda x, n: 0.5 ** np.abs(n) * x(n),
    "avg2": lambda x, n: 0.5 * x(n) + 0.5 * x(n - 1),
    "absx": lambda x, n: np.abs(x(n)) + 0j,
    "expx": lambda x, n: np.exp(x(n)),
    "n_diff": lambda x, n: n * (x(n) - x(n - 1)),
    "max0": lambda x, n: np.maximum(0, r(x(n))) + 0j if np.allclose(np.imag(x(n)), 0) else np.full(n.shape, np.nan),
    "x0_prod": lambda x, n: x(n) * x(np.zeros_like(n)),
    "diff": lambda x, n: x(n) - x(n - 1),
    "down2": lambda x, n: x(2 * n),
    "fold": lambda x, n: x(np.abs(n)),
    "avg3": lambda x, n: (x(n) + x(n - 1) + x(n - 2)) / 3,
    "idx_absn_minus_n": lambda x, n: x(np.abs(n) - n),
    "np1_diff": lambda x, n: (n + 1) * (x(n) - x(n - 1)),
    "sin2": lambda x, n: np.sin(x(n)) ** 2,
    "ln_abs": lambda x, n: np.log(np.abs(x(n))) + 0j,
    "median": lambda x, n: (np.median(np.stack([r(x(n)), r(x(n - 1)), r(x(n - 2))]), axis=0) + 0j
                            if np.allclose(np.imag(x(n)), 0) else np.full(n.shape, np.nan)),
    "square": lambda x, n: x(n) ** 2,
    "n_x3n": lambda x, n: n * x(3 * n),
    "fir3": lambda x, n: x(n) + x(n - 2) - x(n - 4),
    "x10_exp": lambda x, n: x(n) ** 10 + np.exp(x(n)),
    "conv_half_np2": lambda x, n: conv(x, n, lambda k: 0.5 ** k * _u(k + 2), -2, None),
    "x_plus_xrev": lambda x, n: x(n) + x(-n),
    "accum": lambda x, n: x.cum(n),
    "cos_times_prev": lambda x, n: np.cos(x(n)) * x(n - 1),
    "half_n_prev": lambda x, n: x(n) + 0.5 ** n * x(n - 1),
    "future_sum": lambda x, n: conv(x, n, lambda k: 3.0 ** np.minimum(k, 2) * _u(2 - k), None, 2),
}
REAL_ONLY = {"clip", "max0", "median"}

# ============================================================================ generated systems
GAIN = {
    "one": lambda n: np.ones(n.shape),
    "n": lambda n: n.astype(float),
    "half_n": lambda n: 0.5 ** n.astype(float),
    "m1n": lambda n: np.where(n % 2 == 0, 1.0, -1.0),
    "cos3": lambda n: np.cos(np.pi * n / 3),
    "un": lambda n: (n >= 0).astype(float),
    "half_abs": lambda n: 0.5 ** np.abs(n).astype(float),
}
FUNC = {"id": lambda v: v, "sq": lambda v: v ** 2, "cube": lambda v: v ** 3, "exp": np.exp, "cos": np.cos}
MAP = {"rev": lambda n: -n, "down2": lambda n: 2 * n, "fold": lambda n: np.abs(n)}


def _term(t):
    if t["type"] == "basic":
        g, f = GAIN[t["g"]], FUNC[t["f"]]
        kind, k = t["phi"]
        phi = (lambda n, k=k: n - k) if kind == "shift" else MAP[kind]
        return lambda x, n: g(n) * f(x(phi(n)))
    kind, a, b = t["kind"], t["a"], t["b"]
    if kind == "mavg":
        return lambda x, n: sum(x(n + j) for j in range(-a, b + 1))
    if kind == "accum":
        return lambda x, n: x.cum(n - a)
    if kind == "geo":
        return lambda x, n: conv(x, n, lambda k: 0.5 ** np.maximum(k, 0) * _u(k), 0, None)
    if kind == "geo_future":
        return lambda x, n: conv(x, n, lambda k: 0.5 ** np.maximum(-k, 0) * _u(-k), None, 0)
    if kind == "accum_future":
        return lambda x, n: x.total() - x.cum(n - 1)
    raise ValueError(kind)


def build(origin):
    if origin["kind"] == "bank":
        return BANK[origin["id"]], origin["id"] in REAL_ONLY
    st = origin["struct"]
    terms = [(t["c"][0] / t["c"][1], _term(t)) for t in st["terms"]]
    d = st["d"][0] / st["d"][1]
    return (lambda x, n: sum(c * T(x, n) for c, T in terms) + d), False


# ============================================================================ property probes
NS = np.arange(-20, 21)


def _close(a, b):
    a, b = np.asarray(a, complex), np.asarray(b, complex)
    fa, fb = np.isfinite(a), np.isfinite(b)
    if not np.array_equal(fa, fb):
        return False
    if not np.all(fa):                    # identical non-finite samples (e.g. ln 0 = -inf) count as equal
        if not np.array_equal(a[~fa], b[~fb], equal_nan=True):
            return False
        a, b = a[fa], b[fb]
        if a.size == 0:
            return True
    scale = max(1.0, float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    return bool(np.allclose(a, b, rtol=1e-7, atol=1e-9 * scale))


def is_linear(T, real_only, rng):
    for cplx in ([False] if real_only else [False, True]):
        for _ in range(3):
            def rnd(m):
                v = rng.uniform(-2, 2, m)
                return v + 1j * rng.uniform(-2, 2, m) if cplx else v
            v1, v2 = rnd(17), rnd(17)
            a, b = (rnd(1)[0], rnd(1)[0])
            with np.errstate(all="ignore"):
                lhs = T(Sig(a * v1 + b * v2, -8), NS)
                rhs = a * T(Sig(v1, -8), NS) + b * T(Sig(v2, -8), NS)
            if not _close(lhs, rhs):
                return False
    return True


def is_ti(T, rng):
    for _ in range(3):
        v = rng.uniform(-2, 2, 17)
        for n0 in (1, 2, 3, -2, 5):
            with np.errstate(all="ignore"):
                if not _close(T(Sig(v, -8 + n0), NS), T(Sig(v, -8), NS - n0)):
                    return False
    return True


def is_causal(T, rng):
    ks = np.arange(-20, 21)
    for n0 in range(-8, 9):
        for _ in range(2):
            v = rng.uniform(0.5, 2, 41) * rng.choice([-1, 1], 41)
            w = v.copy()
            fut = ks > n0
            w[fut] = rng.uniform(0.5, 2, fut.sum()) * rng.choice([-1, 1], fut.sum())
            with np.errstate(all="ignore"):
                if not _close(T(Sig(v, -20), np.array([n0])), T(Sig(w, -20), np.array([n0]))):
                    return False
    return True


def _bounded_inputs(W, real_only):
    k = np.arange(-W, W + 1)
    ins = [Sig([0.0], 0), Sig([1.0], 0), Sig(np.ones(2 * W + 1), -W), Sig(np.where(k % 2 == 0, 1.0, -1.0), -W),
           Sig((k >= 0).astype(float), -W), Sig((k <= 0).astype(float), -W),
           Sig(np.cos(np.pi * k / 3), -W)]
    if not real_only:
        ins += [Sig((1j) ** (k % 4), -W), Sig((-1j) ** (k % 4), -W)]
    return ins


def _peak(T, W, real_only):
    n = np.arange(-2 * W, 2 * W + 1)
    peak = 0.0
    for x in _bounded_inputs(W, real_only):
        with np.errstate(all="ignore"):
            y = T(x, n)
        if not np.all(np.isfinite(y)):
            return np.inf
        peak = max(peak, float(np.max(np.abs(y))))
    return peak


def is_stable(T, real_only):
    """Bounded outputs for every bounded test input, and no growth with the window: the small-window peak
    is taken over W = 25..36 (every alignment of the window edges with gains of period 1, 2, 3, 4, 6, 12)."""
    p_small = max(_peak(T, W, real_only) for W in range(25, 37))
    p_large = max(_peak(T, W, real_only) for W in (400, 401))
    if not (np.isfinite(p_small) and np.isfinite(p_large)):
        return False
    return p_large <= 1.05 * p_small + 1e-9


_CACHE = {}


def measure(origin):
    key = json.dumps(origin, sort_keys=True)
    if key not in _CACHE:
        T, real_only = build(origin)
        rng = np.random.default_rng(12345)
        _CACHE[key] = [is_linear(T, real_only, rng), is_ti(T, rng), is_causal(T, rng), is_stable(T, real_only)]
    return _CACHE[key]


NAMES = ["lin", "ti", "caus", "stab"]


def check(params, correct):
    probs = []
    got = measure(params["origin"])
    for name, flag, row in zip(NAMES, got, params["rows"]):
        want = correct[name]["html"].strip()
        if want != ("Yes" if flag else "No"):
            probs.append(f"{name}: graded answer {want} but numerical probe says {flag} for y[n] = {params['sys_tex']}")
        if row["ans"] != want:
            probs.append(f"{name}: answer panel says {row['ans']} but grading uses {want}")
    return probs


def _key(params, name, text):
    for opt in params[name]:
        if opt["html"].strip() == text:
            return opt["key"]
    raise KeyError(f"{name}: no option {text}")


def submissions(params, correct):
    cases = []
    right = {nm: correct[nm]["html"].strip() for nm in NAMES}
    for nm in NAMES:                                       # each property flipped -> 0
        wrong = "No" if right[nm] == "Yes" else "Yes"
        cases.append(({nm: _key(params, nm, wrong)}, {nm: 0}))
    cases.append(({nm: _key(params, nm, right[nm]) for nm in NAMES}, {nm: 1 for nm in NAMES}))   # explicit keys
    cases.append(({"lin": correct["lin"]["key"], "stab": correct["stab"]["key"]}, {"lin": 1, "stab": 1}))
    cases.append(({nm: _key(params, nm, "Yes") for nm in NAMES},                # "all Yes" guess
                  {nm: (1 if right[nm] == "Yes" else 0) for nm in NAMES}))
    return cases
