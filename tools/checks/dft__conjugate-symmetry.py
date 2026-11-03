"""Independent check of questions/dft/conjugate-symmetry (numpy only; the displayed LaTeX is parsed).

fa23 (FA2023 MT2 #7): Z is parsed from the displayed list, z = ifft(Z) is split into x = Re z and
     y = Im z, and X = fft(x), Y = fft(y) give the answers at k1 and k2 (k2 values must be real; for
     N = 4, x and y are integers). (a): for a random complex b, fft(conj(b)) is compared with every
     displayed option, evaluated from its LaTeX; exactly the marked one must match. Probes: X = Re Z,
     Y = Im Z, no conjugate, no 1/2, plain flip, 1/(2j) -> j/2 or 1/2, Z^*[N] read as 0 or Z^*[N-1].
sp21 (SP2021 MT2 #5): every displayed list is inverted with ifft; the marked options must be exactly
     the real signals; the shortcuts "X[0], X[N/2] real" and "plus |X[k]| = |X[N-k]|" must fail and a
     plain-flip list must be among the options. (b) x[0] from ifft; (c) the marked justification is
     the true one (real / x[0] complex / x[0] real but not real).
fa25 (FA2025 MT2 #6(a)): each displayed relation X[k] = X_d(w) or X_d^*(w) is evaluated for three
     random real signals (fft against a direct DTFT sum); true iff it holds for all of them. (b): brute
     force over k for X[k] = conj(X_d(w)). Shortcuts and giveaways between (a) and (b) are rejected."""

import re
from fractions import Fraction

import numpy as np

RNG = np.random.default_rng(3101)

LINT = [(re.compile(r"e\^\{\s*[+-]?\s*j\s*[-+]"), "sign after j in an exponent"),
        (re.compile(r"(?:\+|(?<![<>])-)\s*-(?![->])"), "double sign"),
        (re.compile(r"\}1\s*[a-z]\b"), "factor 1 glued to a symbol")]


def _strings(obj):
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _strings(v)


def _lint(params):
    out = []
    for s in _strings(params):
        for rx, what in LINT:
            if rx.search(s):
                out.append(f"display lint ({what}): {s[:120]}")
    return out


def _cplx(item):
    s = item.replace(" ", "").replace(r"\,", "")
    m = re.fullmatch(r"(-?)\\tfrac\{(\d+)\}\{(\d+)\}", s)
    if m:
        return (-1 if m.group(1) else 1) * int(m.group(2)) / int(m.group(3))
    s = re.sub(r"(?<![\d.])j", "1j", s)
    return complex(s)


def _parse_list(tex):
    body = tex.strip().strip("$")
    body = body.split(r"\{", 1)[1].rsplit(r"\}", 1)[0]
    return np.array([_cplx(t) for t in body.split(r",\ ")])


def _pi_mult(tex):
    """'-\\tfrac{3\\pi}{8}' -> -0.375, '2\\pi' -> 2, '\\pi' -> 1, '0' -> 0 (multiples of pi)."""
    s = tex.strip()
    if s == "0":
        return 0.0
    sgn = -1.0 if s.startswith("-") else 1.0
    s = s.lstrip("-")
    m = re.fullmatch(r"\\tfrac\{(\d*)\\pi\}\{(\d+)\}", s)
    if m:
        return sgn * int(m.group(1) or 1) / int(m.group(2))
    m = re.fullmatch(r"(\d*)\\pi", s)
    return sgn * int(m.group(1) or 1)


def _dtft(x, w):
    return complex(np.sum(x * np.exp(-1j * w * np.arange(len(x)))))


def _num(v):
    v = float(v)
    return str(int(round(v))) if abs(v - round(v)) < 1e-12 else repr(v)


def _differs(wrong, right):
    return abs(float(wrong) - float(right)) > 2 * (1e-6 + 1e-3 * abs(float(right))) and abs(float(wrong) - float(right)) > 0.02 * abs(float(right))


# ============================================================================= fa23
CONJ_EVAL = {
    r"B^*[\langle -k\rangle_N]": lambda B, k, N: np.conj(B[(-k) % N]),
    r"B^*[\langle N-k\rangle_N]": lambda B, k, N: np.conj(B[(N - k) % N]),
    r"B^*[k]": lambda B, k, N: np.conj(B),
    r"B[\langle -k\rangle_N]": lambda B, k, N: B[(-k) % N],
    r"B^*[N+k]": lambda B, k, N: np.conj(B[(N + k) % N]),
    r"-B[k]": lambda B, k, N: -B,
    r"B^*[N-1-k]": lambda B, k, N: np.conj(B[N - 1 - k]),
}


def _conj_option(html):
    m = re.fullmatch(r"\$\\\{(.+)\\\}_\{k=0\}\^\{N-1\}\$", html.strip())
    return m.group(1)


def _fa23_truth(params):
    Z = _parse_list(params["Z_tex"])
    N = len(Z)
    z = np.fft.ifft(Z)
    x, y = z.real, z.imag
    return Z, N, x, y, np.fft.fft(x), np.fft.fft(y)


def _fa23_probes(params):
    Z, N, x, y, X, Y = _fa23_truth(params)
    k, k2 = params["k1"], params["k2"]
    km = (-k) % N
    Zk, Zm, Zf = Z[k], Z[km], Z[N - 1 - k]
    P = {
        "X": {"reZ": complex(Zk.real, 0), "noconj": (Zk + Zm) / 2, "nohalf": Zk + np.conj(Zm), "flip": (Zk + np.conj(Zf)) / 2},
        "Y": {"imZ": complex(Zk.imag, 0), "noj": (Zk - np.conj(Zm)) / 2, "sign": -Y[k], "noconj": (Zk - Zm) / 2j,
              "nohalf": (Zk - np.conj(Zm)) / 1j},
    }
    Z2 = Z[k2]
    P["Xs"] = {"zero": Z2.real / 2, "swap": Z2.imag, "prev": ((Z2 + np.conj(Z[(k2 - 1) % N])) / 2).real}
    P["Ys"] = {"sign": -Z2.imag, "swap": Z2.real, "half": Z2.imag / 2}
    return P


def _check_fa23(params, correct):
    probs = []
    Z, N, x, y, X, Y = _fa23_truth(params)
    k1, k2 = params["k1"], params["k2"]
    if N not in (4, 6) or k1 in (0, N // 2) or k2 not in (0, N // 2):
        probs.append(f"N = {N}, k1 = {k1}, k2 = {k2}")
    if not (np.allclose(Z.real, np.round(Z.real)) and np.allclose(Z.imag, np.round(Z.imag))) or np.max(np.abs(Z)) > 17:
        probs.append(f"Z not small Gaussian integers: {Z}")
    if N == 4 and not (np.allclose(x, np.round(x)) and np.allclose(y, np.round(y))):
        probs.append(f"N = 4 but x = {x}, y = {y} are not integers")
    want = {"reX": X[k1].real, "imX": X[k1].imag, "reY": Y[k1].real, "imY": Y[k1].imag, "Xs": X[k2].real, "Ys": Y[k2].real}
    for name, v in want.items():
        if abs(correct[name] - v) > 1e-9:
            probs.append(f"{name}: answer {correct[name]}, numpy {v}")
    if abs(X[k2].imag) > 1e-9 or abs(Y[k2].imag) > 1e-9:
        probs.append("X[k2], Y[k2] not real")
    P = _fa23_probes(params)
    for name, part, src in (("reX", "real", "X"), ("imX", "imag", "X"), ("reY", "real", "Y"), ("imY", "imag", "Y")):
        n = sum(_differs(getattr(v, part), correct[name]) for v in P[src].values())
        if n < 2:
            probs.append(f"{name}: only {n} classic mistakes differ from the answer")
    for name in ("Xs", "Ys"):
        if sum(_differs(v, correct[name]) for v in P[name].values()) < 2:
            probs.append(f"{name}: fewer than two classic mistakes differ")
    # (a) conjugation property, numerically
    b = RNG.normal(size=7) + 1j * RNG.normal(size=7)
    B, target, kk = np.fft.fft(b), np.fft.fft(np.conj(b)), np.arange(7)
    hits = []
    for o in params["conj"]:
        f = CONJ_EVAL[_conj_option(o["html"])]
        if np.allclose(f(B, kk, 7), target):
            hits.append(o["key"])
    if hits != [correct["conj"]["key"]]:
        probs.append(f"(a) options matching fft(conj b): {hits}, marked {correct['conj']['key']}")
    if len(params["conj"]) != 5:
        probs.append("(a) not 5 options")
    return probs


def _subs_fa23(params, correct):
    cases = []
    for name in ("reX", "imX", "reY", "imY", "Xs", "Ys"):
        v = correct[name]
        cases += [({name: f" {_num(v)} "}, {name: 1}), ({name: f"{2 * v:g}/2"}, {name: 1}),
                  ({name: f"{v:.4f}"}, {name: 1})]
    cases.append(({"reX": "1+2j"}, {"reX": "invalid"}))
    P = _fa23_probes(params)
    for name, part, src in (("reX", "real", "X"), ("imX", "imag", "X"), ("reY", "real", "Y"), ("imY", "imag", "Y")):
        for v in P[src].values():
            w = getattr(v, part)
            if _differs(w, correct[name]):
                cases.append(({name: _num(round(w, 6))}, {name: 0}))
    for name in ("Xs", "Ys"):
        for w in P[name].values():
            if _differs(w, correct[name]):
                cases.append(({name: _num(round(w, 6))}, {name: 0}))
    for o in params["conj"]:
        cases.append(({"conj": o["key"]}, {"conj": 1 if o["key"] == correct["conj"]["key"] else 0}))
    return cases


# ============================================================================= sp21
def _is_real(L):
    return bool(np.max(np.abs(np.fft.ifft(L).imag)) < 1e-9)


def _check_sp21(params, correct):
    probs = []
    N = params["N"]
    opts = params["real"]
    lists = [_parse_list(o["html"]) for o in opts]
    real = [_is_real(L) for L in lists]
    marked = {o["key"] for o in correct["real"]}
    for o, r in zip(opts, real):
        if r != (o["key"] in marked):
            probs.append(f"(a) {o['html']}: real = {r}, marked {o['key'] in marked}")
    if len(opts) != 5 or not 1 <= sum(real) <= 3 or any(len(L) != N for L in lists):
        probs.append(f"(a) {len(opts)} options, {sum(real)} real")
    if len({tuple(L) for L in lists}) != 5:
        probs.append("(a) repeated list")
    h = N // 2
    selfp = [abs(L[0].imag) < 1e-9 and (N % 2 == 1 or abs(L[h].imag) < 1e-9) for L in lists]
    mags = [s and np.allclose(np.abs(L), np.abs(L[(-np.arange(N)) % N])) for s, L in zip(selfp, lists)]
    if sum(selfp) <= sum(real) or sum(mags) <= sum(real):
        probs.append("(a) checking only X[0], X[N/2] (or also the magnitudes) finds the answer")
    flip = [not r and np.allclose(L[1:N - 1], np.conj(L[1:N - 1][::-1])) for L, r in zip(lists, real)]   # inner pairs k <-> N-1-k
    if not any(flip):
        probs.append("(a) no plain-flip near miss (X[k] = X^*[N-1-k]) among the options")
    for L in lists:
        if np.any(np.abs(L) < 1e-9) or np.max(np.abs(L.real)) > 7 or np.max(np.abs(L.imag)) > 4:
            probs.append(f"(a) entries not exam-nice: {L}")
    # (b), (c)
    L = _parse_list(params["XL_tex"])
    x = np.fft.ifft(L)
    if len(L) != N:
        probs.append("(b) list length")
    if abs(correct["rex0"] - x[0].real) > 1e-9 or abs(correct["imx0"] - x[0].imag) > 1e-9:
        probs.append(f"(b) x[0] = {x[0]}, answer {correct['rex0']} + j{correct['imx0']}")
    if any(np.allclose(L, M) for M in lists):
        probs.append("(b) list repeats an option of (a)")
    if _is_real(L):
        want = "Yes:"
    elif abs(x[0].imag) > 1e-9:
        want = r"No: $x[0]$ is not real"
    else:
        want = r"No: $x[0]$ is real"
    hit = [o["key"] for o in params["isreal"] if o["html"].strip().startswith(want)]
    if hit != [correct["isreal"]["key"]]:
        probs.append(f"(c) expected '{want}', marked {correct['isreal']['html'][:40]}")
    return probs


def _subs_sp21(params, correct):
    lists = [_parse_list(o["html"]) for o in params["real"]]
    keys = [o["key"] for o in params["real"]]
    T = [k for k, L in zip(keys, lists) if _is_real(L)]
    F = [k for k, L in zip(keys, lists) if not _is_real(L)]
    cases = [({"real": T}, {"real": 1}), ({"real": T + F[:1]}, {"real": (len(T) - 1) / len(T)}), ({"real": F}, {"real": 0})]
    cases += [({"real": [k]}, {"real": 0}) for k in F]
    if len(T) > 1:
        cases.append(({"real": T[:1]}, {"real": 1 / len(T)}))
    L = _parse_list(params["XL_tex"])
    N = len(L)
    re0, im0 = correct["rex0"], correct["imx0"]
    for name, v in (("rex0", re0), ("imx0", im0)):
        cases += [({name: f"{N * v:g}/{N}"}, {name: 1}), ({name: f" {v:.6f} "}, {name: 1})]
    for w in (L.real.sum(), L[0].real):                              # no 1/N; x[0] confused with X[0]
        if _differs(w, re0):
            cases.append(({"rex0": _num(round(w, 6))}, {"rex0": 0}))
    for w in (L.imag.sum(), -im0, L[1:(N + 1) // 2].imag.sum() / N):  # no 1/N; sign; one half of the pairs
        if _differs(w, im0):
            cases.append(({"imx0": _num(round(w, 6))}, {"imx0": 0}))
    for o in params["isreal"]:
        cases.append(({"isreal": o["key"]}, {"isreal": 1 if o["key"] == correct["isreal"]["key"] else 0}))
    return cases


# ============================================================================= fa25
REL = re.compile(r"\$X\[(\d+)\] = X_d(\^\*)?\\!\\left\((.+)\\right\)\$")


def _relations(params):
    out = []
    for o in params["rel"]:
        m = REL.fullmatch(o["html"].strip())
        out.append((o["key"], int(m.group(1)), bool(m.group(2)), _pi_mult(m.group(3))))
    return out


def _holds(N, k, conj, t):
    gaps = []
    for _ in range(3):
        x = RNG.normal(size=N)
        X = np.fft.fft(x)
        v = _dtft(x, np.pi * t)
        gaps.append(abs(X[k] - (np.conj(v) if conj else v)))
    if max(gaps) < 1e-8:
        return True
    if min(gaps) > 1e-3:
        return False
    return None


def _bin_conj(N, t):
    x = RNG.normal(size=N)
    X, v = np.fft.fft(x), np.conj(_dtft(x, np.pi * t))
    return [k for k in range(N) if abs(X[k] - v) < 1e-8]


def _check_fa25(params, correct):
    probs = []
    N = params["N"]
    rels = _relations(params)
    marked = {o["key"] for o in correct["rel"]}
    truth = []
    for key, k, conj, t in rels:
        h = _holds(N, k, conj, t)
        if h is None:
            probs.append(f"(a) X[{k}] vs X_d{'*' if conj else ''}({t} pi): undecided")
        elif h != (key in marked):
            probs.append(f"(a) X[{k}] = X_d{'*' if conj else ''}({t} pi): numpy {h}, marked {key in marked}")
        truth.append(bool(h))
    if len(rels) != 5 or not 2 <= sum(truth) <= 4:
        probs.append(f"(a) {len(rels)} relations, {sum(truth)} true")
    ks = [k for _, k, _, _ in rels]
    if len(set(ks)) != 5 or ks != sorted(ks):
        probs.append(f"(a) bins {ks} not distinct and sorted")
    shortcuts = {"conjugate <=> true": [c for _, _, c, _ in rels],
                 "no conjugate <=> true": [not c for _, _, c, _ in rels],
                 "ignore the conjugate": [abs(((t - 2 * k / N) / 2) - round((t - 2 * k / N) / 2)) < 1e-9 for _, k, _, t in rels],
                 "only the literal w_k": [(not c) and abs(t - 2 * k / N) < 1e-9 for _, k, c, t in rels]}
    for name, guess in shortcuts.items():
        if guess == truth:
            probs.append(f"(a) the shortcut '{name}' answers every relation")
    for nm, tex in (("ka", params["wa_tex"]), ("kb", params["wb_tex"])):
        hits = _bin_conj(N, _pi_mult(tex))
        if hits != [correct[nm]]:
            probs.append(f"(b) {nm}: bins with X[k] = X_d^*({tex}) are {hits}, answer {correct[nm]}")
    touched = set(ks) | {(N - k) % N for k in ks}
    for nm in ("ka", "kb"):
        if correct[nm] in touched or (N - correct[nm]) % N in touched:
            probs.append(f"(b) {nm} = {correct[nm]} is a bin used in (a)")
    return probs


def _subs_fa25(params, correct):
    N = params["N"]
    rels = _relations(params)
    marked = [o["key"] for o in correct["rel"]]
    T = [key for key, *_ in rels if key in marked]
    F = [key for key, *_ in rels if key not in marked]
    cases = [({"rel": T}, {"rel": 1}), ({"rel": T + F[:1]}, {"rel": (len(T) - 1) / len(T)}), ({"rel": T[:1]}, {"rel": 1 / len(T)}),
             ({"rel": F}, {"rel": 0})]
    cases += [({"rel": [k]}, {"rel": 0}) for k in F]
    ka, kb = correct["ka"], correct["kb"]
    cases += [({"ka": f"+{ka}"}, {"ka": 1}), ({"ka": f" {ka} "}, {"ka": 1}), ({"kb": f"0{kb}"}, {"kb": 1}),
              ({"kb": f"{kb}.0"}, {"kb": "invalid"})]
    ta, tb = _pi_mult(params["wa_tex"]), _pi_mult(params["wb_tex"])
    for wrong in {round(N * ta / 2) % N, round(N * ta) % N, ka - N}:          # no conjugate; pi k/N grid; negative bin
        if wrong != ka:
            cases.append(({"ka": str(wrong)}, {"ka": 0}))
    for wrong in {round(N * tb / 2) % N, round(-N * tb) % N, N - kb}:         # no conjugate; pi grid; mirror
        if wrong != kb:
            cases.append(({"kb": str(wrong)}, {"kb": 0}))
    return cases


# ============================================================================= entry points
def check(params, correct):
    kind = params["kind"]
    probs = {"fa23": _check_fa23, "sp21": _check_sp21, "fa25": _check_fa25}[kind](params, correct)
    return probs + _lint(params)


def submissions(params, correct):
    return {"fa23": _subs_fa23, "sp21": _subs_sp21, "fa25": _subs_fa25}[params["kind"]](params, correct)
