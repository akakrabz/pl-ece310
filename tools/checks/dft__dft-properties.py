"""Independent check of questions/dft/dft-properties.

(a) Y[k] is built with numpy from X = fft(x) exactly as the question defines it (the explicit list of
+-X_i is parsed from the displayed LaTeX), and y = ifft(Y) is compared with the correct answer; the
probes (wrong shift direction, linear instead of circular shift, missing reversal, truncated linear
convolution, ...) are computed here with numpy from the family parameters.
(b) every displayed candidate list is parsed and inverted with ifft: exactly one gives a real signal,
and that one is marked correct."""

import re

import checklib  # noqa: F401  (puts PrairieLearn's python library on sys.path)
import numpy as np
import prairielearn as pl


def _y_correct(correct):
    return np.array(pl.from_json(correct["y"]), dtype=complex).reshape(-1)


def _parse_list(tex):
    """'\\{X_{6},\\ -X_{7}\\}' -> [(1, 6), (-1, 7)]"""
    body = tex.split(r"\{", 2)[-1].rsplit(r"\}", 1)[0]
    out = []
    for item in body.split(r",\ "):
        m = re.fullmatch(r"\s*(-?)X_\{(\d+)\}\s*", item)
        out.append((-1 if m.group(1) else 1, int(m.group(2))))
    return out


def _parse_cplx_list(html):
    body = html.strip().strip("$")
    body = body.split(r"\{", 1)[1].rsplit(r"\}", 1)[0]
    vals = []
    for item in body.split(r",\ "):
        s = item.replace(" ", "")
        s = re.sub(r"(?<![\d.])j", "1j", s)
        vals.append(complex(s))
    return np.array(vals)


def _Y(params):
    x = np.array(params["x"], float)
    N = len(x)
    X = np.fft.fft(x)
    k = np.arange(N)
    fam, d = params["family"], params.get("disp", {})
    if fam == "shift":
        t = d["t"][0] / d["t"][1]
        sgn = 1 if d["advance"] else -1
        return np.exp(sgn * 1j * (np.pi * t * k + (np.pi if d["neg"] else 0))) * X
    if fam == "reverse":
        return np.conj(X) if d["conj"] else X[(-k) % N]
    if fam == "revshift":
        t = d["t"][0] / d["t"][1]
        return np.exp(-1j * np.pi * t * k) * X[(-k) % N]
    if fam == "half":
        return (-1.0) ** k * X if d["alt"] else X[(k - N // 2) % N]
    if fam in ("list", "listf"):
        pairs = _parse_list(params["Y_def_tex"].split("=", 1)[1])
        return np.array([s * X[i] for s, i in pairs])
    if fam == "square":
        return X ** 2
    if fam == "evenodd":
        return X.real + 0j if d["even"] else 1j * X.imag
    raise KeyError(fam)


def _probes(params):
    """Classic wrong answers (complex arrays), computed from the family parameters with numpy."""
    x = np.array(params["x"], float)
    N = len(x)
    n = np.arange(N)
    fam, d = params["family"], params.get("disp", {})
    out = {}

    def lin_delay(v, m):
        return np.concatenate([np.zeros(m), v[:N - m]])

    if fam == "shift":
        t = d["t"][0] / d["t"][1]
        q = round(t * N / 2)
        m = (N - q) if d["advance"] else q                       # delay
        s = -1.0 if d["neg"] else 1.0
        out["wrong_dir"] = s * np.roll(x, -m)
        out["linear"] = s * lin_delay(x, m)
        if d["neg"]:
            out["dropped_pi"] = np.roll(x, m)
    elif fam == "reverse":
        out["plain_flip"] = x[::-1]
        out["identity"] = x.copy()
    elif fam == "revshift":
        m = round(d["t"][0] / d["t"][1] * N / 2)
        out["wrong_dir"] = x[(-n - m) % N]
        out["no_reversal"] = np.roll(x, m)
        out["plain_flip_shift"] = np.roll(x[::-1], m)
    elif fam == "half":
        alt = (-1.0) ** n * x
        sh = np.roll(x, N // 2)
        out["swapped_property"] = alt if d["alt"] else sh
        out["negated"] = -x
        if d["alt"]:
            out["linear"] = lin_delay(x, N // 2)
    elif fam == "list":
        m, s, sigma = d["m"], d["s"], d["sigma"]
        phase = np.exp(2j * np.pi * m * n / N)
        if sigma == -1:
            pre = s * (-1) ** m
            out["wrong_dir"] = pre * np.conj(phase) * np.roll(x, N // 2)
            out["linear"] = pre * phase * lin_delay(x, N // 2)
            out["no_shift"] = pre * phase * x
            out["sign"] = -pre * phase * np.roll(x, N // 2)
        else:
            out["wrong_dir"] = np.conj(phase) * x
            out["shift_instead"] = np.roll(x, m)
    elif fam == "square":
        c = np.convolve(x, x)
        out["truncated"] = c[:N]
        out["elementwise"] = x ** 2
        out["linear_full"] = c
    elif fam == "evenodd":
        r = x[(-n) % N]
        sg = 1 if d["even"] else -1
        out["no_half"] = x + sg * r
        out["plain_flip"] = (x + sg * x[::-1]) / 2
        out["identity"] = x.copy()
    return {k: np.asarray(v, dtype=complex) for k, v in out.items()}


def _differs(a, b):
    return a.shape != b.shape or np.max(np.abs(a - b)) > 1e-3


def _formula_y(x, c, ell, n0):
    N = len(x)
    n = np.arange(N)
    return c * np.exp(2j * np.pi * ell * n / N) * x[(n - n0) % N]


def _majority(lists):
    out = []
    for i in range(len(lists[0])):
        vals = [complex(l[i]) for l in lists]
        counts = sorted(((sum(abs(v - w) < 1e-9 for w in vals), v) for v in vals), key=lambda t: -t[0])
        tied = len(counts) > 1 and counts[0][0] == counts[1][0] and abs(counts[0][1] - counts[1][1]) > 1e-9
        out.append(None if tied else counts[0][1])
    return out


def check(params, correct):
    probs = []
    x = np.array(params["x"], float)
    N = len(x)
    if N not in (6, 8):
        probs.append(f"N = {N}")
    y = np.fft.ifft(_Y(params))
    if params["family"] == "listf":
        c, ell, n0 = correct["coef"], correct["ell"], correct["n0"]
        if c not in (1, -1) or not 0 <= ell < N or not 0 <= n0 < N:
            probs.append(f"(c, l, n0) = {(c, ell, n0)} out of range")
        elif not np.allclose(_formula_y(x, c, ell, n0), y, atol=1e-9):
            probs.append(f"listf: c e^(j2pi l n/N) x[<n-n0>] with {(c, ell, n0)} != ifft of the list {np.round(y, 4)}")
        hits = [(cc, ll, nn) for cc in (1, -1) for ll in range(N) for nn in range(N)
                if np.allclose(_formula_y(x, cc, ll, nn), y, atol=1e-9)]
        if len(hits) != 1:
            probs.append(f"listf: {len(hits)} triples (c, l, n0) reproduce y: {hits}")
        if len(_formula_mistakes(params, correct)) < 2:
            probs.append("listf: fewer than two wrong (c, l, n0) probes")
        shown = _parse_list(params["Y_def_tex"].split("=", 1)[1])
        if [i for _, i in shown] != params["list_idx"] or [s for s, _ in shown] != params["list_sgn"]:
            probs.append("displayed list does not match list_idx/list_sgn")
        probs += _check_b(params, correct)
        return probs
    got = _y_correct(correct)
    if got.shape != y.shape or not np.allclose(got, y, atol=1e-9):
        probs.append(f"{params['family']}: y = {got} but ifft of the defined Y gives {np.round(y, 6)}")
    if not np.allclose(y.real, np.round(y.real)) or not np.allclose(y.imag, np.round(y.imag)):
        probs.append(f"y = {y} has non-integer parts")
    if params["y_complex"] != bool(np.max(np.abs(y.imag)) > 1e-9):
        probs.append("y_complex flag is wrong")
    wrong = [k for k, v in _probes(params).items() if _differs(v, y)]
    if len(wrong) < 2:
        probs.append(f"{params['family']}: only {wrong} differ from the answer")
    if params["family"] == "list":
        shown = _parse_list(params["Y_def_tex"].split("=", 1)[1])
        if [i for _, i in shown] != params["list_idx"] or [s for s, _ in shown] != params["list_sgn"]:
            probs.append("displayed list does not match list_idx/list_sgn")
    probs += _check_b(params, correct)
    return probs


def _check_b(params, correct):
    probs = []
    opts = params["real_list"]
    real = []
    for o in opts:
        W = _parse_cplx_list(o["html"])
        w = np.fft.ifft(W)
        real.append(bool(np.max(np.abs(w.imag)) < 1e-9))
        if len(W) != params["M"]:
            probs.append(f"option {o['html']} has length {len(W)}")
    if sum(real) != 1:
        probs.append(f"{sum(real)} real lists among the options")
    elif opts[real.index(True)]["key"] != correct["real_list"]["key"]:
        probs.append("the real list is not the one marked correct")
    if len({o["html"] for o in opts}) != len(opts) or len(opts) != 4:
        probs.append("options are not 4 distinct lists")
    lists = [_parse_cplx_list(o["html"]) for o in opts]
    good = real.index(True) if sum(real) == 1 else 0
    M = params["M"]
    rev = (-np.arange(M)) % M
    odd_im = [bool(np.allclose(l.imag, -l.imag[rev])) for l in lists]
    even_re = [bool(np.allclose(l.real, l.real[rev])) for l in lists]
    if sum(odd_im) < 2:
        probs.append("(b) checking only that Im Z[k] is odd singles out the answer")
    if sum(even_re) < 2:
        probs.append("(b) checking only that Re Z[k] is even singles out the answer")
    maj = _majority(lists)
    if all(m is not None and abs(m - v) < 1e-9 for m, v in zip(maj, lists[good])):
        probs.append("(b) the per-entry majority vote over the options returns the correct list")
    agree = [sum(abs(a - b) < 1e-9 for j, l in enumerate(lists) if j != i for a, b in zip(lists[i], l)) for i in range(4)]
    if agree[good] > min(a for i, a in enumerate(agree) if i != good):
        probs.append(f"(b) the correct list shares the most entries with the others ({agree})")
    return probs


def _fmt(z):
    def num(v):
        return str(int(round(v))) if abs(v - round(v)) < 1e-9 else f"{v:.6g}"
    re_, im_ = float(np.real(z)), float(np.imag(z))
    if abs(im_) < 1e-9:
        return num(re_)
    jt = ("" if abs(abs(im_) - 1) < 1e-12 else num(abs(im_))) + "j"
    if abs(re_) < 1e-9:
        return ("-" if im_ < 0 else "") + jt
    return f"{num(re_)}{'-' if im_ < 0 else '+'}{jt}"


def _row(v, sep=", "):
    return "[" + sep.join(_fmt(z) for z in v) + "]"


def _formula_mistakes(params, correct):
    """Wrong (c, l, n0) triples: modulation the wrong way, sign, forgotten shift, list rotation read as a time shift."""
    N = params["N"]
    x = np.array(params["x"], float)
    c, ell, n0 = correct["coef"], correct["ell"], correct["n0"]
    right = _formula_y(x, c, ell, n0)
    cands = [(c, (N - ell) % N, n0), (-c, ell, n0), (c, ell, (n0 + N // 2) % N), (c, 0, ell)]
    out = []
    for t in cands:
        if t not in out and not np.allclose(_formula_y(x, *t), right, atol=1e-6):
            out.append(t)
    return out


def submissions(params, correct):
    if params["family"] == "listf":
        c, ell, n0 = correct["coef"], correct["ell"], correct["n0"]
        cases = [
            ({"coef": f" {c} ", "ell": f"+{ell}", "n0": f"0{n0}"}, {"coef": 1, "ell": 1, "n0": 1}),
            ({"coef": f"{c:+d}", "ell": f" {ell}", "n0": f"+{n0}"}, {"coef": 1, "ell": 1, "n0": 1}),
            ({"ell": f"{ell}.0"}, {"ell": "invalid"}),
        ]
        for t in _formula_mistakes(params, correct):
            exp = {k: (1 if v == w else 0) for k, v, w in zip(("coef", "ell", "n0"), t, (c, ell, n0))}
            cases.append(({"coef": str(t[0]), "ell": str(t[1]), "n0": str(t[2])}, exp))
        for o in params["real_list"]:
            cases.append(({"real_list": o["key"]}, {"real_list": 1 if o["key"] == correct["real_list"]["key"] else 0}))
        return cases
    y = _y_correct(correct)
    cases = [
        ({"y": _row(y, " ")}, {"y": 1}),                                         # MATLAB style
        ({"y": "[[" + ", ".join(_fmt(z) for z in y) + "]]"}, {"y": 1}),          # Python style
        ({"y": _row(y).replace("j", "i")}, {"y": 1}),                            # i for the imaginary unit
        ({"y": "[" + ", ".join(f"{z.real:.1f}" if abs(z.imag) < 1e-9 else _fmt(z) for z in y) + "]"}, {"y": 1}),
        ({"y": "[1, 2,, 3]"}, {"y": "invalid"}),
    ]
    for name, v in _probes(params).items():
        if _differs(v, y):
            cases.append(({"y": _row(v)}, {"y": 0}))
    for o in params["real_list"]:
        cases.append(({"real_list": o["key"]}, {"real_list": 1 if o["key"] == correct["real_list"]["key"] else 0}))
    return cases
