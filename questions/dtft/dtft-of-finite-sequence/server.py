"""DTFT X_d(w) of a short finite sequence from the definition, plus |X_d(w0)| at w0 in {0, pi/2, pi}
(HW5 Problem 1(a), Fall 2026; Lecture 13 DTFT definition, Lecture 14 examples {1, 2, 3}, {1, 0, 0, -1}).

Three families:
  sym   odd-length sequence symmetric about n_c != 0:  X_d = e^{-j w n_c} (b + 2a cos w [+ 2a' cos 2w])
  anti  {a, 0, -a} centred at n_c:                       X_d = 2 j a sin(w) e^{-j w n_c}
  gen   2-4 arbitrary integer samples, not even about n = 0.
No family is even about n = 0, so the sign error e^{+jwn} (= X_d(-w) = X_d^*(w)) is always wrong.

Grading of X_d: SymPy's built-in equality test cannot prove e.g. e^{-jw}(1 + 4cos w) =
2 + e^{-jw} + 2e^{-2jw}, so `pl-symbolic-input` is used as a parser only (no correct_answers entry
for it; the element then defers testing to test() below) and grade() compares the parsed
expression numerically with the true DTFT at several frequencies."""

import cmath
import random

import prairielearn as pl
import prairielearn.sympy_utils as psu
from ece310 import fmt

W_TEST = [0.37, 1.23, 2.11, -0.81, 2.93, -2.47, 0.05]
W0_TEX = {"0": "0", "pi/2": r"\tfrac{\pi}{2}", "pi": r"\pi"}


# ----------------------------------------------------------------------------- formatting helpers
def _phase_tex(n):
    """e^{-j w n} in LaTeX ('' for n = 0)."""
    if n == 0:
        return ""
    k = "" if abs(n) == 1 else str(abs(n))
    return f"e^{{{'-' if n > 0 else ''}j{k}\\omega}}"


def _phase_sym(n):
    return "1" if n == 0 else f"exp(-j*({n})*w)"


def _cplx_tex(re, im):
    if im == 0:
        return fmt.tex_num(re)
    jt = ("" if abs(im) == 1 else str(abs(im))) + "j"
    if re == 0:
        return ("-" if im < 0 else "") + jt
    return f"{fmt.tex_num(re)} {'-' if im < 0 else '+'} {jt}"


def _unit(w0, n):
    """e^{-j w0 n} as a Gaussian integer (re, im) for w0 in {0, pi/2, pi}."""
    if w0 == "0":
        return (1, 0)
    if w0 == "pi":
        return (1 if n % 2 == 0 else -1, 0)
    return {0: (1, 0), 1: (0, -1), 2: (-1, 0), 3: (0, 1)}[n % 4]


def _isqrt_exact(m):
    r = int(round(m ** 0.5))
    return r if r * r == m else None


# ----------------------------------------------------------------------------- sequences
def _draw():
    fam = random.choices(["sym", "anti", "gen"], weights=[38, 14, 48])[0]
    if fam == "sym":
        nc = random.choice([-3, -2, -1, 1, 2, 3])
        if random.random() < 0.6:
            a, b = random.choice([-3, -2, -1, 1, 2, 3]), random.randint(-4, 4)
            x, amp = [a, b, a], [(b, 0), (2 * a, 1)]
        else:
            a, b, c = random.choice([-2, -1, 1, 2]), random.choice([-3, -2, -1, 1, 2, 3]), random.randint(-3, 3)
            x, amp = [a, b, c, b, a], [(c, 0), (2 * b, 1), (2 * a, 2)]
        start = nc - (len(x) - 1) // 2
        return fam, x, start, nc, amp
    if fam == "anti":
        nc = random.randint(-2, 2)
        if random.random() < 0.3:
            a = random.choice([-3, -2, -1, 1, 2, 3])
            x, amp = [a, 0, -a], [(2 * a, 1)]                      # 2j a sin(w)
        else:
            a, b = random.choice([-3, -2, -1, 1, 2, 3]), random.choice([-3, -2, -1, 1, 2, 3])
            x, amp = [a, b, 0, -b, -a], [(2 * b, 1), (2 * a, 2)]   # 2j (b sin w + a sin 2w)
        return fam, x, nc - (len(x) - 1) // 2, nc, amp
    while True:
        L = random.choice([2, 3, 3, 4, 4])
        x = [random.randint(-3, 3) for _ in range(L)]
        x[0], x[-1] = random.choice([-3, -2, -1, 1, 2, 3]), random.choice([-3, -2, -1, 1, 2, 3])
        start = random.randint(-3, 2)
        if x == x[::-1] or x == [-v for v in x[::-1]]:        # (anti)symmetric: those are the other families
            continue
        d = {start + i: v for i, v in enumerate(x)}
        if any(d.get(n, 0) != d.get(-n, 0) for n in d):      # not even about n = 0
            return fam, x, start, None, None


def _value(x, start, w0):
    re = im = 0
    for i, v in enumerate(x):
        ur, ui = _unit(w0, start + i)
        re, im = re + v * ur, im + v * ui
    return re, im


def _mistakes(fam, x, start, amp, w0):
    """Classic wrong values for |X_d(w0)| (the checker recomputes these independently)."""
    re, im = _value(x, start, w0)
    mag = _isqrt_exact(re * re + im * im)
    out = []
    if im == 0 and re < 0:
        out.append(re)                                   # magnitude without the absolute value
    if w0 == "pi":
        out.append(sum(x))                               # forgot the (-1)^n
    if w0 == "pi/2":
        out.append(re * re + im * im)                    # forgot the square root
        if re != 0 and im != 0:
            out.append(abs(re) + abs(im))
    if w0 == "0":
        out.append(sum(abs(v) for v in x))               # added magnitudes
    if amp is not None and fam == "sym":                 # forgot the 2 in 2a cos(kw)
        cosk = {"0": lambda k: 1, "pi": lambda k: (-1) ** k, "pi/2": lambda k: [1, 0, -1, 0][k % 4]}[w0]
        out.append(abs(sum((c if k == 0 else c // 2) * cosk(k) for c, k in amp)))
    if fam == "anti" and w0 == "pi/2":
        out.append(abs(x[(len(x) - 1) // 2 - 1]))        # forgot the 2 in 2ja sin w
    uniq = []
    for m in out:
        if m != mag and m not in uniq:
            uniq.append(m)
    return uniq


def generate(data):
    while True:
        fam, x, start, nc, amp = _draw()
        options = []
        for w0 in ("0", "pi/2", "pi"):
            re, im = _value(x, start, w0)
            mag = _isqrt_exact(re * re + im * im)
            if mag and mag > 0 and len(_mistakes(fam, x, start, amp, w0)) >= 2:
                options.append(w0)
        if options:
            w0 = random.choice(options)
            break

    ns = list(range(start, start + len(x)))
    re, im = _value(x, start, w0)
    mag = _isqrt_exact(re * re + im * im)

    exp_tex = fmt.tex_sum((v, _phase_tex(n)) for n, v in zip(ns, x))
    exp_sym = " + ".join(f"{fmt.sym_num(v)}*{_phase_sym(n)}" for n, v in zip(ns, x) if v != 0)
    p = data["params"]
    p["x"], p["start"], p["family"], p["w0"] = x, start, fam, w0
    p["x_tex"] = fmt.tex_seq(x, start)
    p["w0_tex"] = W0_TEX[w0]
    p["def_tex"] = " + ".join(
        f"{fmt.tex_num(v, paren=True)}\\,e^{{-j\\omega({n})}}" if n < 0 else f"{fmt.tex_num(v, paren=True)}\\,e^{{-j\\omega\\cdot {n}}}"
        for n, v in zip(ns, x) if v != 0)
    p["exp_tex"] = exp_tex
    p["X_sym"] = exp_sym
    p["X_sym_conj"] = " + ".join(f"{fmt.sym_num(v)}*{_phase_sym(-n)}" for n, v in zip(ns, x) if v != 0)
    p["closed"] = fam in ("sym", "anti")
    if fam == "sym":
        cos_body = {0: "", 1: r"\cos\omega", 2: r"\cos 2\omega"}
        amp_tex = fmt.tex_sum((c, cos_body[k]) for c, k in amp)
        half = (len(x) - 1) // 2
        pair_terms = []
        if x[half] != 0:
            pair_terms.append((x[half], ""))
        for k in range(1, half + 1):
            body = r"\left(e^{j\omega}+e^{-j\omega}\right)" if k == 1 else r"\left(e^{j2\omega}+e^{-j2\omega}\right)"
            pair_terms.append((x[half - k], body))
        p["group_tex"] = fmt.tex_sum(pair_terms)
        p["closed_tex"] = (_phase_tex(nc) + r"\left(" + amp_tex + r"\right)") if nc != 0 else amp_tex
        p["amp_tex"] = amp_tex
        p["euler_tex"] = r"e^{jk\omega}+e^{-jk\omega} = 2\cos k\omega"
        p["closed_sym"] = f"{_phase_sym(nc)}*(" + " + ".join(
            f"{fmt.sym_num(c)}*cos({k}*w)" if k else fmt.sym_num(c) for c, k in amp if c != 0) + ")"
    elif fam == "anti":
        half = (len(x) - 1) // 2
        sin_body = {1: r"j\sin\omega", 2: r"j\sin 2\omega"}
        diff_body = {1: r"\left(e^{j\omega}-e^{-j\omega}\right)", 2: r"\left(e^{j2\omega}-e^{-j2\omega}\right)"}
        p["group_tex"] = fmt.tex_sum((x[half - k], diff_body[k]) for k in range(1, half + 1))
        amp_tex = fmt.tex_sum((c, sin_body[k]) for c, k in amp)
        p["amp_tex"] = amp_tex
        p["closed_tex"] = (_phase_tex(nc) + r"\left(" + amp_tex + r"\right)") if nc != 0 else amp_tex
        p["euler_tex"] = r"e^{jk\omega}-e^{-jk\omega} = 2j\sin k\omega"
        p["closed_sym"] = f"{_phase_sym(nc)}*(" + " + ".join(
            f"{fmt.sym_num(c)}*j*sin({k}*w)" for c, k in amp) + ")"
    else:
        p["closed_sym"] = " + ".join(
            f"{fmt.sym_num(v)}*(cos({n}*w) - j*sin({n}*w))" for n, v in zip(ns, x) if v != 0)
    if p["closed"]:
        p["phase_nc_tex"] = _phase_tex(nc)
        p["center_tex"] = str(nc)
        p["symmetry_word"] = "symmetric" if fam == "sym" else "antisymmetric"

    # magnitude at w0
    unit_tex = {"0": "1", "pi": "(-1)^n", "pi/2": "(-j)^n"}[w0]
    p["unit_tex"] = unit_tex
    p["terms_w0_tex"] = fmt.tex_sum((v * _unit(w0, n)[0], "") for n, v in zip(ns, x) if _unit(w0, n)[1] == 0) \
        if w0 != "pi/2" else " + ".join(
            f"{fmt.tex_num(v, paren=True)}\\cdot\\left({_cplx_tex(*_unit(w0, n))}\\right)" for n, v in zip(ns, x) if v != 0)
    p["Xw0_tex"] = _cplx_tex(re, im)
    p["mag_tex"] = fmt.tex_num(mag)
    p["mag_work_tex"] = (rf"\left|{_cplx_tex(re, im)}\right| = {mag}" if im == 0
                         else rf"\sqrt{{{fmt.tex_num(re, paren=True)}^2 + {fmt.tex_num(im, paren=True)}^2}} = \sqrt{{{re * re + im * im}}} = {mag}")
    if p["closed"]:
        cosv = {"0": {1: 1, 2: 1}, "pi": {1: -1, 2: 1}, "pi/2": {1: 0, 2: -1}}[w0]
        sinv = {"0": 0, "pi": 0, "pi/2": 1}[w0]
        if fam == "sym":
            amp_val = sum(c * (cosv[k] if k else 1) for c, k in amp)
        else:
            amp_val = 2 * x[(len(x) - 1) // 2 - 1] * sinv     # sin(2 w0) = 0 at w0 = pi/2
        p["amp_val_tex"] = fmt.tex_num(amp_val) if fam == "sym" else _cplx_tex(0, amp_val)
    data["correct_answers"]["mag"] = float(mag)


# ----------------------------------------------------------------------------- grading of X_d
def _true_X(params, w):
    return sum(v * cmath.exp(-1j * w * (params["start"] + i)) for i, v in enumerate(params["x"]))


def _matches(expr, params):
    syms = {s.name: s for s in expr.free_symbols}
    if set(syms) - {"w"}:
        return False
    for w in W_TEST:
        try:
            val = complex(expr.subs(syms["w"], w).evalf(30)) if syms else complex(expr.evalf(30))
        except (TypeError, ValueError, ZeroDivisionError, OverflowError):
            return False
        tru = _true_X(params, w)
        if not abs(val - tru) <= 1e-8 * max(1.0, abs(tru)):
            return False
    return True


def grade(data):
    sub = data["submitted_answers"].get("X")
    if sub is None or "X" in data["format_errors"]:
        return
    try:
        expr = (psu.json_to_sympy(sub, allow_complex=True) if isinstance(sub, dict)
                else psu.convert_string_to_sympy(sub, ["w"], allow_complex=True))
        ok = _matches(expr, data["params"])
    except Exception:
        ok = False
    data["partial_scores"]["X"] = {"score": 1.0 if ok else 0.0, "weight": 1}
    pl.set_weighted_score_data(data)


def test(data):
    tt = data["test_type"]
    if tt == "correct":
        data["raw_submitted_answers"]["X"] = data["params"]["X_sym"]
        data["partial_scores"]["X"] = {"score": 1.0, "weight": 1}
    elif tt == "incorrect":
        data["raw_submitted_answers"]["X"] = data["params"]["X_sym_conj"]
        data["partial_scores"]["X"] = {"score": 0.0, "weight": 1}
    if tt in ("correct", "incorrect"):
        pl.set_weighted_score_data(data)
