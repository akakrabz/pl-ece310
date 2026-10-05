"""Independent check of questions/lccde/lccde-to-transfer-function.

to_H : the correct H(z) string is evaluated at test points and compared with sum_n h[n] z^{-n},
       h from scipy.signal.lfilter run on the (b, a) of the displayed equation.
to_ba: the correct vectors must reproduce the displayed H(z) (sympy parse of the given LaTeX is
       not available, so the generator's (b, a) are compared with numpy polynomial evaluation of the
       vectors at test points, and a[0] == 1, last entries nonzero)."""

from fractions import Fraction

import numpy as np
from checklib import S, Z, eval_z, impulse_response

TEST_Z = [1.7 * np.exp(0.4j), 2.9 * np.exp(-2.2j), 4.1 * np.exp(1.3j)]


def _ba(params):
    return [Fraction(v) for v in params["b"]], [Fraction(v) for v in params["a"]]


def _poly_inv(c, z):
    return sum(complex(v) * z ** (-k) for k, v in enumerate(c))


def check(params, correct):
    probs = []
    b, a = _ba(params)
    if params["to_H"]:
        h = impulse_response(b, a, 400)
        for z in TEST_Z:
            if abs(z) <= max([abs(r) for r in np.roots([float(v) for v in a])] + [0]) * 1.2 + 0.1:
                continue
            want = sum(h[n] * z ** (-n) for n in range(400))
            got = eval_z(correct["H"], z)
            if abs(got - want) > 1e-7 * max(1, abs(want)):
                probs.append(f"H({z:.3f}) = {got} but sum h[n] z^-n = {want}")
    else:
        bv = np.array(correct["bvec"]["_value"][0], dtype=float)
        av = np.array(correct["avec"]["_value"][0], dtype=float)
        if abs(av[0] - 1) > 1e-12:
            probs.append(f"a[0] = {av[0]} (must be normalised to 1)")
        if bv[-1] == 0 or av[-1] == 0:
            probs.append("a trailing zero in b or a")
        for z in TEST_Z:
            want = _poly_inv(b, z) / _poly_inv(a, z)
            got = _poly_inv(bv, z) / _poly_inv(av, z)
            if abs(got - want) > 1e-9 * max(1, abs(want)):
                probs.append(f"vectors give H({z:.3f}) = {got}, generator H = {want}")
        # decimals must be exact (dyadic) so that typed decimals are exact
        for v in list(bv) + list(av):
            if abs(v * 64 - round(v * 64)) > 1e-12:
                probs.append(f"coefficient {v} has no short exact decimal")
        # the displayed H(z) must be the same function: impulse response from lfilter of the vectors
        h1 = impulse_response(bv, av, 30)
        h2 = impulse_response(b, a, 30)
        if not np.allclose(h1, h2):
            probs.append("impulse responses differ")
    return probs


def _fmt_vec(v, sep=", "):
    return "[" + sep.join(repr(float(x)) if float(x) != int(x) else str(int(x)) for x in v) + "]"


def _sym_poly(c, top=None):
    """sum c_k z^{-k} (top=None) or sum c_k z^{top-k} as a SymPy string."""
    terms = []
    for k, v in enumerate(c):
        f = Fraction(v)
        if f == 0:
            continue
        e = -k if top is None else top - k
        terms.append(f"({f.numerator}/{f.denominator})*z^({e})")
    return "(" + " + ".join(terms) + ")"


def submissions(params, correct):
    b, a = _ba(params)
    cases = []
    if params["to_H"]:
        top = max(len(b), len(a)) - 1
        cases.append(({"H": f"{_sym_poly(b, top)}/{_sym_poly(a, top)}"}, {"H": 1}))            # positive powers of z
        cases.append(({"H": f"3*{_sym_poly(b)}/(3*{_sym_poly(a)})"}, {"H": 1}))                  # scaled num and den
        a_bad = [a[0]] + [-v for v in a[1:]]
        cases.append(({"H": f"{_sym_poly(b)}/{_sym_poly(a_bad)}"}, {"H": 0}))                    # sign error on a_k
        cases.append(({"H": f"{_sym_poly(a)}/{_sym_poly(b)}"}, {"H": 0}))                        # upside down
        cases.append(({"H": "1/(1 - 0.5*z^(-1))"}, {"H": "invalid"}))                            # decimals refused
        # sanity: the two equivalent strings really are equal to the correct answer
        assert S(f"{_sym_poly(b, top)}/{_sym_poly(a, top)}".replace("^", "**")).equals(S(correct["H"]))
    else:
        bv = [float(v) for v in b]
        av = [float(v) for v in a]
        cases.append(({"bvec": _fmt_vec(bv, " "), "avec": _fmt_vec(av, " ")}, {"bvec": 1, "avec": 1}))   # MATLAB style
        cases.append(({"bvec": _fmt_vec(bv), "avec": _fmt_vec(av)}, {"bvec": 1, "avec": 1}))             # comma style
        sig4 = lambda v: "[" + ", ".join(f"{x:.4g}" for x in v) + "]"
        cases.append(({"bvec": sig4(bv), "avec": sig4(av)}, {"bvec": 1, "avec": 1}))                     # 4 significant digits
        if len(av) > 1:
            a_bad = [av[0]] + [-v for v in av[1:]]
            cases.append(({"avec": _fmt_vec(a_bad)}, {"avec": 0, "bvec": 1}))                            # sign error on a_k
        c = params["c"] if params["c"] != 1 else 2
        cases.append(({"bvec": _fmt_vec([c * v for v in bv]), "avec": _fmt_vec([c * v for v in av])},
                      {"bvec": 0, "avec": 0}))                                                           # not normalised
        if params["has_lead_zeros"]:
            cases.append(({"bvec": _fmt_vec(bv[params["lead_zeros"]:])}, {"bvec": 0}))                  # dropped b0 = 0
        cases.append(({"bvec": _fmt_vec(bv + [0.0])}, {"bvec": 0}))                                      # trailing zero (wrong length)
        cases.append(({"avec": "[1, 1/2]"}, {"avec": "invalid"}))                                        # fractions refused
    return cases
