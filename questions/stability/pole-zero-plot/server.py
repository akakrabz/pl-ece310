"""Read H(z) off a pole-zero plot, then causality, stability and the ROC (Lecture 6 pole-zero plots,
Lecture 9 section 1.2 poles/zeros/gain, Lecture 11 Table 1 and slide 10; HW3 #4 run backwards;
FA2019 Midterm 1 #10, T/F FA2023 #1(f), SP2025 #1(e), FA2025 #1(c)).

The plot (ece310.svg.pole_zero) shows ALL finite poles and zeros, origin included: 2 or 3 poles (two
real, a conjugate pair with or without a real pole, a double pole with or without another real
pole, three real poles) at nice locations (+-1/4, +-1/2, +-1, +-3/2, +-2, +-j/2, (+-1 + j)/2, +-j,
+-3j/2, +-1 + +-j), and Z = P or P - 1 zeros (0, +-1, +-2, +-1/2, or the pair +-j). Then
    H(z) = K prod(z - z_k) / prod(z - p_k) = K z^{Z-P} prod_{z_k != 0}(1 - z_k z^-1) / prod(1 - p_k z^-1),
proper in z, so the causal ROC |z| > r_max includes infinity. The gain K (never 1) is fixed by b_0 = K
(only when Z = P), by H(1) = c or by H(-1) = c (no pole/zero at that point, c a nice rational).

Asked: (a) H(z), pl-symbolic-input in z with j allowed, graded numerically here (grade/test below:
any factored/expanded form, complex factors included; NaN/zoo/inf are wrong); (b) the ROC if causal;
(c) is the causal system stable (the solution calls it marginally stable only when the unit-circle
poles are simple; a double pole there makes h[n] grow like n); (d) the stable ROC or "none" (a pole ON
the unit circle: an ROC that only touches |z| = 1 does not contain it, Lecture 11); (e) causal /
anti-causal / non-causal (two-sided) / no stable system, the Lecture 11 Table 1 terms. When every pole
is outside the unit circle (stable ROC |z| < r_min) a zero at z = 0 is required, so h[0] = H(0) = 0 and
"anti-causal" is exact (the notes' n0 = 0 corner never arises). Stability cases are drawn 35/40/25
(all poles inside / some outside / one on |z| = 1). Layouts whose labels or markers would overlap are
widened or redrawn (_figure)."""

import cmath
import random
from fractions import Fraction as F

import prairielearn as pl
import prairielearn.sympy_utils as psu
from ece310 import fmt, svg

REAL = [F(1, 4), F(-1, 4), F(1, 2), F(-1, 2), F(1), F(-1), F(3, 2), F(-3, 2), F(2), F(-2)]
PAIRS = [(F(0), F(1, 2)), (F(1, 2), F(1, 2)), (F(-1, 2), F(1, 2)), (F(0), F(1)), (F(0), F(3, 2)),
         (F(1), F(1)), (F(-1), F(1))]
DOUBLE = [F(1, 2), F(-1, 2), F(-1), F(2)]
ZREAL = [F(0), F(0), F(1), F(-1), F(2), F(-2), F(1, 2), F(-1, 2)]
GAINS = [F(2), F(3), F(-1), F(-2), F(1, 2), F(3, 2), F(-1, 2), F(4)]
WIDTH = 400
TEST_Z = [1.37 * cmath.exp(0.41j), 0.83 * cmath.exp(2.2j), 2.71 * cmath.exp(-1.3j), 0.29 * cmath.exp(1.1j),
          1.9 * cmath.exp(2.9j), 0.61 * cmath.exp(-2.6j), 3.3 * cmath.exp(0.9j)]


# ----------------------------------------------------------------------------- drawing the system
def _r2(pt):
    return pt[0] * pt[0] + pt[1] * pt[1]


def _draw_poles(kind):
    while True:
        st = random.choice(["two_real", "two_real", "pair", "pair_real", "pair_real", "double", "double_real",
                            "three_real"])
        if st == "two_real":
            poles = [(r, F(0)) for r in random.sample(REAL, 2)]
        elif st == "pair":
            a, b = random.choice(PAIRS)
            poles = [(a, b), (a, -b)]
        elif st == "pair_real":
            a, b = random.choice(PAIRS)
            poles = [(a, b), (a, -b), (random.choice(REAL), F(0))]
        elif st == "double":
            r = random.choice(DOUBLE)
            poles = [(r, F(0)), (r, F(0))]
        elif st == "double_real":
            r = random.choice(DOUBLE)
            poles = [(r, F(0)), (r, F(0)), (random.choice([q for q in REAL if q != r]), F(0))]
        else:
            poles = [(r, F(0)) for r in random.sample(REAL, 3)]
        r2s = [_r2(p) for p in poles]
        got = "on" if any(v == 1 for v in r2s) else ("inside" if all(v < 1 for v in r2s) else "outside")
        if got == kind:
            return st, poles


def _draw_zeros(poles):
    P = len(poles)
    Z = P if random.random() < 0.65 else P - 1
    while True:
        zeros = []
        if Z >= 2 and random.random() < 0.25:
            zeros += [(F(0), F(1)), (F(0), F(-1))]
        while len(zeros) < Z:
            zeros.append((random.choice(ZREAL), F(0)))
        if any(zz in poles for zz in zeros):
            continue
        return zeros


def _H0(z0, poles, zeros):
    """prod(z0 - z_k) / prod(z0 - p_k) for real z0, exact (conjugate pairs give real products)."""
    val = F(1)
    for pts, sgn in ((zeros, 1), (poles, -1)):
        done = set()
        for i, (a, b) in enumerate(pts):
            if i in done:
                continue
            if b == 0:
                f = z0 - a
            else:
                j = next(k for k in range(len(pts)) if k not in done and k != i and pts[k] == (a, -b))
                done.add(j)
                f = (z0 - a) ** 2 + b * b
            done.add(i)
            val = val * f if sgn == 1 else val / f
    return val


def _nice(q):
    return q != 0 and q.denominator <= 16 and abs(q.numerator) <= 60


def _jlab(v):
    a = abs(v)
    s = "−" if v < 0 else ""
    if a == 1:
        return s + "j"
    if a.denominator == 1:
        return f"{s}{a.numerator}j"
    return s + ("j" if a.numerator == 1 else f"{a.numerator}j") + f"/{a.denominator}"


def _overlap(a, b):
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def _figure_at(poles, zeros, lim):
    """(svg, ok) for one plot range. Real-axis ticks come from svg.pole_zero; the imaginary-axis ticks
    and labels are added here, each label on the side of the axis that has no marker at that height
    (the library writes them on the left, where a pole at -1/2 + j/2 would collide with them).
    ok = False when any two of {markers, multiplicity labels, tick labels, "Re z", "Im z"} overlap
    (estimated text boxes: 7 px per character at 12 px)."""
    pts = poles + zeros
    re_vals = sorted({p[0] for p in pts if p[0] != 0} | {F(1), F(-1)})
    im_vals = sorted({p[1] for p in pts if p[1] != 0} | ({F(1), F(-1)} if any(p[1] != 0 for p in pts) else set()))
    scale = (WIDTH - 52) / (2 * lim)

    def X(v):
        return 26 + (float(v) + lim) * scale

    def Y(v):
        return WIDTH - 26 - (float(v) + lim) * scale
    xz, yz = X(0), Y(0)
    owned = []                                    # (box, owner): boxes of the same owner may touch
    mult = {}
    for q in pts:
        mult[q] = mult.get(q, 0) + 1
    for q, m in mult.items():
        cx, cy = X(q[0]), Y(q[1])
        owned.append(((cx - 9, cy - 9, cx + 9, cy + 9), q))
        if m > 1:
            owned.append(((cx + 8, cy - 17, cx + 17, cy - 5), q))
    owned.append(((X(lim) - 37, yz - 18, X(lim) - 4, yz - 5), "Re"))
    owned.append(((xz + 7, Y(lim) - 7, xz + 40, Y(lim) + 6), "Im"))
    for v in re_vals:
        w = 7 * len(svg.fmt_num(v))
        owned.append(((X(v) - w / 2 - 1, yz + 11, X(v) + w / 2 + 1, yz + 23), ("re", v)))
    extra = []
    ok = True
    for v in im_vals:
        lab = _jlab(v)
        y = Y(v)
        at = [float(a) for a, b in pts if b == v]
        off = 14 if any(a == 0 for a in at) else 6
        need = off + 7 * len(lab) + 10            # label extent plus a marker half-width and a gap
        left_free = all(not (a < 0 and -a * scale < need) for a in at)
        right_free = all(not (a > 0 and a * scale < need) for a in at)
        ok = ok and (left_free or right_free)
        w = 7 * len(lab)
        if left_free:
            side_x, anchor, box = xz - off, "end", (xz - off - w, y - 6, xz - off, y + 6)
        else:
            side_x, anchor, box = xz + off, "start", (xz + off, y - 6, xz + off + w, y + 6)
        owned.append((box, ("im", v)))
        extra.append(f'<line x1="{xz - 3:.1f}" y1="{y:.1f}" x2="{xz + 3:.1f}" y2="{y:.1f}" stroke="currentColor"/>')
        extra.append(f'<text x="{side_x:.1f}" y="{y + 4:.1f}" text-anchor="{anchor}" font-size="12">{lab}</text>')
    for i, (b1, o1) in enumerate(owned):
        for b2, o2 in owned[i + 1:]:
            if o1 != o2 and _overlap(b1, b2):
                ok = False
    fig = svg.pole_zero([complex(float(a), float(b)) for a, b in poles],
                        [complex(float(a), float(b)) for a, b in zeros],
                        lim=lim, ticks=[(float(v), svg.fmt_num(v)) for v in re_vals], iticks=[],
                        width=WIDTH, max_width=WIDTH + 20, title="pole-zero plot of H(z)")
    return fig.replace("</svg>", "".join(extra) + "</svg>"), ok


def _figure(poles, zeros):
    """Plot range = largest coordinate + 0.8 (room for the axis labels), widened in two steps if
    labels would still collide (e.g. the multiplicity "2" of a double root at z = 2 against "Re z")."""
    base = float(max(F(1), max(max(abs(a), abs(b)) for a, b in poles + zeros))) + 0.8
    for lim in (base, base + 0.25, base + 0.5):
        fig, ok = _figure_at(poles, zeros, lim)
        if ok:
            return fig, True
    return fig, False


def _draw():
    kind = random.choices(["inside", "outside", "on"], weights=[35, 40, 25])[0]
    while True:
        st, poles = _draw_poles(kind)
        zeros = _draw_zeros(poles)
        norms = []
        if len(zeros) == len(poles):
            norms += ["b0", "b0"]
        if (F(1), F(0)) not in poles + zeros:
            norms.append("dc")
        if (F(-1), F(0)) not in poles + zeros:
            norms.append("ny")
        if not norms:
            continue
        norm = random.choice(norms)
        gains = random.sample(GAINS, len(GAINS))
        K = c = None
        for g in gains:
            if norm == "b0":
                K, c = g, g
                break
            val = g * _H0(F(1) if norm == "dc" else F(-1), poles, zeros)
            if _nice(val):
                K, c = g, val
                break
        if K is None:
            continue
        if all(_r2(q) > 1 for q in poles) and (F(0), F(0)) not in zeros:
            continue        # the stable ROC is then |z| < r_min: a zero at z = 0 makes H(0) = h[0] = 0 (anti-causal)
        fig, ok = _figure(poles, zeros)
        if ok:
            return kind, st, poles, zeros, norm, K, c, fig


# ----------------------------------------------------------------------------- formatting
def _groups(pts):
    """[(re, im, multiplicity)] with conjugate pairs listed once (im > 0) and flagged."""
    out = []
    seen = []
    for a, b in pts:
        if b < 0:
            continue
        key = (a, b)
        if key in seen:
            out[seen.index(key)][2] += 1
        else:
            seen.append(key)
            out.append([a, b, 1])
    return out


def _loc_tex(a, b):
    if b == 0:
        return fmt.tex_num(a)
    jb = "j" if b == 1 else "j" + fmt.tex_num(b)
    return (r"\pm " + jb) if a == 0 else f"{fmt.tex_num(a)} \\pm {jb}"


def _zinv_factors(pts):
    """[(coeffs in z^-1, multiplicity)] for the nonzero roots: [1, -r] or [1, -2a, a^2 + b^2]."""
    out = []
    for a, b, m in _groups(pts):
        if b == 0 and a == 0:
            continue
        out.append(([F(1), -a] if b == 0 else [F(1), -2 * a, a * a + b * b], m))
    return out


def _z_factors(pts):
    """[(coeffs in descending powers of z, multiplicity)] for every root (z for a root at 0)."""
    out = []
    for a, b, m in _groups(pts):
        out.append(([F(1), -a] if b == 0 else [F(1), -2 * a, a * a + b * b], m))
    return out


def _expand(factors):
    """Product of the factors with multiplicity. No trimming: in descending powers of z a trailing
    zero coefficient is a real factor z (poly.pmul trims, which is right only for z^-1 polynomials)."""
    c = [F(1)]
    for f, m in factors:
        for _ in range(m):
            out = [F(0)] * (len(c) + len(f) - 1)
            for i, u in enumerate(c):
                for j, v in enumerate(f):
                    out[i + j] += u * v
            c = out
    return c


def _tex_zinv_factor(f, m):
    t = r"\left(" + fmt.tex_sum((v, fmt.zpow_tex(i)) for i, v in enumerate(f)) + r"\right)"
    return t + (f"^{{{m}}}" if m > 1 else "")


def _tex_z_factor(f, m):
    deg = len(f) - 1
    body = fmt.tex_sum((v, "z" if deg - i == 1 else (f"z^{{{deg - i}}}" if deg - i > 1 else "")) for i, v in enumerate(f))
    t = "z" if body == "z" else r"\left(" + body + r"\right)"
    return t + (f"^{{{m}}}" if m > 1 else "")


def _sym_poly(f, var_power):
    """SymPy string of sum f[i] * var_power(i)."""
    terms = []
    for i, v in enumerate(f):
        if v == 0:
            continue
        pw = var_power(i)
        terms.append(fmt.sym_num(v) if pw == "1" else f"{fmt.sym_num(v)}*{pw}")
    return "(" + (" + ".join(terms) if terms else "0") + ")"


def _zinv_pow(i):
    return "1" if i == 0 else f"z**(-{i})"


def _sym_factored_zinv(K, d, numf, denf):
    num = "*".join(_sym_poly(f, _zinv_pow) + (f"**{m}" if m > 1 else "") for f, m in numf) or "1"
    den = "*".join(_sym_poly(f, _zinv_pow) + (f"**{m}" if m > 1 else "") for f, m in denf) or "1"
    shift = f"*z**(-{d})" if d else ""
    return f"{fmt.sym_num(K)}{shift}*{num}/({den})"


def _sym_desc(c):
    deg = len(c) - 1
    return _sym_poly(c, lambda i: "1" if deg - i == 0 else f"z**({deg - i})")


def _cplx_root_str(a, b):
    return f"({fmt.sym_num(a)} + {fmt.sym_num(b)}*j)"


def _sym_cplx(K, poles, zeros):
    """K * prod(z - z_k) / prod(z - p_k) with one linear factor per root (j for complex roots)."""
    def lin(pts):
        fs = []
        for a, b in pts:
            fs.append(f"(z - {_cplx_root_str(a, b)})" if b != 0 else ("z" if a == 0 else f"(z - {fmt.sym_num(a)})"))
        return "*".join(fs) or "1"
    return f"{fmt.sym_num(K)}*{lin(zeros)}/({lin(poles)})"


def _H_value(params, z):
    K = float(F(params["K"]))
    val = complex(K)
    for a, b in params["zeros"]:
        val *= z - complex(float(F(a)), float(F(b)))
    for a, b in params["poles"]:
        val /= z - complex(float(F(a)), float(F(b)))
    return val


def _r_tex(r2):
    """Radius sqrt(r2) for r2 in {q^2, 2 q^2}."""
    for mult, suffix in ((1, ""), (2, r"\sqrt{2}")):
        q = r2 / mult
        n, d = q.numerator, q.denominator
        rn, rd = round(n ** 0.5), round(d ** 0.5)
        if rn * rn == n and rd * rd == d:
            root = F(rn, rd)
            if not suffix:
                return fmt.tex_num(root)
            if root.denominator == 1:
                return suffix if root == 1 else f"{root.numerator}{suffix}"
            top = suffix if root.numerator == 1 else f"{root.numerator}{suffix}"
            return r"\tfrac{" + top + "}{" + str(root.denominator) + "}"
    raise ValueError(r2)


def generate(data):
    kind, st, poles, zeros, norm, K, c, fig = _draw()
    P, Z = len(poles), len(zeros)
    d = P - Z
    numf, denf = _zinv_factors(zeros), _zinv_factors(poles)
    N = _expand(numf)
    D = _expand(denf)
    B = [F(0)] * d + [K * v for v in N]
    p = data["params"]
    p["poles"] = [[str(a), str(b)] for a, b in poles]
    p["zeros"] = [[str(a), str(b)] for a, b in zeros]
    p["K"] = str(K)
    p["norm"], p["structure"], p["kind"] = norm, st, kind
    p["fig_html"] = fig
    p["is_b0"], p["is_dc"], p["is_ny"] = norm == "b0", norm == "dc", norm == "ny"
    p["K_tex"] = fmt.tex_num(K)
    p["c_tex"] = fmt.tex_num(c)
    z0 = F(1) if norm == "dc" else F(-1)
    p["z0_tex"] = "1" if norm == "dc" else "-1"

    # ---- solution pieces
    pg, zg = _groups(poles), _groups(zeros)

    def listing(groups):
        items = []
        for a, b, m in groups:
            t = _loc_tex(a, b)
            items.append(t + (r"\ (\text{double})" if m == 2 else ""))
        return r",\quad ".join(items)
    p["poles_list_tex"] = listing(pg)
    p["zeros_list_tex"] = listing(zg)
    p["P"], p["Z"] = P, Z
    p["count_text"] = f"{Z} zero{'' if Z == 1 else 's'} and {P} poles, counted with multiplicity"
    p["has_delay"] = d > 0
    p["d_tex"] = fmt.zpow_tex(d)
    num_z = "".join(_tex_z_factor(f, m) for f, m in _z_factors(zeros)) or "1"
    den_z = "".join(_tex_z_factor(f, m) for f, m in _z_factors(poles))
    p["Hz_K_tex"] = r"K\," + fmt.tex_frac(num_z, den_z)
    num_zi = "".join(_tex_zinv_factor(f, m) for f, m in numf) or "1"
    den_zi = "".join(_tex_zinv_factor(f, m) for f, m in denf)
    p["Hzi_K_tex"] = r"K\," + fmt.zpow_tex(d) + fmt.tex_frac(num_zi, den_zi)
    p["has_pairs"] = any(b != 0 for a, b in poles + zeros)
    pair_lines = []
    for a, b, m in pg + zg:
        if b != 0:
            lhs = (r"\left(" + fmt.tex_sum([(1, "z"), (-a, ""), (-b, "j")]) + r"\right)\left("
                   + fmt.tex_sum([(1, "z"), (-a, ""), (b, "j")]) + r"\right)")
            pair_lines.append(lhs + " = " + _tex_z_factor([F(1), -2 * a, a * a + b * b], 1).replace(r"\left(", "").replace(r"\right)", ""))
    p["pair_lines"] = [{"line": t} for t in dict.fromkeys(pair_lines)]
    if norm == "b0":
        p["gain_tex"] = rf"b_0 = K = {fmt.tex_num(K)}"
    else:
        h0 = _H0(z0, poles, zeros)
        p["H0_tex"] = fmt.tex_num(h0, paren=True)
        p["gain_tex"] = (rf"H({p['z0_tex']}) = K\cdot {fmt.tex_num(h0, paren=True)} = {fmt.tex_num(c)}"
                         rf"\ \Rightarrow\ K = {fmt.tex_num(K)}")
    lead = "" if K == 1 else ("-" if K == -1 else fmt.tex_num(K))
    p["H_fact_tex"] = lead + fmt.zpow_tex(d) + fmt.tex_frac(num_zi, den_zi)
    p["H_exp_tex"] = fmt.tex_frac(fmt.tex_poly_zinv(B), fmt.tex_poly_zinv(D))
    p["H_pos_tex"] = (lead + r"\," if lead else "") + fmt.tex_frac(num_z, den_z)

    # ---- symbolic strings (Test button, checker)
    p["H_sym"] = _sym_factored_zinv(K, d, numf, denf)
    Nz, Dz = _expand(_z_factors(zeros)), _expand(_z_factors(poles))
    p["H_pos_sym"] = f"{fmt.sym_num(K)}*{_sym_desc(Nz)}/{_sym_desc(Dz)}"
    p["H_cplx_sym"] = _sym_cplx(K, poles, zeros)
    p["H_swap_sym"] = f"{fmt.sym_num(K)}*{_sym_desc(Dz)}/{_sym_desc(Nz)}"

    # ---- ROCs
    rs = sorted({_r2(q) for q in poles})
    rocs = [(None, rs[0])] + list(zip(rs, rs[1:])) + [(rs[-1], None)]
    Z_ = r"\lvert z\rvert"

    def roc_tex(lo, hi):
        if lo is None:
            return f"{Z_} \\lt {_r_tex(hi)}"
        if hi is None:
            return f"{Z_} \\gt {_r_tex(lo)}"
        return f"{_r_tex(lo)} \\lt {Z_} \\lt {_r_tex(hi)}"

    def has_unit(lo, hi):
        return (lo is None or lo < 1) and (hi is None or hi > 1)
    stable = next(((lo, hi) for lo, hi in rocs if has_unit(lo, hi)), None)
    causal_stable = rs[-1] < 1
    p["roc_c_choices"] = [{"text": f"${roc_tex(lo, hi)}$", "correct": hi is None} for lo, hi in rocs]
    p["stab_c_choices"] = [{"text": "Yes", "correct": causal_stable}, {"text": "No", "correct": not causal_stable}]
    p["roc_s_choices"] = [{"text": f"${roc_tex(lo, hi)}$", "correct": (lo, hi) == stable} for lo, hi in rocs] + \
        [{"text": "none: no ROC of this $H(z)$ gives a BIBO stable system", "correct": stable is None}]
    skind = "none" if stable is None else ("causal" if stable[1] is None else ("left" if stable[0] is None else "two"))
    p["kind_s_choices"] = [
        {"text": r"causal: $h[n] = 0$ for $n \lt 0$ (right-sided)", "correct": skind == "causal"},
        {"text": r"anti-causal: $h[n] = 0$ for $n \ge 0$ (left-sided)", "correct": skind == "left"},
        {"text": r"non-causal: $h[n]$ is two-sided", "correct": skind == "two"},
        {"text": "there is no BIBO stable system with this $H(z)$", "correct": skind == "none"},
    ]
    p["mags_tex"] = r",\ ".join(_r_tex(v) for v in rs)
    p["nrocs"] = len(rocs)
    p["roc_list"] = [{"tex": roc_tex(lo, hi), "note": ("contains" if has_unit(lo, hi) and not (lo == 1 or hi == 1)
                                                       else "does not contain")} for lo, hi in rocs]
    p["causal_roc_tex"] = roc_tex(rs[-1], None)
    p["rmax_tex"] = _r_tex(rs[-1])
    p["causal_stable"] = causal_stable
    p["causal_unstable"] = not causal_stable
    p["rmax_on"] = rs[-1] == 1
    p["stable_none"] = stable is None
    p["stable_exists"] = stable is not None
    p["stable_roc_tex"] = roc_tex(*stable) if stable else ""
    p["skind_text"] = {
        "causal": "the exterior of the outermost pole circle, including $z = \\infty$: $h[n]$ is right-sided and causal",
        "left": ("the disc inside the innermost pole circle, including $z = 0$, so $h[n] = 0$ for $n \\gt 0$, and "
                 "$h[0] = H(0) = 0$ because of the zero at $z = 0$. Hence $h[n] = 0$ for $n \\ge 0$: the system is anti-causal "
                 "(Lecture 11, Table 1)"),
        "two": "a ring between two pole circles: $h[n]$ is two-sided, so the system is non-causal",
        "none": ""}[skind]
    on = [g for g in pg if g[0] * g[0] + g[1] * g[1] == 1]
    p["on_tex"] = r",\ ".join(_loc_tex(a, b) for a, b, m in on)
    if len(on) == 1 and on[0][1] == 0:
        p["on_words"] = "The double pole at" if on[0][2] == 2 else "The pole at"
        p["on_verb"] = "lies"
    else:
        p["on_words"], p["on_verb"] = "The poles at", "lie"
    # marginal stability needs SIMPLE poles on the unit circle; a repeated one makes h[n] grow like n
    p["marginal"] = rs[-1] == 1 and all(m == 1 for a, b, m in on)
    p["on_repeated"] = rs[-1] == 1 and any(m > 1 for a, b, m in on)


# ----------------------------------------------------------------------------- grading of H(z)
def _matches(expr, params):
    syms = {s.name: s for s in expr.free_symbols}
    if set(syms) - {"z"}:
        return False
    for zt in TEST_Z:
        try:
            val = complex(expr.subs(syms["z"], zt).evalf(30)) if syms else complex(expr.evalf(30))
        except (TypeError, ValueError, ZeroDivisionError, OverflowError):
            return False
        tru = _H_value(params, zt)
        if not abs(val - tru) <= 1e-8 * max(1.0, abs(tru)):
            return False
    return True


def grade(data):
    sub = data["submitted_answers"].get("H")
    if sub is None or "H" in data["format_errors"]:
        return
    try:
        expr = (psu.json_to_sympy(sub, allow_complex=True) if isinstance(sub, dict)
                else psu.convert_string_to_sympy(sub, ["z"], allow_complex=True))
        ok = _matches(expr, data["params"])
    except Exception:
        ok = False
    data["partial_scores"]["H"] = {"score": 1.0 if ok else 0.0, "weight": 1}
    pl.set_weighted_score_data(data)


def test(data):
    tt = data["test_type"]
    if tt == "correct":
        data["raw_submitted_answers"]["H"] = data["params"]["H_sym"]
        data["partial_scores"]["H"] = {"score": 1.0, "weight": 1}
    elif tt == "incorrect":
        data["raw_submitted_answers"]["H"] = data["params"]["H_swap_sym"]
        data["partial_scores"]["H"] = {"score": 0.0, "weight": 1}
    if tt in ("correct", "incorrect"):
        pl.set_weighted_score_data(data)
