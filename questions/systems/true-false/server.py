"""Six True/False statements (exam family "True/False", problem 1 on 7 of 7 past Midterm 1 exams).

Statements come from the study site's True/False bank (official answers). Where a statement has a
parameter (a pole, a gain, an index), it is randomized and the answer recomputed: e.g. FA2023 #1(f)
"a stable H = 1/(1 - p z^-1) must be non-causal" is True iff |p| > 1. At most one statement per topic
group, six groups out of eight."""

import random
from fractions import Fraction as Fr

from ece310 import fmt, mc


def T(x):
    return fmt.tex_num(x)


def zfac(p):
    """1 - p z^{-1} in LaTeX."""
    return fmt.tex_sum([(1, ""), (-Fr(p), "z^{-1}")])


# ----------------------------------------------------------------------------- fixed statements
# sid: (source, statement, answer, reason)
FIXED = {
    "h_any_system": ("FA2024 #1(a)", r"For a discrete-time system with impulse response $h[n]$, the system output $y[n]$ to any input signal $x[n]$ is always determined using $h[n]$.",
                     False, r"Only an LTI system is determined by $h[n]$: $y[n] = n\,x[n]$ has $h[n] = n\,\delta[n] = 0$, yet it maps $\delta[n-1]$ to $\delta[n-1]$."),
    "h_conv_any": ("FA2023 #1(d)", r"For a discrete-time system with impulse response $h[n]$, the output to any input $x[n]$ is always given by $y[n] = x[n] * h[n]$.",
                   False, r"$y = x*h$ needs linearity and time-invariance: $y[n] = x^2[n]$ has $h[n] = \delta[n]$, so $x*h = x \ne x^2$."),
    "fully_described_lti": ("SP2023 #1(a)", r"If the response $y[n]$ of a discrete-time system to any possible input $x[n]$ is fully described by its unit pulse response, then the system must be LTI.",
                            True, r"Then $y = x*h$ for every input, and every convolution system is linear and time-invariant (HW2 #2)."),
    "abs_sum_any_system": ("SP2025 #1(b)", r"If the unit pulse response $h[n]$ of a system $S$ is absolutely summable, then $S$ must be BIBO stable regardless of whether $S$ is LTI or not.",
                           False, r"The $\sum\lvert h\rvert$ test holds for LTI systems only: $y[n] = n\,x[n]$ has $h = 0$, yet $x = u[n]$ gives $y = n\,u[n]$."),
    "causal_conv_causal": ("FA2025 #1(b)", r"The convolution of two causal signals (signals that are zero for all $n \lt 0$) is always a causal signal.",
                           True, r"$y[n] = \sum_k x[k]h[n-k]$ needs $k \ge 0$ and $n - k \ge 0$; for $n \lt 0$ no $k$ qualifies, so $y[n] = 0$."),
    "conv_sum_causal": ("SP2025 #1(a)", r"If $y[n] = \sum_{\ell=-\infty}^{\infty} h[\ell]\,x[n-\ell]$ for some well-defined $h[n]$, the system must be causal.",
                        False, r"The convolution sum makes the system LTI, not causal: $h[n] = \delta[n+1]$ gives $y[n] = x[n+1]$."),
    "right_sided_causal": ("SP2023 #1(c)", r"If a system has a right-sided unit pulse response, then it must be causal.",
                           False, r"Right-sided allows a start at a negative index: $h[n] = u[n+1]$ has $h[-1] = 1$."),
    "causal_implies_ti": ("FA2023 #1(a)", r"Any causal discrete-time system must also be time-invariant.",
                          False, r"The properties are independent: $y[n] = n\,x[n]$ is causal and time-varying."),
    "tv_not_causal": ("FA2019 #1(f)", r"A time-varying system cannot be causal.",
                      False, r"$y[n] = n\,x[n]$ is time-varying and causal (so is the HW1 #6 window)."),
    "bounded_h_stable": ("FA2025 #1(a)", r"The impulse response of a given LTI system is known to be a bounded sequence. This system must be BIBO stable.",
                         False, r"Bounded is not absolutely summable: $h[n] = u[n]$ is bounded, but the bounded input $u[n]$ gives $(n+1)u[n]$."),
    "unstable_h_unbounded": ("FA2024 #1(b)", r"The impulse response $h[n]$ of an unstable LTI system will always be an unbounded sequence.",
                             False, r"$h[n] = u[n]$ is bounded by $1$ and unstable ($\sum\lvert h\rvert = \infty$)."),
    "two_sided_never_stable": ("SP2025 #1(f)", r"An LTI system with a two-sided impulse response is never BIBO stable.",
                               False, r"$h[n] = \left(\tfrac12\right)^{\lvert n\rvert}$ is two-sided with $\sum\lvert h\rvert = 3$."),
    "unstable_any_input": ("FA2019 #1(j)", r"The response $y[n]$ of a BIBO unstable LTI system to any nonzero input $x[n]$ is always unbounded.",
                           False, r"Unstable means some bounded input blows up, not every input: $h = u[n]$ with $x = \delta[n] - \delta[n-1]$ gives $y = \delta[n]$."),
    "stable_unbounded_in_out": ("FA2025 #1(e)", r"An unbounded signal is passed as input to a BIBO stable LTI system. The resulting output is always an unbounded signal.",
                                False, r"BIBO says nothing about unbounded inputs: $h = \delta[n] - 2\delta[n-1]$ maps $2^n u[n]$ to $\delta[n]$."),
    "parallel_stable": ("SP2025 #1(d)", r"The parallel connection of two BIBO stable LTI systems always results in another BIBO stable LTI system.",
                        True, r"$\sum\lvert h_1 + h_2\rvert \le \sum\lvert h_1\rvert + \sum\lvert h_2\rvert \lt \infty$."),
    "parallel_stable_any": ("FA2023 #1(c)", r"Two BIBO stable systems connected in parallel always form a BIBO stable system.",
                            True, r"Even without LTI: if $\lvert y_1\rvert \le B_1$ and $\lvert y_2\rvert \le B_2$, then $\lvert y_1 + y_2\rvert \le B_1 + B_2$."),
    "parallel_unstable": ("FA2025 #1(f)", r"The parallel connection of two unstable LTI systems always forms another unstable LTI system.",
                          False, r"$h_1 = u[n]$ and $h_2 = \delta[n] - u[n]$ are unstable, but $h_1 + h_2 = \delta[n]$."),
    "cascade_unstable": ("SP2021 #1(c)", r"A cascade of two BIBO unstable LTI systems cannot be stable.",
                         False, r"Pole-zero cancellation: the causal $H_1 = \tfrac{1-2z^{-1}}{1-3z^{-1}}$ and $H_2 = \tfrac{1-3z^{-1}}{1-2z^{-1}}$ are unstable, but $H_1H_2 = 1$."),
    "series_order": ("FA2023 #1(e)", r"Two LTI systems with impulse responses $h_1[n]$ and $h_2[n]$ are connected in series. For the same input $x[n]$, the output is the same for either ordering of the two systems.",
                     True, r"Convolution is commutative and associative: $(x*h_1)*h_2 = x*(h_2*h_1) = (x*h_2)*h_1$."),
    "roc_no_zeros": ("FA2023 #1(b)", r"The ROC of a given z-transform cannot contain any poles or zeros.",
                     False, r"Poles never, zeros yes: $X(z) = 1 - z^{-1}$ has ROC $z \ne 0$, which contains its zero $z = 1$."),
    "poles_of_sum": ("SP2021 #1(b)", r"Let $X_1(z)$, $X_2(z)$ be the rational z-transforms of $x_1[n]$, $x_2[n]$. Then the poles of $X_1(z)$ and $X_2(z)$ must be poles of the z-transform of $x_1[n] + x_2[n]$.",
                     False, r"Poles can cancel: $x_1 = (\tfrac12)^n u[n]$ and $x_2 = \delta[n] - (\tfrac12)^n u[n]$ add up to $\delta[n]$, with $X = 1$."),
    "lccde_finite_poles": ("FA2024 #1(d)", r"An LCCDE system can only have a finite number of poles.",
                           True, r"$H(z) = B(z)/A(z)$ is a ratio of polynomials of finite degree."),
    "fir_stable_or_unstable": ("FA2019 #1(e)", r"An LTI system with a finite-length impulse response can be BIBO stable or unstable.",
                               False, r"For an FIR system $\sum\lvert h\rvert$ is a finite sum of finite numbers: always stable."),
}

GROUPS = {
    "h": ["h_any_system", "h_conv_any", "fully_described_lti", "abs_sum_any_system"],
    "causal": ["causal_conv_causal", "conv_sum_causal", "right_sided_causal", "causal_implies_ti", "tv_not_causal"],
    "stab": ["bounded_h_stable", "unstable_h_unbounded", "two_sided_never_stable", "h_an_stable", "unstable_any_input"],
    "unbounded": ["stable_unbounded_in_out", "unbounded_input_unstable", "causal_ramp"],
    "connect": ["parallel_stable", "parallel_stable_any", "parallel_unstable", "cascade_unstable", "series_order"],
    "z": ["roc_no_zeros", "poles_of_sum", "step_pole", "ztrans_exp_dne", "stable_pole_causality", "two_poles_two_sided",
          "no_roc_not_stable"],
    "lccde": ["lccde_finite_poles", "fir_stable_or_unstable", "lccde_poles_distinct", "lccde_causal_or_anticausal"],
    "delta": ["delta_pick", "delta_all"],
}


def pj(x):
    x = Fr(x)
    return [x.numerator, x.denominator]


# ----------------------------------------------------------------------------- parameterized statements
def stable_pole_causality():
    p = random.choice([Fr(3), Fr(2), Fr(3, 2), Fr(-2), Fr(1, 2), Fr(1, 3), Fr(-1, 2), Fr(3, 4)])
    word = random.choice(["non-causal", "anti-causal"]) if abs(p) > 1 else random.choice(["non-causal", "anti-causal"])
    ans = abs(p) > 1
    text = rf"A BIBO stable LTI system with transfer function $H(z) = \dfrac{{1}}{{{zfac(p)}}}$ must be {word}."
    if ans:
        why = (rf"The pole $z = {T(p)}$ is outside the unit circle, so stability forces the ROC $\lvert z\rvert \lt {T(abs(p))}$ and "
               rf"$h[n] = -\left({T(p)}\right)^n u[-n-1]$, which is zero for $n \ge 0$: {word}.")
    else:
        why = (rf"The pole $z = {T(p)}$ is inside the unit circle, so stability forces the ROC $\lvert z\rvert \gt {T(abs(p))}$ and "
               rf"$h[n] = \left({T(p)}\right)^n u[n]$, which is causal.")
    return "FA2023 #1(f), SP2025 #1(e)", text, ans, why, {"p": pj(p), "word": word}


def two_poles_two_sided():
    pool = [Fr(1, 4), Fr(1, 2), Fr(2, 3), Fr(3), Fr(2), Fr(4), Fr(-1, 2), Fr(-3), Fr(3, 4), Fr(5, 4)]
    while True:
        p1, p2 = sorted(random.sample(pool, 2), key=abs)
        if abs(p1) != abs(p2):
            break
    ans = abs(p1) < 1 < abs(p2)
    text = (rf"A BIBO stable LTI system has a transfer function with two poles, at $z = {T(p1)}$ and $z = {T(p2)}$. "
            r"The impulse response of this system must be two-sided.")
    if ans:
        why = rf"Stability needs the ROC to contain $\lvert z\rvert = 1$: that is the ring ${T(abs(p1))} \lt \lvert z\rvert \lt {T(abs(p2))}$, so $h$ is two-sided."
    elif abs(p2) < 1:
        why = rf"Both poles are inside the unit circle; the stable ROC is $\lvert z\rvert \gt {T(abs(p2))}$, so $h$ is right-sided (causal)."
    else:
        why = rf"Both poles are outside the unit circle; the stable ROC is $\lvert z\rvert \lt {T(abs(p1))}$, so $h$ is left-sided."
    return "FA2025 #1(c)", text, ans, why, {"p1": pj(p1), "p2": pj(p2)}


def no_roc_not_stable():
    pool = [Fr(1, 2), Fr(1, 3), Fr(-1, 2), Fr(2), Fr(3), Fr(-2), Fr(1), Fr(-1)]
    while True:
        a, b = random.sample(pool, 2)
        if abs(a) != abs(b):
            break
    ans = abs(a) == 1 or abs(b) == 1
    text = (rf"An LTI system with transfer function $H(z) = \dfrac{{1}}{{{zfac(a)}}} + \dfrac{{1}}{{{zfac(b)}}}$ "
            r"must <b>not</b> be BIBO stable.")
    lo, hi = sorted([abs(a), abs(b)])
    if ans:
        on = a if abs(a) == 1 else b
        why = rf"The pole $z = {T(on)}$ lies on the unit circle, and an ROC never contains a pole, so no ROC contains $\lvert z\rvert = 1$: never stable."
    elif lo < 1 < hi:
        why = rf"No ROC is given: the ring ${T(lo)} \lt \lvert z\rvert \lt {T(hi)}$ contains the unit circle and gives a stable (two-sided) system."
    elif hi < 1:
        why = rf"No ROC is given: $\lvert z\rvert \gt {T(hi)}$ contains the unit circle and gives a stable (causal) system."
    else:
        why = rf"No ROC is given: $\lvert z\rvert \lt {T(lo)}$ contains the unit circle and gives a stable (anti-causal) system."
    return "FA2024 #1(f), SP2021 #1(a)", text, ans, why, {"a": pj(a), "b": pj(b)}


def causal_ramp():
    p = random.choice([Fr(1), Fr(1), Fr(1, 2), Fr(-1, 2), Fr(2), Fr(1, 3), Fr(3, 2), Fr(-1)])
    ans = p == 1 or abs(p) > 1
    text = (rf"A causal LTI system with transfer function $H(z) = \dfrac{{z^{{-1}}}}{{{zfac(p)}}}$ produces an unbounded output "
            r"for the input $x[n] = u[n]$.")
    if p == 1:
        why = r"$Y(z) = \dfrac{z^{-1}}{(1-z^{-1})^2}$, so $y[n] = n\,u[n]$: the input pole lands on the system pole at $z = 1$ (double pole)."
    elif p == -1:
        why = (r"$Y(z) = \dfrac{z^{-1}}{(1-z^{-1})(1+z^{-1})}$ has two distinct simple poles on the unit circle: "
               r"$y[n] = \tfrac12 u[n] - \tfrac12(-1)^n u[n]$ stays bounded.")
    elif abs(p) > 1:
        why = rf"The causal system has a pole at $z = {T(p)}$ outside the unit circle; $y[n]$ contains a nonzero multiple of $\left({T(p)}\right)^n u[n]$."
    else:
        why = rf"The pole $z = {T(p)}$ is inside the unit circle, so the causal system is stable, and a bounded input gives a bounded output."
    return "SP2021 #1(d)", text, ans, why, {"p": pj(p)}


def delta_pick():
    m = random.choice([2, 3, 4])
    v = random.choice([3, 4, 5, -2, 6])
    kind = random.choice(["true", "value", "index"])
    mc_, vc = m, v
    if kind == "value":
        vc = random.choice([w for w in (2, 3, 4, 5, 8, -1) if w != v])
    elif kind == "index":
        mc_ = random.choice([w for w in (1, 2, 3, 4, 5) if w != m])
    ans = (mc_ == m and vc == v)
    text = (rf"Suppose $\sum_{{n=-\infty}}^{{\infty}} x[n]\,\delta\!\left[2^{{n}}u[n] - {2 ** m}\right] = {v}$. "
            rf"Then $x[{mc_}] = {vc}$.")
    why = (rf"$2^n u[n] - {2 ** m}$ equals $-{2 ** m}$ for $n \lt 0$ and vanishes only at $n = {m}$, so the sum is just $x[{m}] = {v}$"
           + ("." if ans else rf", which does not give $x[{mc_}] = {vc}$."))
    return "FA2019 #1(g)", text, ans, why, {"m": m, "v": v, "mc": mc_, "vc": vc}


def delta_all():
    S = random.choice([2, 4, 5, 7])
    Sc = random.choice([S, S - 1, S + 1])
    ans = Sc == S
    text = (rf"Suppose $\sum_{{n=-\infty}}^{{\infty}} x[n]\,\delta\!\left[4\cos\!\left(2n\pi + \tfrac{{\pi}}{{2}}\right) - 6\sin(n\pi)\right] = {S}$. "
            rf"Then $\sum_{{n=-\infty}}^{{\infty}} x[n] = {Sc}$.")
    why = (r"For every integer $n$, $\cos(2\pi n + \tfrac\pi2) = 0$ and $\sin(\pi n) = 0$, so the delta is $1$ for all $n$ and "
           rf"the given sum is $\sum_n x[n] = {S}$" + ("." if ans else rf", not ${Sc}$."))
    return "SP2021 #1(e)", text, ans, why, {"S": S, "Sc": Sc}


def step_pole():
    p = random.choice([Fr(1, 2), Fr(1, 3), Fr(-1, 2), Fr(2), Fr(1), Fr(1)])
    ans = p != 1
    text = (rf"Suppose the step response $g[n]$ of an LTI system (the output for $x[n] = u[n]$) has a z-transform with a pole at "
            rf"$z = {T(p)}$. Then $H(z)$ also has a pole at $z = {T(p)}$.")
    if ans:
        why = (rf"$G(z) = H(z)\cdot\dfrac{{1}}{{1-z^{{-1}}}}$, and the step's factor is finite and nonzero at $z = {T(p)}$, "
               "so it can neither create nor cancel that pole: it must come from $H$.")
    else:
        why = r"The step itself contributes the pole at $z = 1$: $H(z) = 1$ gives $G(z) = \dfrac{1}{1-z^{-1}}$ with a pole at $1$, but $H$ has none."
    return "FA2019 #1(h)", text, ans, why, {"p": pj(p)}


def ztrans_exp_dne():
    w = random.choice([(1, 4), (1, 3), (1, 2), (2, 3), (1, 6)])
    wt = rf"\frac{{{'' if w[0] == 1 else w[0]}\pi}}{{{w[1]}}}"
    text = rf"The z-transform of $x[n] = e^{{j{wt}n}}$ (for all $n$) does not exist because the ROC is empty."
    why = r"$x$ is two-sided with $\lvert x[n]\rvert = 1$: the $n \ge 0$ half needs $\lvert z\rvert \gt 1$, the $n \lt 0$ half needs $\lvert z\rvert \lt 1$."
    return "SP2025 #1(c)", text, True, why, {"w": list(w)}


def h_an_stable():
    a = random.choice([Fr(1, 2), Fr(-1, 2), Fr(3, 4), Fr(1, 3), Fr(2), Fr(-2), Fr(1), Fr(-1), Fr(5, 4)])
    ans = abs(a) < 1
    base = "1" if a == 1 else rf"\left({T(a)}\right)^{{n}}" if (a < 0 or a.denominator != 1) else rf"{T(a)}^{{n}}"
    htex = "u[n]" if a == 1 else rf"{base}u[n]"
    text = rf"The LTI system with impulse response $h[n] = {htex}$ is BIBO stable."
    if ans:
        why = rf"$\sum_n \lvert h[n]\rvert = \sum_{{n\ge0}} {fmt.tex_pow(abs(a), 'n')} = \dfrac{{1}}{{1 - {T(abs(a))}}} \lt \infty$."
    elif abs(a) == 1:
        why = r"$\lvert h[n]\rvert = 1$ for every $n \ge 0$: bounded, but $\sum\lvert h\rvert = \infty$ (same trap as FA2025 #1(a))."
    else:
        why = rf"$\lvert h[n]\rvert = {fmt.tex_pow(abs(a), 'n')}$ grows without bound, so $\sum\lvert h\rvert = \infty$."
    return "Lecture 4 / FA2025 #1(a)", text, ans, why, {"a": pj(a)}


def unbounded_input_unstable():
    a = random.choice([Fr(3), Fr(2), Fr(-2), Fr(3, 2), Fr(1, 2), Fr(-1, 2)])
    ans = abs(a) <= 1
    xt = rf"{T(a)}^{{n}}u[n]" if (a > 0 and a.denominator == 1) else rf"\left({T(a)}\right)^{{n}}u[n]"
    text = rf"If the response $y[n]$ of an LTI system to the input $x[n] = {xt}$ is unbounded, the system must be BIBO unstable."
    if ans:
        why = r"This input is bounded, and a bounded input with an unbounded output is exactly what BIBO instability means."
    else:
        why = r"This input is itself unbounded, so it proves nothing: the stable identity system $h = \delta[n]$ returns it unchanged."
    return "FA2019 #1(i)", text, ans, why, {"a": pj(a)}


def lccde_poles_distinct():
    N = random.choice([2, 3, 4, 5])
    text = rf"An LTI system given by $y[n] = y[n-{N}] + x[n]$ has ${N}$ distinct poles in its transfer function."
    kind = {2: "square", 3: "cube", 4: "fourth", 5: "fifth"}[N]
    why = (rf"$H(z) = \dfrac{{1}}{{1 - z^{{-{N}}}}} = \dfrac{{z^{{{N}}}}}{{z^{{{N}}} - 1}}$: the poles are the ${N}$ distinct "
           rf"{kind} roots of unity $e^{{j2\pi k/{N}}}$, $k = {', '.join(str(k) for k in range(N))}$.")
    return "FA2024 #1(e)", text, True, why, {"N": N}


def lccde_causal_or_anticausal():
    a = random.choice([Fr(1, 2), Fr(1, 3), Fr(2), Fr(-1, 2), Fr(3, 4)])
    text = rf"An LTI system specified by the difference equation $y[n] {'-' if a > 0 else '+'} {T(abs(a))}\,y[n-1] = x[n]$ can be causal or anti-causal."
    why = (rf"The equation fixes $H(z) = \dfrac{{1}}{{{zfac(a)}}}$ but not its ROC: $\left({T(a)}\right)^n u[n]$ and "
           rf"$-\left({T(a)}\right)^n u[-n-1]$ both satisfy it.")
    return "FA2019 #1(a)", text, True, why, {"a": pj(a)}


PARAM = {f.__name__: f for f in (stable_pole_causality, two_poles_two_sided, no_roc_not_stable, causal_ramp, delta_pick,
                                 delta_all, step_pole, ztrans_exp_dne, h_an_stable, unbounded_input_unstable,
                                 lccde_poles_distinct, lccde_causal_or_anticausal)}


def make(sid):
    if sid in FIXED:
        src, text, ans, why = FIXED[sid]
        return src, text, ans, why, {}
    return PARAM[sid]()


def generate(data):
    groups = random.sample(list(GROUPS), 6)
    items = []
    for i, g in enumerate(groups):
        sid = random.choice(GROUPS[g])
        src, text, ans, why, prm = make(sid)
        items.append({"name": f"tf{i + 1}", "label": "abcdef"[i], "sid": sid, "stmt_html": text, "src": src,
                      "param": prm, "ans": "True" if ans else "False", "why": why,
                      "choices": mc.true_false(ans)})
    data["params"]["statements"] = items
