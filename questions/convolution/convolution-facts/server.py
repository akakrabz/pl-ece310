"""Convolution shortcuts that need no full convolution (Lecture 4 §2.1: start point, end point,
length, eqs. 15-17; identity and shift properties; HW2 #5 sum checks).

x occupies a1..b1 and h occupies a2..b2 (short integer sequences with nonzero end samples, shown
with the n = 0 arrow). Asked, for y = x * h:
  (a) first index a1 + a2, last index b1 + b2, length Lx + Lh - 1        (integers)
  (b) sum_n y[n] = (sum_n x[n]) (sum_n h[n])                              (integer)
  (c) y at the first and last index = products of the end samples         (integers)
  (d) y[n0] for one interior n0 with 2-3 overlapping terms                (integer)
  (e) one conceptual multiple choice: x * delta[n-k], x * (delta[n+s] - delta[n+s-1]),
      x * u[n-k] (running sum), and x[n-a] * h[n-b] = y[n-a-b].
The full y is never shown (convolution/finite-convolution asks for it)."""

import random

from ece310 import convlib, fmt, mc

NZ = [-3, -2, -1, 1, 2, 3]


def _rand_seq(length):
    vals = [random.randint(-3, 3) for _ in range(length)]
    vals[0], vals[-1] = random.choice(NZ), random.choice(NZ)
    return vals


def _arg(k):
    """index n - k as LaTeX: n, n-3, n+2."""
    return convlib.tex_arg(k).replace("-", " - ").replace("+", " + ")


def _sig(name, k):
    return f"{name}[{_arg(k)}]"


def _conv(x, h):
    y = [0] * (len(x) + len(h) - 1)
    for i, a in enumerate(x):
        for j, b in enumerate(h):
            y[i + j] += a * b
    return y


def _mc_item():
    kind = random.choice(["delta", "diff", "step", "shift2"])
    if kind == "delta":
        k = random.choice([-2, -1, 1, 2, 3])
        stem = rf"Which expression equals $v[n] = x[n] * \delta[{_arg(k)}]$?"
        correct = f"${_sig('x', k)}$"
        wrong = [f"${_sig('x', -k)}$", rf"$x[{k}]\,\delta[{_arg(k)}]$", f"$x[{-k} - n]$" if k else "$x[-n]$",
                 rf"$x[n]\,\delta[{_arg(k)}]$"]
        prm = {"k": k}
        why = (rf"Convolving with a shifted impulse shifts: $x[n] * \delta[{_arg(k)}] = \sum_m x[m]\,\delta[n - {k} - m]"
               rf" = {_sig('x', k)}$ (only $m = n - {k}$ survives)." if k > 0 else
               rf"Convolving with a shifted impulse shifts: $x[n] * \delta[{_arg(k)}] = \sum_m x[m]\,\delta[n + {-k} - m]"
               rf" = {_sig('x', k)}$ (only $m = n + {-k}$ survives).")
    elif kind == "diff":
        s = random.choice([0, 1])
        prm = {"s": s}
        h_tex = rf"\delta[{_arg(-s)}] - \delta[{_arg(1 - s)}]"
        stem = rf"Which expression equals $v[n] = x[n] * \big({h_tex}\big)$?"
        correct = f"${_sig('x', -s)} - {_sig('x', 1 - s)}$"
        wrong = [f"${_sig('x', -s)} - {_sig('x', -s - 1)}$", f"${_sig('x', 1 - s)} - {_sig('x', -s)}$",
                 rf"$\sum_{{m=-\infty}}^{{{_arg(-s)}}} x[m]$", f"${_sig('x', -s)}\\,{_sig('x', 1 - s)}$"]
        why = (rf"Distribute over the two impulses: $x[n] * \delta[{_arg(-s)}] - x[n] * \delta[{_arg(1 - s)}]"
               rf" = {_sig('x', -s)} - {_sig('x', 1 - s)}$, a first difference (it undoes a running sum).")
    elif kind == "step":
        k = random.choice([0, 1, 2, -1])
        top = _arg(k)
        prm = {"k": k}
        stem = rf"Which expression equals $v[n] = x[n] * u[{_arg(k)}]$?"
        correct = rf"$\sum_{{m=-\infty}}^{{{top}}} x[m]$"
        wrong = [rf"$\sum_{{m={top}}}^{{\infty}} x[m]$", rf"$\sum_{{m=-\infty}}^{{{_arg(-k)}}} x[m]$",
                 rf"${_sig('x', k)}\,u[{_arg(k)}]$", r"$\sum_{m=-\infty}^{\infty} x[m]$"]
        why = (rf"$x[n] * u[{_arg(k)}] = \sum_m x[m]\,u[{_arg(k)} - m]$, and $u[{_arg(k)} - m] = 1$ exactly when "
               rf"$m \le {top}$: a running sum of $x$ up to ${top}$.")
    else:
        while True:
            a, b = random.choice([-2, -1, 1, 2, 3]), random.choice([-2, -1, 1, 2, 3])
            if a + b != 0 and a - b != 0:
                break
        prm = {"a": a, "b": b}
        stem = (rf"Let $y[n] = x[n] * h[n]$. Which expression equals $w[n] = {_sig('x', a)} * {_sig('h', b)}$?")
        correct = f"${_sig('y', a + b)}$"
        wrong = [f"${_sig('y', a - b)}$", f"${_sig('y', a)}$", f"${_sig('y', -a - b)}$", f"${_sig('y', a * b)}$"]
        why = (rf"Shifting either factor shifts the convolution: $x[n - a] * h[n - b] = y[n - a - b]$ "
               rf"(commutativity and associativity with $\delta[n-a] * \delta[n-b] = \delta[n-a-b]$). "
               rf"Here $a = {a}$ and $b = {b}$, so $w[n] = {_sig('y', a + b)}$.")
    return kind, stem, correct, wrong, why, prm


def generate(data):
    while True:
        lx, lh = random.choice([(3, 2), (3, 3), (4, 2), (4, 3), (5, 2), (3, 4)])
        x, h = _rand_seq(lx), _rand_seq(lh)
        a1, a2 = random.randint(-3, 2), random.randint(-2, 2)
        if sum(x) == 0 or sum(h) == 0 or (a1, a2) == (0, 0):
            continue
        y = _conv(x, h)
        b1, b2 = a1 + lx - 1, a2 + lh - 1
        nf, nl = a1 + a2, b1 + b2
        # interior index with 2 or 3 overlapping sample pairs
        cand = []
        for n0 in range(nf + 1, nl):
            ks = [k for k in range(a1, b1 + 1) if a2 <= n0 - k <= b2]
            if 2 <= len(ks) <= 3:
                cand.append((n0, ks))
        if cand:
            break
    n0, ks = random.choice(cand)
    yn0 = y[n0 - nf]

    p = data["params"]
    p["x"], p["a1"], p["h"], p["a2"] = x, a1, h, a2
    p["x_tex"], p["h_tex"] = fmt.tex_seq(x, a1), fmt.tex_seq(h, a2)
    p["b1"], p["b2"], p["lx"], p["lh"] = b1, b2, lx, lh
    p["n0"] = n0
    p["first_tex"] = f"{a1} {'-' if a2 < 0 else '+'} {abs(a2)}"
    p["last_tex"] = f"{b1} {'-' if b2 < 0 else '+'} {abs(b2)}"
    p["len_tex"] = f"{lx} + {lh} - 1"
    p["sx_tex"] = " + ".join(fmt.tex_num(v, paren=True) for v in x) + f" = {sum(x)}"
    p["sh_tex"] = " + ".join(fmt.tex_num(v, paren=True) for v in h) + f" = {sum(h)}"
    p["sy_tex"] = f"{fmt.tex_num(sum(x), paren=True)}\\cdot{fmt.tex_num(sum(h), paren=True)} = {sum(x) * sum(h)}"
    p["yfirst_tex"] = (f"x[{a1}]\\,h[{a2}] = {fmt.tex_num(x[0], paren=True)}\\cdot{fmt.tex_num(h[0], paren=True)}"
                       f" = {x[0] * h[0]}")
    p["ylast_tex"] = (f"x[{b1}]\\,h[{b2}] = {fmt.tex_num(x[-1], paren=True)}\\cdot{fmt.tex_num(h[-1], paren=True)}"
                      f" = {x[-1] * h[-1]}")
    p["kmin"], p["kmax"] = ks[0], ks[-1]
    p["terms_tex"] = " + ".join(f"x[{k}]\\,h[{n0 - k}]" for k in ks)
    p["vals_tex"] = " + ".join(f"{fmt.tex_num(x[k - a1], paren=True)}\\cdot{fmt.tex_num(h[n0 - k - a2], paren=True)}"
                               for k in ks)
    p["yn0_val"] = yn0
    p["y_tex"] = fmt.tex_seq(y, nf)

    kind, stem, correct, wrong, why, prm = _mc_item()
    p["mc_kind"], p["mc_param"], p["mc_stem_html"], p["mc_why_html"] = kind, prm, stem, why
    p["mc_correct_html"] = correct
    p["prop_choices"] = mc.choices(correct, wrong, n=4)

    c = data["correct_answers"]
    c["nfirst"], c["nlast"], c["leny"] = nf, nl, lx + lh - 1
    c["sumy"] = sum(x) * sum(h)
    c["yfirst"], c["ylast"] = x[0] * h[0], x[-1] * h[-1]
    c["yn0"] = yn0
