"""Index transformation y[n] = x[a n + b] of a finite sequence (HW1 #2: x[-n + 3], x[2n + 1]).

a is 1 (shift), -1 (reversal + shift), 2 or -2 (decimation + shift); b is a small nonzero integer.
The student enters the first index where y is nonzero and the samples of y from its first to its
last nonzero sample. For |a| = 2 the variant is redrawn until decimation really drops at least one
nonzero sample of x and at least two nonzero samples survive."""

import random

import numpy as np
import prairielearn as pl
from ece310 import fmt

NONZERO = [-3, -2, -1, 1, 2, 3, 4, 5]


def _arg_tex(a, b):
    """a n + b written the way HW1 writes it: n + 2, n - 3, -n + 3, 2n + 1, -2n - 1."""
    an = {1: "n", -1: "-n", 2: "2n", -2: "-2n"}[a]
    if b == 0:
        return an
    return f"{an} {'+' if b > 0 else '-'} {abs(b)}"


def _xval(x, nx, m):
    i = m - nx
    return x[i] if 0 <= i < len(x) else 0


def _transform(x, nx, a, b):
    """(values, start) of y[n] = x[a n + b], trimmed to first..last nonzero sample."""
    ns = list(range(-40, 41))
    vals = [_xval(x, nx, a * n + b) for n in ns]
    nz = [i for i, v in enumerate(vals) if v != 0]
    if not nz:
        return None, None
    return vals[nz[0]:nz[-1] + 1], ns[nz[0]]


def _draw():
    L = random.randint(4, 6)
    x = [random.randint(-3, 5) for _ in range(L)]
    x[0], x[-1] = random.choice(NONZERO), random.choice(NONZERO)
    nx = random.randint(-3, 1)
    a = random.choice([1, -1, -1, 2, 2, -2, -2])
    b = random.choice([-3, -2, -1, 1, 2, 3])
    return x, nx, a, b


def _table_html(ns, a, b, x, nx):
    head = "".join(f"<th>${n}$</th>" for n in ns)
    ms = "".join(f"<td>${a * n + b}$</td>" for n in ns)
    vals = "".join(f"<td><b>${_xval(x, nx, a * n + b)}$</b></td>" for n in ns)
    return ('<table class="table table-sm table-bordered text-center" style="width:auto">'
            f"<tr><th>$n$</th>{head}</tr>"
            f"<tr><th>$m = {_arg_tex(a, b)}$</th>{ms}</tr>"
            f"<tr><th>$y[n] = x[m]$</th>{vals}</tr></table>")


def generate(data):
    while True:
        x, nx, a, b = _draw()
        y, ny = _transform(x, nx, a, b)
        if y is None:
            continue
        if abs(a) == 2:
            kept = [m for m in range(nx, nx + len(x)) if (m - b) % 2 == 0 and _xval(x, nx, m) != 0]
            dropped = [m for m in range(nx, nx + len(x)) if (m - b) % 2 != 0 and _xval(x, nx, m) != 0]
            if len(kept) < 2 or len(dropped) < 1:
                continue
        if y != y[::-1] or abs(a) == 2:      # a palindrome makes "reversed vector" undetectable
            break

    nx_last = nx + len(x) - 1
    ny_last = ny + len(y) - 1
    p = data["params"]
    p["x"], p["nx"], p["a"], p["b"] = x, nx, a, b
    p["x_tex"] = fmt.tex_seq(x, nx)
    p["arg_tex"] = _arg_tex(a, b)
    p["nx_last"] = nx_last
    p["ny_last"] = ny_last
    p["y_tex"] = fmt.tex_seq(y, ny)

    # support of the argument: nx <= a n + b <= nx_last  <=>  n in [lo, hi]
    lo_r, hi_r = sorted([(nx - b) / a, (nx_last - b) / a])
    lo, hi = int(np.ceil(lo_r - 1e-12)), int(np.floor(hi_r + 1e-12))
    p["support_tex"] = (fr"{nx} \le {_arg_tex(a, b)} \le {nx_last} \iff {lo} \le n \le {hi}")
    p["trimmed"] = (lo != ny or hi != ny_last)
    p["table_html"] = _table_html(list(range(lo, hi + 1)), a, b, x, nx)

    # two-step reading: shift first, then scale/reverse the CURRENT n
    w_tex = "x[n]" if b == 0 else f"x[n {'+' if b > 0 else '-'} {abs(b)}]"
    p["w_tex"] = w_tex
    p["shift_words"] = (f"advance $x$ by ${b}$ (move it {b} to the left)" if b > 0
                        else f"delay $x$ by ${-b}$ (move it {-b} to the right)")
    p["scale_tex"] = {1: "w[n]", -1: "w[-n]", 2: "w[2n]", -2: "w[-2n]"}[a]
    p["scale_words"] = {1: "nothing else to do (a pure shift)",
                        -1: "reverse $w$ about $n = 0$",
                        2: "keep only the even-index samples of $w$ ($w[2n]$) and close the gaps",
                        -2: "keep only the even-index samples of $w$, close the gaps, and reverse"}[a]
    p["decimate"] = abs(a) == 2
    if abs(a) == 2:
        dropped = [m for m in range(nx, nx_last + 1) if (m - b) % 2 != 0 and _xval(x, nx, m) != 0]
        p["dropped_tex"] = ",\\ ".join(f"x[{m}] = {_xval(x, nx, m)}" for m in dropped)
        p["parity"] = "odd" if b % 2 else "even"
    a_str = {1: "", -1: "-", 2: "2", -2: "-2"}[a]
    p["wrong_order_tex"] = (f"x[{a_str}(n {'+' if b > 0 else '-'} {abs(b)})] = x[{_arg_tex(a, a * b)}]"
                            if a != 1 else "")
    p["show_wrong_order"] = a != 1

    data["correct_answers"]["ny"] = ny
    data["correct_answers"]["y"] = pl.to_json(np.array([[float(v) for v in y]]))
