"""Lecture 1 beyond index maps: sequences written with delta[n] and u[n], sampling a sinusoid, and
the continuous-/discrete-time, continuous-/discrete-valued classification. Every variant asks:

  (a) a sequence written with shifted/scaled steps and impulses (HW1 #1, #3 style): the first
      nonzero index and the samples from the first to the last nonzero one. Families:
        steps   c1 u[n-k1] + c2 u[n-k2] + c3 u[n-k3] with c1 + c2 + c3 = 0 (HW1 #3(b) form)
        window  C (-1)^n (u[n-a] - u[n-b]) + c delta[n-m]   (or a constant C instead of (-1)^n)
        ramp    C (n-s) (u[n-a] - u[n-b])  or  C (n-s) u[n-a] u[l-n]  (HW1 #1(c), FA2025 #3(b));
                the factor n - s may vanish at an end sample, which then is not "nonzero";
  (b) x_a(t) = A cos(2 pi f0 t + theta) (or A cos(Omega0 t + theta)) sampled at f_s (or with period
      T), f0 < f_s/2: the digital frequency w0 = 2 pi f0 / f_s = Omega0 T as a multiple of pi and
      one sample x[k], chosen so that x[k] is a nonzero rational number;
  (c) one described signal to classify with Lecture 1's two yes/no questions (time axis, values).
"""

import math
import random
from fractions import Fraction as Fr

import numpy as np
import prairielearn as pl
from ece310 import fmt
from ece310.w1_signals import cos_pi, plain_frac, principal, tex_phase_term, tex_pi

NZ = [-3, -2, -1, 1, 2, 3]


# ----------------------------------------------------------------------------- part (a) helpers
def _arg(k):
    """n - k as an index: n, n - 3, n + 2."""
    return "n" if k == 0 else (f"n - {k}" if k > 0 else f"n + {-k}")


def _arg_rev(m):
    """m - n as an index: -n, 3 - n, -n - 2."""
    return "-n" if m == 0 else (f"{m} - n" if m > 0 else f"-n - {-m}")


def _lin(s):
    """n - s as a factor: n, (n - 2), (n + 1)."""
    return "n" if s == 0 else f"\\left({_arg(s)}\\right)"


def _trim(vals):
    """vals: dict n -> value. (first index, list of samples first..last nonzero) or (None, [])."""
    nz = sorted(n for n, v in vals.items() if v != 0)
    if not nz:
        return None, []
    return nz[0], [vals.get(n, 0) for n in range(nz[0], nz[-1] + 1)]


def _fam_steps():
    while True:
        k1 = random.randint(-3, 1)
        L1, L2 = random.randint(1, 3), random.randint(1, 3)
        if L1 + L2 < 3:
            continue
        k2, k3 = k1 + L1, k1 + L1 + L2
        c1, lvl2 = random.choice(NZ), random.choice(NZ)
        if lvl2 == c1:
            continue
        terms = [(c1, k1), (lvl2 - c1, k2), (-lvl2, k3)]
        x_tex = fmt.tex_sum((c, f"u[{_arg(k)}]") for c, k in terms)
        vals = {n: sum(c * (n >= k) for c, k in terms) for n in range(-12, 13)}
        ext = dict(vals)                       # off by one: the last level kept one sample too long
        ext[k3] = lvl2
        expl = (rf"Read it from left to right: nothing is on before $n = {k1}$; $u[{_arg(k1)}]$ switches on the level "
                rf"${c1}$ at $n = {k1}$; $u[{_arg(k2)}]$ adds ${lvl2 - c1}$ at $n = {k2}$ (level ${lvl2}$); "
                rf"$u[{_arg(k3)}]$ adds ${-lvl2}$ at $n = {k3}$, and from there on the levels cancel: "
                rf"${c1} {'+' if lvl2 - c1 > 0 else '-'} {abs(lvl2 - c1)} {'+' if -lvl2 > 0 else '-'} {abs(lvl2)} = 0$. "
                rf"Each step switches on where its argument is zero, so the level ${lvl2}$ lasts for "
                rf"${k2} \le n \le {k3 - 1}$ (not up to $n = {k3}$).")
        return {"fam": "steps", "x_tex": x_tex, "vals": vals, "ext": ext, "expl": expl,
                "lo": k1 - 1, "hi": k3 + 1}


def _fam_window():
    while True:
        a = random.randint(-3, 1)
        L = random.randint(3, 5)
        b = a + L
        alt = random.random() < 0.65
        C = random.choice([1, 1, -1, 2]) if alt else random.choice([2, 3, -1, -2])
        m = random.randint(a - 2, b + 1)
        c = random.choice(NZ)

        def wv(n, bb=b):
            return (C * (-1) ** (n % 2) if alt else C) if a <= n < bb else 0

        vals = {n: wv(n) + (c if n == m else 0) for n in range(-12, 13)}
        first, ys = _trim(vals)
        if first is None or len(ys) < 3 or len(ys) > 8 or sum(v != 0 for v in ys) < 3:
            continue
        ext = {n: wv(n, b + 1) + (c if n == m else 0) for n in range(-12, 13)}
        coef = {1: "", -1: "-"}.get(C, str(C))
        win = rf"\left(u[{_arg(a)}] - u[{_arg(b)}]\right)"
        w_tex = (coef + r"(-1)^{n}" + win) if alt else (coef + win)
        d_tex = fmt.tex_sum([(c, rf"\delta[{_arg(m)}]")])
        if random.random() < 0.7:
            x_tex = w_tex + (" - " if c < 0 else " + ") + d_tex.lstrip("-")
        else:
            x_tex = d_tex + (" - " + w_tex[1:] if w_tex.startswith("-") else " + " + w_tex)
        inside = a <= m < b
        where = (rf"inside the window, so it adds ${c}$ to the window sample there: $x[{m}] = {wv(m)} {'+' if c > 0 else '-'} {abs(c)} = {vals[m]}$"
                 if inside else rf"outside the window, so $x[{m}] = {c}$")
        expl = (rf"$u[{_arg(a)}] - u[{_arg(b)}]$ is $1$ for ${a} \le n \le {b - 1}$ ($u[{_arg(b)}]$ is already on at "
                rf"$n = {b}$, so that sample is cancelled): {b - a} samples. "
                + (rf"On them the factor ${coef}(-1)^n$ gives alternating signs, "
                   rf"starting with ${wv(a)}$ at $n = {a}$. " if alt else rf"On them the value is ${C}$. ")
                + rf"The impulse $\delta[{_arg(m)}]$ sits where its argument is zero, at $n = {m}$: {where}.")
        return {"fam": "window", "x_tex": x_tex, "vals": vals, "ext": ext, "expl": expl,
                "lo": min(a, m) - 1, "hi": max(b, m + 1)}


def _fam_ramp():
    while True:
        a = random.randint(-3, 1)
        L = random.randint(3, 5)
        last = a + L - 1
        s = random.randint(a - 1, last)
        C = random.choice([1, 1, -1, 2])
        prod = random.random() < 0.45
        lin = _lin(s)
        coef = {1: "", -1: "-"}.get(C, str(C))
        if prod:
            x_tex = coef + lin + rf"\,u[{_arg(a)}]\,u[{_arg_rev(last)}]"
        else:
            x_tex = coef + lin + rf"\left(u[{_arg(a)}] - u[{_arg(last + 1)}]\right)"

        def v(n, lo=a, hi=last):
            return C * (n - s) if lo <= n <= hi else 0

        vals = {n: v(n) for n in range(-12, 13)}
        first, ys = _trim(vals)
        if first is None or len(ys) < 2 or sum(x != 0 for x in ys) < 2:
            continue
        # off by one at the right edge: the window read one sample too long (difference form) or
        # u[l - n] read as "n < l" (product form)
        ext = {n: (v(n, a, last + 1) if not prod else v(n, a, last - 1)) for n in range(-12, 13)}
        win = (rf"u[{_arg(a)}]\,u[{_arg_rev(last)}]$ is $1$ where both steps are on: $n \ge {a}$ and "
               rf"${_arg_rev(last)} \ge 0$, i.e. $n \le {last}" if prod else
               rf"u[{_arg(a)}] - u[{_arg(last + 1)}]$ is $1$ for ${a} \le n \le {last}")
        zero = ""
        if a <= s <= last:
            zero = (rf" The factor ${_arg(s)}$ vanishes at $n = {s}$, so $x[{s}] = 0$"
                    + (" and that sample is not part of the nonzero range." if s in (a, last) else
                       " (an interior zero, which stays in the list)."))
        expl = (rf"${win}$. On those samples $x[n] = {fmt.tex_sum([(C, 'n'), (-C * s, '')])}$." + zero)
        return {"fam": "ramp", "x_tex": x_tex, "vals": vals, "ext": ext, "expl": expl,
                "lo": a - 1, "hi": last + 1, "lead_zero": s == a, "a": a}


def _table_html(vals, lo, hi):
    ns = list(range(lo, hi + 1))
    head = "".join(f"<th>${n}$</th>" for n in ns)
    row = "".join(f"<td><b>${vals.get(n, 0)}$</b></td>" for n in ns)
    return ('<table class="table table-sm table-bordered text-center" style="width:auto">'
            f"<tr><th>$n$</th>{head}</tr><tr><th>$x[n]$</th>{row}</tr></table>")


def _gen_seq(p):
    fam = random.choice(["steps", "window", "window", "ramp", "ramp"])
    d = {"steps": _fam_steps, "window": _fam_window, "ramp": _fam_ramp}[fam]()
    first, ys = _trim(d["vals"])
    ef, eys = _trim(d["ext"])
    p["seq_family"] = d["fam"]
    p["seq_tex"] = d["x_tex"]
    p["seq_expl_html"] = d["expl"]
    p["seq_table_html"] = _table_html(d["vals"], min(d["lo"], first - 1), max(d["hi"], first + len(ys)))
    last = first + len(ys) - 1
    p["seq_samples_tex"] = (fmt.tex_seq(ys, first, pad_to_zero=False)
                            + rf",\qquad {first} \le n \le {last}\ \text{{(zero elsewhere)}}")
    p["seq_last"] = first + len(ys) - 1
    p["seq_ext"] = [ef, eys]                      # the off-by-one reading (for the checker's probes)
    p["seq_lead_zero"] = bool(d.get("lead_zero", False))
    p["seq_a"] = d.get("a", first)
    return first, ys


# ----------------------------------------------------------------------------- part (b)
W0 = [Fr(1, 2), Fr(1, 3), Fr(2, 3), Fr(1, 4), Fr(3, 4), Fr(1, 6), Fr(5, 6)]
THETA = [Fr(0), Fr(1, 6), Fr(-1, 6), Fr(1, 4), Fr(-1, 4), Fr(1, 3), Fr(-1, 3), Fr(1, 2), Fr(-1, 2),
         Fr(2, 3), Fr(-2, 3), Fr(3, 4), Fr(-3, 4), Fr(5, 6), Fr(-5, 6)]
FS = [1000, 2000, 4000, 5000, 6000, 8000, 9000, 10000, 12000, 16000, 20000, 24000, 40000, 48000]
RATIONAL = {Fr(0), Fr(1, 3), Fr(2, 3), Fr(1), Fr(4, 3), Fr(5, 3)}     # cos(t pi) in {+-1, +-1/2}


def _hz_tex(f):
    """8000 -> 8\\text{ kHz}; 1500 -> 1500\\text{ Hz}; 12500 -> 12.5\\text{ kHz}."""
    if f % 1000 == 0:
        return f"{f // 1000}" + r"\text{ kHz}"
    if f % 100 == 0 and f > 1000:
        return f"{f / 1000:g}" + r"\text{ kHz}"
    return f"{f}" + r"\text{ Hz}"


def _period_tex(fs):
    """T = 1/fs written in microseconds / milliseconds when exact, else as a fraction of a second."""
    us = Fr(10 ** 6, fs)
    if us.denominator == 1:
        us = us.numerator
        if us % 1000 == 0:
            return f"{us // 1000}" + r"\text{ ms}"
        return f"{us}" + r"\ \mu\text{s}"
    return r"\tfrac{1}{" + str(fs) + r"}\text{ s}"


def _xk_mistakes(A, w, th, k, val):
    """Distinct wrong values of x[k]: phase sign flipped, x_a(k) = A cos(theta) (t = k seconds),
    w0 without the factor 2, phase dropped."""
    out = []
    for m in (A * math.cos(math.pi * (w * k - th)), A * math.cos(math.pi * th),
              A * math.cos(math.pi * (w * k / 2 + th)), A * math.cos(math.pi * w * k)):
        if abs(m - float(val)) > 0.02 * max(1, abs(val)) and all(abs(m - o) > 1e-6 for o in out):
            out.append(m)
    return out


def _gen_sampling(p):
    while True:
        w, th = random.choice(W0), random.choice(THETA)
        if th == 0 and random.random() < 0.6:
            continue
        A = random.choice([2, 3, 4, 5, 6, 8, 10])
        k = random.randint(1, 5)
        ang = (w * k + th) % 2
        if ang not in RATIONAL:
            continue
        fs = random.choice(FS)
        f0 = w * fs / 2
        if f0.denominator != 1 or f0 % 50 != 0:
            continue
        f0 = int(f0)
        val = A * cos_pi(ang).q
        if abs(A * math.cos(math.pi * th) - val) < 0.02 * abs(val) or len(_xk_mistakes(A, w, th, k, val)) < 2:
            continue                  # x_a(k) (t = k s) or the usual slips would give the right number
        break
    mode = random.choice(["fs", "T", "Omega"])
    p["mode_fs"], p["mode_T"], p["mode_Omega"] = mode == "fs", mode == "T", mode == "Omega"
    ph = tex_phase_term(th)
    if mode == "Omega":
        p["xa_tex"] = rf"{A}\cos\!\left({2 * f0}\pi t{ph}\right)"
        p["rate_tex"] = rf"f_s = {_hz_tex(fs)}"
        p["w0_steps_tex"] = (rf"\omega_0 = \Omega_0 T = \dfrac{{\Omega_0}}{{f_s}} = \dfrac{{{2 * f0}\pi}}{{{fs}}}"
                             rf" = {tex_pi(w)}")
    else:
        p["xa_tex"] = rf"{A}\cos\!\left(2\pi({f0})t{ph}\right)"
        if mode == "fs":
            p["rate_tex"] = rf"f_s = {_hz_tex(fs)}"
            p["w0_steps_tex"] = (rf"\omega_0 = 2\pi f_0 T = \dfrac{{2\pi f_0}}{{f_s}} = \dfrac{{2\pi\cdot {f0}}}{{{fs}}}"
                                 rf" = {tex_pi(w)}")
        else:
            p["rate_tex"] = rf"T = {_period_tex(fs)}"
            p["w0_steps_tex"] = (rf"\omega_0 = 2\pi f_0 T = 2\pi\cdot {f0}\cdot \tfrac{{1}}{{{fs}}}"
                                 rf" = {tex_pi(w)}")
    p["fs"], p["f0"], p["A"], p["k"] = fs, f0, A, k
    half_tex = _hz_tex(fs // 2)
    p["nyq_tex"] = {"fs": rf"f_0 = {_hz_tex(f0)} \lt f_s/2 = {half_tex}",
                    "T": rf"f_0 = {_hz_tex(f0)} \lt \tfrac{{1}}{{2T}} = {half_tex}",
                    "Omega": rf"f_0 = \tfrac{{\Omega_0}}{{2\pi}} = {_hz_tex(f0)} \lt f_s/2 = {half_tex}"}[mode]
    p["w0_pq"], p["theta"] = [w.numerator, w.denominator], [th.numerator, th.denominator]
    p["fs_words"] = f"{fs:,}".replace(",", "{,}")
    p["T_tex"] = _period_tex(fs)
    p["w0_tex"] = tex_pi(w)
    p["w0_frac"] = plain_frac(w)
    p["xn_tex"] = rf"{A}\cos\!\left({tex_pi(w)}n{ph}\right)"
    p["xk_steps_tex"] = (rf"x[{k}] = {A}\cos\!\left({tex_pi(w)}\cdot {k}{ph}\right) = {A}\cos\!\left({tex_pi(w * k + th)}\right)"
                         + (rf" = {A}\cos\!\left({tex_pi(principal(ang))}\right)" if principal(ang) != w * k + th else "")
                         + rf" = {A}\cdot\left({fmt.tex_num(cos_pi(ang).q)}\right) = {fmt.tex_num(val)}")
    p["xk_frac"] = plain_frac(val)
    p["ang"] = [ang.numerator, ang.denominator]
    return float(w), float(val)


# ----------------------------------------------------------------------------- part (c)
CLASSES = {
    "AN": "continuous-time and continuous-valued: an analog signal",
    "CQ": "continuous-time and discrete-valued (quantized)",
    "DT": "discrete-time and continuous-valued",
    "DG": "discrete-time and discrete-valued (quantized): a digital signal",
}
ORDER = ["AN", "CQ", "DT", "DG"]
WHY_T = {True: "it exists at every instant $t$ (continuous time)",
         False: "it exists only at integer indices $n$ (discrete time)"}

# (description, class, reason for the time axis, reason for the values)
BANK = [
    ("The voltage $v(t)$ produced by a microphone while someone speaks, at every instant $t$.", "AN",
     "the voltage can take any real value"),
    ("The air pressure at your eardrum while a violin plays, at every instant $t$.", "AN",
     "pressure varies continuously and can take any real value"),
    ("The height of the tide at a pier, as it varies over a day, at every instant $t$.", "AN",
     "a height can take any real value"),
    ("The temperature of a cooling cup of coffee at every instant $t$, as an exact real number.", "AN",
     "the temperature is any real number"),
    ("The voltage on a digital logic line, which at every instant $t$ is either $0$ V or $5$ V.", "CQ",
     "only two values are possible"),
    ("The number of people inside a store at each instant $t$ of the day.", "CQ",
     "a head count is an integer"),
    ("The floor number of an elevator car, recorded at every instant $t$ (the floor it is at or last passed).", "CQ",
     "floor numbers are integers"),
    ("The number of cars waiting at a red light, at every instant $t$.", "CQ",
     "a count of cars is an integer"),
    ("The samples $x[n] = x_a(nT)$ of a microphone voltage taken every $T = 1/44100$ s, each kept as an exact real number.", "DT",
     "each sample is an exact real number (no rounding)"),
    ("The samples of a seismometer voltage taken every $T = 10$ ms, each kept as an exact real number (no rounding).", "DT",
     "each sample can be any real number"),
    ("The exact air temperature at noon on day $n$ of the year (any real value).", "DT",
     "each value is an exact real number"),
    ("The exact average speed of a runner during lap $n$ of a race (any positive real value).", "DT",
     "an exact average can be any real number"),
    ("Hourly temperature readings, each stored to the nearest $0.1\\,^{\\circ}$C.", "DG",
     "rounding to the nearest $0.1\\,^{\\circ}$C quantizes the values"),
    ("A song on a CD: $44{,}100$ samples per second, each stored as a 16-bit integer.", "DG",
     "a 16-bit integer has only $2^{16}$ possible values"),
    ("The daily closing price of a stock, recorded in whole cents.", "DG",
     "prices in whole cents are quantized"),
    ("The number of cars passing a sensor in each one-minute interval $n$.", "DG",
     "a count is an integer"),
]


def _gen_classify(p):
    desc, cls, why_v = random.choice(BANK)
    p["sig_desc_html"] = desc
    p["sig_class"] = cls
    cont_t = cls in ("AN", "CQ")
    quant = cls in ("CQ", "DG")
    p["sig_why_html"] = (f"Time axis: {WHY_T[cont_t]}. Values: {why_v}, so they are "
                         + ("<b>discrete-valued (quantized)</b>" if quant else "<b>continuous-valued</b>")
                         + f". Lecture 1: {CLASSES[cls]}.")
    p["sig_choices"] = [{"text": CLASSES[c], "correct": c == cls} for c in ORDER]


# ----------------------------------------------------------------------------- generate
def generate(data):
    p = data["params"]
    first, ys = _gen_seq(p)
    w, val = _gen_sampling(p)
    _gen_classify(p)
    c = data["correct_answers"]
    c["nfirst"] = first
    c["xs"] = pl.to_json(np.array([[float(v) for v in ys]]))
    c["w0"] = w
    c["xk"] = val
