"""Series/parallel networks of three FIR blocks and a delay (Lecture 10 section 2 "system algebra":
series multiplies, parallel adds, Table 1 and the slide-11 five-block example; Lecture 4 convolution
properties; FA2024 Midterm 1 #3(b) and SP2025 #3(b): the overall h[n] of a series / parallel pair).

Shapes (drawn as an inline SVG block diagram; no feedback, as in Lecture 10):
    T1  (h1 * h2) +- h3        T2  h1 * (h2 +- h3)        T3  (h1 +- h2) * h3
with a delay block z^{-k} (k = 1 or 2) in series inside one of the two parallel branches, and a minus
sign on the second adder input in about a quarter of the variants. Blocks are 2-3 integer samples
starting at n = 0; in a fifth of the variants one block starts at n = -1 (non-causal, as FA2024 #3).

Asked: (a) the overall h[n] (first index + row vector), (b) the overall H(z) (pl-symbolic-input in z,
any equivalent form: SymPy proves polynomial identities in z^-1), (c) one output value: either
y[n0] for x[n] = delta[n - m], or the constant output for x[n] = c (all n), which is c H(1).
Every variant is checked to make the classic mistakes visible: forgetting the delay, adding where
the figure convolves (series read as parallel), ignoring the minus sign."""

import random

import numpy as np
import prairielearn as pl
from ece310 import fmt

ENDS = [-2, -1, 1, 2]
SUB = {1: "₁", 2: "₂", 3: "₃"}
SUP = {1: "⁻¹", 2: "⁻²"}
BW, BH, GAP, VGAP, R = 70, 36, 30, 18, 11


# ----------------------------------------------------------------------------- Laurent polynomials {n: coef}
def lmul(a, b):
    out = {}
    for i, u in a.items():
        for j, v in b.items():
            out[i + j] = out.get(i + j, 0) + u * v
    return {k: v for k, v in out.items() if v != 0}


def ladd(a, b, sb=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + sb * v
    return {k: v for k, v in out.items() if v != 0}


def lseq(d):
    """{n: c} -> (values from first to last nonzero, first index)."""
    lo, hi = min(d), max(d)
    return [d.get(n, 0) for n in range(lo, hi + 1)], lo


def ltex(d):
    return fmt.tex_sum((d[n], fmt.zpow_tex(n)) for n in sorted(d)) if d else "0"


def lsym(d):
    return " + ".join(f"{fmt.sym_num(d[n])}*{fmt.zpow_sym(n)}" for n in sorted(d)) if d else "0"


# ----------------------------------------------------------------------------- networks
def leaf(i):
    return {"type": "h", "i": i}


def evaluate(node, blocks, drop_delay=False, series_as_parallel=False, ignore_minus=False):
    t = node["type"]
    if t == "h":
        vals, start = blocks[node["i"] - 1]
        return {start + j: v for j, v in enumerate(vals) if v != 0}
    if t == "d":
        return {0: 1} if drop_delay else {node["k"]: 1}
    kids = [evaluate(c, blocks, drop_delay, series_as_parallel, ignore_minus) for c in node["kids"]]
    if t == "ser":
        if series_as_parallel and node.get("main"):          # blocks added instead of convolved (delays kept)
            out, dl = {}, {0: 1}
            for c, k in zip(node["kids"], kids):
                if c["type"] == "d":
                    dl = lmul(dl, k)
                else:
                    out = ladd(out, k)
            return lmul(out, dl)
        out = {0: 1}
        for k in kids:
            out = lmul(out, k)
        return out
    signs = [1, 1] if ignore_minus else node["signs"]
    return ladd(ladd({}, kids[0], signs[0]), kids[1], signs[1])


def _with_delay(branch, k):
    d = {"type": "d", "k": k}
    if branch["type"] == "h":
        kids = [branch, d] if random.random() < 0.5 else [d, branch]
        return {"type": "ser", "kids": kids}
    kids = list(branch["kids"])
    kids.insert(random.randint(0, len(kids)), d)
    return {"type": "ser", "kids": kids, "main": branch.get("main", False)}


def _network():
    shape = random.choice(["T1", "T2", "T3"])
    signs = [1, -1] if random.random() < 0.25 else [1, 1]
    k = random.choice([1, 1, 2])
    which = random.randint(0, 1)
    if shape == "T1":
        branches = [{"type": "ser", "kids": [leaf(1), leaf(2)], "main": True}, leaf(3)]
        branches[which] = _with_delay(branches[which], k)
        return shape, {"type": "par", "signs": signs, "kids": branches}
    if shape == "T2":
        branches = [leaf(2), leaf(3)]
        branches[which] = _with_delay(branches[which], k)
        return shape, {"type": "ser", "main": True, "kids": [leaf(1), {"type": "par", "signs": signs, "kids": branches}]}
    branches = [leaf(1), leaf(2)]
    branches[which] = _with_delay(branches[which], k)
    return shape, {"type": "ser", "main": True, "kids": [{"type": "par", "signs": signs, "kids": branches}, leaf(3)]}


def _block():
    L = random.choice([2, 2, 3])
    vals = [random.choice(ENDS)] + [random.randint(-2, 2) for _ in range(L - 2)] + [random.choice(ENDS)]
    return vals


def _has_minus(node):
    if node["type"] == "par" and -1 in node["signs"]:
        return True
    return any(_has_minus(c) for c in node.get("kids", []))


# ----------------------------------------------------------------------------- symbolic / LaTeX of the network
def net_tex(node, blocks, steps, zform=True):
    """LaTeX of the network in terms of H_i(z) (zform) and the list of evaluation steps."""
    t = node["type"]
    if t == "h":
        return f"H_{node['i']}(z)"
    if t == "d":
        return fmt.zpow_tex(node["k"])
    parts = [net_tex(c, blocks, steps) for c in node["kids"]]
    if t == "ser":
        out = ""
        for c, s in zip(node["kids"], parts):
            out += (r"\big[" + s + r"\big]") if c["type"] == "par" else s
        return out
    s0, s1 = node["signs"]
    return parts[0] + (" - " if s1 < 0 else " + ") + parts[1]


def collect_steps(node, blocks):
    """[(tex of the sub-network, tex of the factors with values, resulting polynomial tex)] for every
    series node with two or more children and every parallel node, innermost first."""
    steps = []

    def rec(nd):
        t = nd["type"]
        if t in ("h", "d"):
            return evaluate(nd, blocks)
        vals = [rec(c) for c in nd["kids"]]
        res = evaluate(nd, blocks)
        if t == "ser":
            fac = "".join(r"\left(" + ltex(v) + r"\right)" if len(v) > 1 else
                          (ltex(v) if list(v.values()) != [1] or list(v.keys()) != [0] else "")
                          for v in vals)
        else:
            s1 = nd["signs"][1]
            fac = r"\left(" + ltex(vals[0]) + r"\right)" + (" - " if s1 < 0 else " + ") + r"\left(" + ltex(vals[1]) + r"\right)"
        steps.append((net_tex(nd, blocks, []), fac, ltex(res)))
        return res
    rec(node)
    return steps


# ----------------------------------------------------------------------------- SVG block diagram
def _layout(node):
    t = node["type"]
    if t in ("h", "d"):
        return {"w": BW, "h": BH, "ym": BH / 2}
    ls = [_layout(c) for c in node["kids"]]
    if t == "ser":
        ym = max(l["ym"] for l in ls)
        below = max(l["h"] - l["ym"] for l in ls)
        return {"w": sum(l["w"] for l in ls) + GAP * (len(ls) - 1), "h": ym + below, "ym": ym, "kids": ls}
    wb = max(l["w"] for l in ls)
    yt = ls[0]["ym"]
    yb = ls[0]["h"] + VGAP + ls[1]["ym"]
    return {"w": 24 + wb + 28 + 2 * R, "h": ls[0]["h"] + VGAP + ls[1]["h"], "ym": (yt + yb) / 2,
            "kids": ls, "wb": wb, "yt": yt, "yb": yb}


def _line(x1, y1, x2, y2, arrow=True):
    m = ' marker-end="url(#w3sysarrow)"' if arrow else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="currentColor" '
            f'stroke-width="1.5"{m}/>')


def _draw(node, lay, x, y, out):
    """Draw node with its layout box's top-left corner at (x, y)."""
    t = node["type"]
    if t in ("h", "d"):
        label = f"h{SUB[node['i']]}[n]" if t == "h" else f"z{SUP[node['k']]}"
        out.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{BW}" height="{BH}" rx="3" fill="none" '
                   'stroke="currentColor" stroke-width="1.5"/>')
        out.append(f'<text x="{x + BW / 2:.1f}" y="{y + BH / 2 + 5:.1f}" text-anchor="middle" font-size="15">{label}</text>')
        return
    if t == "ser":
        ym = y + lay["ym"]
        cx = x
        for i, (c, cl) in enumerate(zip(node["kids"], lay["kids"])):
            if i:
                out.append(_line(cx - GAP, ym, cx - 1, ym))
            _draw(c, cl, cx, ym - cl["ym"], out)
            cx += cl["w"] + GAP
        return
    ym, yt, yb = y + lay["ym"], y + lay["yt"], y + lay["yb"]
    xs, xb = x + 10, x + 24
    xa = x + 24 + lay["wb"] + 28 + R
    out.append(_line(x, ym, xs, ym, arrow=False))
    out.append(f'<circle cx="{xs:.1f}" cy="{ym:.1f}" r="3" fill="currentColor"/>')
    out.append(_line(xs, yt, xs, yb, arrow=False))
    tops = [y, y + lay["kids"][0]["h"] + VGAP]
    for c, cl, top, yl in zip(node["kids"], lay["kids"], tops, (yt, yb)):
        out.append(_line(xs, yl, xb - 1, yl))
        _draw(c, cl, xb, top, out)
        xe = xb + cl["w"]
        yend = ym - R if yl < ym else ym + R
        out.append(f'<polyline points="{xe:.1f},{yl:.1f} {xa:.1f},{yl:.1f} {xa:.1f},{yend + (-2 if yl < ym else 2):.1f}" '
                   'fill="none" stroke="currentColor" stroke-width="1.5" marker-end="url(#w3sysarrow)"/>')
    out.append(f'<circle cx="{xa:.1f}" cy="{ym:.1f}" r="{R}" fill="none" stroke="currentColor" stroke-width="1.5"/>')
    out.append(_line(xa - 6, ym, xa + 6, ym, arrow=False))
    out.append(_line(xa, ym - 6, xa, ym + 6, arrow=False))
    for s, yl, dy in zip(node["signs"], (yt, yb), (-R - 5, R + 15)):
        out.append(f'<text x="{xa + 8:.1f}" y="{ym + dy:.1f}" font-size="15">{"+" if s > 0 else "−"}</text>')


def diagram(node):
    lay = _layout(node)
    left, right, pad = 56, 70, 14
    W = left + lay["w"] + right
    H = lay["h"] + 2 * pad
    ym = pad + lay["ym"]
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.0f} {H:.0f}" width="100%" '
           f'style="max-width:{W:.0f}px;display:block;margin:0.5em auto;color:inherit" font-family="sans-serif" '
           'role="img" aria-label="block diagram of the system">',
           '<defs><marker id="w3sysarrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
           'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>']
    out.append(f'<text x="4" y="{ym + 5:.1f}" font-size="15" font-style="italic">x[n]</text>')
    out.append(_line(40, ym, left - 1, ym))
    _draw(node, lay, left, pad, out)
    xr = left + lay["w"]
    out.append(_line(xr, ym, xr + 28, ym))
    out.append(f'<text x="{xr + 32:.1f}" y="{ym + 5:.1f}" font-size="15" font-style="italic">y[n]</text>')
    out.append("</svg>")
    return "".join(out)


# ----------------------------------------------------------------------------- generate
def _draw_variant():
    while True:
        shape, net = _network()
        blocks = [(_block(), 0) for _ in range(3)]
        if random.random() < 0.2:
            i = random.randint(0, 2)
            blocks[i] = (blocks[i][0], -1)
        h = evaluate(net, blocks)
        if not h:
            continue
        vals, start = lseq(h)
        if len(vals) < 3 or max(abs(v) for v in vals) > 12:
            continue
        wrong = {"delay": evaluate(net, blocks, drop_delay=True),
                 "add": evaluate(net, blocks, series_as_parallel=True)}
        if _has_minus(net):
            wrong["minus"] = evaluate(net, blocks, ignore_minus=True)
        if any(w == h or not w for w in wrong.values()) or wrong["delay"] == wrong["add"]:
            continue
        H1 = sum(h.values())
        if H1 == 0:
            continue
        if random.random() < 0.5:
            m = random.choice([1, 2, 3])
            cands = []
            for n0 in range(start + m, start + m + len(vals)):
                right = h.get(n0 - m, 0)
                if right == 0:
                    continue
                mist = {h.get(n0 + m, 0), h.get(n0, 0), wrong["delay"].get(n0 - m, 0), wrong["add"].get(n0 - m, 0)}
                mist.discard(right)
                if len(mist) >= 2 and n0 - m != start:
                    cands.append((n0, right))
            if not cands:
                continue
            n0, yval = random.choice(cands)
            out = {"kind": "delta", "m": m, "n0": n0, "y": yval}
        else:
            c = random.choice([2, 3, -1, -2])
            mist = {H1, c * sum(wrong["add"].values())}
            if "minus" in wrong:
                mist.add(c * sum(wrong["minus"].values()))
            mist.discard(c * H1)
            if len(mist) < 2:
                continue
            out = {"kind": "const", "c": c, "y": c * H1}
        return shape, net, blocks, h, wrong, out


def generate(data):
    shape, net, blocks, h, wrong, out = _draw_variant()
    vals, start = lseq(h)
    p = data["params"]
    p["shape"], p["net"] = shape, net
    p["blocks"] = [{"vals": v, "start": s} for v, s in blocks]
    p["fig_html"] = diagram(net)
    p["blocks_tex"] = r",\qquad ".join(f"h_{i + 1}[n] = {fmt.tex_seq(v, s)}" for i, (v, s) in enumerate(blocks))
    p["Hi_list"] = [{"tex": f"H_{i + 1}(z) = {ltex({s + j: x for j, x in enumerate(v) if x})}"}
                    for i, (v, s) in enumerate(blocks)]
    p["net_tex"] = net_tex(net, blocks, [])
    p["steps"] = [{"what": a, "fac": b, "res": c} for a, b, c in collect_steps(net, blocks)]
    p["H_tex"] = ltex(h)
    p["h_tex"] = fmt.tex_seq(vals, start)
    p["h_last"] = start + len(vals) - 1
    p["h_len"] = len(vals)
    p["h_plain"] = "[" + ", ".join(map(str, vals)) + "]"
    p["H_sym_alt"] = lsym(h)
    p["is_delta"] = out["kind"] == "delta"
    p["is_const"] = out["kind"] == "const"
    p["out"] = out
    if out["kind"] == "delta":
        m, n0 = out["m"], out["n0"]
        p["m_tex"] = f"n-{m}"
        p["n0"] = n0
        p["shift_tex"] = f"y[n] = h[n-{m}]"
        p["y_work_tex"] = f"y[{n0}] = h[{n0} - {m}] = h[{n0 - m}] = {out['y']}"
        p["y_label"] = f"$y[{n0}] =$"
    else:
        c = out["c"]
        p["c_tex"] = str(c)
        H1 = sum(h.values())
        p["H1_tex"] = " + ".join(fmt.tex_num(v, paren=True) for v in vals if v) + f" = {H1}"
        p["y_work_tex"] = rf"y[n] = {c}\cdot H(1) = {c}\cdot({H1}) = {out['y']}"
        p["y_label"] = r"$y[n] =$"
    c = data["correct_answers"]
    c["nh"] = start
    c["h"] = pl.to_json(np.array([[float(v) for v in vals]]))
    c["H"] = lsym(h)
    c["y0"] = int(out["y"])
