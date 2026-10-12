"""Independent check of questions/convolution/system-interconnection.

The block sequences are parsed back from the displayed LaTeX; the network (params["net"]) is
evaluated with numpy.convolve on (samples, start index) pairs; the SVG is checked to show every block
label, one adder per parallel connection and a minus sign exactly when the network has one. H(z) is
compared with SymPy, the output value recomputed from the numpy h[n]."""

import re

import numpy as np
import sympy
from checklib import S, Z

SUB = {1: "₁", 2: "₂", 3: "₃"}
SUP = {1: "⁻¹", 2: "⁻²"}


def parse_seq(tex):
    items = tex.strip()[2:-2].split(r",\ ")
    vals, zero_at = [], None
    for i, it in enumerate(items):
        m = re.fullmatch(r"\\underset\{\\uparrow\}\{(-?\d+)\}", it)
        if m:
            zero_at = i
        vals.append(int(m.group(1)) if m else int(it))
    start = -zero_at
    while vals[0] == 0:
        vals.pop(0)
        start += 1
    while vals[-1] == 0:
        vals.pop()
    return vals, start


def blocks_from_tex(params):
    out = []
    for part in params["blocks_tex"].split(r",\qquad "):
        lhs, rhs = part.split(" = ")
        out.append(parse_seq(rhs))
    return out


def trim(v, s):
    v = list(v)
    while v and v[0] == 0:
        v.pop(0)
        s += 1
    while v and v[-1] == 0:
        v.pop()
    return v, s


def add(a, b, sign=1):
    (va, sa), (vb, sb) = a, b
    if not va:
        return trim([sign * x for x in vb], sb)
    if not vb:
        return a
    lo = min(sa, sb)
    hi = max(sa + len(va), sb + len(vb))
    out = np.zeros(hi - lo, dtype=int)
    out[sa - lo:sa - lo + len(va)] += np.array(va, dtype=int)
    out[sb - lo:sb - lo + len(vb)] += sign * np.array(vb, dtype=int)
    return trim(out.tolist(), lo)


def conv(a, b):
    (va, sa), (vb, sb) = a, b
    return trim(np.convolve(va, vb).astype(int).tolist(), sa + sb)


def ev(node, blocks, mode=None):
    t = node["type"]
    if t == "h":
        return blocks[node["i"] - 1]
    if t == "d":
        return ([1], 0) if mode == "delay" else ([1], node["k"])
    kids = [ev(c, blocks, mode) for c in node["kids"]]
    if t == "ser":
        if mode == "add" and node.get("main"):
            acc, dl = ([], 0), ([1], 0)
            for c, k in zip(node["kids"], kids):
                if c["type"] == "d":
                    dl = conv(dl, k)
                else:
                    acc = add(acc, k)
            return conv(acc, dl) if acc[0] else acc
        acc = ([1], 0)
        for k in kids:
            acc = conv(acc, k)
        return acc
    signs = [1, 1] if mode == "minus" else node["signs"]
    return add(add(([], 0), kids[0], signs[0]), kids[1], signs[1])


def count(node, t):
    return (node["type"] == t) + sum(count(c, t) for c in node.get("kids", []))


def has_minus(node):
    return (node["type"] == "par" and -1 in node["signs"]) or any(has_minus(c) for c in node.get("kids", []))


def delays(node):
    return ([node["k"]] if node["type"] == "d" else []) + [k for c in node.get("kids", []) for k in delays(c)]


def at(seq, n):
    v, s = seq
    return v[n - s] if 0 <= n - s < len(v) else 0


def zpoly(seq):
    v, s = seq
    return sum(c * Z ** (-(s + i)) for i, c in enumerate(v))


def mistakes(params):
    blocks = blocks_from_tex(params)
    net = params["net"]
    out = {"delay": ev(net, blocks, "delay"), "add": ev(net, blocks, "add")}
    if has_minus(net):
        out["minus"] = ev(net, blocks, "minus")
    return out


def y_mistakes(params, correct):
    blocks = blocks_from_tex(params)
    h = ev(params["net"], blocks)
    w = mistakes(params)
    o = params["out"]
    right = correct["y0"]
    if o["kind"] == "delta":
        m, n0 = o["m"], o["n0"]
        c = [at(h, n0 + m), at(h, n0), at(w["delay"], n0 - m), at(w["add"], n0 - m)]
    else:
        c = [sum(h[0]), o["c"] * sum(w["add"][0])] + ([o["c"] * sum(w["minus"][0])] if "minus" in w else [])
    return [v for i, v in enumerate(c) if v != right and v not in c[:i]]


def check(params, correct):
    probs = []
    blocks = blocks_from_tex(params)
    if [list(b) for b in blocks] != [[b["vals"], b["start"]] for b in params["blocks"]]:
        probs.append(f"blocks_tex {blocks} != params {params['blocks']}")
    net = params["net"]
    h = ev(net, blocks)
    if correct["nh"] != h[1]:
        probs.append(f"first index {correct['nh']} != {h[1]}")
    got = [int(round(v)) for v in correct["h"]["_value"][0]]
    if got != h[0]:
        probs.append(f"h {got} != numpy {h[0]}")
    if sympy.expand(S(correct["H"]) - zpoly(h)) != 0:
        probs.append(f"H = {correct['H']} != {zpoly(h)}")
    o = params["out"]
    want = at(h, o["n0"] - o["m"]) if o["kind"] == "delta" else o["c"] * sum(h[0])
    if correct["y0"] != want:
        probs.append(f"y0 = {correct['y0']} != {want}")
    w = mistakes(params)
    if any(v == h for v in w.values()) or w["delay"] == w["add"]:
        probs.append("a classic mistake gives the correct h (probe meaningless)")
    if len(y_mistakes(params, correct)) < 2:
        probs.append("fewer than two classic mistakes for the output value")
    fig = params["fig_html"]
    for i in (1, 2, 3):
        if fig.count(f"h{SUB[i]}[n]") != 1:
            probs.append(f"figure shows h{i} {fig.count(f'h{SUB[i]}[n]')} times")
    for k in delays(net):
        if f"z{SUP[k]}" not in fig:
            probs.append(f"figure lacks the delay z^-{k}")
    if fig.count('r="11"') != count(net, "par"):
        probs.append("number of adders in the figure != parallel connections")
    if ("−</text>" in fig) != has_minus(net):
        probs.append("minus sign in the figure does not match the network")
    if count(net, "d") != 1 or count(net, "h") != 3:
        probs.append("network must have three blocks and one delay")
    return probs


def _vec(v):
    return "[" + ", ".join(map(str, v)) + "]"


def _poly_str(seq):
    v, s = seq
    return " + ".join(f"({c})*z^({-(s + i)})" for i, c in enumerate(v) if c) or "0"


def _product_form(node, blocks):
    t = node["type"]
    if t == "h":
        return "(" + _poly_str(blocks[node["i"] - 1]) + ")"
    if t == "d":
        return f"z^(-{node['k']})"
    parts = [_product_form(c, blocks) for c in node["kids"]]
    if t == "ser":
        return "*".join(parts)
    return f"({parts[0]} {'-' if node['signs'][1] < 0 else '+'} {parts[1]})"


def submissions(params, correct):
    blocks = blocks_from_tex(params)
    h = ev(params["net"], blocks)
    v, s = h
    top = s + len(v) - 1
    over = " + ".join(f"({c})*z^({top - (s + i)})" for i, c in enumerate(v) if c)     # all over z^top
    cases = [
        ({"h": "[" + " ".join(map(str, v)) + "]"}, {"h": 1}),
        ({"h": _vec([f"{c}.0" for c in v])}, {"h": 1}),
        ({"H": _product_form(params["net"], blocks)}, {"H": 1}),                     # the diagram as a formula
        ({"H": f"({over})/z^({top})"}, {"H": 1}),
        ({"H": "x + 1"}, {"H": "invalid"}),
        ({"nh": "0.5"}, {"nh": "invalid"}),
    ]
    for name, wseq in mistakes(params).items():
        wv, ws = wseq
        cases.append(({"nh": str(ws), "h": _vec(wv)}, {"nh": 0 if ws != s else 1, "h": 0 if wv != v else 1}))
        cases.append(({"H": _poly_str(wseq)}, {"H": 0}))
    y = correct["y0"]
    cases += [({"y0": f"+{y}" if y >= 0 else f"-{abs(y)}"}, {"y0": 1}), ({"y0": f" {y} "}, {"y0": 1})]
    cases += [({"y0": str(m)}, {"y0": 0}) for m in y_mistakes(params, correct)]
    return cases
