"""Independent check of questions/sampling/sketch-sampled-spectrum.

X_a(Omega) is rebuilt in numpy (rad/s) and X_d(omega) = (1/T) sum_k X_a((omega - 2 pi k)/T) is
summed over |k| <= 60 on a dense grid. Checked against it: the solution polyline (vertex values and
the straight segments between them), the pl-sketch `function` string (evaluated with Python's own
abs), the `match` key points (= the corners), the corner count and the canvas range. Probes mark
corners with the point tool via checklib.sketch_raw, computed here (own breakpoints, numpy values),
never from the server's vertex list; the helper polyline is drawn as the browser sends it."""

import math
from fractions import Fraction

import numpy as np
from checklib import sketch_raw


def _xa(spec):
    A, B = float(Fraction(spec["A"])), float(Fraction(spec["B"]))
    WB = 2 * math.pi * B
    shape = spec["shape"]

    def f(W):
        W = np.abs(np.asarray(W, float))
        if shape == "tri":
            return np.where(W < WB, A * (1 - W / WB), 0.0)
        W1 = float(Fraction(spec["r"])) * WB
        return np.where(W <= W1, A, np.where(W < WB, A * (WB - W) / (WB - W1), 0.0))
    return f


def xd(spec, u, ks=range(-60, 61), fs=None):
    """X_d at omega = pi*u (array); ks selects the copies that are added."""
    fs = float(Fraction(spec["fs"])) if fs is None else fs
    u = np.asarray(u, float)
    return sum(_xa(spec)((np.pi * u - 2 * np.pi * k) * fs) for k in ks) * fs


def breakpoints(spec, fs=None):
    fs = Fraction(spec["fs"]) if fs is None else fs
    b = 2 * Fraction(spec["B"]) / fs
    edges = [b, -b, Fraction(0)]
    if spec["shape"] == "trap":
        c = Fraction(spec["r"]) * b
        edges += [c, -c]
    pts = {Fraction(-1), Fraction(1)}
    for k in range(-4, 5):
        for e in edges:
            if -1 < e + 2 * k < 1:
                pts.add(e + 2 * k)
    return sorted(pts)


def polyline(spec, fs=None, ks=range(-60, 61), scale=1.0):
    xs = breakpoints(spec, fs)
    ys = xd(spec, [float(x) for x in xs], ks, None if fs is None else float(fs)) * scale
    return [(float(x), float(y)) for x, y in zip(xs, ys)]


def check(params, correct):
    probs = []
    spec = params["spec"]
    H = float(Fraction(spec["A"]) * Fraction(spec["fs"]))
    grid = np.linspace(-1, 1, 2001)
    exact = xd(spec, grid)
    verts = [(float(Fraction(x)), float(Fraction(y))) for x, y in params["vertices"]]
    vx, vy = zip(*verts)
    if abs(vx[0] + 1) > 1e-12 or abs(vx[-1] - 1) > 1e-12 or list(vx) != sorted(vx):
        probs.append(f"solution polyline does not run from -1 to 1: {vx}")
    if np.max(np.abs(np.interp(grid, vx, vy) - exact)) > 1e-6 * H:
        probs.append("solution polyline differs from the copy sum")
    for x, _ in verts:
        if abs(x * 60 - round(x * 60)) > 1e-9:
            probs.append(f"corner at x = {x} is not a multiple of 1/60")
    fn = params["fn"]
    vals = np.array([eval(fn, {"__builtins__": {}}, {"abs": abs, "x": float(x)}) for x in grid])
    if np.max(np.abs(vals - exact)) > 1e-6 * H:
        probs.append("pl-sketch function string differs from the copy sum")
    for kp in params["key_points"]:
        x, y = float(kp["x"]), float(kp["y"])
        if abs(float(xd(spec, x)) - y) > 1e-6 * H:
            probs.append(f"key point ({x}, {y}) is not on X_d")
    sol = [tuple(map(float, d["c"].strip("()").split(","))) for d in params["sol_pts"]]
    if len(sol) != len(verts) or any(abs(a - c) > 1e-9 or abs(b - d) > 1e-6 * H for (a, b), (c, d) in zip(sol, verts)):
        probs.append("solution points differ from the vertex list")
    if params["ncorners"] != len(verts) or len(params["key_points"]) != len(verts):
        probs.append("corner count / key points do not match the corners")
    mine = _simplify(polyline(spec))
    if len(mine) != len(verts) or any(abs(a[0] - c[0]) > 1e-9 or abs(a[1] - c[1]) > 1e-6 * H for a, c in zip(mine, verts)):
        probs.append(f"corners {verts} differ from the independent corners {mine}")
    if min(b[0] - a[0] for a, b in zip(verts, verts[1:])) < 0.05 - 1e-12:
        probs.append("two corners closer than 0.05 (point tool merges points within 10 px)")
    ymax = float(params["ymax"])
    if not (1.2 * exact.max() - 1e-9 <= ymax <= 3 * exact.max()):
        probs.append(f"canvas top {ymax} does not fit the peak {exact.max()}")
    if exact.min() < -1e-9:
        probs.append("negative X_d")
    r = params["xd"]["ranges"]
    axis_px = r["height"] - (0 - r["y_end"]) * r["height"] / (r["y_start"] - r["y_end"])   # pixels above the bottom edge
    if axis_px < 25:
        probs.append(f"the axis y = 0 is only {axis_px:.1f} px above the canvas bottom edge")
    if params["xd"]["ranges"]["width"] != 720:
        probs.append("canvas width changed")
    return probs


def _simplify(pts):
    """Drop collinear interior points (corners only)."""
    out = [pts[0]]
    for i in range(1, len(pts) - 1):
        (x0, y0), (x1, y1), (x2, y2) = out[-1], pts[i], pts[i + 1]
        if abs((y1 - y0) * (x2 - x1) - (y2 - y1) * (x1 - x0)) > 1e-9 * max(1.0, abs(y1)):
            out.append(pts[i])
    out.append(pts[-1])
    return out


def _pts(corners):
    return [[(x, y)] for x, y in corners]


def submissions(params, correct):
    spec = params["spec"]
    fs = Fraction(spec["fs"])
    ymax = float(params["ymax"])
    r = params["xd"]["ranges"]
    px_x = (r["x_end"] - r["x_start"]) / r["width"]            # graph units per pixel
    px_y = (r["y_end"] - r["y_start"]) / r["height"]
    corners = _simplify(polyline(spec))
    raw = lambda pts, **extra: sketch_raw(params, "xd", {"pt": _pts(pts), **extra})  # noqa: E731
    cases = [(raw(corners), {"xd": 1})]
    # every corner jittered by 3-8 px, the two ends 8 px inwards
    jit = []
    for i, (x, y) in enumerate(corners):
        dx = (3 + i % 6) * px_x * (1 if i % 2 else -1)
        dy = (8 - i % 6) * px_y * (1 if i % 3 else -1)
        if i == 0:
            dx = 8 * px_x
        if i == len(corners) - 1:
            dx = -8 * px_x
        jit.append((x + dx, max(y + dy, 0.0) if y == 0 else y + dy))
    cases.append((raw(jit), {"xd": 1}))
    # zero corners 9 px below the axis: 1; dragged onto the bottom edge of the canvas (clamped there,
    # ~30 px below the axis): 0
    if any(y == 0 for _, y in corners):
        cases.append((raw([(x, -9 * px_y if y == 0 else y) for x, y in corners]), {"xd": 1}))
        cases.append((raw([(x, r["y_start"] if y == 0 else y) for x, y in corners]), {"xd": 0}))
    # extra points on X_d (midpoints of the edges) and the helper line drawn through the corners: still 1
    mids = [((a[0] + b[0]) / 2, (a[1] + b[1]) / 2) for a, b in zip(corners, corners[1:])]
    cases.append((raw(sorted(corners + mids)), {"xd": 1}))
    cases.append((raw(corners, pl=[corners]), {"xd": 1}))
    # wrong height (3/4 of it), the 1/T forgotten
    cases.append((raw([(x, 0.75 * y) for x, y in corners]), {"xd": 0}))
    cases.append((raw([(x, y / float(fs)) for x, y in corners]), {"xd": 0}))
    # wrong band edge: Omega_B T computed with 3/4 of the true T (edge at 3/4), same height
    wrong_fs = fs * Fraction(4, 3)
    edge = _simplify([(x, y * float(fs / wrong_fs)) for x, y in polyline(spec, fs=wrong_fs)])
    cases.append((raw(edge), {"xd": 0}))
    # a corner left out (an interior one when there is one, else an end)
    drop = len(corners) // 2 if len(corners) > 2 else 0
    cases.append((raw(corners[:drop] + corners[drop + 1:]), {"xd": 0}))
    # overlap not added: only the copy centred at 0 (also the upper envelope of the copies on [-1, 1])
    single = polyline(spec, ks=[0])
    if max(abs(a[1] - b[1]) for a, b in zip(single, polyline(spec))) > 0.1 * ymax:
        cases.append((raw(_simplify(single)), {"xd": 0}))
        # ... and its corners plus the neighbouring copy's corners marked separately (not added)
        b = 2 * Fraction(spec["B"]) / fs
        right = [(x, y) for x, y in polyline(spec, ks=[1]) if x >= float(2 - b) - 1e-12]
        left = [(-x, y) for x, y in right]
        cases.append((raw(sorted(set(_simplify(single) + right + left))), {"xd": 0}))
    return cases
