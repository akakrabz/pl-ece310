#!/usr/bin/env python3
"""Generate the inline-SVG figures for the ECE 329 notes.

Each figure is written to /home/claude/work/figs/<name>.html as a single
<figure class="ece-fig"> block with NO blank lines (so CommonMark keeps it as
one raw-HTML block when it is inlined into the markdown).

Conventions: strokes/text use currentColor (theme-aware); accents use the CSS
variables --accent (secondary), --accent2 (tertiary), --hi (hot red), --muted.
"""
import math, os, textwrap
import numpy as np

import sys
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
os.makedirs(OUT, exist_ok=True)

FONT = "font-family:'Source Sans 3','Source Sans Pro',system-ui,sans-serif"
MATHF = "font-family:'STIX Two Math','Latin Modern Math','Cambria Math','Times New Roman',serif;font-style:italic"

# --------------------------------------------------------------------------- helpers
def svg_open(w, h, vb=None):
    vb = vb or f"0 0 {w} {h}"
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb}" width="{w}" height="{h}" '
            f'role="img" style="{FONT};font-size:14px">'
            '<defs>'
            '<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            '<path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>'
            '<marker id="ahs" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">'
            '<path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>'
            '</defs>')

def svg_close():
    return "</svg>"

def line(x1, y1, x2, y2, stroke="currentColor", w=1.5, dash=None, arrow=False, small=False, opacity=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = ' marker-end="url(#ahs)"' if (arrow and small) else (' marker-end="url(#ah)"' if arrow else "")
    o = f' opacity="{opacity}"' if opacity is not None else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{w}"{d}{m}{o} stroke-linecap="round"/>'

def arrow(x1, y1, x2, y2, stroke="var(--accent)", w=2, small=False, opacity=None):
    return line(x1, y1, x2, y2, stroke=stroke, w=w, arrow=True, small=small, opacity=opacity)

def text(x, y, s, size=14, anchor="middle", fill="currentColor", math=False, weight=None, dy=None, rotate=None):
    st = f"font-size:{size}px;" + (MATHF + ";" if math else "") + (f"font-weight:{weight};" if weight else "")
    r = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    d = f' dy="{dy}"' if dy is not None else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" fill="{fill}" style="{st}"{r}{d}>{s}</text>'

def circle(cx, cy, r, fill="none", stroke="currentColor", w=1.5, dash=None, opacity=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' opacity="{opacity}"' if opacity is not None else ""
    return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}" stroke="{stroke}" stroke-width="{w}"{d}{o}/>'

def path(d, fill="none", stroke="currentColor", w=1.5, dash=None, opacity=None, arrow=False):
    dd = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' opacity="{opacity}"' if opacity is not None else ""
    m = ' marker-end="url(#ah)"' if arrow else ""
    return f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{w}"{dd}{o}{m} stroke-linejoin="round" stroke-linecap="round"/>'

def rect(x, y, w, h, fill="none", stroke="currentColor", sw=1.5, dash=None, opacity=None, rx=0):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' opacity="{opacity}"' if opacity is not None else ""
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}{o}/>'

def charge(cx, cy, sign="+", r=11, fill="var(--hi)"):
    s = [circle(cx, cy, r, fill=fill, stroke="none")]
    s.append(line(cx-5, cy, cx+5, cy, stroke="white", w=2.2))
    if sign == "+":
        s.append(line(cx, cy-5, cx, cy+5, stroke="white", w=2.2))
    return "".join(s)

def figure(name, svg, caption):
    html = f'<figure class="ece-fig">{svg}<figcaption>{caption}</figcaption></figure>'
    html = "\n".join(l for l in html.splitlines() if l.strip())  # no blank lines
    with open(f"{OUT}/{name}.html", "w") as f:
        f.write(html + "\n")
    print("wrote", name)

# 2-D "isometric-ish" projection for the 3-D sketches: x toward viewer-left, y right, z up
def P(x, y, z, ox=0, oy=0, s=1.0):
    X = ox + s * (y * 1.0 - x * 0.55)
    Y = oy - s * (z * 1.0 - x * 0.35)
    return X, Y

# --------------------------------------------------------------------------- 1. dl / dS on a cube
def fig_dl_ds_cube():
    W, H = 640, 300
    o = [svg_open(W, H)]
    # --- left: a curve with tangent dl elements
    ox, oy = 40, 60
    o.append(text(ox+90, oy-25, "dl is always tangent to the path", size=13, anchor="middle"))
    d = f"M{ox},{oy+170} C{ox+40},{oy+40} {ox+120},{oy+200} {ox+200},{oy+60}"
    o.append(path(d, w=2))
    # tangent arrows at a few parameter values (approximate derivative of the cubic)
    def bez(t):
        p0=(ox,oy+170); p1=(ox+40,oy+40); p2=(ox+120,oy+200); p3=(ox+200,oy+60)
        x=(1-t)**3*p0[0]+3*(1-t)**2*t*p1[0]+3*(1-t)*t**2*p2[0]+t**3*p3[0]
        y=(1-t)**3*p0[1]+3*(1-t)**2*t*p1[1]+3*(1-t)*t**2*p2[1]+t**3*p3[1]
        dx=3*(1-t)**2*(p1[0]-p0[0])+6*(1-t)*t*(p2[0]-p1[0])+3*t**2*(p3[0]-p2[0])
        dy=3*(1-t)**2*(p1[1]-p0[1])+6*(1-t)*t*(p2[1]-p1[1])+3*t**2*(p3[1]-p2[1])
        n=math.hypot(dx,dy); return x,y,dx/n,dy/n
    for t in (0.18, 0.5, 0.82):
        x,y,ux,uy = bez(t)
        o.append(arrow(x, y, x+34*ux, y+34*uy, stroke="var(--accent)", w=2.2))
    x,y,ux,uy = bez(0.5)
    o.append(text(x+18, y+30, "dl", size=15, math=True, fill="var(--accent)"))
    o.append(text(ox-6, oy+178, "P₁", size=13, anchor="end"))
    o.append(text(ox+206, oy+58, "P₂", size=13, anchor="start"))
    o.append(text(ox+100, oy+215, "Cartesian: dl = x̂ dx + ŷ dy + ẑ dz", size=13))
    # --- right: cube with outward dS on three visible faces
    cx, cy, s = 430, 190, 95
    def Q(x,y,z): return P(x,y,z,cx,cy,s)
    verts = {k: Q(*v) for k, v in {
        "000":(0,0,0),"100":(1,0,0),"010":(0,1,0),"110":(1,1,0),
        "001":(0,0,1),"101":(1,0,1),"011":(0,1,1),"111":(1,1,1)}.items()}
    def poly(keys, fill, op=0.18):
        pts = " ".join(f"{verts[k][0]:.1f},{verts[k][1]:.1f}" for k in keys)
        return f'<polygon points="{pts}" fill="{fill}" fill-opacity="{op}" stroke="currentColor" stroke-width="1.4" stroke-linejoin="round"/>'
    # hidden edges dashed
    for a,b in (("000","010"),("000","001"),("000","100")):
        o.append(line(*verts[a], *verts[b], dash="4 4", w=1.1, opacity=0.6))
    o.append(poly(["100","110","111","101"], "var(--accent)"))   # front face (x=1)
    o.append(poly(["010","110","111","011"], "var(--accent2)"))  # right face (y=1)
    o.append(poly(["001","101","111","011"], "var(--muted)"))    # top face (z=1)
    # face-centre normals
    fx = Q(1,0.5,0.5); fx2 = Q(1.55,0.5,0.5)
    fy = Q(0.5,1,0.5); fy2 = Q(0.5,1.55,0.5)
    fz = Q(0.5,0.5,1); fz2 = Q(0.5,0.5,1.55)
    o.append(arrow(*fx, *fx2, stroke="var(--accent)", w=2.4))
    o.append(arrow(*fy, *fy2, stroke="var(--accent2)", w=2.4))
    o.append(arrow(*fz, *fz2, stroke="currentColor", w=2.4))
    o.append(text(fx2[0]-10, fx2[1]+22, "dS = x̂ dy dz", size=12, anchor="middle"))
    o.append(text(fy2[0]+4, fy2[1]+4, "dS = ŷ dz dx", size=12, anchor="start"))
    o.append(text(fz2[0], fz2[1]-6, "dS = ẑ dx dy", size=12))
    # edge labels
    e1a, e1b = verts["100"], verts["110"]
    o.append(text((e1a[0]+e1b[0])/2, (e1a[1]+e1b[1])/2+16, "dy", size=12, math=True))
    e2a, e2b = verts["110"], verts["111"]
    o.append(text(e2a[0]+12, (e2a[1]+e2b[1])/2+4, "dz", size=12, math=True, anchor="start"))
    e3a, e3b = verts["100"], verts["000"]
    o.append(text((e3a[0]+e3b[0])/2-14, (e3a[1]+e3b[1])/2+2, "dx", size=12, math=True))
    o.append(text(cx+10, 292, "closed surface: every dS points outward", size=13))
    o.append(svg_close())
    figure("dl-ds-cube", "".join(o),
           "<strong>The two new objects of Lecture 1.</strong> Left: the differential length vector d<b>l</b> is tangent to the path and changes direction along it. Right: a differential surface vector d<b>S</b> is normal to the surface; for a closed surface the convention is <em>outward</em>. Its direction comes from the order of the cross product of two tangent vectors (swap the order and it flips).")

# --------------------------------------------------------------------------- 2. Coulomb pair
def fig_coulomb_pair():
    W, H = 640, 200
    o = [svg_open(W, H)]
    x1, x2, y = 150, 470, 95
    o.append(line(x1, y, x2, y, dash="6 5", w=1.2, opacity=0.6))
    o.append(text((x1+x2)/2, y-10, "R", size=16, math=True))
    o.append(charge(x1, y, "+"))
    o.append(charge(x2, y, "+"))
    o.append(text(x1, y+32, "Q₁", size=15, math=True))
    o.append(text(x2, y+32, "Q₂", size=15, math=True))
    # unit vector r̂ from 1 to 2, drawn at 2
    o.append(arrow(x2+16, y-40, x2+60, y-40, stroke="currentColor", w=1.8))
    o.append(text(x2+38, y-48, "R̂ (from 1 to 2)", size=12))
    # forces
    o.append(arrow(x2+16, y, x2+120, y, stroke="var(--hi)", w=3))
    o.append(text(x2+70, y+20, "F₂ on Q₂", size=13, fill="var(--hi)"))
    o.append(arrow(x1-16, y, x1-120, y, stroke="var(--hi)", w=3))
    o.append(text(x1-70, y+20, "F₁ on Q₁", size=13, fill="var(--hi)"))
    o.append(text(W/2, 160, "F₂ = Q₁Q₂ R̂ / (4πε₀R²) ;  like charges: F₂ ∥ +R̂ (repel), unlike: F₂ ∥ −R̂ (attract);  F₁ = −F₂", size=13))
    o.append(svg_close())
    figure("coulomb-pair", "".join(o),
           "<strong>Coulomb's law, with the bookkeeping the course uses.</strong> Label the <em>source</em> 1 and the point where you want the field 2; the unit vector R̂ always runs from the source to the field point. The force on Q₂ is along +R̂ for like charges and along −R̂ for unlike charges. Newton's third law gives F₁ = −F₂.")

# --------------------------------------------------------------------------- 3. Finite line charge geometry
def fig_line_charge_side():
    W, H = 640, 330
    o = [svg_open(W, H)]
    ax, top, bot = 150, 30, 300
    zc = (top+bot)/2
    # z axis with charged segment
    o.append(line(ax, bot+10, ax, top-10, w=1.2, arrow=True, small=True, opacity=0.7))
    o.append(text(ax-4, top-14, "z", size=15, math=True, anchor="end"))
    o.append(line(ax, zc+95, ax, zc-95, stroke="var(--hi)", w=5))
    for k in range(-4, 5):
        o.append(text(ax-11, zc-k*22+4, "+", size=12, fill="var(--hi)"))
    o.append(text(ax-30, zc-95, "z = a", size=12, anchor="end"))
    o.append(text(ax-30, zc+100, "z = −a", size=12, anchor="end"))
    o.append(text(ax-40, zc+4, "ρₗ [C/m]", size=13, anchor="end", fill="var(--hi)"))
    # field point
    px, py = 470, zc
    o.append(line(ax, zc, px, py, dash="5 5", w=1.2, opacity=0.6))
    o.append(text((ax+px)/2, zc+18, "r", size=15, math=True))
    o.append(circle(px, py, 4, fill="currentColor"))
    o.append(text(px+8, py+18, "P", size=15, math=True, anchor="start"))
    # element dQ at +z and mirror at -z
    zq = zc-62
    o.append(rect(ax-5, zq-6, 10, 12, fill="var(--accent)", stroke="none"))
    o.append(text(ax+12, zq-8, "dQ = ρₗ dz", size=13, anchor="start", fill="var(--accent)"))
    o.append(rect(ax-5, zc+62-6, 10, 12, fill="var(--accent2)", stroke="none"))
    o.append(text(ax+12, zc+62+18, "mirror element", size=12, anchor="start", fill="var(--accent2)"))
    # R vectors from elements to P
    o.append(line(ax, zq, px, py, stroke="var(--accent)", w=1.6))
    o.append(line(ax, zc+62, px, py, stroke="var(--accent2)", w=1.6, dash="4 3"))
    o.append(text((ax+px)/2-20, (zq+py)/2-8, "R", size=15, math=True, fill="var(--accent)"))
    # angle alpha at P between R and horizontal
    o.append(path(f"M{px-40},{py} A40,40 0 0 0 {px-40*math.cos(math.atan2(zc-zq, px-ax)):.1f},{py-40*math.sin(math.atan2(zc-zq, px-ax)):.1f}", w=1.2))
    o.append(text(px-52, py-12, "α", size=15, math=True))
    # dE at P: along R direction away from element (down-right), components
    ux, uy = (px-ax), (py-zq); n = math.hypot(ux,uy); ux/=n; uy/=n
    L = 90
    o.append(arrow(px, py, px+L*ux, py+L*uy, stroke="var(--accent)", w=2.6))
    o.append(text(px+L*ux+10, py+L*uy+6, "dE", size=15, math=True, anchor="start", fill="var(--accent)"))
    o.append(arrow(px, py, px+L*ux, py, stroke="var(--accent)", w=1.6, opacity=0.8))
    o.append(text(px+L*ux/2+6, py+16, "dEᵣ = dE cos α", size=12, fill="var(--accent)"))
    o.append(line(px+L*ux, py, px+L*ux, py+L*uy, stroke="var(--accent)", w=1.2, dash="3 3", opacity=0.8))
    # mirror element contribution (up-right)
    o.append(arrow(px, py, px+L*ux, py-L*uy, stroke="var(--accent2)", w=2.2, opacity=0.9))
    o.append(text(px+L*ux+10, py-L*uy-2, "dE′ (mirror)", size=13, anchor="start", fill="var(--accent2)"))
    o.append(text(W/2+40, H-12, "z-components cancel in pairs → only the radial (r̂) part survives the integral", size=13))
    o.append(svg_close())
    figure("line-charge-side", "".join(o),
           "<strong>Field of a finite line charge — the symmetry step of the 5-step recipe.</strong> Every element dQ at height z has a mirror element at −z. Their fields at P have equal magnitude; the z-components cancel and the radial components add, so only dE cos α needs integrating.")

# --------------------------------------------------------------------------- 4. Dipole-like pair field map (computed)
def fig_dipole_map():
    W, H = 640, 420
    o = [svg_open(W, H)]
    # world coords: x in [-3,3], y in [-2,2] mapped to canvas
    def M(x, y): return 40 + (x+3)/6*560, 20 + (2-y)/4*380
    q = [(+2.0, -1.0, 0.0), (-1.0, 1.0, 0.0)]  # (charge, x, y)
    def E(x, y):
        ex = ey = 0.0
        for qq, qx, qy in q:
            dx, dy = x-qx, y-qy
            r3 = (dx*dx+dy*dy)**1.5 + 1e-12
            ex += qq*dx/r3; ey += qq*dy/r3
        return ex, ey
    # unit-vector samples
    for gx in np.linspace(-2.8, 2.8, 15):
        for gy in np.linspace(-1.8, 1.8, 10):
            if min(math.hypot(gx-qx, gy-qy) for _,qx,qy in q) < 0.25: continue
            ex, ey = E(gx, gy); n = math.hypot(ex, ey)
            ux, uy = ex/n, ey/n
            X, Y = M(gx, gy); L = 13
            o.append(line(X, Y, X+L*ux, Y-L*uy, stroke="currentColor", w=1.3, arrow=True, small=True, opacity=0.75))
            o.append(circle(X, Y, 1.6, fill="var(--muted)", stroke="none"))
    # a few field lines from the +2Q charge (RK2), stop near -Q or off-canvas
    for k in range(12):
        th = 2*math.pi*k/12 + 0.13
        x, y = q[0][1]+0.22*math.cos(th), q[0][2]+0.22*math.sin(th)
        pts = [M(x, y)]
        for _ in range(900):
            ex, ey = E(x, y); n = math.hypot(ex, ey); h = 0.02
            xm, ym = x + 0.5*h*ex/n, y + 0.5*h*ey/n
            ex2, ey2 = E(xm, ym); n2 = math.hypot(ex2, ey2)
            x, y = x + h*ex2/n2, y + h*ey2/n2
            pts.append(M(x, y))
            if math.hypot(x-q[1][1], y-q[1][2]) < 0.2 or abs(x) > 3.1 or abs(y) > 2.1: break
        d = "M" + " L".join(f"{X:.1f},{Y:.1f}" for X, Y in pts)
        o.append(path(d, stroke="var(--accent)", w=1.4, opacity=0.85))
    # charges
    X, Y = M(q[0][1], q[0][2]); o.append(charge(X, Y, "+", r=13)); o.append(text(X, Y+30, "+2Q", size=14, weight=600))
    X, Y = M(q[1][1], q[1][2]); o.append(charge(X, Y, "-", r=11, fill="var(--accent)")); o.append(text(X, Y+28, "−Q", size=14, weight=600))
    # null point: on the x axis beyond -Q where 2/(x+1)^2 = 1/(x-1)^2 -> x = 3+2√2 ≈ 5.83 (off canvas) -> annotate
    o.append(text(W-14, H-8, "field null lies on the axis at x ≈ 5.8 (off the map)", size=11, anchor="end", fill="var(--muted)"))
    o.append(svg_close())
    figure("dipole-map", "".join(o),
           "<strong>Superposition, drawn.</strong> Unit vectors of <b>E</b> (grey samples) and a few field lines for a +2Q charge at (−1, 0) and a −Q charge at (+1, 0). Lines start on positive charge and end on negative charge or run off to infinity; because the charges are unequal, only some of the lines leaving +2Q end on −Q. Field lines never cross: <b>E</b> has one direction at each point.")

# --------------------------------------------------------------------------- 5. Flux through a patch
def fig_flux_patch():
    W, H = 640, 230
    o = [svg_open(W, H)]
    # uniform field arrows across the whole width
    for k in range(6):
        y = 45 + k*28
        o.append(arrow(20, y, 300, y, stroke="var(--accent)", w=1.6, opacity=0.55))
        o.append(arrow(340, y, 620, y, stroke="var(--accent)", w=1.6, opacity=0.55))
    # left: patch facing the field (normal parallel)
    cx, cy = 160, 115
    o.append(f'<ellipse cx="{cx}" cy="{cy}" rx="10" ry="70" fill="var(--accent2)" fill-opacity="0.25" stroke="currentColor" stroke-width="1.6"/>')
    o.append(arrow(cx, cy, cx+70, cy, stroke="currentColor", w=2.4))
    o.append(text(cx+62, cy-10, "n̂", size=15, math=True))
    o.append(text(cx, 210, "ψ = E · A  (max)", size=14))
    # right: tilted patch
    cx, cy, ang = 480, 115, 40
    o.append(f'<ellipse cx="{cx}" cy="{cy}" rx="10" ry="70" transform="rotate({-ang} {cx} {cy})" fill="var(--accent2)" fill-opacity="0.25" stroke="currentColor" stroke-width="1.6"/>')
    a = math.radians(ang)
    o.append(arrow(cx, cy, cx+70*math.cos(a), cy-70*math.sin(a), stroke="currentColor", w=2.4))
    o.append(text(cx+72*math.cos(a)+8, cy-72*math.sin(a), "n̂", size=15, math=True, anchor="start"))
    o.append(path(f"M{cx+38},{cy} A38,38 0 0 0 {cx+38*math.cos(a):.1f},{cy-38*math.sin(a):.1f}", w=1.2))
    o.append(text(cx+50, cy-14, "α", size=15, math=True))
    o.append(text(cx, 210, "ψ = E A cos α = Eₙ A", size=14))
    o.append(text(W/2, 22, "uniform E", size=13, fill="var(--accent)"))
    o.append(svg_close())
    figure("flux-patch", "".join(o),
           "<strong>Flux = how many field lines pierce the patch.</strong> Only the component of the field along the patch normal n̂ counts: Δψ = <b>E</b>·Δ<b>S</b> = (E cos α) ΔS. Tilting the patch by α reduces the flux by cos α; edge-on (α = 90°) nothing passes through.")

# --------------------------------------------------------------------------- 6. Gaussian cylinder about a line charge
def fig_gauss_cylinder():
    W, H = 640, 340
    o = [svg_open(W, H)]
    cx, top, bot, rx, ry = 250, 70, 250, 110, 28
    # line charge along the axis
    o.append(line(cx, 330, cx, 12, stroke="var(--hi)", w=4, arrow=True))
    o.append(line(cx, 330, cx, 12, stroke="var(--hi)", w=4))
    o.append(text(cx+8, 26, "ρₗ [C/m]", size=13, anchor="start", fill="var(--hi)"))
    # cylinder body
    o.append(f'<path d="M{cx-rx},{top} L{cx-rx},{bot} A{rx},{ry} 0 0 0 {cx+rx},{bot} L{cx+rx},{top}" fill="var(--accent)" fill-opacity="0.10" stroke="currentColor" stroke-width="1.5"/>')
    o.append(f'<ellipse cx="{cx}" cy="{top}" rx="{rx}" ry="{ry}" fill="var(--accent)" fill-opacity="0.18" stroke="currentColor" stroke-width="1.5"/>')
    o.append(f'<path d="M{cx-rx},{bot} A{rx},{ry} 0 0 1 {cx+rx},{bot}" fill="none" stroke="currentColor" stroke-width="1.2" stroke-dasharray="4 4"/>')
    # radius and height labels
    o.append(line(cx, top, cx+rx, top, dash="4 3", w=1.2))
    o.append(text(cx+rx/2, top-8, "r", size=15, math=True))
    o.append(line(cx+rx+22, top, cx+rx+22, bot, w=1.2))
    o.append(line(cx+rx+16, top, cx+rx+28, top, w=1.2)); o.append(line(cx+rx+16, bot, cx+rx+28, bot, w=1.2))
    o.append(text(cx+rx+34, (top+bot)/2+5, "L", size=15, math=True, anchor="start"))
    # dS arrows: top cap up, side radial
    o.append(arrow(cx-55, top-10, cx-55, top-45, stroke="currentColor", w=2))
    o.append(text(cx-55, top-52, "dS", size=13, math=True))
    for yy in (140, 200):
        o.append(arrow(cx+rx, yy, cx+rx+40, yy, stroke="currentColor", w=2))
        o.append(arrow(cx-rx, yy, cx-rx-40, yy, stroke="currentColor", w=2))
    o.append(text(cx+rx+46, 145, "dS", size=13, math=True, anchor="start"))
    # D / E arrows (radial, outside surface too)
    for yy in (115, 170, 225):
        o.append(arrow(cx+12, yy, cx+rx-8, yy, stroke="var(--accent)", w=2.2, opacity=0.9))
        o.append(arrow(cx-12, yy, cx-rx+8, yy, stroke="var(--accent)", w=2.2, opacity=0.9))
    o.append(text(cx+40, 108, "D, E", size=13, fill="var(--accent)"))
    # bookkeeping on the right
    X = 410
    o.append(text(X, 90, "∮ D·dS = top + bottom + side", size=13, anchor="start", weight=600))
    o.append(text(X, 118, "top:  D ⊥ dS  →  0", size=14, anchor="start"))
    o.append(text(X, 143, "bottom:  0", size=14, anchor="start"))
    o.append(text(X, 168, "side:  D · 2πrL", size=14, anchor="start"))
    o.append(text(X, 205, "Q_enc = ρₗ L", size=14, anchor="start"))
    o.append(text(X, 240, "⇒  D = ρₗ / (2πr)   [C/m²]", size=14, anchor="start", weight=600))
    o.append(text(X, 265, "⇒  E = ρₗ / (2πε₀r) r̂   [V/m]", size=14, anchor="start", weight=600))
    o.append(svg_close())
    figure("gauss-cylinder", "".join(o),
           "<strong>Gauss's law with cylindrical symmetry.</strong> Choose a coaxial cylinder of radius r and length L. By symmetry <b>D</b> is radial and depends only on r, so the caps contribute nothing (field parallel to the caps) and the side contributes D·(2πrL). The L cancels against the enclosed charge ρₗL.")

# --------------------------------------------------------------------------- 7. Sheet pillbox + slab with E_x plot
def fig_sheet_slab():
    W, H = 640, 300
    o = [svg_open(W, H)]
    # ---- left: sheet at z=0 with a pillbox
    cx, cy = 150, 150
    o.append(line(30, cy, 270, cy, stroke="var(--hi)", w=4))
    for k in range(8):
        o.append(text(45+k*30, cy-8, "+", size=12, fill="var(--hi)"))
    o.append(text(272, cy+5, "ρₛ", size=14, anchor="start", fill="var(--hi)"))
    # pillbox (two caps area A)
    o.append(f'<rect x="{cx-45}" y="{cy-55}" width="90" height="110" fill="var(--accent)" fill-opacity="0.10" stroke="currentColor" stroke-width="1.5"/>')
    o.append(text(cx+52, cy-60, "cap area A", size=12, anchor="start"))
    o.append(arrow(cx-25, cy-55, cx-25, cy-90, stroke="currentColor", w=2)); o.append(text(cx-25, cy-96, "dS", size=12, math=True))
    o.append(arrow(cx-25, cy+55, cx-25, cy+90, stroke="currentColor", w=2))
    o.append(arrow(cx+15, cy-12, cx+15, cy-48, stroke="var(--accent)", w=2.4)); o.append(text(cx+32, cy-30, "E", size=15, math=True, fill="var(--accent)"))
    o.append(arrow(cx+15, cy+12, cx+15, cy+48, stroke="var(--accent)", w=2.4))
    o.append(text(cx, 275, "2·(D·A) = ρₛ A  ⇒  E = ẑ sgn(z) ρₛ/(2ε₀)", size=13))
    # ---- right: slab and E_x(x) plot
    ox, oy = 360, 200
    xw = 60  # half-width in px
    o.append(f'<rect x="{ox+120-xw}" y="{oy-150}" width="{2*xw}" height="150" fill="var(--hi)" fill-opacity="0.15" stroke="none"/>')
    o.append(line(ox+120-xw, oy-150, ox+120-xw, oy, dash="4 4", w=1.1, opacity=0.6))
    o.append(line(ox+120+xw, oy-150, ox+120+xw, oy, dash="4 4", w=1.1, opacity=0.6))
    o.append(text(ox+120-xw+6, oy-150+14, "ρ", size=13, anchor="start", fill="var(--hi)"))
    o.append(text(ox+120-xw, oy+16, "−W/2", size=12)); o.append(text(ox+120+xw, oy+16, "W/2", size=12))
    # axes
    o.append(line(ox, oy, ox+250, oy, w=1.3, arrow=True, small=True)); o.append(text(ox+258, oy+4, "x", size=15, math=True, anchor="start"))
    o.append(line(ox+120, oy+10, ox+120, oy-140, w=1.3, arrow=True, small=True)); o.append(text(ox+128, oy-128, "Eₓ", size=14, anchor="start"))
    # E_x plot: linear inside, flat outside
    a = 45
    pts = [(ox+10, oy+a), (ox+120-xw, oy+a), (ox+120+xw, oy-a), (ox+240, oy-a)]
    o.append(path("M" + " L".join(f"{x},{y}" for x, y in pts), stroke="var(--accent)", w=2.6))
    o.append(line(ox+120, oy-a, ox+120+xw, oy-a, dash="3 3", w=1, opacity=0.6))
    o.append(text(ox+126, oy-a-6, "ρW/(2ε₀)", size=12, anchor="start"))
    o.append(text(ox+126, oy+a+14, "−ρW/(2ε₀)", size=12, anchor="start"))
    o.append(text(ox+125, 275, "inside: Eₓ = ρx/ε₀ ; outside: ±ρW/(2ε₀)", size=13))
    o.append(svg_close())
    figure("sheet-slab", "".join(o),
           "<strong>Planar symmetry: sheet and slab.</strong> Left: a pillbox straddling a charged sheet; both caps carry equal flux D·A and the sides carry none, so 2DA = ρₛA. Right: for a uniform slab the same pillbox argument gives a field that grows linearly inside (more enclosed charge as the caps move out) and saturates outside at the sheet value with ρₛ = ρW.")

# --------------------------------------------------------------------------- 8. Half-space flux argument
def fig_flux_plane():
    W, H = 640, 260
    o = [svg_open(W, H)]
    cx, cy = 200, 150
    # plane (perspective parallelogram)
    o.append(f'<polygon points="{cx-150},{cy+30} {cx+110},{cy+30} {cx+150},{cy-10} {cx-110},{cy-10}" fill="var(--accent2)" fill-opacity="0.18" stroke="currentColor" stroke-width="1.4"/>')
    o.append(text(cx-146, cy+44, "z = 0 plane", size=12, anchor="start"))
    o.append(arrow(cx-120, cy+10, cx-120, cy-40, stroke="currentColor", w=2)); o.append(text(cx-120, cy-48, "n̂ = ẑ", size=13))
    # charge above, radial lines
    qx, qy = cx, cy-80
    for k in range(12):
        th = 2*math.pi*k/12 + math.pi/12
        ex, ey = math.cos(th), math.sin(th)
        L = 62
        col = "var(--hi)" if ey > 0.05 else "var(--accent)"
        o.append(arrow(qx+14*ex, qy+14*ey, qx+L*ex, qy+L*ey, stroke=col, w=1.8, small=True, opacity=0.9))
    o.append(charge(qx, qy, "+"))
    o.append(text(qx+22, qy-12, "+Q", size=14, anchor="start", weight=600))
    o.append(text(cx+20, cy+58, "half the lines go up (never cross the plane),", size=12))
    o.append(text(cx+20, cy+74, "half go down through it, against n̂", size=12))
    o.append(text(cx+20, cy+94, "⇒  ψ(through plane) = −Q/2", size=13, weight=600))
    # right: bookkeeping with two charges
    X = 400
    o.append(text(X, 60, "Two charges: add the fluxes", size=13, anchor="start", weight=600))
    o.append(text(X, 88, "+3Q at z = +h:  −3Q/2", size=13, anchor="start"))
    o.append(text(X, 106, "(half its lines go down, against n̂)", size=11, anchor="start", fill="var(--muted)"))
    o.append(text(X, 130, "−Q at z = −h:   −Q/2", size=13, anchor="start"))
    o.append(text(X, 148, "(lines come down into it)", size=11, anchor="start", fill="var(--muted)"))
    o.append(text(X, 176, "total  ψ = −2Q", size=14, anchor="start", weight=600))
    o.append(text(X, 206, "Flip n̂ → every sign flips.", size=12, anchor="start"))
    o.append(text(X, 224, "A closed surface would count", size=12, anchor="start"))
    o.append(text(X, 240, "only the charge inside it.", size=12, anchor="start"))
    o.append(svg_close())
    figure("flux-plane", "".join(o),
           "<strong>The half-space argument.</strong> A point charge sends its field lines out uniformly in all directions, so exactly half of them cross any infinite plane that does not contain the charge: the flux through the plane has magnitude Q/2. Whether it counts as + or − depends only on which way you chose the normal. With several charges, add the contributions (superposition).")

# --------------------------------------------------------------------------- 9. Curl vs divergence twin panels
def fig_curl_div_panels():
    W, H = 640, 300
    o = [svg_open(W, H)]
    def panel(x0, y0, title, mode):
        o.append(f'<rect x="{x0}" y="{y0}" width="280" height="110" rx="6" fill="var(--muted)" fill-opacity="0.12" stroke="currentColor" stroke-width="1.2"/>')
        o.append(text(x0+140, y0-8, title, size=13, weight=600))
        for i in range(6):
            for j in range(2):
                X = x0 + 30 + i*44; Y = y0 + 35 + j*40
                L = 8 + 5*i
                if mode == "x":
                    o.append(arrow(X-L/2, Y, X+L/2, Y, stroke="var(--accent)", w=2.2, small=True))
                else:
                    o.append(arrow(X, Y+L/2, X, Y-L/2, stroke="var(--accent)", w=2.2, small=True))
        o.append(circle(x0+140, y0+55, 3, fill="currentColor"))
    panel(40, 40, "E = Eₓ(x) x̂ — arrows along x, growing along x", "x")
    panel(40, 185, "E = Eᵧ(x) ŷ — arrows along y, growing along x", "y")
    # annotations
    X = 340
    o.append(text(X, 70, "variation ALONG the flow", size=13, anchor="start"))
    o.append(text(X, 92, "→ ∇·E ≠ 0  (a source)", size=14, anchor="start", weight=600))
    o.append(text(X, 118, "no variation ACROSS the flow", size=13, anchor="start"))
    o.append(text(X, 140, "→ ∇×E = 0", size=14, anchor="start", weight=600))
    o.append(text(X, 215, "no variation ALONG the flow", size=13, anchor="start"))
    o.append(text(X, 237, "→ ∇·E = 0", size=14, anchor="start", weight=600))
    o.append(text(X, 263, "variation ACROSS the flow", size=13, anchor="start"))
    o.append(text(X, 285, "→ ∇×E ≠ 0  (a paddlewheel spins)", size=14, anchor="start", weight=600))
    o.append(svg_close())
    figure("curl-div-panels", "".join(o),
           "<strong>The instructor's rule of thumb.</strong> Divergence measures how a field changes <em>along</em> its own direction (a scalar: net outflow); curl measures how it changes <em>across</em> its direction (a vector: net circulation). The same two pictures answer both questions — with opposite answers.")

# --------------------------------------------------------------------------- 10. Paddlewheel in a river
def fig_paddlewheel():
    W, H = 640, 262
    o = [svg_open(W, H)]
    def river(x0, label, wheel_x, spin):
        # banks
        o.append(line(x0, 30, x0, 220, w=2))
        o.append(line(x0+150, 30, x0+150, 220, w=2))
        # velocity profile: parabolic
        for k in range(7):
            X = x0 + 12 + k*21
            u = (X - x0)/150.0
            L = 90*4*u*(1-u) + 6
            o.append(arrow(X, 200, X, 200-L, stroke="var(--accent)", w=1.8, small=True, opacity=0.8))
        # wheel
        wx, wy = wheel_x, 120
        o.append(circle(wx, wy, 16, fill="var(--paper)", stroke="currentColor", w=1.8))
        for a in (0, 45, 90, 135):
            r = math.radians(a)
            o.append(line(wx-16*math.cos(r), wy-16*math.sin(r), wx+16*math.cos(r), wy+16*math.sin(r), w=1.6))
        if spin:
            sweep = 1 if spin == "cw" else 0
            o.append(path(f"M{wx-26},{wy-4} A26,26 0 0 {sweep} {wx+26},{wy-4}", stroke="var(--hi)", w=2.2, arrow=True))
        l1, l2 = label
        o.append(text(x0+75, 234, l1, size=12)); o.append(text(x0+75, 248, l2, size=12, weight=600))
    river(40, ("mid-stream: same speed both sides", "no spin: curl = 0"), 115, None)
    river(245, ("near left bank: right side faster", "spins CCW: curl ≠ 0"), 285, "ccw")
    river(450, ("near right bank: left side faster", "spins CW: curl ≠ 0"), 560, "cw")
    o.append(svg_close())
    figure("paddlewheel", "".join(o),
           "<strong>Curl as a paddlewheel test.</strong> Drop a tiny paddlewheel into the flow. For a straight flow it spins where the speed differs from one side of the wheel to the other — where the field varies <em>across</em> its own direction (curved flow adds a geometric part: see the text). The spin axis (right-hand rule) is the direction of ∇×<b>v</b>; a uniform stream, however fast, has zero curl.")

# --------------------------------------------------------------------------- 11. Stokes tiling
def fig_stokes_tiling():
    W, H = 640, 220
    o = [svg_open(W, H)]
    def blob(cx, cy, s=1.0):
        return (f"M{cx-70*s},{cy} C{cx-70*s},{cy-60*s} {cx-20*s},{cy-70*s} {cx+10*s},{cy-55*s} "
                f"S{cx+75*s},{cy-40*s} {cx+70*s},{cy+5*s} S{cx+20*s},{cy+70*s} {cx-20*s},{cy+60*s} S{cx-70*s},{cy+50*s} {cx-70*s},{cy}")
    # panel 1: boundary circulation
    cx, cy = 110, 110
    o.append(f'<path d="{blob(cx, cy)}" fill="var(--accent)" fill-opacity="0.10" stroke="currentColor" stroke-width="2"/>')
    o.append(path(f"M{cx-40},{cy+38} A45,45 0 1 1 {cx+42},{cy+30}", stroke="var(--hi)", w=2.4, arrow=True))
    o.append(text(cx, cy+8, "C", size=16, math=True))
    o.append(text(cx, 200, "∮ v·dl around C", size=14))
    o.append(text(230, 115, "=", size=26))
    # panel 2: tiled with small loops; interior edges cancel
    cx, cy = 330, 110
    o.append(f'<path d="{blob(cx, cy)}" fill="var(--accent)" fill-opacity="0.10" stroke="currentColor" stroke-width="2"/>')
    for i in range(-2, 3):
        for j in range(-2, 3):
            X, Y = cx + i*28, cy + j*26
            if math.hypot(i*28, j*26) > 62: continue
            r = 8; a = math.radians(70)
            o.append(path(f"M{X+r},{Y} A{r},{r} 0 1 1 {X+r*math.cos(a):.1f},{Y+r*math.sin(a):.1f}", stroke="var(--hi)", w=1.4, arrow=True, opacity=0.9))
    o.append(text(cx, 200, "Σ small circulations", size=14))
    o.append(text(455, 115, "=", size=26))
    # panel 3: curl dots
    cx, cy = 550, 110
    o.append(f'<path d="{blob(cx, cy)}" fill="var(--accent)" fill-opacity="0.10" stroke="currentColor" stroke-width="2"/>')
    for i in range(-2, 3):
        for j in range(-2, 3):
            X, Y = cx + i*28, cy + j*26
            if math.hypot(i*28, j*26) > 62: continue
            o.append(circle(X, Y, 6, fill="none", stroke="var(--hi)", w=1.4))
            o.append(circle(X, Y, 1.8, fill="var(--hi)", stroke="none"))
    o.append(text(cx, 200, "∬ (∇×v)·dS over S", size=14))
    o.append(svg_close())
    figure("stokes-tiling", "".join(o),
           "<strong>Stokes' theorem by tiling.</strong> Cut the surface into small cells. Adjacent cells share an edge traversed in opposite directions, so all interior contributions cancel and only the outer boundary survives. Each cell's circulation is (∇×<b>v</b>)·d<b>S</b>: curl is <em>circulation per unit area</em>. The divergence theorem is the same idea with flux through small volumes.")

# --------------------------------------------------------------------------- 12. pn-junction by superposition of slabs
def fig_pn_superposition():
    W, H = 640, 400
    o = [svg_open(W, H)]
    ox, oy = 90, 120            # origin of the E plots
    sx = 70                     # px per unit width (W1 = W2 = 1 unit)
    def X(x): return ox + 160 + x*sx
    # ---- top: charge profile
    y0 = 55
    o.append(line(60, y0, 560, y0, w=1.2, arrow=True, small=True)); o.append(text(568, y0+4, "x", size=14, math=True, anchor="start"))
    o.append(f'<rect x="{X(-1)}" y="{y0}" width="{sx}" height="24" fill="var(--accent2)" fill-opacity="0.35" stroke="currentColor" stroke-width="1"/>')
    o.append(f'<rect x="{X(0)}" y="{y0-24}" width="{sx}" height="24" fill="var(--hi)" fill-opacity="0.35" stroke="currentColor" stroke-width="1"/>')
    o.append(text(X(-0.5), y0+38, "−ρ₁ (acceptors)", size=12)); o.append(text(X(0.5), y0-30, "+ρ₂ (donors)", size=12))
    o.append(text(X(-1)-4, y0-6, "−W₁", size=11, anchor="end")); o.append(text(X(1)+4, y0-6, "W₂", size=11, anchor="start"))
    o.append(text(60, y0-8, "ρ(x)", size=13, anchor="start"))
    o.append(text(560, y0-30, "neutral: ρ₁W₁ = ρ₂W₂", size=12, anchor="end"))
    # ---- middle: each slab's field (shifted slab result), amplitude a px
    a = 34
    ym = 185
    o.append(line(60, ym, 560, ym, w=1.2, arrow=True, small=True)); o.append(text(60, ym-a-22, "Eₓ of each slab", size=13, anchor="start"))
    # slab 1 (−ρ1 on (−W1,0)): plateau +ρ1W1/2ε0 for x<−W1, ramp down to −ρ1W1/2ε0 at 0, flat after
    p1 = [(X(-2.3), ym-a), (X(-1), ym-a), (X(0), ym+a), (X(2.3), ym+a)]
    o.append(path("M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in p1), stroke="var(--accent2)", w=2.4))
    # slab 2 (+ρ2 on (0,W2)): plateau −ρ2W2/2ε0 for x<0, ramp up to +ρ2W2/2ε0 at W2, flat after
    p2 = [(X(-2.3), ym+a), (X(0), ym+a), (X(1), ym-a), (X(2.3), ym-a)]
    o.append(path("M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in p2), stroke="var(--hi)", w=2.4, dash="6 4"))
    o.append(text(X(-2.2), ym-a-6, "+ρ₁W₁/2ε₀", size=11, anchor="start", fill="var(--accent2)"))
    o.append(text(X(1.3), ym-a-6, "+ρ₂W₂/2ε₀", size=11, anchor="start", fill="var(--hi)"))
    # ---- bottom: total
    yb = 300
    o.append(line(60, yb, 560, yb, w=1.2, arrow=True, small=True)); o.append(text(60, yb-14, "total Eₓ", size=13, anchor="start"))
    pt = [(X(-2.3), yb), (X(-1), yb), (X(0), yb+2*a), (X(1), yb), (X(2.3), yb)]
    o.append(path("M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pt), stroke="var(--accent)", w=2.8))
    o.append(text(X(0)+6, yb+2*a+4, "−ρ₁W₁/ε₀ = −ρ₂W₂/ε₀", size=12, anchor="start"))
    o.append(text(X(-1.6), yb+18, "E = 0 outside", size=12)); o.append(text(X(1.7), yb+18, "E = 0 outside", size=12))
    o.append(svg_close())
    figure("pn-superposition", "".join(o),
           "<strong>A pn junction's depletion field, by superposition of two slabs.</strong> Each charged layer is a slab whose field you already know (linear inside, constant outside). Because the layers carry equal and opposite charge per area, the two outside plateaus cancel and the two inside ramps add: the field is a triangle, strongest at the metallurgical junction and zero outside the depletion region.")

# =========================================================================== Lectures 5–10
def _tint(x, y, w, h, color, op=0.12, rx=0):
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{color}" fill-opacity="{op}" stroke="none"/>'

# --------------------------------------------------------------------------- 13. Equipotentials and the gradient
def fig_equipotentials_gradient():
    W, H = 640, 360
    o = [svg_open(W, H)]
    # ---- left: V = x^2 - 6y on [-5,5]^2, contours are parabolas y = (x^2 - c)/6
    x0, y0, S = 30, 20, 290                    # panel origin and size (px)
    def M(x, y): return x0 + (x+5)/10*S, y0 + (5-y)/10*S
    o.append(f'<rect x="{x0}" y="{y0}" width="{S}" height="{S}" fill="var(--muted)" fill-opacity="0.08" stroke="currentColor" stroke-width="1.2"/>')
    for c in range(-30, 51, 10):
        pts = []
        for x in np.linspace(-5, 5, 81):
            y = (x*x - c)/6
            if -5 <= y <= 5: pts.append(M(x, y))
            else:
                if len(pts) > 1:
                    o.append(path("M" + " L".join(f"{X:.1f},{Y:.1f}" for X, Y in pts), stroke="var(--accent2)", w=1.4, opacity=0.9)); pts = []
                elif pts: pts = []
        if len(pts) > 1:
            o.append(path("M" + " L".join(f"{X:.1f},{Y:.1f}" for X, Y in pts), stroke="var(--accent2)", w=1.4, opacity=0.9))
    # contour labels along x = 0 (y = -c/6) where they fit
    for c in (-30, -10, 10, 30):
        X, Y = M(0.0, -c/6.0)
        o.append(text(X+4, Y-3, f"{c} V", size=10, anchor="start", fill="var(--accent2)"))
    # E = -grad V = (-2x, 6): unit arrows on a coarse grid
    for gx in np.linspace(-4, 4, 5):
        for gy in np.linspace(-4, 4, 5):
            ex, ey = -2*gx, 6.0; n = math.hypot(ex, ey)
            X, Y = M(gx, gy); L = 16
            o.append(line(X, Y, X+L*ex/n, Y-L*ey/n, stroke="var(--accent)", w=1.6, arrow=True, small=True, opacity=0.9))
    o.append(text(x0+S+12, y0+S+4, "x", size=14, math=True, anchor="start")); o.append(text(x0-12, y0+S/2, "y", size=14, math=True))
    o.append(text(W/2, H-8, "V = x² − 6y :  E = −∇V = −2x x̂ + 6 ŷ,  everywhere perpendicular to the equipotentials", size=12))
    # ---- right: crowded vs sparse equipotentials
    rx0, ry0 = 350, 40
    # (a) uniform field: equally spaced vertical lines
    for k in range(6):
        X = rx0 + 10 + k*24
        o.append(line(X, ry0+8, X, ry0+92, stroke="var(--accent2)", w=1.4))
    for Y in (ry0+30, ry0+70):
        o.append(arrow(rx0+2, Y, rx0+140, Y, stroke="var(--accent)", w=2))
    o.append(text(rx0+70, ry0+112, "uniform: even spacing", size=11))
    # (b) diverging field: arcs whose spacing grows to the right
    cx, cy = rx0+156, ry0+52
    for k, r in enumerate((14, 26, 42, 62, 86)):
        d = f"M{cx+r*math.cos(-0.9):.1f},{cy+r*math.sin(-0.9):.1f} A{r},{r} 0 0 1 {cx+r*math.cos(0.9):.1f},{cy+r*math.sin(0.9):.1f}"
        o.append(path(d, stroke="var(--accent2)", w=1.4))
    for th in (-0.55, 0.0, 0.55):
        o.append(arrow(cx+12*math.cos(th), cy+12*math.sin(th), cx+92*math.cos(th), cy+92*math.sin(th), stroke="var(--accent)", w=2))
    o.append(text(rx0+215, ry0+112, "spreading: gaps widen", size=11))
    o.append(text(rx0+140, ry0+132, "(equipotentials drawn 10 V apart)", size=11, fill="var(--muted)"))
    o.append(text(rx0+140, ry0+152, "|E| = |dV/dl|: crowded lines = strong field", size=12, weight=600))
    o.append(text(rx0+140, ry0+172, "E points downhill (toward lower V),", size=12))
    o.append(text(rx0+140, ry0+190, "∇V points uphill; they are antiparallel.", size=12))
    o.append(text(rx0+140, ry0+228, "Contour-map analogy: V is altitude,", size=12, fill="var(--muted)"))
    o.append(text(rx0+140, ry0+246, "E is the steepest way down.", size=12, fill="var(--muted)"))
    o.append(svg_close())
    figure("equipotentials-gradient", "".join(o),
           "<strong>Potential landscape.</strong> Left: equipotentials of V = x² − 6y (parabolas) with the direction of <b>E</b> = −∇V at sample points. The arrows cross every contour at right angles and point toward lower V. Right: equal potential steps drawn as contours — where the contours crowd together the field is strong, where they spread it is weak, because |E| is the potential drop per unit distance.")

# --------------------------------------------------------------------------- 14. Path independence on the unit square
def fig_path_independence():
    W, H = 640, 320
    o = [svg_open(W, H)]
    def panel(px, sign, title, res_u, res_l, verdict):
        S = 200; x0, y0 = px, 60
        def M(x, y): return x0 + x*S, y0 + (1-y)*S
        o.append(f'<rect x="{x0}" y="{y0}" width="{S}" height="{S}" fill="var(--muted)" fill-opacity="0.08" stroke="currentColor" stroke-width="1.2"/>')
        o.append(text(x0+S/2, y0-30, title, size=14, weight=600))
        for gx in np.linspace(0.1, 0.9, 5):
            for gy in np.linspace(0.1, 0.9, 5):
                ex, ey = gy, sign*gx; n = math.hypot(ex, ey)
                if n < 1e-9: continue
                L = 8 + 18*n
                X, Y = M(gx, gy)
                o.append(line(X-L*ex/n/2, Y+L*ey/n/2, X+L*ex/n/2, Y-L*ey/n/2, stroke="currentColor", w=1.3, arrow=True, small=True, opacity=0.55))
        # paths: C_u = (0,0)->(0,1)->(1,1)  ;  C_l = (0,0)->(1,0)->(1,1)
        o.append(path(f"M{M(0,0)[0]},{M(0,0)[1]} L{M(0,1)[0]},{M(0,1)[1]} L{M(1,1)[0]},{M(1,1)[1]}", stroke="var(--hi)", w=3, arrow=True))
        o.append(path(f"M{M(0,0)[0]},{M(0,0)[1]} L{M(1,0)[0]},{M(1,0)[1]} L{M(1,1)[0]},{M(1,1)[1]}", stroke="var(--accent2)", w=3, arrow=True))
        o.append(circle(*M(0,0), 4, fill="currentColor")); o.append(circle(*M(1,1), 4, fill="currentColor"))
        o.append(text(M(0,0)[0]-10, M(0,0)[1]+16, "o", size=14, math=True))
        o.append(text(M(1,1)[0]+10, M(1,1)[1]-8, "p", size=14, math=True))
        o.append(text(M(0.5,1)[0], M(0.5,1)[1]-8, "Cᵤ", size=13, fill="var(--hi)", weight=600))
        o.append(text(M(1,0.5)[0]+14, M(1,0.5)[1]+4, "Cₗ", size=13, anchor="start", fill="var(--accent2)", weight=600))
        o.append(text(x0+S/2, y0+S+22, f"along Cᵤ: {res_u}    along Cₗ: {res_l}", size=13))
        o.append(text(x0+S/2, y0+S+42, verdict, size=13, weight=600))
    panel(40, +1, "E = y x̂ + x ŷ   (∇×E = 0)", "1", "1", "same answer: path independent")
    panel(380, -1, "E = y x̂ − x ŷ   (∇×E = −2ẑ)", "+1", "−1", "different: no potential exists")
    o.append(svg_close())
    figure("path-independence", "".join(o),
           "<strong>Same endpoints, two paths.</strong> For the curl-free field (left) the line integral ∫<sub>o</sub><sup>p</sup> <b>E</b>·d<b>l</b> is 1 along either path, so V(o) − V(p) is well defined. For the swirling field (right) the two paths give +1 and −1; their difference, 2, is the flux of ∇×<b>E</b> = −2ẑ through the square taken with the loop's (clockwise, −ẑ) orientation — Stokes' theorem in action.")

# --------------------------------------------------------------------------- 15. Boundary conditions: pillbox and loop
def fig_bc_pillbox_loop():
    W, H = 640, 330
    o = [svg_open(W, H)]
    yI = 175   # interface height
    o.append(_tint(20, 40, 600, yI-40, "var(--accent)", 0.07))
    o.append(_tint(20, yI, 600, 290-yI, "var(--muted)", 0.12))
    o.append(line(20, yI, 620, yI, w=2))
    o.append(text(30, 58, "medium 1  (ε₁, σ₁, μ₁)", size=13, anchor="start"))
    o.append(text(30, 282, "medium 2  (ε₂, σ₂, μ₂)", size=13, anchor="start"))
    # normal
    o.append(arrow(320, yI, 320, yI-60, stroke="currentColor", w=2.4)); o.append(text(330, yI-44, "n̂", size=16, anchor="start", weight=600))
    o.append(text(330, yI-26, "from 2 into 1", size=11, anchor="start", fill="var(--muted)"))
    # ---- pillbox (left)
    cx = 150; hw, hh = 70, 26
    o.append(f'<rect x="{cx-hw}" y="{yI-hh}" width="{2*hw}" height="{2*hh}" fill="var(--accent)" fill-opacity="0.10" stroke="currentColor" stroke-width="1.5" stroke-dasharray="5 3"/>')
    for k in range(6): o.append(text(cx-55+k*22, yI-5, "+", size=13, fill="var(--hi)"))
    o.append(text(cx-hw-6, yI+4, "ρₛ", size=14, anchor="end", fill="var(--hi)"))
    o.append(arrow(cx-40, yI-hh, cx-40, yI-hh-28, stroke="currentColor", w=1.8)); o.append(text(cx-40, yI-hh-34, "dS", size=12, math=True))
    o.append(arrow(cx-40, yI+hh, cx-40, yI+hh+28, stroke="currentColor", w=1.8)); o.append(text(cx-40, yI+hh+42, "dS", size=12, math=True))
    o.append(text(cx-hw+4, yI-hh-6, "cap area A", size=11, anchor="start"))
    o.append(line(cx+hw+18, yI-hh, cx+hw+18, yI+hh, w=1)); o.append(text(cx+hw+24, yI+4-14, "h→0", size=11, anchor="start"))
    # D1 above (steep, long), D2 below (shallower, shorter)
    o.append(arrow(cx+20, yI-8, cx+50, yI-70, stroke="var(--accent)", w=2.6)); o.append(text(cx+56, yI-62, "D₁", size=15, anchor="start", fill="var(--accent)", weight=600))
    o.append(arrow(cx+20, yI+70, cx+44, yI+12, stroke="var(--accent)", w=2.6, opacity=0.7)); o.append(text(cx+50, yI+50, "D₂", size=15, anchor="start", fill="var(--accent)", weight=600))
    o.append(text(cx, 312, "D₁ₙA − D₂ₙA = ρₛA  ⇒  n̂·(D₁ − D₂) = ρₛ", size=13, weight=600))
    # ---- loop (right)
    lx = 470; L, hh2 = 80, 22
    o.append(f'<rect x="{lx-L/2}" y="{yI-hh2}" width="{L}" height="{2*hh2}" fill="var(--accent2)" fill-opacity="0.10" stroke="currentColor" stroke-width="1.5" stroke-dasharray="5 3"/>')
    o.append(arrow(lx-L/2+10, yI-hh2, lx+L/2-10, yI-hh2, stroke="currentColor", w=1.8)); o.append(text(lx, yI-hh2-8, "dl  (length L)", size=11, math=False))
    o.append(arrow(lx+L/2-10, yI+hh2, lx-L/2+10, yI+hh2, stroke="currentColor", w=1.8)); o.append(text(lx, yI+hh2+16, "dl", size=12, math=True))
    o.append(line(lx+L/2+14, yI-hh2, lx+L/2+14, yI+hh2, w=1)); o.append(text(lx+L/2+20, yI+4-14, "h→0", size=11, anchor="start"))
    # E1 and E2 tangential-ish arrows (same tangential component)
    o.append(arrow(lx-60, yI-75, lx-10, yI-45, stroke="var(--accent2)", w=2.6)); o.append(text(lx-4, yI-48, "E₁", size=15, anchor="start", fill="var(--accent2)", weight=600))
    o.append(arrow(lx-60, yI+45, lx-10, yI+62, stroke="var(--accent2)", w=2.6, opacity=0.7)); o.append(text(lx-4, yI+66, "E₂", size=15, anchor="start", fill="var(--accent2)", weight=600))
    o.append(text(lx, 312, "E₁ₜL − E₂ₜL = 0  ⇒  n̂×(E₁ − E₂) = 0", size=13, weight=600))
    o.append(svg_close())
    figure("bc-pillbox-loop", "".join(o),
           "<strong>Two integral laws, squeezed onto an interface.</strong> Left: Gauss's law on a pillbox whose height shrinks to zero — the side wall contributes no flux, the two caps give D₁ₙA − D₂ₙA, and the enclosed charge is the surface charge ρₛA. Right: ∮<b>E</b>·d<b>l</b> = 0 around a loop whose short sides shrink — only the two long edges survive, traversed in opposite directions, so the tangential components must match. The unit normal n̂ always points from medium 2 into medium 1.")

# --------------------------------------------------------------------------- 16. "Remember this drawing": the four boundary conditions
def fig_bc_summary_panels():
    W, H = 640, 300
    o = [svg_open(W, H)]
    yI = 150; x0 = 60; cw = 140
    o.append(_tint(x0, 50, 4*cw, yI-50, "var(--accent)", 0.07))
    o.append(_tint(x0, yI, 4*cw, 240-yI, "var(--muted)", 0.12))
    o.append(line(x0, yI, x0+4*cw, yI, w=2))
    for k in range(1, 4): o.append(line(x0+k*cw, 50, x0+k*cw, 240, w=1, dash="3 3", opacity=0.5))
    o.append(arrow(30, yI, 30, yI-50, stroke="currentColor", w=2.2)); o.append(text(30, yI-58, "n̂", size=15, weight=600))
    o.append(text(30, 70, "1", size=13, weight=600)); o.append(text(30, 232, "2", size=13, weight=600))
    titles = ["tangential E", "tangential H", "normal D", "normal B"]
    for k, t in enumerate(titles): o.append(text(x0+k*cw+cw/2, 40, t, size=13, weight=600))
    # col 1: E_t equal
    c = x0+cw/2
    o.append(arrow(c-40, yI-40, c+40, yI-40, stroke="var(--accent2)", w=2.6)); o.append(text(c+46, yI-36, "E₁ₜ", size=12, anchor="start"))
    o.append(arrow(c-40, yI+40, c+40, yI+40, stroke="var(--accent2)", w=2.6, opacity=0.7)); o.append(text(c+46, yI+44, "E₂ₜ", size=12, anchor="start"))
    o.append(text(c, 262, "E₁ₜ = E₂ₜ", size=13)); o.append(text(c, 282, "n̂×(E₁−E₂) = 0", size=12, fill="var(--muted)"))
    # col 2: H_t jumps by J_s (⊗ dots on the interface)
    c = x0+cw+cw/2
    o.append(arrow(c-50, yI-40, c+50, yI-40, stroke="var(--accent2)", w=2.6)); o.append(text(c+56, yI-36, "H₁ₜ", size=12, anchor="start"))
    o.append(arrow(c-25, yI+40, c+25, yI+40, stroke="var(--accent2)", w=2.6, opacity=0.7)); o.append(text(c+31, yI+44, "H₂ₜ", size=12, anchor="start"))
    for dx in (-36, -12, 12, 36):
        o.append(circle(c+dx, yI, 5, fill="var(--paper)", stroke="var(--hi)", w=1.5))
        o.append(line(c+dx-3, yI-3, c+dx+3, yI+3, stroke="var(--hi)", w=1.3)); o.append(line(c+dx-3, yI+3, c+dx+3, yI-3, stroke="var(--hi)", w=1.3))
    o.append(text(c+50, yI+4, "Jₛ", size=12, anchor="start", fill="var(--hi)"))
    o.append(text(c, 262, "H₁ₜ − H₂ₜ = Jₛ", size=13)); o.append(text(c, 282, "n̂×(H₁−H₂) = Jₛ", size=12, fill="var(--muted)"))
    # col 3: D_n jumps by rho_s
    c = x0+2*cw+cw/2
    for dx in (-40, -20, 0, 20, 40): o.append(text(c+dx, yI-4, "+", size=13, fill="var(--hi)"))
    o.append(text(c+52, yI+4, "ρₛ", size=12, anchor="start", fill="var(--hi)"))
    o.append(arrow(c-15, yI-10, c-15, yI-85, stroke="var(--accent)", w=2.6)); o.append(text(c-9, yI-70, "D₁ₙ", size=12, anchor="start"))
    o.append(arrow(c-15, yI+60, c-15, yI+15, stroke="var(--accent)", w=2.6, opacity=0.7)); o.append(text(c-9, yI+50, "D₂ₙ", size=12, anchor="start"))
    o.append(text(c, 262, "D₁ₙ − D₂ₙ = ρₛ", size=13)); o.append(text(c, 282, "n̂·(D₁−D₂) = ρₛ", size=12, fill="var(--muted)"))
    # col 4: B_n equal
    c = x0+3*cw+cw/2
    o.append(arrow(c, yI-10, c, yI-70, stroke="var(--accent)", w=2.6)); o.append(text(c+6, yI-56, "B₁ₙ", size=12, anchor="start"))
    o.append(arrow(c, yI+70, c, yI+10, stroke="var(--accent)", w=2.6, opacity=0.7)); o.append(text(c+6, yI+56, "B₂ₙ", size=12, anchor="start"))
    o.append(text(c, 262, "B₁ₙ = B₂ₙ", size=13)); o.append(text(c, 282, "n̂·(B₁−B₂) = 0", size=12, fill="var(--muted)"))
    o.append(svg_close())
    figure("bc-summary-panels", "".join(o),
           "<strong>The four boundary conditions in one drawing</strong> (the instructor's \"remember this drawing!!\"). Tangential <b>E</b> and normal <b>B</b> are continuous; tangential <b>H</b> jumps by the surface current density Jₛ [A/m] and normal <b>D</b> jumps by the surface charge density ρₛ [C/m²]. Each one is a Maxwell equation with ∇ → n̂, fields → (field₁ − field₂), ρ → ρₛ, <b>J</b> → <b>J</b>ₛ, and ∂/∂t → 0 (static, or the h→0 limit kills the flux terms).")

# --------------------------------------------------------------------------- 17. Laplace between plates
def fig_laplace_plates():
    W, H = 640, 316
    o = [svg_open(W, H)]
    # left: side view of the plates
    x0, yb, yt = 40, 230, 70   # bottom plate y, top plate y
    o.append(_tint(x0, yt, 260, yb-yt, "var(--accent)", 0.06))
    o.append(line(x0, yb, x0+260, yb, stroke="currentColor", w=6, opacity=0.85)); o.append(line(x0, yt, x0+260, yt, stroke="currentColor", w=6, opacity=0.85))
    for k in range(8): o.append(text(x0+22+k*32, yb+18, "+", size=13, fill="var(--hi)"))
    for k in range(8): o.append(text(x0+22+k*32, yt-8, "−", size=14, fill="var(--accent)"))
    o.append(text(x0+265, yb+4, "z = 0,  V = 0", size=12, anchor="start")); o.append(text(x0+265, yt+4, "z = d = 2 m,  V = −3 V", size=12, anchor="start"))
    for X in (x0+60, x0+130, x0+200):
        o.append(arrow(X, yb-14, X, yt+14, stroke="var(--accent)", w=2.4))
    o.append(text(x0+130+14, (yb+yt)/2, "E = +1.5 ẑ V/m", size=12, anchor="start"))
    o.append(text(x0+130, yt-30, "E = 0 above", size=11, fill="var(--muted)")); o.append(text(x0+130, yb+38, "E = 0 below", size=11, fill="var(--muted)"))
    # normals into the gap
    o.append(arrow(x0+20, yb-2, x0+20, yb-30, stroke="currentColor", w=1.6)); o.append(text(x0+8, yb-20, "n̂", size=12, anchor="end"))
    o.append(arrow(x0+20, yt+2, x0+20, yt+30, stroke="currentColor", w=1.6)); o.append(text(x0+8, yt+22, "n̂", size=12, anchor="end"))
    o.append(text(x0+130, 300, "ρₛ(0) = +1.5 ε₀ ,   ρₛ(d) = −1.5 ε₀  [C/m²]", size=12))
    # right: V(z) plot (z horizontal)
    px, py = 380, 200
    o.append(line(px, py, px+230, py, w=1.3, arrow=True, small=True)); o.append(text(px+238, py+4, "z", size=14, math=True, anchor="start"))
    o.append(line(px, py+70, px, py-90, w=1.3, arrow=True, small=True)); o.append(text(px-6, py-80, "V", size=14, math=True, anchor="end"))
    o.append(line(px, py, px+180, py+60, stroke="var(--accent2)", w=2.6))
    o.append(line(px+180, py, px+180, py+60, dash="3 3", w=1, opacity=0.6)); o.append(text(px+186, py-8, "d = 2 m", size=11, anchor="start"))
    o.append(line(px, py+60, px+180, py+60, dash="3 3", w=1, opacity=0.6)); o.append(text(px-6, py+64, "−3 V", size=11, anchor="end"))
    o.append(text(px+128, py-4, "V(z) = −1.5 z", size=13, fill="var(--accent2)", weight=600))
    o.append(text(px+115, 300, "slope dV/dz = −1.5 V/m  ⇒  E_z = +1.5 V/m", size=12))
    o.append(svg_close())
    figure("laplace-plates", "".join(o),
           "<strong>Laplace's equation between two plates.</strong> With nothing but the plates to depend on, V = V(z), so ∇²V = d²V/dz² = 0 forces a straight line V = Az + B; the plate potentials fix A and B. The field is minus the slope, and points from the 0 V plate toward the −3 V plate. The plate charges follow from the boundary condition ρₛ = n̂·<b>D</b> with n̂ pointing out of each conductor into the gap.")

# --------------------------------------------------------------------------- 18. pn junction: E and V
def fig_pn_potential():
    W, H = 640, 330
    o = [svg_open(W, H)]
    W1, W2, r1, r2 = 1.0, 2.0, 2.0, 1.0        # rho1*W1 = rho2*W2 (neutral); eps0 = 1 units
    xs = np.linspace(-2.0, 3.0, 301)
    def Ex(x):
        if -W1 < x < 0: return -r1*(x+W1)
        if 0 <= x < W2: return r2*(x-W2)
        return 0.0
    def V(x):
        if x <= -W1: return -r1*W1**2/2
        if x < 0: return r1*(x+W1)**2/2 - r1*W1**2/2
        if x < W2: return -r2*(x-W2)**2/2 + r2*W2**2/2
        return r2*W2**2/2
    def axes(py, ylab, sy, y_lo, y_hi):
        px = 60
        def M(x, y): return px + (x+2)/5*520, py - y*sy
        o.append(line(px, py, px+530, py, w=1.2, arrow=True, small=True)); o.append(text(px+540, py+4, "x", size=14, math=True, anchor="start"))
        o.append(line(px+208, py-y_hi*sy-8, px+208, py-y_lo*sy+8, w=1.0, opacity=0.5))
        o.append(text(px+200, py-y_hi*sy-12, ylab, size=13, anchor="end"))
        return M
    # shaded slabs (shared)
    M0 = lambda x: 60 + (x+2)/5*520
    o.append(_tint(M0(-W1), 30, M0(0)-M0(-W1), 270, "var(--hi)", 0.10))
    o.append(_tint(M0(0), 30, M0(W2)-M0(0), 270, "var(--accent)", 0.10))
    o.append(text((M0(-W1)+M0(0))/2, 24, "−ρ₁  (width W₁)", size=12, fill="var(--hi)")); o.append(text((M0(0)+M0(W2))/2, 24, "+ρ₂  (width W₂)", size=12, fill="var(--accent)"))
    # E panel
    M = axes(108, "Eₓ", 26, -2.2, 0.4)
    pts = [M(x, Ex(x)) for x in xs]
    o.append(path("M" + " L".join(f"{X:.1f},{Y:.1f}" for X, Y in pts), stroke="var(--accent)", w=2.4))
    o.append(text(M(0, -2.0)[0]+8, M(0, -2.0)[1]+4, "−ρ₁W₁/ε₀", size=11, anchor="start"))
    # V panel
    M = axes(250, "V", 24, -1.2, 2.2)
    pts = [M(x, V(x)) for x in xs]
    o.append(path("M" + " L".join(f"{X:.1f},{Y:.1f}" for X, Y in pts), stroke="var(--accent2)", w=2.4))
    o.append(text(M(-1.9, -1.0)[0], M(-1.9, -1.0)[1]-8, "V₁ = −ρ₁W₁²/(2ε₀)", size=11, anchor="start"))
    o.append(text(M(2.2, 2.0)[0], M(2.2, 2.0)[1]-8, "V₂ = +ρ₂W₂²/(2ε₀)", size=11, anchor="start"))
    o.append(text(M(0.5, 0.0)[0]+150, M(0.5,0.0)[1]+34, "concave up where ρ&lt;0, concave down where ρ&gt;0", size=11, fill="var(--muted)"))
    o.append(text(320, 318, "V(0) = 0 chosen as reference;  V₂₁ = V₂ − V₁ = (ρ₁W₁² + ρ₂W₂²)/(2ε₀)", size=12))
    o.append(svg_close())
    figure("pn-potential", "".join(o),
           "<strong>Poisson's equation across a pn junction.</strong> Top: the depletion-region field, a triangle pointing in −x (from the positive layer toward the negative one). Bottom: the potential, obtained by integrating dV/dx = −Eₓ twice; it is quadratic in each slab, matched by continuity at x = 0, and settles to plateaus V₁ and V₂ outside. Its curvature is V″ = −ρ/ε₀: bowl-shaped in the negative slab, dome-shaped in the positive one. The step V₂₁ is the junction's built-in potential.")

# --------------------------------------------------------------------------- 19. Conductor in a field
def fig_conductor_in_field():
    W, H = 640, 300
    o = [svg_open(W, H)]
    def sheets(x0, label):
        yt, yb = 60, 240
        o.append(line(x0, yt, x0+250, yt, stroke="var(--accent)", w=3)); o.append(line(x0, yb, x0+250, yb, stroke="var(--hi)", w=3))
        for k in range(8): o.append(text(x0+18+k*31, yt-8, "−", size=14, fill="var(--accent)"))
        for k in range(8): o.append(text(x0+18+k*31, yb+18, "+", size=13, fill="var(--hi)"))
        o.append(text(x0+255, yt+4, "−ρₛ", size=12, anchor="start")); o.append(text(x0+255, yb+4, "+ρₛ", size=12, anchor="start"))
        o.append(text(x0+125, 285, label, size=13, weight=600))
        return yt, yb
    # (a) empty gap
    yt, yb = sheets(30, "(a) E₀ = ẑ ρₛ/ε₀ between the sheets")
    for X in (80, 160, 240): o.append(arrow(X, yb-12, X, yt+12, stroke="var(--accent)", w=2.4))
    o.append(text(175, 150, "E₀", size=15, anchor="start", fill="var(--accent)", weight=600))
    # (b) conducting slab inserted
    x0 = 350; yt, yb = sheets(x0, "(b) slab inserted: E = 0 inside")
    st, sb = 115, 185
    o.append(f'<rect x="{x0+20}" y="{st}" width="210" height="{sb-st}" fill="var(--muted)" fill-opacity="0.35" stroke="currentColor" stroke-width="1.4"/>')
    for k in range(7): o.append(text(x0+38+k*30, st+13, "+", size=13, fill="var(--hi)"))
    for k in range(7): o.append(text(x0+38+k*30, sb-4, "−", size=14, fill="var(--accent)"))
    o.append(text(x0+125, 155, "σ > 0 :  E = 0,  V = const", size=12))
    o.append(text(x0+236, st+8, "+ρₛ", size=11, anchor="start")); o.append(text(x0+236, sb-2, "−ρₛ", size=11, anchor="start"))
    for X in (x0+60, x0+190):
        o.append(arrow(X, yb-10, X, sb+8, stroke="var(--accent)", w=2.2)); o.append(arrow(X, st-8, X, yt+10, stroke="var(--accent)", w=2.2))
    o.append(text(x0+125, yt+38, "E₀", size=13, fill="var(--accent)", weight=600)); o.append(text(x0+125, yb-24, "E₀", size=13, fill="var(--accent)", weight=600))
    o.append(svg_close())
    figure("conductor-in-field", "".join(o),
           "<strong>How a conductor kills the field inside itself.</strong> (a) Two charged sheets set up a uniform field E₀. (b) A conducting slab is inserted: free electrons drift against E₀ toward the bottom face, leaving the top face positive. The induced pair ±ρₛ produces −E₀ inside the slab, so the total interior field is zero and current stops. Outside the slab nothing changes. The induced charge sits only on the surfaces; the interior stays neutral and at one potential.")

# --------------------------------------------------------------------------- 20. Polarized dielectric slab
def fig_dielectric_slab_polarization():
    W, H = 640, 340
    o = [svg_open(W, H)]
    # left: sheets, slab with molecular dipoles, fields
    x0 = 30; yt, yb = 50, 270
    o.append(line(x0, yt, x0+300, yt, stroke="var(--accent)", w=3)); o.append(line(x0, yb, x0+300, yb, stroke="var(--hi)", w=3))
    for k in range(9): o.append(text(x0+22+k*32, yt-8, "−", size=14, fill="var(--accent)"))
    for k in range(9): o.append(text(x0+22+k*32, yb+18, "+", size=13, fill="var(--hi)"))
    o.append(text(x0+305, yt+4, "−ρₛ₀", size=12, anchor="start")); o.append(text(x0+305, yb+4, "+ρₛ₀", size=12, anchor="start"))
    st, sb = 110, 210
    o.append(f'<rect x="{x0+20}" y="{st}" width="{260}" height="{sb-st}" fill="var(--accent2)" fill-opacity="0.12" stroke="currentColor" stroke-width="1.3"/>')
    # molecules: ellipses with - at bottom and + at top? Field points up (from + sheet at bottom to - sheet on top),
    # so each dipole has + displaced upward: + on top, - on bottom.
    for i in range(5):
        for j in range(2):
            cx, cy = x0+72+i*40, st+28+j*44
            o.append(f'<ellipse cx="{cx}" cy="{cy}" rx="10" ry="15" fill="var(--paper)" stroke="currentColor" stroke-width="1"/>')
            o.append(text(cx, cy-4, "+", size=11, fill="var(--hi)")); o.append(text(cx, cy+11, "−", size=12, fill="var(--accent)"))
    # bound surface charge rows
    for i in range(5): o.append(text(x0+72+i*40, st+8, "+", size=11, fill="var(--hi)"))
    for i in range(5): o.append(text(x0+72+i*40, sb-3, "−", size=12, fill="var(--accent)"))
    o.append(text(x0+285, st+8, "+ρₛb", size=11, anchor="start")); o.append(text(x0+285, sb-2, "−ρₛb", size=11, anchor="start"))
    # field arrows: outside long, inside short
    o.append(arrow(x0+300-8, yb-10, x0+300-8, sb+6, stroke="var(--accent)", w=2.4)); o.append(arrow(x0+300-8, st-6, x0+300-8, yt+10, stroke="var(--accent)", w=2.4))
    o.append(text(x0+300-8-6, yb-30, "E₀", size=12, anchor="end", fill="var(--accent)", weight=600))
    o.append(arrow(x0+36, sb-12, x0+36, st+12, stroke="var(--accent)", w=2.4)); o.append(text(x0+44, (st+sb)/2+4, "E", size=13, anchor="start", fill="var(--accent)", weight=600))
    o.append(text(x0+30, sb+16, "E = E₀ − P/ε₀", size=11, anchor="start"))
    o.append(text(x0+150, 322, "dipole moment per volume: P = N p  [C/m²],  P ∥ E", size=12))
    # right: bar diagram D = eps0 E + P
    bx, by = 400, 250; bw = 54; scale = 150
    o.append(text(bx+110, 40, "D is the same inside and out (field ⟂ faces)", size=13, weight=600))
    # outside bar: eps0 E0
    o.append(f'<rect x="{bx}" y="{by-scale}" width="{bw}" height="{scale}" fill="var(--accent)" fill-opacity="0.35" stroke="currentColor" stroke-width="1.2"/>')
    o.append(text(bx+bw/2, by-scale/2+4, "ε₀E₀", size=12)); o.append(text(bx+bw/2, by+18, "outside", size=12))
    # inside bar: eps0 E (lower part) + P (upper part), E = E0/3 say
    fE = 0.35
    o.append(f'<rect x="{bx+120}" y="{by-fE*scale}" width="{bw}" height="{fE*scale}" fill="var(--accent)" fill-opacity="0.35" stroke="currentColor" stroke-width="1.2"/>')
    o.append(f'<rect x="{bx+120}" y="{by-scale}" width="{bw}" height="{(1-fE)*scale}" fill="var(--accent2)" fill-opacity="0.35" stroke="currentColor" stroke-width="1.2"/>')
    o.append(text(bx+120+bw/2, by-fE*scale/2+4, "ε₀E", size=12)); o.append(text(bx+120+bw/2, by-scale+(1-fE)*scale/2+4, "P", size=13, weight=600))
    o.append(text(bx+120+bw/2, by+18, "inside", size=12))
    o.append(line(bx-10, by-scale, bx+120+bw+10, by-scale, dash="4 3", w=1, opacity=0.6)); o.append(text(bx+120+bw+14, by-scale+4, "D", size=14, anchor="start", weight=600))
    o.append(text(bx+110, 290, "D ≡ ε₀E + P = ε₀E₀  (this slab)", size=12))
    o.append(text(bx+110, 310, "general law: ∇·D = ρ_free", size=12, fill="var(--muted)"))
    o.append(svg_close())
    figure("dielectric-slab-polarization", "".join(o),
           "<strong>A polarized dielectric weakens the field but does not kill it.</strong> Left: the applied field stretches each neutral molecule into a small dipole aligned with <b>E</b>; inside the slab neighbouring + and − ends cancel, leaving bound surface charge ±ρₛb on the faces, whose field opposes E₀. Right: the bookkeeping. Outside, ε₀E₀; inside, the smaller ε₀E plus the polarization <b>P</b> add up to the same total. That total is <b>D</b> — the quantity whose flux counts free charge only.")

# --------------------------------------------------------------------------- 21. Two-layer parallel plates
def fig_two_layer_plates():
    W, H = 640, 300
    o = [svg_open(W, H)]
    x0, yb, yt, ym = 40, 240, 60, 150
    o.append(_tint(x0, yt, 260, ym-yt, "var(--accent2)", 0.14)); o.append(_tint(x0, ym, 260, yb-ym, "var(--muted)", 0.08))
    o.append(line(x0, yb, x0+260, yb, stroke="currentColor", w=6, opacity=0.85)); o.append(line(x0, yt, x0+260, yt, stroke="currentColor", w=6, opacity=0.85))
    o.append(line(x0, ym, x0+260, ym, dash="5 4", w=1.2))
    for k in range(8): o.append(text(x0+22+k*32, yb+18, "−", size=14, fill="var(--accent)"))
    for k in range(8): o.append(text(x0+22+k*32, yt-8, "+", size=13, fill="var(--hi)"))
    o.append(text(x0+265, yb+4, "z = 0: −2ε₀,  V = 0", size=11, anchor="start"))
    o.append(text(x0+265, yt+4, "z = 2 m: +2ε₀,  V = ?", size=11, anchor="start"))
    o.append(text(x0+265, ym+4, "z = 1 m", size=11, anchor="start"))
    o.append(text(x0+10, ym-12, "ε = 2ε₀", size=12, anchor="start")); o.append(text(x0+10, yb-12, "ε = ε₀", size=12, anchor="start"))
    # D arrows (equal) and E arrows (different)
    for X in (x0+60, x0+140):
        o.append(arrow(X, yt+12, X, ym-8, stroke="var(--accent)", w=2.4)); o.append(arrow(X, ym+8, X, yb-12, stroke="var(--accent)", w=2.4))
    o.append(text(x0+60+8, ym-40, "D", size=13, anchor="start", fill="var(--accent)", weight=600))
    # E: length proportional to 1 (top) and 2 (bottom)
    X = x0+228
    o.append(arrow(X, ym-20-14, X, ym-20+14, stroke="var(--accent2)", w=3)); o.append(text(X-8, ym-16, "E = −1 ẑ", size=11, anchor="end", fill="var(--accent2)"))
    o.append(arrow(X, ym+45-28, X, ym+45+28, stroke="var(--accent2)", w=3)); o.append(text(X-8, ym+49, "E = −2 ẑ", size=11, anchor="end", fill="var(--accent2)"))
    o.append(text(x0+130, 278, "D = −2ε₀ ẑ in both layers (no free charge at z = 1)", size=12))
    # right: V(z) polyline
    px, py = 440, 230
    o.append(line(px, py, px+200, py, w=1.3, arrow=True, small=True)); o.append(text(px+208, py+4, "z", size=14, math=True, anchor="start"))
    o.append(line(px, py+10, px, py-170, w=1.3, arrow=True, small=True)); o.append(text(px-6, py-160, "V", size=14, math=True, anchor="end"))
    sz, sv = 80, 45
    pts = [(px, py), (px+sz, py-2*sv), (px+2*sz, py-3*sv)]
    o.append(path("M" + " L".join(f"{X},{Y}" for X, Y in pts), stroke="var(--accent2)", w=2.6))
    for (X, Y), lab in zip(pts[1:], ("2 V", "3 V")):
        o.append(line(px, Y, X, Y, dash="3 3", w=1, opacity=0.6)); o.append(text(px-6, Y+4, lab, size=11, anchor="end"))
        o.append(line(X, py, X, Y, dash="3 3", w=1, opacity=0.6))
    o.append(text(px+sz, py+16, "1", size=11)); o.append(text(px+2*sz, py+16, "2", size=11))
    o.append(text(px+40, py-2*sv-30, "slope 2", size=11, fill="var(--accent2)")); o.append(text(px+120, py-3*sv-14, "slope 1", size=11, fill="var(--accent2)"))
    o.append(text(px+100, 278, "V(2) = 2·1 + 1·1 = 3 V", size=12, weight=600))
    o.append(svg_close())
    figure("two-layer-plates", "".join(o),
           "<strong>Two dielectric layers between charged plates.</strong> With no free charge at the interface, the normal component of <b>D</b> is the same in both layers, so <b>E</b> = <b>D</b>/ε is twice as large in the ε₀ layer as in the 2ε₀ layer. The potential is piecewise linear with a kink at the interface — Laplace's equation holds in each homogeneous layer but not across the boundary between them.")

# --------------------------------------------------------------------------- 22. Refraction of E at a dielectric interface
def fig_field_refraction():
    W, H = 640, 320
    o = [svg_open(W, H)]
    yI = 160; cx = 215
    o.append(_tint(20, 30, 380, yI-30, "var(--accent)", 0.06)); o.append(_tint(20, yI, 380, 290-yI, "var(--muted)", 0.12))
    o.append(line(20, yI, 400, yI, w=2))
    o.append(text(392, 48, "medium 1:  ε₁ = 4ε₀", size=13, anchor="end")); o.append(text(392, 282, "medium 2:  ε₂ = 6ε₀", size=13, anchor="end"))
    o.append(arrow(cx, yI, cx, yI-110, stroke="currentColor", w=1.6)); o.append(text(cx+6, yI-100, "n̂", size=14, anchor="start", weight=600))
    o.append(line(cx, yI, cx, yI+110, w=1.2, dash="4 3", opacity=0.6))
    s = 16  # px per unit
    # E1 = (3, -6): drawn in medium 1 ending at the interface point
    o.append(arrow(cx-3*s, yI-6*s, cx, yI, stroke="var(--accent)", w=3)); o.append(text(cx-3*s-8, yI-6*s-6, "E₁ = 3x̂ − 6ẑ", size=13, anchor="end", fill="var(--accent)", weight=600))
    # E2 = (3, -4): drawn in medium 2 starting at the interface point
    o.append(arrow(cx, yI, cx+3*s, yI+4*s, stroke="var(--accent)", w=3)); o.append(text(cx+3*s+8, yI+4*s+6, "E₂ = 3x̂ − 4ẑ", size=13, anchor="start", fill="var(--accent)", weight=600))
    # component construction (dashed): tangential legs along the top/bottom, normal legs along n̂
    o.append(line(cx-3*s, yI-6*s, cx, yI-6*s, dash="3 3", w=1.2, opacity=0.7)); o.append(line(cx, yI+4*s, cx+3*s, yI+4*s, dash="3 3", w=1.2, opacity=0.7))
    o.append(arrow(cx-3*s, yI-6*s-10, cx, yI-6*s-10, stroke="var(--accent2)", w=2, small=True)); o.append(text(cx-1.5*s, yI-6*s-16, "Eₜ = 3", size=11, fill="var(--accent2)"))
    o.append(arrow(cx, yI+4*s+10, cx+3*s, yI+4*s+10, stroke="var(--accent2)", w=2, small=True)); o.append(text(cx+1.5*s, yI+4*s+26, "Eₜ = 3", size=11, fill="var(--accent2)"))
    o.append(text(cx+6, yI-3*s, "Eₙ = −6", size=11, anchor="start", fill="var(--muted)")); o.append(text(cx-6, yI+2*s+4, "Eₙ = −4", size=11, anchor="end", fill="var(--muted)"))
    # angles from the normal
    o.append(text(cx-6, yI-56, "θ₁", size=13, anchor="end")); o.append(text(cx+38, yI+34, "θ₂", size=13, anchor="start"))
    # right: rules
    rx = 418
    o.append(text(rx, 56, "tangential E continuous:", size=12, anchor="start", weight=600))
    o.append(text(rx, 76, "E₁ₜ = E₂ₜ = 3", size=12, anchor="start"))
    o.append(text(rx, 112, "normal D continuous", size=12, anchor="start", weight=600))
    o.append(text(rx, 130, "(no free surface charge):", size=11, anchor="start", fill="var(--muted)"))
    o.append(text(rx, 150, "ε₁E₁ₙ = ε₂E₂ₙ", size=12, anchor="start"))
    o.append(text(rx, 168, "E₂ₙ = (4/6)(−6) = −4", size=12, anchor="start"))
    o.append(text(rx, 204, "Entering the higher-ε side", size=12, anchor="start"))
    o.append(text(rx, 222, "the normal part shrinks, so", size=12, anchor="start"))
    o.append(text(rx, 240, "the field tilts toward the surface:", size=12, anchor="start"))
    o.append(text(rx, 266, "tan θ₁ / tan θ₂ = ε₁ / ε₂", size=13, anchor="start", weight=600))
    o.append(svg_close())
    figure("field-refraction", "".join(o),
           "<strong>How E bends at a dielectric interface.</strong> Split the field into a tangential part (parallel to the surface) and a normal part. The tangential part is the same on both sides; the normal part scales by ε₁/ε₂ so that <b>D</b>ₙ = εEₙ is the same on both sides. Here the field enters the higher-permittivity medium, its normal component shrinks from 6 to 4, and the arrow tilts toward the surface. Angles are measured from the normal.")

# --------------------------------------------------------------------------- 23. Coaxial capacitor
def fig_coax_capacitor():
    W, H = 640, 340
    o = [svg_open(W, H)]
    cx, cy = 170, 165; a, b, rg = 42, 120, 82
    o.append(circle(cx, cy, b+12, fill="var(--muted)", stroke="currentColor", w=1.2, opacity=0.9))
    o.append(circle(cx, cy, b, fill="var(--accent2)", stroke="currentColor", w=1.2, opacity=0.9))
    o.append(f'<circle cx="{cx}" cy="{cy}" r="{b}" fill="var(--paper)" fill-opacity="0.75" stroke="none"/>')
    o.append(circle(cx, cy, a, fill="var(--muted)", stroke="currentColor", w=1.2))
    o.append(text(cx, cy+5, "+λ", size=14, weight=600))
    for k in range(12):
        th = 2*math.pi*k/12
        o.append(text(cx+(b-9)*math.cos(th), cy+(b-9)*math.sin(th)+4, "−", size=13, fill="var(--accent)"))
    o.append(text(20, cy+b+30, "−λ induced on the inner face of the outer conductor", size=11, anchor="start"))
    for k in range(8):
        th = 2*math.pi*k/8 + 0.2
        o.append(arrow(cx+(a+6)*math.cos(th), cy+(a+6)*math.sin(th), cx+(b-16)*math.cos(th), cy+(b-16)*math.sin(th), stroke="var(--accent)", w=1.8, small=True, opacity=0.85))
    o.append(circle(cx, cy, rg, fill="none", stroke="var(--hi)", w=1.6, dash="6 4"))
    o.append(arrow(cx, cy, cx+a*math.cos(-2.4), cy+a*math.sin(-2.4), stroke="currentColor", w=1.4, small=True)); o.append(text(cx-22, cy-30, "a", size=13, math=True))
    o.append(arrow(cx, cy, cx+rg*math.cos(0.6), cy+rg*math.sin(0.6), stroke="var(--hi)", w=1.4, small=True)); o.append(text(cx+42, cy+56, "r", size=13, math=True, fill="var(--hi)"))
    o.append(arrow(cx, cy, cx+b*math.cos(-0.7), cy+b*math.sin(-0.7), stroke="currentColor", w=1.4, small=True)); o.append(text(cx+66, cy-100, "b", size=13, math=True))
    o.append(text(20, cy+b+48, "ε (and σ, if lossy) fills a &lt; r &lt; b", size=12, anchor="start"))
    o.append(text(20, 28, "Gauss on the dashed cylinder (length L):  D_r · 2πrL = λL", size=12, anchor="start", weight=600))
    # right: E(r) and V(r)
    px, py = 380, 250; w = 220; h = 170
    o.append(line(px, py, px+w+10, py, w=1.3, arrow=True, small=True)); o.append(text(px+w+18, py+4, "r", size=14, math=True, anchor="start"))
    o.append(line(px, py+8, px, py-h-10, w=1.3, arrow=True, small=True))
    ra, rb = 1.0, 3.0
    def X(r): return px + (r-0.4)/(3.4-0.4)*w
    o.append(line(X(ra), py, X(ra), py-h, dash="3 3", w=1, opacity=0.6)); o.append(text(X(ra), py+16, "a", size=12, math=True))
    o.append(line(X(rb), py, X(rb), py-h, dash="3 3", w=1, opacity=0.6)); o.append(text(X(rb), py+16, "b", size=12, math=True))
    rs = np.linspace(ra, rb, 60)
    E = [1.0/r for r in rs]; Emax = 1.0/ra
    o.append(path("M" + " L".join(f"{X(r):.1f},{py-0.9*h*e/Emax:.1f}" for r, e in zip(rs, E)), stroke="var(--accent)", w=2.4))
    o.append(text(px+112, py-150, "E_r ∝ 1/r", size=12, anchor="start", fill="var(--accent)", weight=600))
    Vv = [math.log(rb/r) for r in rs]; Vmax = math.log(rb/ra)
    o.append(path("M" + " L".join(f"{X(r):.1f},{py-0.9*h*v/Vmax:.1f}" for r, v in zip(rs, Vv)), stroke="var(--accent2)", w=2.4))
    o.append(text(px+112, py-132, "V(r) ∝ ln(b/r)", size=12, anchor="start", fill="var(--accent2)", weight=600))
    o.append(text(px+w/2, 46, "V(a) − V(b) = (λ/2πε) ln(b/a)", size=12, weight=600))
    o.append(text(px+w/2, 66, "𝒞 = λ / [V(a) − V(b)] = 2πε / ln(b/a)  [F/m]", size=12))
    o.append(svg_close())
    figure("coax-capacitor", "".join(o),
           "<strong>The coaxial capacitor.</strong> Left: cross-section with +λ per unit length on the inner conductor and −λ induced on the inner face of the outer conductor; the field is radial and a dashed Gaussian cylinder of radius r encloses λL. Right: E falls as 1/r and the potential falls logarithmically from the inner to the outer conductor. The capacitance per unit length depends only on ε and the radius ratio b/a.")

# --------------------------------------------------------------------------- 24. Lossy capacitor and its RC model
def fig_lossy_capacitor_rc():
    W, H = 640, 300
    o = [svg_open(W, H)]
    # left: plates with lossy medium
    x0, yt, yb = 40, 60, 230
    o.append(_tint(x0, yt, 240, yb-yt, "var(--accent2)", 0.14))
    o.append(line(x0, yt, x0+240, yt, stroke="currentColor", w=6, opacity=0.85)); o.append(line(x0, yb, x0+240, yb, stroke="currentColor", w=6, opacity=0.85))
    for k in range(7): o.append(text(x0+25+k*32, yt+16, "+", size=13, fill="var(--hi)"))
    for k in range(7): o.append(text(x0+25+k*32, yb-6, "−", size=14, fill="var(--accent)"))
    o.append(text(x0+245, yt+4, "+Q,  V", size=12, anchor="start")); o.append(text(x0+245, yb+4, "−Q,  0", size=12, anchor="start"))
    o.append(text(x0+120, (yt+yb)/2-30, "ε,  σ", size=14, weight=600))
    o.append(arrow(x0+60, yt+26, x0+60, yb-26, stroke="var(--accent)", w=2.4)); o.append(text(x0+68, (yt+yb)/2+4, "E", size=14, anchor="start", fill="var(--accent)", weight=600))
    o.append(arrow(x0+170, yt+26, x0+170, yb-26, stroke="var(--hi)", w=2.4)); o.append(text(x0+178, (yt+yb)/2+4, "J = σE", size=13, anchor="start", fill="var(--hi)", weight=600))
    o.append(text(x0+120, 268, "leakage current I = ∫J·dS = GV  drains the charge", size=12))
    # right: parallel RC
    px, py = 370, 60; wd = 230; ht = 170
    # top and bottom rails
    o.append(line(px, py, px+wd, py, w=1.6)); o.append(line(px, py+ht, px+wd, py+ht, w=1.6))
    # capacitor branch (left)
    xc = px+60
    o.append(line(xc, py, xc, py+ht/2-10, w=1.6)); o.append(line(xc, py+ht/2+10, xc, py+ht, w=1.6))
    o.append(line(xc-22, py+ht/2-10, xc+22, py+ht/2-10, stroke="var(--accent)", w=3)); o.append(line(xc-22, py+ht/2+10, xc+22, py+ht/2+10, stroke="var(--accent)", w=3))
    o.append(text(xc-30, py+ht/2+4, "C", size=14, anchor="end", math=True))
    # resistor branch (right): zigzag
    xr = px+170
    o.append(line(xr, py, xr, py+ht/2-30, w=1.6)); o.append(line(xr, py+ht/2+30, xr, py+ht, w=1.6))
    zz = [(xr, py+ht/2-30)]
    for i in range(6): zz.append((xr + (10 if i % 2 == 0 else -10), py+ht/2-30+(i+0.5)*10))
    zz.append((xr, py+ht/2+30))
    o.append(path("M" + " L".join(f"{X},{Y}" for X, Y in zz), stroke="var(--hi)", w=2.2))
    o.append(text(xr+18, py+ht/2+4, "R = 1/G", size=13, anchor="start"))
    # terminals and labels
    o.append(circle(px, py, 3.5, fill="currentColor")); o.append(circle(px, py+ht, 3.5, fill="currentColor"))
    o.append(text(px-8, py+4, "+", size=14, anchor="end")); o.append(text(px-8, py+ht+4, "−", size=14, anchor="end"))
    o.append(text(px-8, py+ht/2+4, "V(t)", size=13, anchor="end"))
    o.append(arrow(px-40, py-18, px-6, py-18, stroke="currentColor", w=1.6, small=True)); o.append(text(px-44, py-14, "I", size=14, anchor="end", math=True))
    o.append(text(px+wd/2, py+ht+38, "I = C dV/dt + G V", size=13, weight=600))
    o.append(text(px+wd/2, py+ht+58, "τ = RC = C/G = ε/σ  (independent of geometry)", size=12))
    o.append(svg_close())
    figure("lossy-capacitor-rc", "".join(o),
           "<strong>A capacitor with a leaky dielectric.</strong> If the filling has conductivity σ, the same field that stores charge also drives a current density <b>J</b> = σ<b>E</b> from plate to plate. The structure behaves as a capacitance C in parallel with a conductance G = (σ/ε)C: the two share one field pattern, so their ratio is a material property. Left alone, the charge decays as e<sup>−t/τ</sup> with τ = ε/σ — the same relaxation time that makes conductors field-free inside.")

# --------------------------------------------------------------------------- 25. Capacitance recipes as flowcharts
def fig_capacitance_recipes():
    W, H = 640, 250
    o = [svg_open(W, H)]
    def chain(y, title, nodes, labels, color):
        o.append(text(30, y-26, title, size=13, anchor="start", weight=600, fill=color))
        n = len(nodes); x0 = 40; step = (W-80)/(n-1) if n > 1 else 0
        for i, nd in enumerate(nodes):
            X = x0 + i*step
            o.append(f'<rect x="{X-24}" y="{y-16}" width="48" height="32" rx="7" fill="{color}" fill-opacity="0.16" stroke="currentColor" stroke-width="1.2"/>')
            o.append(text(X, y+5, nd, size=14, weight=600))
            if i < n-1:
                o.append(arrow(X+26, y, X+step-26, y, stroke="currentColor", w=1.5, small=True))
                o.append(text(X+step/2, y-10, labels[i], size=10.5, fill="var(--muted)"))
    chain(60, "given V (potentials on the conductors): Laplace first", ["V", "E", "D", "ρₛ", "Q", "C"],
          ["E = −∇V", "D = εE", "ρₛ = n̂·D", "Q = ∫ρₛ dA", "C = Q/V"], "var(--accent)")
    chain(140, "given Q (charge on the conductors): Gauss first", ["Q", "D", "E", "V", "C"],
          ["Gauss's law", "E = D/ε", "V = −∫E·dl", "C = Q/V"], "var(--accent2)")
    chain(220, "conductance of the same geometry (lossy filling)", ["V", "E", "J", "I", "G"],
          ["E = −∇V", "J = σE", "I = ∫J·dS", "G = I/V"], "var(--hi)")
    o.append(svg_close())
    figure("capacitance-recipes", "".join(o),
           "<strong>Three chains, one habit.</strong> Capacitance and conductance are never computed from a formula you remember; they fall out of a chain of steps whose direction depends on what the problem hands you. Given the voltage, start from Laplace's equation and walk to the charge; given the charge, start from Gauss's law and walk to the voltage. For conductance, replace ε by σ and D by J: the geometry factor is the same, which is why G/C = σ/ε whenever one homogeneous medium fills the field region.")

if __name__ == "__main__":
    for f in (fig_dl_ds_cube, fig_coulomb_pair, fig_line_charge_side, fig_dipole_map, fig_flux_patch,
              fig_gauss_cylinder, fig_sheet_slab, fig_flux_plane, fig_curl_div_panels, fig_paddlewheel,
              fig_stokes_tiling, fig_pn_superposition,
              fig_equipotentials_gradient, fig_path_independence, fig_bc_pillbox_loop, fig_bc_summary_panels,
              fig_laplace_plates, fig_pn_potential, fig_conductor_in_field, fig_dielectric_slab_polarization,
              fig_two_layer_plates, fig_field_refraction, fig_coax_capacitor, fig_lossy_capacitor_rc,
              fig_capacitance_recipes):
        f()
