#!/usr/bin/env python3
"""Offline checks for the ECE 329 PrairieLearn questions (no PrairieLearn needed).

For every question directory:
  * import server.py and run generate() on N seeds (random seeded like PL does);
  * data must be JSON-serializable; every {{params.x}} in question.html must exist;
    every answers-name must have a correct_answers entry (or be a choice element);
  * multiple-choice: exactly one correct choice, all choice texts distinct;
  * physics: closed-form answers re-checked numerically (finite differences / quadrature).
Run:  python3 tools/test_questions.py [N]
"""
import importlib.util, json, math, os, random, re, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QDIR = os.path.join(ROOT, "questions")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
EPS0 = 8.8541878128e-12


def load(qpath):
    spec = importlib.util.spec_from_file_location("server", os.path.join(qpath, "server.py"))
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def run(mod, seed):
    random.seed(seed); np.random.seed(seed % (2**32 - 1))
    data = {"params": {}, "correct_answers": {}, "variant_seed": seed, "options": {}}
    mod.generate(data)
    json.dumps(data)   # must be serializable
    return data


def mustache_refs(html):
    refs = set(re.findall(r"\{\{\{?\s*params\.([A-Za-z0-9_]+)\s*\}?\}\}", html))
    loops = set(re.findall(r"\{\{#params\.([A-Za-z0-9_]+)\}\}", html))
    names = re.findall(r'answers-name="([^"]+)"', html)
    choice_elems = re.findall(r"<pl-(multiple-choice|checkbox)[^>]*answers-name=\"([^\"]+)\"", html)
    return refs, loops, names, {n for _, n in choice_elems}


def sym_eval(expr, **vals):
    """Evaluate a sympy-style string with numpy (only +,-,*,/,** and names)."""
    if not re.fullmatch(r"[0-9A-Za-z_+\-*/().\s]+", expr):
        raise ValueError("unexpected characters in expression: " + expr)
    return eval(expr, {"__builtins__": {}}, vals)


errors = []
def check(cond, msg):
    if not cond:
        errors.append(msg)


# ------------------------------------------------------------------ per-question physics checks
def phys_conservative(d):
    p = d["params"]; alpha, beta, m, n, pp = p["alpha"], p["beta"], p["m"], p["n"], p["p"]
    E0, eps = 1.7, 2.3
    def E(x, y, z):
        return np.array([alpha*m*x**(m-1)*y**n, alpha*n*x**m*y**(n-1), beta*pp*z**(pp-1)]) * E0
    Vs, rhos = d["correct_answers"]["V"], d["correct_answers"]["rho"]
    h = 1e-5
    for _ in range(5):
        x, y, z = np.random.uniform(0.5, 2, 3)
        V = lambda x, y, z: sym_eval(Vs, x=x, y=y, z=z, E0=E0)
        grad = np.array([(V(x+h,y,z)-V(x-h,y,z)), (V(x,y+h,z)-V(x,y-h,z)), (V(x,y,z+h)-V(x,y,z-h))])/(2*h)
        check(np.allclose(-grad, E(x,y,z), rtol=1e-4, atol=1e-6), f"conservative: -grad V != E at {(x,y,z)}")
        divE = ((E(x+h,y,z)-E(x-h,y,z))[0] + (E(x,y+h,z)-E(x,y-h,z))[1] + (E(x,y,z+h)-E(x,y,z-h))[2])/(2*h)
        rho = sym_eval(rhos, x=x, y=y, z=z, E0=E0, epsilon=eps)
        check(abs(rho - eps*divE) < 1e-4*max(1, abs(rho)), f"conservative: rho != eps*div E ({rho} vs {eps*divE})")
        curlz = ((E(x+h,y,z)-E(x-h,y,z))[1] - (E(x,y+h,z)-E(x,y-h,z))[0])/(2*h)
        check(abs(curlz) < 1e-6, "conservative: curl_z != 0")
    # line integral: V(0)-V(B) in units of E0
    bx, by, bz = [int(v) for v in re.findall(r"-?\d+", p["B_tex"])]
    val = (sym_eval(Vs, x=0, y=0, z=0, E0=1.0) - sym_eval(Vs, x=bx, y=by, z=bz, E0=1.0))
    check(abs(val - d["correct_answers"]["lineint"]) < 1e-9, "conservative: line integral mismatch")
    check(sum(c["correct"] == "true" for c in p["curl_choices"]) == 1, "conservative: MC needs exactly one correct")
    check(len({c["text"] for c in p["curl_choices"]}) == 4, f"conservative: duplicate MC choices {[c['text'] for c in p['curl_choices']]}")


def phys_slab(d):
    p = d["params"]; k = p["k"]; rho0 = p["rho0_nC"]*1e-9; a = p["a_mm"]*1e-3; b = p["b_over_a"]*a
    xin = a*p["xin_num"]/p["xin_den"]
    xs = np.linspace(-a, a, 200001); rho = rho0*np.abs(xs/a)**k
    sigma = np.trapezoid(rho, xs)
    check(abs(sigma*1e9 - d["correct_answers"]["sigma"]) < 1e-6*sigma*1e9, "slab: sigma")
    def Ex(x):   # Gauss with symmetric pillbox: 2 eps0 E = int_{-x}^{x} rho
        xx = np.linspace(-x, x, 20001); return np.trapezoid(rho0*np.abs(xx/a)**k*(np.abs(xx) <= a), xx)/(2*EPS0)
    check(abs(Ex(xin) - d["correct_answers"]["E_in"]) < 1e-5*abs(Ex(xin)), "slab: E_in")
    E_out = sigma/(2*EPS0)   # pillbox enclosing the whole slab
    check(abs(E_out - d["correct_answers"]["E_out"]) < 1e-5*abs(E_out), "slab: E_out")
    xs2 = np.linspace(0, b, 4001); Vb = -np.trapezoid([Ex(x) for x in xs2], xs2)
    check(abs(Vb - d["correct_answers"]["V_b"]) < 2e-4*abs(Vb), f"slab: V_b {Vb} vs {d['correct_answers']['V_b']}")


def phys_interface(d):
    p = d["params"]; er1, er2 = p["er1"], p["er2"]; c = d["correct_answers"]
    check(c["E2x"] == p["E1x"], "interface: E2x")
    check(abs(er2*c["E2z"] - er1*p["E1z"]) < 1e-9, "interface: normal D continuity")
    check(abs(c["D1z"] - er1*p["E1z"]) < 1e-9 and abs(c["D2x"] - er2*p["E1x"]) < 1e-9, "interface: D")
    check(abs(c["P1x"] - (er1-1)*p["E1x"]) < 1e-9 and abs(c["P2z"] - (er2-1)*c["E2z"]) < 1e-9, "interface: P")
    check(abs(c["rho_sb"] - (p["E1z"] - c["E2z"])) < 1e-9, "interface: bound charge = E1z - E2z")
    check(float(c["E2z"]).is_integer(), "interface: E2z not integer")


def phys_coax(d):
    p = d["params"]; a, c, b = p["a_mm"]*1e-3, p["c_mm"]*1e-3, p["b_mm"]*1e-3; er1, er2 = p["er1"], p["er2"]; lam = p["lam_nC"]*1e-9
    r1, r2 = 0.5*(a+c), 0.5*(c+b)
    check(abs(d["correct_answers"]["D_r1"] - lam/(2*math.pi*r1)*1e9) < 1e-9, "coax: D")
    check(abs(d["correct_answers"]["E_r2"] - lam/(2*math.pi*EPS0*er2*r2)) < 1e-6, "coax: E")
    rs = np.linspace(a, b, 200001); E = np.where(rs < c, lam/(2*math.pi*EPS0*er1*rs), lam/(2*math.pi*EPS0*er2*rs))
    Vab = np.trapezoid(E, rs)
    check(abs(Vab - d["correct_answers"]["Vab"]) < 1e-4*Vab, f"coax: Vab {Vab} vs {d['correct_answers']['Vab']}")
    Cseries = 1/(math.log(c/a)/(2*math.pi*EPS0*er1) + math.log(b/c)/(2*math.pi*EPS0*er2))
    check(abs(Cseries*1e12 - d["correct_answers"]["Cp"]) < 1e-6*Cseries*1e12, "coax: C' vs series formula")
    check(d["correct_answers"]["q_outer"] == -p["lam_nC"], "coax: outer charge")


def phys_flux(d):
    p = d["params"]; ch = p["choices"]
    check(sum(c["correct"] == "true" for c in ch) == 1, "flux: exactly one correct")
    check(len({c["text"] for c in ch}) == len(ch) and len(ch) >= 4, f"flux: choices not distinct or fewer than 4: {[c['text'] for c in ch]}")
    up = p["normal_tex"].startswith("+")
    psi = (p["q_below"] - p["q_above"])/2 * (1 if up else -1)
    correct = [c["text"] for c in ch if c["correct"] == "true"][0]
    # parse the correct text back to a number
    m = re.search(r"\$(-?)(?:\\dfrac\{(\d*)Q\}\{(\d+)\}|(\d*)Q|0)\$", correct)
    val = 0.0
    if m and m.group(3): val = int(m.group(2) or 1)/int(m.group(3))
    elif m and m.group(4) is not None: val = float(m.group(4) or 1)
    if m and m.group(1) == "-": val = -val
    check(abs(val - psi) < 1e-9, f"flux: correct choice {correct} != {psi}")


def _one_correct(ch, name):
    check(sum(c["correct"] == "true" for c in ch) == 1, f"{name}: MC needs exactly one correct")
    check(len({c["text"] for c in ch}) == len(ch), f"{name}: duplicate MC choices")


def phys_two_layer(d):
    p = d["params"]; c = d["correct_answers"]; er1, er2, dd, t2 = p["er1"], p["er2"], p["d"], p["t2"]
    check(abs(er1*p["E1z"] - er2*c["E2z"]) < 1e-9, "two-layer: D_z continuity")
    check(float(c["E2z"]).is_integer() and p["E1z"] < 0, "two-layer: E2z integer / E1z sign")
    n = (dd + t2)*10000; h = (dd + t2)/n; zm = (np.arange(n) + 0.5)*h; Ez = np.where(zm < dd, p["E1z"], c["E2z"])   # midpoint rule, cells aligned to z = d
    check(abs(-np.sum(Ez)*h - c["Vp"]) < 1e-9*abs(c["Vp"]), "two-layer: Vp = -int E dz")
    check(abs(c["Vd"] + p["E1z"]*dd) < 1e-9 and c["Vp"] > c["Vd"] > 0, "two-layer: V(d), monotone potential")
    check(abs(c["rho_top"] + er2*c["E2z"]) < 1e-9 and c["rho_top"] > 0, "two-layer: top plate charge")
    C_series = 1/(dd/(er1*EPS0) + t2/(er2*EPS0))
    check(abs(C_series*1e12 - c["C_per_A"]) < 1e-6*C_series*1e12, "two-layer: C/A series")
    check(abs(c["rho_top"]*EPS0/c["Vp"] - C_series) < 1e-9*C_series, "two-layer: C/A = rho_s/Vp")


def phys_sheet(d):
    p = d["params"]; c = d["correct_answers"]; er1, er2, dd, t2, s = p["er1"], p["er2"], p["d"], p["t2"], p["s"]
    check(abs(er2*c["E2z"] - er1*c["E1z"] - s) < 1e-9, "sheet: D jump = rho_s")
    check(abs(c["E1z"]*dd + c["E2z"]*t2) < 1e-9, "sheet: V(z0) = 0 (int E over both regions)")
    check(abs(-c["E1z"]*dd - c["V0"]) < 1e-9, "sheet: V0 = -E1z d")
    check(abs(c["rho_0"] + c["rho_z0"] + s) < 1e-9, "sheet: plate charges sum to -rho_s")
    check(abs(c["rho_0"] - er1*c["E1z"]) < 1e-9 and abs(c["rho_z0"] + er2*c["E2z"]) < 1e-9, "sheet: plate BCs")
    check(np.sign(c["V0"]) == np.sign(s), "sheet: sign of V0")


def phys_insertion(d):
    p = d["params"]; c = d["correct_answers"]; er, E0 = p["er"], p["E0"]
    if p["mode"] == "Q":
        check(c["D"] == E0 and abs(c["E"] - E0/er) < 1e-12 and abs(c["Wratio"] - 1/er) < 1e-12, "insertion Q: D fixed")
    else:
        check(c["E"] == E0 and c["D"] == er*E0 and c["Wratio"] == er, "insertion V: E fixed")
    check(abs(c["P"] - (c["D"] - c["E"])) < 1e-9, "insertion: P = D - eps0 E")
    check(abs(c["D"] - er*c["E"]) < 1e-9, "insertion: D = eps E")
    check(abs(c["Wratio"] - (c["D"]*c["E"])/(E0*E0)) < 1e-9, "insertion: energy ratio = DE/D0E0")
    _one_correct(p["choices"], "insertion")
    correct = [x["text"] for x in p["choices"] if x["correct"] == "true"][0]
    check(("= 0$" in correct) == (p["mode"] == "Q"), "insertion: MC correct choice does not match the mode")


def phys_lossy(d):
    p = d["params"]; c = d["correct_answers"]; eps = p["er"]*EPS0; sig = p["sigma"]
    g = p["geometry"]
    if g == "plates":
        area = lambda x: p["A_cm2"]*1e-4 + 0*x; lo, hi = 0.0, p["d_mm"]*1e-3
    elif g == "coax":
        area = lambda r: 2*math.pi*r*p["L_m"]; lo, hi = p["a_mm"]*1e-3, p["b_mm"]*1e-3
    else:
        area = lambda r: 4*math.pi*r*r; lo, hi = p["a_cm"]*1e-2, p["b_cm"]*1e-2
    xs = np.linspace(lo, hi, 400001)
    inv_area = 1/area(xs)
    R = np.trapezoid(inv_area/sig, xs)          # sum of shell resistances dR = dr/(sigma area)
    Vunit = np.trapezoid(inv_area/eps, xs)      # V for unit charge, from D = Q/area
    check(abs(1/Vunit*1e12 - c["C"]) < 1e-4*c["C"], f"lossy {g}: C by integration {1/Vunit*1e12} vs {c['C']}")
    check(abs(1/R*1e9 - c["G"]) < 1e-4*c["G"], f"lossy {g}: G by shell resistances")
    check(abs(c["tau"] - eps/sig*1e6) < 1e-9*c["tau"] and abs(c["C"]*1e-12/(c["G"]*1e-9) - eps/sig) < 1e-9*eps/sig, "lossy: tau = C/G = eps/sigma")
    check(abs(c["Qt"] - p["Q0"]*math.exp(-p["t1"]/c["tau"])) < 1e-9, "lossy: Q(t1)")
    check(abs(c["I0"] - p["Q0"]*1e-9/(c["tau"]*1e-6)*1e6) < 1e-9*c["I0"], "lossy: I0 = Q0/tau")


def phys_ampere_coax(d):
    p = d["params"]; c = d["correct_answers"]; a, b, cc = p["a_mm"]*1e-3, p["b_mm"]*1e-3, p["c_mm"]*1e-3; I = p["I"]
    Ja, Jb = I/(math.pi*a*a), -I/(math.pi*(cc*cc - b*b))
    def Ienc(r):        # integrate J(r') 2 pi r' dr' by the midpoint rule; cells (1 um) align with a, b, c (integer mm)
        n = int(round(r/1e-6)); h = r/n; rm = (np.arange(n) + 0.5)*h; J = np.where(rm < a, Ja, np.where((rm > b) & (rm < cc), Jb, 0.0))
        return np.sum(J*2*math.pi*rm)*h
    for key, r in (("H1", a/2), ("H2", (a+b)/2), ("H3", (b+cc)/2), ("H4", 2*cc)):
        H = Ienc(r)/(2*math.pi*r)
        check(abs(H - c[key]) < 1e-6*I/(2*math.pi*r), f"ampere-coax: {key} {H} vs {c[key]}")
    check(abs(c["H4same"] - 2*I/(2*math.pi*2*cc)) < 1e-9, "ampere-coax: H4same")
    check(c["H1"] > 0 and c["H2"] > c["H3"] > 0 and c["H4"] == 0, "ampere-coax: ordering")


def phys_sheets(d):
    p = d["params"]; c = d["correct_answers"]; J1, J2 = p["Js1"], p["Js2"]; MU0 = 4e-7*math.pi
    for key, x in (("Hleft", -1.0), ("Hmid", p["d"]/2), ("Hright", p["d"] + 1.0)):
        H = 0.5*J1*np.sign(x) + 0.5*J2*np.sign(x - p["d"])      # H = 1/2 Js x n per sheet
        check(abs(H - c[key]) < 1e-9, f"sheets: {key}")
    f = np.cross(np.array([0, 0, J2]), np.array([0, MU0*J1/2, 0]))   # Js2 x B1(x=d)
    check(abs(f[0]*1e6 - c["fx"]) < 1e-9 and abs(f[1]) + abs(f[2]) < 1e-15, "sheets: force")
    check((c["fx"] < 0) == (J1*J2 > 0), "sheets: parallel currents attract")
    check(("+" in p["dir_correct"]) == (J1 > 0), "sheets: direction MC")


def phys_solenoid(d):
    p = d["params"]; c = d["correct_answers"]; MU0 = 4e-7*math.pi
    n = p["N"]/(p["L_cm"]*1e-2); A = math.pi*(p["a_mm"]*1e-3)**2
    check(abs(c["H"] - n*p["I"]) < 1e-9, "solenoid: H")
    check(abs(c["Psi"] - MU0*n*p["I"]*A*1e6) < 1e-9, "solenoid: Psi")
    check(abs(c["Lind"] - p["N"]*c["Psi"]*1e-6/p["I"]*1e6) < 1e-9, "solenoid: L = N Psi / I")
    check(abs(c["W"] - 0.5*MU0*c["H"]**2*A*p["L_cm"]*1e-2*1e6) < 1e-6*c["W"], "solenoid: W = energy density x volume")
    check(abs(c["tau"] - c["Lind"]/p["R"]) < 1e-9, "solenoid: tau = L/R")
    check(p["a_mm"]*1e-3 <= p["L_cm"]*1e-2/5, "solenoid: not long")


def phys_faraday(d):
    p = d["params"]; c = d["correct_answers"]; A = p["a_cm"]*p["b_cm"]*1e-4
    Psi = lambda t: (p["B0"] + p["k"]*t)*A
    h = 1e-4
    check(abs(c["Psi"] - Psi(p["t1"])*1e3) < 1e-9, "faraday: Psi(t1)")
    check(abs(c["emf"] + p["N"]*(Psi(0.5+h) - Psi(0.5-h))/(2*h)) < 1e-6, "faraday: emf = -N dPsi/dt")
    check(abs(c["I"] - c["emf"]/p["R"]*1e3) < 1e-9, "faraday: I = emf/R")
    _one_correct(p["sense_choices"], "faraday sense"); _one_correct(p["lenz_choices"], "faraday lenz")
    sense = [x["text"] for x in p["sense_choices"] if x["correct"] == "true"][0]
    lenz = [x["text"] for x in p["lenz_choices"] if x["correct"] == "true"][0]
    check(sense.startswith("counterclockwise") == (c["emf"] > 0), "faraday: sense MC vs sign of emf")
    check(("+" in lenz) == (p["k"] < 0), "faraday: Lenz — induced field +z iff flux decreasing")


PHYS = {"conservative-field-analysis": phys_conservative, "nonuniform-slab": phys_slab,
        "dielectric-interface": phys_interface, "coax-two-layer": phys_coax, "flux-through-plane": phys_flux,
        "two-layer-parallel-plates": phys_two_layer, "sheet-between-grounded-plates": phys_sheet,
        "dielectric-insertion": phys_insertion, "lossy-capacitor-relaxation": phys_lossy,
        "ampere-coax-fields": phys_ampere_coax, "current-sheets-superposition": phys_sheets,
        "solenoid-inductance-energy": phys_solenoid, "faraday-loop-emf": phys_faraday}

for topic in sorted(os.listdir(QDIR)):
    for q in sorted(os.listdir(os.path.join(QDIR, topic))):
        qpath = os.path.join(QDIR, topic, q)
        info = json.load(open(os.path.join(qpath, "info.json")))
        check(info.get("type") == "v3" and "uuid" in info and "topic" in info, f"{q}: info.json")
        html = open(os.path.join(qpath, "question.html")).read()
        stripped = re.sub(r"\{\{\{\s*[A-Za-z0-9_.]+\s*\}\}\}|\{\{[#/^]?\s*[A-Za-z0-9_.]+\s*\}\}", "", html)
        check("{{" not in stripped, f"{q}: stray Mustache braces after removing valid tags (LaTeX brace collision?)")
        for block in re.findall(r"<pl-multiple-choice.*?</pl-multiple-choice>", html, re.S):
            if "{{#params" not in block:      # static choices: exactly one must be marked correct
                check(block.count('correct="true"') == 1, f"{q}: static multiple-choice block needs exactly one correct answer")
        refs, loops, names, choice_names = mustache_refs(html)
        mod = load(qpath)
        variants = set()
        for seed in range(1, N + 1):
            d = run(mod, seed)
            for r in refs | loops:
                check(r in d["params"], f"{q}: question.html references params.{r} which generate() did not set (seed {seed})")
            for nme in names:
                if nme not in choice_names:
                    check(nme in d["correct_answers"], f"{q}: answers-name {nme} has no correct_answers entry")
            PHYS[q](d)
            variants.add(json.dumps(d["params"], sort_keys=True))
            if errors and len(errors) > 20: break
        print(f"{q:32s} {N} seeds, {len(variants)} distinct variants, {'OK' if not errors else 'ERRORS'}")

if errors:
    print("\n".join(errors[:30])); sys.exit(1)
print("all checks passed")
