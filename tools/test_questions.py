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


PHYS = {"conservative-field-analysis": phys_conservative, "nonuniform-slab": phys_slab,
        "dielectric-interface": phys_interface, "coax-two-layer": phys_coax, "flux-through-plane": phys_flux}

for topic in sorted(os.listdir(QDIR)):
    for q in sorted(os.listdir(os.path.join(QDIR, topic))):
        qpath = os.path.join(QDIR, topic, q)
        info = json.load(open(os.path.join(qpath, "info.json")))
        check(info.get("type") == "v3" and "uuid" in info and "topic" in info, f"{q}: info.json")
        html = open(os.path.join(qpath, "question.html")).read()
        stripped = re.sub(r"\{\{\{\s*[A-Za-z0-9_.]+\s*\}\}\}|\{\{[#/^]?\s*[A-Za-z0-9_.]+\s*\}\}", "", html)
        check("{{" not in stripped, f"{q}: stray Mustache braces after removing valid tags (LaTeX brace collision?)")
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
