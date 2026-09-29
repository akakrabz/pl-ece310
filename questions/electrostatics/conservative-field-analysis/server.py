"""Conservative polynomial field: curl, potential, charge density, line integral.

Family of FA26 Exam 1, problem 1. The field is generated FROM a potential
    V = -E0 * (alpha * x^m * y^n + beta * z^p),   p = m + n,
so it is curl-free by construction and every part has a closed form.
Symbolic answers are provided as sympy-parsable strings (accepted by
pl-symbolic-input); no sympy import is needed at generate time.
"""
import random

# ----------------------------------------------------------------- string helpers
def mono(coef, xp=0, yp=0, zp=0, tex=False):
    """coef * x^xp * y^yp * z^zp, as a sympy string (tex=False) or LaTeX (tex=True). '' if coef == 0."""
    if coef == 0:
        return ""
    parts = []
    for v, pw in (("x", xp), ("y", yp), ("z", zp)):
        if pw == 1:
            parts.append(v)
        elif pw > 1:
            parts.append(f"{v}^{{{pw}}}" if tex else f"{v}**{pw}")
    body = ("" if tex else "*").join(parts)
    if not body:
        return str(coef)
    if coef == 1:
        return body
    if coef == -1:
        return "-" + body
    return f"{coef}{body}" if tex else f"{coef}*{body}"


def join(terms):
    """Sum of non-empty terms with explicit signs; '0' if nothing survives."""
    terms = [t for t in terms if t]
    if not terms:
        return "0"
    out = terms[0]
    for t in terms[1:]:
        out += (" - " + t[1:]) if t.startswith("-") else (" + " + t)
    return out


def generate(data):
    alpha = random.choice([-3, -2, -1, 1, 2, 3])
    beta = random.choice([-2, -1, 1, 2])
    m, n = random.choice([(2, 1), (1, 2), (2, 2), (1, 1)])
    p = m + n

    # E / E0 = -grad(V)/E0 :  Ex = alpha*m x^(m-1) y^n,  Ey = alpha*n x^m y^(n-1),  Ez = beta*p z^(p-1)
    data["params"]["field_tex"] = join([
        mono(alpha * m, m - 1, n, 0, tex=True) + r"\,\hat{x}",
        mono(alpha * n, m, n - 1, 0, tex=True) + r"\,\hat{y}",
        mono(beta * p, 0, 0, p - 1, tex=True) + r"\,\hat{z}",
    ])

    # ---- (a) z-component of the curl: correct answer 0; distractors from real mistakes
    cross = alpha * m * n                       # both cross partials equal cross * x^(m-1) y^(n-1)
    def curl_like(c):                           # c * x^(m-1) y^(n-1) E0, in TeX
        return "$" + mono(c, m - 1, n - 1, 0, tex=True) + r"\,E_0$"
    # a "divergence-looking" distractor: one second-derivative term of the field
    if n >= 2:
        div_like = mono(alpha * n * (n - 1), m, n - 2, 0, tex=True)
    elif m >= 2:
        div_like = mono(alpha * m * (m - 1), m - 2, n, 0, tex=True)
    else:
        div_like = mono(beta * p * (p - 1), 0, 0, p - 2, tex=True)
    correct = {"text": "$0$ (the field is curl-free)", "correct": "true"}
    candidates = [
        curl_like(2 * cross),                 # added the two cross partials instead of subtracting
        curl_like(-2 * cross),                # ... with the sign reversed
        "$" + div_like + r"\,E_0$",           # a divergence term mistaken for the curl
        curl_like(cross),                     # kept only one of the two partials
        curl_like(-cross),
        "$" + mono(beta * p, 0, 0, p - 1, tex=True) + r"\,E_0$",   # E_z itself
    ]
    distractors = []
    for cnd in candidates:
        if cnd != correct["text"] and cnd not in distractors and cnd not in ("$\\,E_0$", "$0\\,E_0$"):
            distractors.append(cnd)
        if len(distractors) == 3:
            break
    choices = [correct] + [{"text": t, "correct": "false"} for t in distractors]
    random.shuffle(choices)
    data["params"]["curl_choices"] = choices

    # ---- (b) potential with V(0,0,0) = 0
    V_inner = join([mono(alpha, m, n, 0), mono(beta, 0, 0, p)])
    data["correct_answers"]["V"] = f"-E0*({V_inner})"
    data["params"]["V_tex"] = r"-E_0\left(" + join([mono(alpha, m, n, 0, tex=True), mono(beta, 0, 0, p, tex=True)]) + r"\right)"

    # ---- (c) rho = epsilon0 * div E
    rho_terms = lambda tex: [
        mono(alpha * m * (m - 1), m - 2, n, 0, tex=tex) if m >= 2 else "",
        mono(alpha * n * (n - 1), m, n - 2, 0, tex=tex) if n >= 2 else "",
        mono(beta * p * (p - 1), 0, 0, p - 2, tex=tex),
    ]
    data["correct_answers"]["rho"] = f"epsilon*E0*({join(rho_terms(False))})"
    data["params"]["rho_tex"] = r"\epsilon_0 E_0\left(" + join(rho_terms(True)) + r"\right)"

    # ---- (d) line integral A=(0,0,0) -> B equals V(A) - V(B) = E0 * (alpha bx^m by^n + beta bz^p)
    bx, by, bz = random.choice([(1, 1, 0), (1, 2, 0), (2, 1, 0), (1, 1, 1), (2, 1, 1), (1, 2, 1)])
    val = alpha * bx**m * by**n + beta * bz**p
    data["params"]["B_tex"] = f"({bx},\\,{by},\\,{bz})"
    data["params"]["lineint_val"] = val
    data["correct_answers"]["lineint"] = float(val)

    # for the answer panel
    data["params"].update({"alpha": alpha, "beta": beta, "m": m, "n": n, "p": p})
