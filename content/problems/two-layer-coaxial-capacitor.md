---
title: "Worked problem — Two-layer coaxial capacitor"
description: "A coax whose dielectric has two concentric layers: D from Gauss's law (one expression for both layers), E layer by layer, the induced charge on the outer conductor, the potential difference, and the capacitance per unit length. Modelled on FA26 Exam 1, problem 4, with the permittivities reordered."
tags: [problem, electrostatics, exam-1]
---

*Problem · style of FA26 Exam 1 #4 (30 pts) · uses [[1-electrostatics/03-gauss-law-at-work|Lecture 3]], [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]], [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]], [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]] · concepts: [[concepts/gauss-law]], [[concepts/conductors]], [[concepts/permittivity]], [[concepts/capacitance]]*

> [!question] Problem
> A long coaxial cable has an inner conductor of radius $a$ and an outer conductor of inner radius $b$. The space between is filled with two perfect dielectrics separated by the cylinder $r = c$ ($a<c<b$):
> $$
> \epsilon_1 = 3\epsilon_0\ \ (a<r<c),\qquad \epsilon_2 = 2\epsilon_0\ \ (c<r<b).
> $$
> The inner conductor carries a charge $+\lambda$ per unit length ($\lambda$ is the exam's symbol for the line density $\rho_l$). Keep the radii and $\lambda$ symbolic. (a) $\mathbf{D}$ for $a<r<b$. (b) $\mathbf{E}$ for $a<r<b$. (c) The charge per unit length on the outer conductor. (d) $V(b)-V(a)$. (e) The capacitance per unit length $\mathcal{C}$. Then put in $a = 1$ mm, $c = 2$ mm, $b = 4$ mm.

## (a) D from the free charge alone

Symmetry: the cable is long and the charge is uniform along it, so $\mathbf{D} = D_r(r)\,\hat{r}$. Take a Gaussian cylinder of radius $r$ ($a<r<b$) and length $L$, coaxial with the cable. No flux crosses the end caps ($\mathbf{D}\perp d\mathbf{S}$ there); the curved wall has area $2\pi rL$ and $\mathbf{D}\parallel d\mathbf{S}$ on it. The enclosed **free** charge is $\lambda L$ (the dielectric's bound charge is not counted — that is the whole point of $\mathbf{D}$):

$$
\oint\mathbf{D}\cdot d\mathbf{S} = D_r\cdot2\pi rL = \lambda L .
$$

> [!key] Answer (a)
> $$
> \mathbf{D} = \frac{\lambda}{2\pi r}\,\hat{r}\quad[\text{C/m}^2],\qquad a<r<b\ \text{— one expression, valid in both layers.}
> $$

Nothing in the derivation mentioned $\epsilon$, so the answer holds on both sides of $r=c$. Equivalently: normal $\mathbf{D}$ is continuous across the dielectric interface (no free charge there), and $\mathbf{D}$ here is entirely normal to it.

> [!trap]
> - Writing $\mathbf{D} = \lambda/(2\pi\epsilon r)$, i.e. giving $\mathbf{D}$ a permittivity — that is $\mathbf{E}$.
> - Using a sphere's $4\pi r^2$ instead of the cylinder's $2\pi rL$.
> - Dropping $\hat{r}$: $\mathbf{D}$ is a vector.

## (b) E, layer by layer

$\mathbf{E} = \mathbf{D}/\epsilon$ with the $\epsilon$ of whichever layer $r$ is in:

> [!key] Answer (b)
> $$
> \mathbf{E} = \begin{cases}\dfrac{\lambda}{2\pi(3\epsilon_0)r}\hat{r} = \dfrac{\lambda}{6\pi\epsilon_0 r}\hat{r}, & a<r<c\\[10pt] \dfrac{\lambda}{2\pi(2\epsilon_0)r}\hat{r} = \dfrac{\lambda}{4\pi\epsilon_0 r}\hat{r}, & c<r<b\end{cases}\qquad[\text{V/m}]
> $$

At $r = c$ the field *jumps up* by the factor $\epsilon_1/\epsilon_2 = 3/2$ on leaving the higher-$\epsilon$ layer: normal $\mathbf{E}$ is discontinuous even though nothing charged sits there (the bound charge on the interface — see the bonus below — is what does it). With the numbers $a = 1$, $c = 2$, $b = 4$ mm and $\lambda = 1$ nC/m: $E(a^+) \approx 5990$ V/m, $E(c^-)\approx3000$ V/m, $E(c^+)\approx4490$ V/m, $E(b^-)\approx2250$ V/m. The largest field is at the inner conductor, as always in a coax.

## (c) Charge on the outer conductor

Inside the metal of the outer conductor the field is zero (a conductor in equilibrium, [[concepts/conductors]]). A Gaussian cylinder whose curved wall lies *inside that metal* therefore has zero flux and must enclose zero net charge:

$$
\lambda + \lambda_{\text{outer, inner surface}} = 0 .
$$

> [!key] Answer (c)
> The outer conductor carries $-\lambda$ [C/m], all of it on its **inner** surface ($r = b$), as surface charge $\rho_s = -\lambda/(2\pi b)$. (Check with the boundary condition at $r=b$, $\hat{n} = -\hat{r}$ pointing out of the metal into the dielectric: $\rho_s = \hat{n}\cdot\mathbf{D}(b) = -\lambda/2\pi b$ ✓.)

Standard coax assumption: the outer conductor is the return path, carrying equal and opposite charge, with no field outside the cable. If the problem said the outer conductor was *isolated and neutral*, it would carry $-\lambda$ on its inner surface and $+\lambda$ on its outer surface, with a field $\lambda/(2\pi\epsilon_0r)$ outside; the answers to (a)–(e) would not change.

> [!trap]
> "Conductors are neutral, so 0" and "$+\lambda$" are the two wrong reflexes. The metal's *interior field* must vanish, and that forces $-\lambda$ onto the surface facing the inner conductor.

## (d) The potential difference

$V(b)-V(a) = -\int_a^b\mathbf{E}\cdot d\mathbf{l}$ along a radial path, $d\mathbf{l} = \hat{r}\,dr$. The integrand changes form at $r = c$, so split there:

$$
V(b)-V(a) = -\int_a^c\frac{\lambda\,dr}{6\pi\epsilon_0r} - \int_c^b\frac{\lambda\,dr}{4\pi\epsilon_0r}
= -\frac{\lambda}{6\pi\epsilon_0}\ln\frac ca - \frac{\lambda}{4\pi\epsilon_0}\ln\frac bc .
$$

> [!key] Answer (d)
> $$
> V(b)-V(a) = -\frac{\lambda}{2\pi\epsilon_0}\left[\frac13\ln\frac ca + \frac12\ln\frac bc\right] = -\frac{\lambda}{12\pi\epsilon_0}\left[2\ln\frac ca + 3\ln\frac bc\right]\quad[\text{V}].
> $$
> **Negative**, as it must be: the inner conductor is positive and field lines run outward, so the potential *falls* from $a$ to $b$. Numerically, with $\lambda = 1$ nC/m and $c/a = b/c = 2$: $V(b)-V(a) = -\dfrac{5\ln2}{12\pi\epsilon_0}\times10^{-9} \approx -10.4$ V.

> [!trap]
> - Getting $V(b)-V(a)>0$: the minus sign in $V = -\int\mathbf{E}\cdot d\mathbf{l}$ was lost, or the limits reversed twice. The FA26 key absorbs the sign by writing $\int_b^c$ and $\int_c^a$; that is correct, but decide on *one* way to handle the sign and stick to it.
> - Using one $\epsilon$ throughout, or forgetting to split the integral at $c$.
> - Inverting the logarithms ($\ln(a/c)$ is negative).

## (e) Capacitance per unit length

$\mathcal{C} = \lambda/V$ where $V$ is the **drop from the positive conductor to the negative one**, $V(a)-V(b) = +\dfrac{\lambda}{2\pi\epsilon_0}\Big[\tfrac13\ln\tfrac ca + \tfrac12\ln\tfrac bc\Big]$:

> [!key] Answer (e)
> $$
> \mathcal{C} = \frac{\lambda}{V(a)-V(b)} = \frac{2\pi\epsilon_0}{\dfrac13\ln\dfrac ca + \dfrac12\ln\dfrac bc} = \frac{12\pi\epsilon_0}{2\ln\dfrac ca + 3\ln\dfrac bc}\quad[\text{F/m}].
> $$
> For $c/a = b/c = 2$: $\mathcal{C} = \dfrac{12\pi\epsilon_0}{5\ln2}\approx 96.3$ pF/m. (The $\lambda$ cancels, as it must — $\mathcal{C}$ is a property of the cable.)

**Series check** — the one-line route to (e). The two layers are crossed one after the other by the same $\mathbf{D}$, so they are two coax capacitors in **series**, each with the single-dielectric formula $\mathcal{C}_i = 2\pi\epsilon_i/\ln(r_{\text{out}}/r_{\text{in}})$:

$$
\frac{1}{\mathcal{C}} = \frac{1}{\mathcal{C}_1}+\frac{1}{\mathcal{C}_2} = \frac{\ln(c/a)}{2\pi\cdot3\epsilon_0} + \frac{\ln(b/c)}{2\pi\cdot2\epsilon_0} \quad\checkmark
$$

**Limit check.** With $\epsilon_1 = \epsilon_2 = \epsilon$ the bracket becomes $\ln(b/a)/\epsilon_r$ and $\mathcal{C}\to2\pi\epsilon/\ln(b/a)$, the single-dielectric coax of [[1-electrostatics/10-capacitance-and-conductance#3-the-coaxial-cable|Lecture 10 §3]] ✓.

> [!trap]
> - A **negative** $\mathcal{C}$ from dividing by $V(b)-V(a)$; the exam key's own intermediate line does this and silently fixes the sign in the box.
> - Adding the layers in parallel ($\mathcal{C}_1+\mathcal{C}_2$): parallel is for dielectrics side by side *along* the field (e.g. each filling half the angle), where $\mathbf{E}$ is common.
> - Units: F/m, not F. If a length $L$ is given, $C = \mathcal{C}L$.

## Bonus: where the bound charge sits, and the stored energy

In each layer $\mathbf{P} = \mathbf{D}-\epsilon_0\mathbf{E} = \Big(1-\dfrac1{\epsilon_r}\Big)\dfrac{\lambda}{2\pi r}\hat{r}$: $\tfrac23\cdot\dfrac{\lambda}{2\pi r}$ in the inner layer and $\tfrac12\cdot\dfrac{\lambda}{2\pi r}$ in the outer. Since $rP_r$ is constant in each layer, $\rho_b = -\nabla\cdot\mathbf{P} = 0$ in the bulk; the bound charge is all on the three cylindrical surfaces, $\rho_{sb} = \mathbf{P}\cdot\hat{n}$ with $\hat{n}$ out of the dielectric. Per unit length:

- at $r = a$ ($\hat{n} = -\hat{r}$): $-\tfrac23\lambda$;
- at $r = c$: two layers meet, and each contributes — the inner layer's outer face ($\hat{n} = +\hat{r}$) carries $+\tfrac23\lambda$, the outer layer's inner face ($\hat{n} = -\hat{r}$) carries $-\tfrac12\lambda$, net $+\tfrac16\lambda$;
- at $r = b$ ($\hat{n} = +\hat{r}$): $+\tfrac12\lambda$.

Total: $-\tfrac23+\tfrac16+\tfrac12 = 0$ ✓. The negative layer hugging the inner conductor is why $E$ in the inner layer is $\lambda/(6\pi\epsilon_0r)$ rather than $\lambda/(2\pi\epsilon_0r)$: Gauss's law for $\epsilon_0\mathbf{E}$ with *all* charge, $(\lambda - \tfrac23\lambda)/(2\pi\epsilon_0r) = \lambda/(6\pi\epsilon_0r)$ ✓. The positive layer at $r=c$ is what makes $E$ jump *up* on crossing into the outer layer.

Stored energy per unit length: $W' = \dfrac{\lambda^2}{2\mathcal{C}} = \dfrac{\lambda^2}{4\pi\epsilon_0}\Big[\tfrac13\ln\tfrac ca + \tfrac12\ln\tfrac bc\Big]\approx5.2$ nJ/m for $\lambda = 1$ nC/m — the same as $\int\tfrac12\epsilon E^2\,dA$ over the cross-section, layer by layer.

## The general pattern

For layers $\epsilon_{r1}$ ($a<r<c$) and $\epsilon_{r2}$ ($c<r<b$):

$$
\mathbf{D} = \frac{\lambda}{2\pi r}\hat{r},\qquad \mathbf{E} = \frac{\lambda}{2\pi\epsilon_0\epsilon_{ri}r}\hat{r},\qquad
V(a)-V(b) = \frac{\lambda}{2\pi\epsilon_0}\left[\frac{\ln(c/a)}{\epsilon_{r1}}+\frac{\ln(b/c)}{\epsilon_{r2}}\right],\qquad
\mathcal{C} = \frac{2\pi\epsilon_0}{\dfrac{\ln(c/a)}{\epsilon_{r1}}+\dfrac{\ln(b/c)}{\epsilon_{r2}}} .
$$

The exam used $\epsilon_{r1} = 2$, $\epsilon_{r2} = 4$ (giving $\mathcal{C} = 8\pi\epsilon_0/[2\ln(c/a)+\ln(b/c)]$); this page used $3$ and $2$. Three or more layers simply add more terms to the sum in the denominator. The "flip the given" variant — a stated voltage $V_0$ instead of $\lambda$ — is solved by the same $\mathcal{C}$: $\lambda = \mathcal{C}V_0$, then everything else follows.
