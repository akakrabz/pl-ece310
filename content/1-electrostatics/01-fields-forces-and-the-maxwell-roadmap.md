---
title: "Lecture 1 — Fields, forces, and the Maxwell roadmap"
description: "What a field is, how E and B are defined through the force on a charge, the four Maxwell equations as the course roadmap, and the two new tools of the course: dl and dS."
tags: [lecture, electrostatics, exam-1]
lecture: 1
---

*Lecture 1 · course notes §1 · slides "Introduction" · assumes MATH 241 vectors and PHYS 212 · next: [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]]*

> [!abstract] In one breath
> Charges push on other charges. Instead of tracking every pair, we describe the push with two **vector fields**, $\mathbf{E}(\mathbf{r},t)$ and $\mathbf{B}(\mathbf{r},t)$, defined by the force they exert on a test charge. Four equations — Maxwell's — say how charges and currents *create* those fields. This whole course is those four equations, first in the easy static case, then with time variation, then as waves on cables. Today: the definitions, the roadmap, and the two small vectors we'll integrate along ($d\mathbf{l}$) and through ($d\mathbf{S}$) for the next fifteen weeks.

## 1. Why fields at all?

Electrons and protons interact "at a distance". The bookkeeping that makes this tractable is a **field**: a quantity defined at every point $(x,y,z)$ and every time $t$. A scalar field assigns a number to each point (temperature, elevation, potential $V$); a vector field assigns a vector (wind velocity, $\mathbf{E}$, $\mathbf{B}$).

The electromagnetic field is a pair of vector fields, and their *operational definition* is the force they produce on a charge $q$ moving with velocity $\mathbf{v}$:

> [!key] Lorentz force — the definition of E and B
> $$
> \mathbf{F} = q\left(\mathbf{E} + \mathbf{v}\times\mathbf{B}\right)
> $$
> - $\mathbf{E}$ is **force per unit stationary charge** (put $\mathbf{v}=0$): units N/C = **V/m**.
> - $\mathbf{B}$ is defined by the *extra* force on a *moving* charge, $\mathbf{F}_M = q\,\mathbf{v}\times\mathbf{B}$. That force is always perpendicular to $\mathbf{v}$ (and to $\mathbf{B}$), so it can do no work; a charge moving *along* $\mathbf{B}$ feels none. Units: **T** (tesla) = Wb/m².
>
> Newton's second law then reads $m\,d\mathbf{v}/dt = q(\mathbf{E}+\mathbf{v}\times\mathbf{B})$, valid for $|\mathbf{v}|\ll c$.

Two things worth noticing now, because they return later:

- **Fields are frame-dependent.** A charge at rest in your frame is a current in mine. Consequently $\mathbf{E}$ and $\mathbf{B}$ mix under a change of frame (to first order, $\mathbf{E}' = \mathbf{E} + \mathbf{v}\times\mathbf{B}$). This is *why* the equations for $\mathbf{E}$ and $\mathbf{B}$ must be coupled — [[1-electrostatics/04-divergence-and-curl#7-the-differential-form-of-maxwells-equations|Faraday's law]] is the first place we'll see it.
- **Charge comes in both signs, mass doesn't.** Ordinary matter is nearly neutral (equal numbers of protons and electrons), which is why gravity — far weaker per particle — wins in the macroscopic world, and why the *net* charge density $\rho$ in a material is usually a tiny imbalance.

See the concept pages: [[concepts/lorentz-force]] · [[concepts/electric-field]].

## 2. The roadmap: Maxwell's equations on day one

The slides show all four equations on day one, in both integral and differential form, and the course counts them off as each is derived ("our first Maxwell equation!" arrives in [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]]). Here they are in free space, differential form:

| | Divergence equations | Curl equations |
|---|---|---|
| **electric** | $\nabla\cdot\mathbf{E} = \rho/\epsilon_0$ | $\nabla\times\mathbf{E} = -\dfrac{\partial\mathbf{B}}{\partial t}$ |
| **magnetic** | $\nabla\cdot\mathbf{B} = 0$ | $\nabla\times\mathbf{B} = \mu_0\mathbf{J} + \mu_0\epsilon_0\dfrac{\partial\mathbf{E}}{\partial t}$ |

with the constants

$$
\mu_0 \approx 4\pi\times 10^{-7}\ \tfrac{\text{H}}{\text{m}}\ (\text{exact before the 2019 SI redefinition}),\qquad
\epsilon_0 = \frac{1}{\mu_0 c^2}\approx \frac{1}{36\pi\times 10^{9}}\ \tfrac{\text{F}}{\text{m}}\approx 8.854\times10^{-12}\ \tfrac{\text{F}}{\text{m}},\qquad
c = \frac{1}{\sqrt{\mu_0\epsilon_0}}\approx 3\times10^{8}\ \tfrac{\text{m}}{\text{s}}.
$$

> [!intuition] How to read the table (this is the whole course in one paragraph)
> The **left column** says what fields *emanate from*: electric field lines start and end on charge ($\rho$); magnetic field lines never start or end anywhere. The **right column** says what makes fields *circulate*: a changing $\mathbf{B}$ curls $\mathbf{E}$ around it (Faraday: generators, transformers), and a current $\mathbf{J}$ or a changing $\mathbf{E}$ curls $\mathbf{B}$ (Ampère–Maxwell: electromagnets, antennas). Sources on the right-hand sides, fields on the left.

Two consequences we will *derive* later but should be *announced* now:

1. **Statics decouple.** If nothing changes in time, both $\partial/\partial t$ terms vanish, and the table splits into two independent problems:
   - *Electrostatics* (Lectures 1–11): $\nabla\cdot\mathbf{E}=\rho/\epsilon_0$ and $\nabla\times\mathbf{E}=0$ — a **curl-free** field fixed by charge.
   - *Magnetostatics* (Lectures 12–13): $\nabla\cdot\mathbf{B}=0$ and $\nabla\times\mathbf{B}=\mu_0\mathbf{J}$ — a **divergence-free** field fixed by current. (Lectures 14–15, Faraday's law and inductance, are where $\partial\mathbf{B}/\partial t$ first re-enters.)
2. **Waves.** With the time derivatives kept, the equations support solutions $\mathbf{E},\mathbf{B}\propto\cos\!\big(2\pi f\,(t - z/c)+\phi\big)$ in empty space — disturbances travelling at $c$ with wavelength $\lambda = c/f$. That is Lecture 18; it is also why the course ends with transmission lines and the Smith chart.

### Quasi-statics: when is a circuit "a circuit"?

A time-varying source of period $T$ can be treated *statically* — using the instantaneous values of $\rho$ and $\mathbf{J}$ — whenever the field has time to "communicate" across the system before the source changes:

$$
T \gg \frac{L}{c}\quad\Longleftrightarrow\quad L \ll \lambda .
$$

That is the **quasi-static approximation**, and it is exactly the assumption behind lumped-circuit theory (ECE 210): a circuit is a circuit only while it is *electrically small*. Reduce $\lambda$ (raise the clock frequency) or enlarge $L$ (a long cable, a PCB trace at GHz) and delays $L/c$ matter — that is the world of Lectures 27–38.

> [!tip] A number to carry around
> $c \approx 3\times10^{8}$ m/s $= 300$ m/µs $= 30$ cm/ns. So 1 MHz ↔ 300 m wavelength; 1 GHz ↔ 30 cm; a 3 GHz clock has $\lambda = 10$ cm, which is *not* large compared with a motherboard.

Concept page: [[concepts/maxwells-equations]].

## 3. Units and the source → field bookkeeping

| Quantity | Symbol | Units | Sourced by |
|---|---|---|---|
| electric field | $\mathbf{E}$ | V/m | charge |
| electric flux density | $\mathbf{D}=\epsilon_0\mathbf{E}$ (+$\mathbf{P}$ in matter) | C/m² | charge |
| magnetic field | $\mathbf{H}$ | A/m | current |
| magnetic flux density | $\mathbf{B}=\mu_0\mathbf{H}$ (in vacuum) | T = Wb/m² | current |
| charge densities | point $Q$; line $\rho_l$; surface $\rho_s$; volume $\rho$ | C; C/m; C/m²; C/m³ | — |
| current densities | wire $I$; sheet $\mathbf{J}_s$; volume $\mathbf{J}$ | A; A/m; A/m² | — |

The instructor's blue brace on the slide is the point: **charges → $\mathbf{E},\mathbf{D}$; currents → $\mathbf{H},\mathbf{B}$.** Keep the "dimension ladder" of the densities in your head — it is how you decide whether an answer's units can possibly be right. Details on [[concepts/charge-density]] and [[0-toolkit/03-units-and-constants]].

## 4. Vector algebra you must be fluent in (review at home)

The slides mark this block "review at home"; the course assumes it. If any line below is not automatic, do the exercises on [[0-toolkit/02-vector-calculus-cheatsheet]] tonight.

- A vector in Cartesian components: $\mathbf{A} = A_x\hat{x} + A_y\hat{y} + A_z\hat{z}$ (the slides also write $\hat{a}_x$ and $\langle A_x, A_y, A_z\rangle$ — same thing).
- Magnitude $|\mathbf{A}| = \sqrt{A_x^2+A_y^2+A_z^2}$; unit vector $\hat{a} = \mathbf{A}/|\mathbf{A}|$.
- **Separation vector** from point 1 to point 2: $\mathbf{R}_{12} = \mathbf{r}_2 - \mathbf{r}_1$ — *final minus initial*. Its unit vector $\hat{R}_{12}$ points from 1 toward 2. This convention runs through every Coulomb's-law calculation.
- Dot product: $\mathbf{A}\cdot\mathbf{B} = |\mathbf{A}||\mathbf{B}|\cos\theta = A_xB_x + A_yB_y + A_zB_z$. Projection of $\mathbf{A}$ on $\hat{b}$ is $\mathbf{A}\cdot\hat{b}$.
- Cross product: $|\mathbf{A}\times\mathbf{B}| = |\mathbf{A}||\mathbf{B}|\sin\theta$, direction by the right-hand rule, computed by the determinant
  $$
  \mathbf{A}\times\mathbf{B} = \begin{vmatrix}\hat{x}&\hat{y}&\hat{z}\\ A_x&A_y&A_z\\ B_x&B_y&B_z\end{vmatrix},
  \qquad \hat{x}\times\hat{y}=\hat{z},\ \hat{y}\times\hat{z}=\hat{x},\ \hat{z}\times\hat{x}=\hat{y}.
  $$

> [!example] In-class discussion problem
> $\mathbf{A} = 3\hat{x}+2\hat{y}+\hat{z}$, $\mathbf{B} = \hat{x}+\hat{y}-\hat{z}$, $\mathbf{C} = \hat{x}+2\hat{y}+3\hat{z}$. Find $\mathbf{B}\times\mathbf{C}$ and $\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C})$.
>
> Expand the determinant along the top row, keeping the alternating sign on $\hat{y}$:
> $$
> \mathbf{B}\times\mathbf{C} = \hat{x}(1\cdot3-(-1)\cdot2) - \hat{y}(1\cdot3-(-1)\cdot1) + \hat{z}(1\cdot2-1\cdot1) = 5\hat{x}-4\hat{y}+\hat{z}.
> $$
> Then $\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C}) = 3\cdot5 + 2\cdot(-4) + 1\cdot1 = 8$. (Geometrically: the signed volume of the parallelepiped spanned by the three vectors.)

## 5. The two new tools: dl and dS

Everything in the integral form of Maxwell's equations is a line integral or a surface integral. Both are built from a small vector element.

<figure class="ece-fig"><svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 300" width="640" height="300" role="img" style="font-family:'Source Sans 3','Source Sans Pro',system-ui,sans-serif;font-size:14px"><defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker><marker id="ahs" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker></defs><text x="130.0" y="35.0" text-anchor="middle" fill="currentColor" style="font-size:13px;">dl is always tangent to the path</text><path d="M40,230 C80,100 160,260 240,120" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><line x1="65.3" y1="184.5" x2="91.3" y2="162.6" stroke="var(--accent)" stroke-width="2.2" marker-end="url(#ah)" stroke-linecap="round"/><line x1="125.0" y1="178.8" x2="158.5" y2="184.7" stroke="var(--accent)" stroke-width="2.2" marker-end="url(#ah)" stroke-linecap="round"/><line x1="197.0" y1="169.9" x2="225.5" y2="151.4" stroke="var(--accent)" stroke-width="2.2" marker-end="url(#ah)" stroke-linecap="round"/><text x="143.0" y="208.8" text-anchor="middle" fill="var(--accent)" style="font-size:15px;font-family:'STIX Two Math','Latin Modern Math','Cambria Math','Times New Roman',serif;font-style:italic;">dl</text><text x="34.0" y="238.0" text-anchor="end" fill="currentColor" style="font-size:13px;">P₁</text><text x="246.0" y="118.0" text-anchor="start" fill="currentColor" style="font-size:13px;">P₂</text><text x="140.0" y="275.0" text-anchor="middle" fill="currentColor" style="font-size:13px;">Cartesian: dl = x̂ dx + ŷ dy + ẑ dz</text><line x1="430.0" y1="190.0" x2="525.0" y2="190.0" stroke="currentColor" stroke-width="1.1" stroke-dasharray="4 4" opacity="0.6" stroke-linecap="round"/><line x1="430.0" y1="190.0" x2="430.0" y2="95.0" stroke="currentColor" stroke-width="1.1" stroke-dasharray="4 4" opacity="0.6" stroke-linecap="round"/><line x1="430.0" y1="190.0" x2="377.8" y2="223.2" stroke="currentColor" stroke-width="1.1" stroke-dasharray="4 4" opacity="0.6" stroke-linecap="round"/><polygon points="377.8,223.2 472.8,223.2 472.8,128.2 377.8,128.2" fill="var(--accent)" fill-opacity="0.18" stroke="currentColor" stroke-width="1.4" stroke-linejoin="round"/><polygon points="525.0,190.0 472.8,223.2 472.8,128.2 525.0,95.0" fill="var(--accent2)" fill-opacity="0.18" stroke="currentColor" stroke-width="1.4" stroke-linejoin="round"/><polygon points="430.0,95.0 377.8,128.2 472.8,128.2 525.0,95.0" fill="var(--muted)" fill-opacity="0.18" stroke="currentColor" stroke-width="1.4" stroke-linejoin="round"/><line x1="425.2" y1="175.8" x2="396.5" y2="194.0" stroke="var(--accent)" stroke-width="2.4" marker-end="url(#ah)" stroke-linecap="round"/><line x1="498.9" y1="159.1" x2="551.1" y2="159.1" stroke="var(--accent2)" stroke-width="2.4" marker-end="url(#ah)" stroke-linecap="round"/><line x1="451.4" y1="111.6" x2="451.4" y2="59.4" stroke="currentColor" stroke-width="2.4" marker-end="url(#ah)" stroke-linecap="round"/><text x="386.5" y="216.0" text-anchor="middle" fill="currentColor" style="font-size:12px;">dS = x̂ dy dz</text><text x="555.1" y="163.1" text-anchor="start" fill="currentColor" style="font-size:12px;">dS = ŷ dz dx</text><text x="451.4" y="53.4" text-anchor="middle" fill="currentColor" style="font-size:12px;">dS = ẑ dx dy</text><text x="425.2" y="239.2" text-anchor="middle" fill="currentColor" style="font-size:12px;font-family:'STIX Two Math','Latin Modern Math','Cambria Math','Times New Roman',serif;font-style:italic;">dy</text><text x="484.8" y="179.8" text-anchor="start" fill="currentColor" style="font-size:12px;font-family:'STIX Two Math','Latin Modern Math','Cambria Math','Times New Roman',serif;font-style:italic;">dz</text><text x="389.9" y="208.6" text-anchor="middle" fill="currentColor" style="font-size:12px;font-family:'STIX Two Math','Latin Modern Math','Cambria Math','Times New Roman',serif;font-style:italic;">dx</text><text x="440.0" y="292.0" text-anchor="middle" fill="currentColor" style="font-size:13px;">closed surface: every dS points outward</text></svg><figcaption><strong>The two new objects of Lecture 1.</strong> Left: the differential length vector d<b>l</b> is tangent to the path and changes direction along it. Right: a differential surface vector d<b>S</b> is normal to the surface; for a closed surface the convention is <em>outward</em>. Its direction comes from the order of the cross product of two tangent vectors (swap the order and it flips).</figcaption></figure>

**Differential length vector** $d\mathbf{l}$: an infinitesimal step along a curve. It is *always tangent* to the curve, and in Cartesian coordinates

$$
d\mathbf{l} = \hat{x}\,dx + \hat{y}\,dy + \hat{z}\,dz .
$$

The line integral $\int_{P_1}^{P_2}\mathbf{F}\cdot d\mathbf{l}$ adds up the *tangential* component of $\mathbf{F}$ along the path — the work done by a force field on something moving from $P_1$ to $P_2$ (the "walking through a river" picture on the slides). Only the part of $\mathbf{F}$ along the step counts; a step perpendicular to the flow costs nothing. This is what we will plug into Faraday's law, and, in [[concepts/conservative-field|Lecture 5]], what defines voltage.

**Differential surface vector** $d\mathbf{S}$: an infinitesimal patch of surface, represented by a vector *normal* to it with magnitude equal to its area. Any patch, zoomed in far enough, looks flat, so we can build it from two tangent steps:

$$
d\mathbf{S} = d\mathbf{l}_1\times d\mathbf{l}_2 .
$$

> [!trap] The order of the cross product is the orientation
> $d\mathbf{l}_1\times d\mathbf{l}_2$ and $d\mathbf{l}_2\times d\mathbf{l}_1$ are opposite normals. For a **closed** surface the convention is always **outward**; for an open surface you must *choose* a normal (and then stick with it). In Cartesian coordinates the six faces of a box have
> $$
> d\mathbf{S} = \pm\hat{x}\,dy\,dz,\quad \pm\hat{y}\,dz\,dx,\quad \pm\hat{z}\,dx\,dy ,
> $$
> with the sign fixed by "outward". A wrong sign here silently flips the sign of a flux — see the [[problems/flux-through-a-plane-from-two-charges|flux-through-a-plane problem]].

The surface integral $\int_S \mathbf{F}\cdot d\mathbf{S}$ adds up the *normal* component of $\mathbf{F}$ over the surface: how much of the field "pokes through". We call it **flux**, and it is the subject of [[1-electrostatics/03-gauss-law-at-work|Lecture 3]]. Cylindrical and spherical versions of $d\mathbf{l}$, $d\mathbf{S}$, $dV$ are tabulated in [[0-toolkit/01-coordinates-and-differential-elements]].

## 6. Where this goes next

- **Lecture 2** defines $\mathbf{E}$ through force, states Coulomb's law, builds fields of charge distributions by superposition, and derives Gauss's law — the first Maxwell equation.
- **Lecture 3** turns Gauss's law into a calculational tool (line, sheet, slab) and introduces the δ-function description of point charges.
- **Lecture 4** converts the integral laws into the differential (divergence and curl) form in the table above.

> [!exam] What Lecture 1 contributes to Exam 1
> Nothing is asked "about Lecture 1" directly, but every problem uses it: vector components, the $\mathbf{R}_{12}$ convention, $d\mathbf{l}$ in a line integral (Exam 1, problem 1d), outward $d\mathbf{S}$ on a Gaussian surface (problems 2 and 4), and units on every answer. Marks are lost for missing units far more often than for missing physics.

### Sources for this page
Kudeki, *ECE 329 Lecture Notes*, Lecture 1 (big picture; frame-dependence of fields; quasi-statics). Shao, Lecture 1 slides (units table, $d\mathbf{l}$/$d\mathbf{S}$ construction, discussion problem). Rao, *Fundamentals of Electromagnetics for Electrical and Computer Engineering*, ch. 1.
