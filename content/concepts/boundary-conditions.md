---
title: "Boundary conditions at an interface"
description: "Normal D jumps by the free surface charge, tangential E is continuous, normal B is continuous, tangential H jumps by the surface current — with n̂ pointing from medium 2 into medium 1. Special cases: conductor surfaces and dielectric interfaces."
tags: [concept, electrostatics, exam-1]
aliases: ["jump conditions", "interface conditions", "normal and tangential components"]
---

> [!key] The four conditions ($\hat{n}$ points from medium 2 into medium 1)
> $$
> \hat{n}\cdot(\mathbf{D}_1-\mathbf{D}_2) = \rho_s,\qquad
> \hat{n}\times(\mathbf{E}_1-\mathbf{E}_2) = 0,\qquad
> \hat{n}\cdot(\mathbf{B}_1-\mathbf{B}_2) = 0,\qquad
> \hat{n}\times(\mathbf{H}_1-\mathbf{H}_2) = \mathbf{J}_s .
> $$
> In components: $D_{1n} - D_{2n} = \rho_s$ [C/m², **free** charge]; $E_{1t} = E_{2t}$; $B_{1n} = B_{2n}$; $H_{1t} - H_{2t} = J_s$ [A/m]. They come from the integral Maxwell equations applied to a pillbox (divergence laws → normal components) and a thin loop (curl laws → tangential components) that straddle the interface and are squashed onto it. They hold for time-varying fields too, because the flux terms vanish with the area. The course notes write the same thing as $\hat{n}\cdot(\mathbf{D}^+-\mathbf{D}^-) = \rho_s$ with "$+$" the side $\hat{n}$ points into.

**What they do not say.** Nothing about tangential $\mathbf{D}$ or normal $\mathbf{E}$ (both jump at a dielectric interface), nor about normal $\mathbf{H}$ or tangential $\mathbf{B}$. And they fix *differences* only; the level on one side must come from elsewhere (a measurement, symmetry, superposition, a conductor's potential).

**Mnemonic.** In each differential Maxwell equation replace $\nabla\to\hat{n}$, field $\to$ (field$_1$ − field$_2$), $\rho\to\rho_s$, $\mathbf{J}\to\mathbf{J}_s$, $\partial/\partial t\to0$.

**Special cases.**
- *Conductor surface* (field zero inside, $\hat{n}$ out of the metal): $\hat{n}\cdot\mathbf{D} = \rho_s$, $\hat{n}\times\mathbf{E} = 0$ — the field leaves perpendicularly and its normal $\mathbf{D}$ *is* the surface charge. Used to read plate charges off a field ([[1-electrostatics/07-poisson-and-laplace|Lecture 7]], [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]]).
- *Two perfect dielectrics, no free charge:* $D_{1n} = D_{2n}$ and $E_{1t} = E_{2t}$, so $E_{2n} = (\epsilon_1/\epsilon_2)E_{1n}$ and field lines refract with $\tan\theta_1/\tan\theta_2 = \epsilon_1/\epsilon_2$ ([[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]]). The bound surface charge left there is $\rho_{sb} = P_{2n} - P_{1n} = \epsilon_0(E_{1n}-E_{2n})$; in general $\hat{n}\cdot(\mathbf{P}_1-\mathbf{P}_2) = -\rho_{sb}$.
- *Planar problems:* a sheet $\rho_s$ makes $D_n$ step by $\rho_s$ as you cross it in the direction of $\hat{n}$; tangential components pass through unchanged (in vacuum on both sides, tangential $\mathbf{D}$ too).

**Examples.** $\mathbf{D} = 0$ for $x<0$, sheets $+2$ at $x=0$ and $-2$ C/m² at $x=5$ m ⇒ $\mathbf{D} = 2\hat{x}$ between, $0$ outside. $\epsilon_1 = \epsilon_0$ over $\epsilon_2 = 2\epsilon_0$ with $\rho_s = 3$ C/m² and $\mathbf{D}_2 = 3\hat{x}+2\hat{y}$ ($\hat{y}$ normal) ⇒ $\mathbf{D}_1 = \tfrac32\hat{x} + 5\hat{y}$: tangential $\mathbf{D}$ scaled by $\epsilon_1/\epsilon_2$, normal $\mathbf{D}$ raised by $\rho_s$.

> [!trap]
> - $E_t$ and $D_n$ are the continuous ones — not $E_n$, not $D_t$.
> - Draw $\hat{n}$ at every interface before writing anything; it flips between the two faces of a slab.
> - $\rho_s$ is free charge. Bound charge is inside $\mathbf{P}$, hence inside $\mathbf{D}$.
> - A stray $\epsilon_0$ in a jump condition means you wrote $\mathbf{E}$ where $\mathbf{D}$ belongs.

Slips in the source materials around these conditions (the missing minus sign in the $\mathbf{P}$ condition, $D_1-D_2=\rho$) are listed on the [[0-toolkit/04-errata-in-the-course-materials|errata page]].

**Where it appears.** [[1-electrostatics/06-circulation-and-boundary-conditions|Lecture 6]] (derived), [[1-electrostatics/07-poisson-and-laplace|Lecture 7]] and [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]] (conductors), [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]] (dielectrics), [[problems/fields-across-a-dielectric-interface]], [[problems/two-layer-coaxial-capacitor]]; later, reflection of waves at interfaces (Unit 3) is these same four conditions at work.

Related: [[concepts/gauss-law]] · [[concepts/electric-flux-density]] · [[concepts/conductors]] · [[concepts/permittivity]].
