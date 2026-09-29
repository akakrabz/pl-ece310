---
title: "Gauss's law"
description: "The flux of D out of any closed surface equals the charge enclosed — always true, and a two-line calculation of E whenever the symmetry is spherical, cylindrical or planar."
tags: [concept, electrostatics, exam-1]
aliases: ["Gauss law", "Gauss's law for E", "Maxwell equation 1"]
---

> [!key] Integral and differential forms
> $$
> \oint_S\mathbf{D}\cdot d\mathbf{S} = Q_{\text{enc}} = \int_V\rho\,dV
> \qquad\Longleftrightarrow\qquad
> \nabla\cdot\mathbf{D} = \rho
> $$
> $S$ any closed surface, $d\mathbf{S}$ outward, $V$ the volume it bounds. Flux out = charge inside. Charges outside $S$ affect $\mathbf{D}$ on $S$ but contribute zero *net* flux.

**Where it comes from.** For a point charge, $|\mathbf{E}|\propto 1/r^2$ while a sphere's area $\propto r^2$: the flux $\epsilon_0 E\cdot 4\pi r^2 = Q$ is independent of $r$. Field lines never stop, so the same count holds for any enclosing surface; superposition extends it to any collection of charges ([[1-electrostatics/02-coulombs-law-superposition-and-gauss#5-from-coulomb-to-gauss|Lecture 2 §5]]). It is derived from Coulomb's law but is *more* general: it holds for moving charges and time-varying fields, where Coulomb's law does not.

**Always true, sometimes useful.** As a *calculation tool* it works only when symmetry makes $|\mathbf{D}|$ constant on a surface you can name:

| symmetry | source | Gaussian surface | result |
|---|---|---|---|
| spherical | point, uniformly charged sphere/shell | concentric sphere | $D = Q_{\text{enc}}/(4\pi r^2)$ |
| cylindrical | infinite line, cylinder, coax | coaxial cylinder (caps give 0) | $D = \rho_{l,\text{enc}}/(2\pi r)$ |
| planar | infinite sheet or slab, charge symmetric about its mid-plane | pillbox straddling the symmetry plane (sides give 0; both caps count) | $2D = \rho_{s,\text{enc}}$ |
| planar, *not* symmetric | parallel plates $\pm\rho_s$, layered sheets | superpose sheet fields, or a one-cap pillbox with its other cap where $E=0$ | plates: $D = \rho_s$ between, $0$ outside |

> [!recipe] The four lines
> 1. State the symmetry: direction of $\mathbf{D}$ and what it depends on.
> 2. Choose the surface so each piece is "$\mathbf{D}\parallel d\mathbf{S}$, constant" or "$\mathbf{D}\perp d\mathbf{S}$".
> 3. Write *top + bottom + side* $= Q_{\text{enc}}$; cancel the arbitrary $L$ or $A$.
> 4. $\mathbf{E} = \mathbf{D}/\epsilon$; attach direction and units.

> [!trap]
> - Both caps of a pillbox carry flux (factor 2) unless one cap sits where $E=0$.
> - $Q_{\text{enc}}$ is an *integral* when $\rho$ varies inside the surface.
> - Gauss's law gives $\mathbf{D}$, which does not care about the material; $\mathbf{E} = \mathbf{D}/\epsilon$ does.
> - A surface with *less* symmetry than the source (a cube around a point charge) still satisfies the law — it just doesn't let you pull $D$ out of the integral.

**Differential form.** Shrink the surface to a point: $\nabla\cdot\mathbf{D} = \rho$ ([[1-electrostatics/04-divergence-and-curl#2-divergence-flux-per-unit-volume|Lecture 4 §2]]). Given a field, differentiate to find the charge that made it (FA26 problem 1c).

**Magnetic twin.** $\oint_S\mathbf{B}\cdot d\mathbf{S} = 0$, $\nabla\cdot\mathbf{B} = 0$: no magnetic charge.

**Where it appears.** [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]] (derived), [[1-electrostatics/03-gauss-law-at-work|Lecture 3]] (used), [[1-electrostatics/04-divergence-and-curl|Lecture 4]] (differential form), Lectures 8–10 (conductors, dielectrics, capacitance), [[problems/charged-slab-with-a-power-law-profile]].

Related: [[concepts/flux]] · [[concepts/electric-flux-density]] · [[concepts/divergence]] · [[concepts/superposition]] · [[demos/point-charges-and-gauss]].
