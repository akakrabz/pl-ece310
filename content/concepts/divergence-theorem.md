---
title: "Divergence theorem"
description: "Flux out of a closed surface equals the volume integral of the divergence inside — the bridge between the integral and differential forms of Gauss's law."
tags: [concept, toolkit, exam-1]
aliases: ["Gauss's theorem", "Gauss–Ostrogradsky theorem"]
---

> [!key] Statement
> $$
> \oint_{S}\mathbf{A}\cdot d\mathbf{S} = \int_V\nabla\cdot\mathbf{A}\,dV,\qquad S = \partial V,\ d\mathbf{S}\ \text{outward}
> $$

**Why it is true.** Tile $V$ with tiny boxes. The flux through a face shared by two neighbours is counted once as "out of" each — with opposite signs — so it cancels. Adding up "flux per unit volume × volume" over all boxes therefore leaves only the outer surface. (Same logic as the fundamental theorem of calculus: interior contributions telescope, the boundary survives.)

**What it is for.**
- Converts Gauss's law between forms: $\oint\mathbf{D}\cdot d\mathbf{S} = Q_{\text{enc}}$ for *every* $V$ ⟺ $\nabla\cdot\mathbf{D} = \rho$ pointwise. Same for $\nabla\cdot\mathbf{B} = 0$.
- Interprets divergence as *flux per unit volume*.
- Gives the continuity equation its integral form: $\oint_S\mathbf{J}\cdot d\mathbf{S} = -\dfrac{d}{dt}\int_V\rho\,dV$ — current out = rate of loss of charge inside.

> [!tip] Check a divergence by a flux
> For the slab field $E_x = \rho x/\epsilon_0$, a pillbox of caps at $\pm x$ has flux $2\epsilon_0E_xA = 2\rho xA$ = enclosed charge. Differentiating instead: $\epsilon_0\,dE_x/dx = \rho$. Same content, two languages.

**Where it appears.** [[1-electrostatics/04-divergence-and-curl#2-divergence-flux-per-unit-volume|Lecture 4 §2]]. Its curl-side twin is [[concepts/stokes-theorem]].
