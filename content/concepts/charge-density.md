---
title: "Charge densities and δ-functions"
description: "ρ, ρ_s, ρ_l and point charges; how to write all of them as a volume density with impulses; and the unit checks that go with each."
tags: [concept, electrostatics, exam-1]
aliases: ["charge density", "delta function charge"]
---

> [!key] The four kinds of charge, one formula for the total
> | distribution | density | total charge |
> |---|---|---|
> | point | $Q$ [C] | $Q$ |
> | line | $\rho_l$ [C/m] (also written $\lambda$) | $\int_L\rho_l\,dl$ |
> | surface | $\rho_s$ [C/m²] | $\int_S\rho_s\,dS$ |
> | volume | $\rho$ [C/m³] | $\int_V\rho\,dV$ |
>
> Each is "charge per unit of the thing it lives on"; the differential element $dQ$ is density × element.

**Impulse (δ) notation** makes every case a volume density, so that $Q_{\text{enc}} = \int_V\rho\,dV$ and $\nabla\cdot\mathbf{D}=\rho$ cover all of them:

$$
\int_{-\infty}^{\infty}\delta(x)\,dx = 1,\qquad [\delta(x)] = \text{m}^{-1},
$$

| charge | $\rho(x,y,z)$ |
|---|---|
| point $Q$ at $(x_0,y_0,z_0)$ | $Q\,\delta(x-x_0)\,\delta(y-y_0)\,\delta(z-z_0)$ |
| line $\rho_l$ along $z$ through $(x_0,y_0)$ | $\rho_l\,\delta(x-x_0)\,\delta(y-y_0)$ |
| sheet $\rho_s(y,z)$ on $x=x_0$ | $\rho_s(y,z)\,\delta(x-x_0)$ |

Unit check for the point charge: C·m⁻³ ✓. A row of point charges $Q$ at $z = n\Delta z$ is $\rho = Q\sum_n\delta(x)\delta(y)\delta(z-n\Delta z)$; on scales coarser than $\Delta z$ it acts like a line charge $\rho_l = Q/\Delta z$ — this is how the course notes motivate the macroscopic line charge.

**Converting between kinds.** A slab of volume density $\rho(x)$ seen from far away is a sheet with $\rho_s = \int\rho(x)\,dx$ (FA26 Exam 1, problem 2a); a thin cylinder of $\rho$ is a line with $\rho_l = \rho\cdot\pi a^2$. Conversely, a "surface charge" is really a thin layer of large $\rho$.

> [!trap] Signs and absolute values
> A density such as $\rho(x) = \rho_0|x|/a$ is *even*: its total is $2\int_0^a$, not zero. Dropping the absolute value (which makes the integral vanish) or the factor 2 is the classic error on slab problems like the [[problems/charged-slab-with-a-power-law-profile|worked one]].

**Where it appears.** [[1-electrostatics/01-fields-forces-and-the-maxwell-roadmap#3-units-and-the-source--field-bookkeeping|Lecture 1]] (units), [[1-electrostatics/03-gauss-law-at-work#5-charge-densities-and-δ-functions|Lecture 3]] (δ formalism), [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]] (bound charge $\rho_b = -\nabla\cdot\mathbf{P}$, surface charge on conductors $\rho_s = \hat{n}\cdot\mathbf{D}$).
