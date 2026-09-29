---
title: "Electrostatic energy"
description: "A capacitor stores W = ½CV² = Q²/2C; the energy lives in the field at a density w = ½εE² = ½D·E joules per cubic metre. Potential energy of a charge is qV."
tags: [concept, electrostatics]
aliases: ["stored energy", "energy density", "field energy"]
---

> [!key] Stored energy
> $$
> W = \int_0^Q\frac{q}{C}\,dq = \frac{Q^2}{2C} = \tfrac12CV^2 = \tfrac12QV\ [\text{J}],\qquad\qquad
> w = \tfrac12\epsilon E^2 = \tfrac12\mathbf{D}\cdot\mathbf{E}\ [\text{J/m}^3].
> $$
> The first is the work done to move charge across a growing voltage; the second follows for a parallel plate from $W = \tfrac12(\epsilon A/d)(Ed)^2 = \tfrac12\epsilon E^2\cdot(Ad)$ and is taken as the field's energy density in general: $W = \int w\,dV$ over all space. The potential energy of a point charge $q$ at potential $V$ is $U = qV$.

**Physics.** Energy is stored in the field, wherever the field is — a coax stores $W' = \tfrac12\mathcal{C}V^2$ per metre in its dielectric, most of it near the inner conductor where $E$ is largest. The same $W$ is also $\int_0^t VI\,dt$, since $VI = VC\,dV/dt = \frac{d}{dt}(\tfrac12CV^2)$. A battery at fixed $V$ supplies $QV = CV^2$ while charging a capacitor from empty; half is stored, half is lost in the charging resistance no matter how small it is.

**Examples.** Capacitor $C$ at $V$: $\tfrac12CV^2$. Two-layer coax with $\pm\lambda$: $W' = \lambda^2/2\mathcal{C}$. Field between plates $E = 3$ V/m in vacuum: $w = \tfrac12\epsilon_0(9)\approx4\times10^{-11}$ J/m³.

> [!trap]
> - Joules per cubic metre, not watts (a slide slips here). $W$ is energy; $P = dW/dt$ is power.
> - With a dielectric, $\epsilon$ replaces $\epsilon_0$ in $w$; $\tfrac12\mathbf{D}\cdot\mathbf{E}$ is the form that needs no thought.
> - $\tfrac12QV$ uses the total charge and the *final* voltage; it is half of $QV$ because the voltage grew from zero during charging.

**Where it appears.** [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]], [[1-electrostatics/05-curl-free-fields-and-the-electrostatic-potential|Lecture 5]] ($U = qV$, the electron-volt); Unit 3 (Poynting's theorem adds the magnetic $\tfrac12\mu H^2$ and the flow of energy).

Related: [[concepts/capacitance]] · [[concepts/electrostatic-potential]].
