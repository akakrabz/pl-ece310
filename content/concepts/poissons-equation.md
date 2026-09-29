---
title: "Poisson's and Laplace's equations"
description: "∇²V = −ρ/ε is Gauss's law written for the potential, valid where ε is constant; ∇²V = 0 in charge-free regions. One-dimensional solutions are linear, logarithmic or 1/r; the general solution is the Coulomb superposition integral."
tags: [concept, electrostatics, exam-1]
aliases: ["Laplace's equation", "Laplacian", "Poisson equation"]
---

> [!key] Definition
> $$
> \nabla^2V = -\frac{\rho}{\epsilon}\ \ (\text{Poisson}),\qquad \nabla^2V = 0\ \ (\text{Laplace, where }\rho = 0),\qquad
> \nabla^2V \equiv \nabla\cdot\nabla V = \frac{\partial^2V}{\partial x^2}+\frac{\partial^2V}{\partial y^2}+\frac{\partial^2V}{\partial z^2}.
> $$
> Obtained from $\nabla\cdot\mathbf{D} = \rho$, $\mathbf{D} = \epsilon\mathbf{E}$, $\mathbf{E} = -\nabla V$, pulling $\epsilon$ out of the divergence — so it requires a **static** field and a **homogeneous** medium ($\epsilon$ independent of position). Otherwise use $\nabla\cdot(\epsilon\nabla V) = -\rho$, or Gauss's law directly. The solution with $V$ prescribed on the boundary of a region is unique, and solutions superpose.

**Meaning.** $\nabla^2V$ is proportional to (average of $V$ on a small sphere) − (value at its centre). Laplace's equation says $V$ has no local maxima or minima in empty space; positive charge makes $V$ curve downward like a dome ($\nabla^2V<0$), negative charge upward like a bowl.

**One-dimensional solutions** (two constants from two boundary values):

| symmetry | equation | solution |
|---|---|---|
| planar $V(z)$ | $V'' = 0$ | $Az+B$ |
| cylindrical $V(r)$ | $\frac1r(rV')' = 0$ | $A\ln r+B$ |
| spherical $V(r)$ | $\frac1{r^2}(r^2V')' = 0$ | $A/r+B$ |

Then $\mathbf{E} = -\nabla V$, $\mathbf{D} = \epsilon\mathbf{E}$, and $\rho_s = \hat{n}\cdot\mathbf{D}$ on each conductor ($\hat{n}$ out of the conductor). With charge present, integrate $V'' = -\rho/\epsilon$ twice and match $V$ (and, where there is no surface charge, $\epsilon V'$) at interfaces: e.g. the pn junction's $V_{21} = \rho_2W_2(W_1+W_2)/2\epsilon_0$.

**General solution** (finite charge distribution, $V\to0$ at infinity): $V(\mathbf{r}) = \displaystyle\int\frac{\rho(\mathbf{r}')\,d^3\mathbf{r}'}{4\pi\epsilon_0\lvert\mathbf{r}-\mathbf{r}'\rvert}$ — the point-charge potential is the impulse response of a linear, shift-invariant system, and this is its superposition integral (Green's function). It also proves the Helmholtz theorem: a field vanishing at infinity is fixed by its divergence and curl.

**Examples.** Plates at $V(0)=0$, $V(2\,\text{m}) = -3$ V: $V = -\tfrac32z$, $\mathbf{E} = \tfrac32\hat{z}$ V/m, plate charges $\pm\tfrac32\epsilon_0$ C/m². Coax with $V(a)=V_0$, $V(b)=0$: $V = V_0\ln(b/r)/\ln(b/a)$.

> [!trap]
> - Not across a dielectric interface (ε jumps) and not in a graded medium: solve each homogeneous layer separately and match, or use $\nabla\cdot\mathbf{D} = 0$ ⇒ $D_n$ constant.
> - Not through a charged sheet in the middle of a gap: two Laplace regions joined by the $D_n$ jump.
> - The exam often runs it backwards: $\rho = \epsilon\nabla\cdot\mathbf{E} = -\epsilon\nabla^2V$. Zero curl says nothing about $\rho$.

**Where it appears.** [[1-electrostatics/07-poisson-and-laplace|Lecture 7]], [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]] (piecewise, and where it fails), [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]] (the given-$V$ route to capacitance), [[problems/curl-potential-and-charge-from-a-field]] (the backwards use).

Related: [[concepts/electrostatic-potential]] · [[concepts/gauss-law]] · [[concepts/divergence]] · [[concepts/boundary-conditions]].
