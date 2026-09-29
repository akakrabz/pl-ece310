---
title: "Conservative (curl-free) fields and the potential"
description: "Four equivalent statements — zero curl, zero circulation, path-independent line integrals, and the existence of a scalar potential — and why every static E field satisfies them."
tags: [concept, electrostatics, exam-1]
aliases: ["curl-free field", "irrotational field", "path independence", "zero circulation"]
---

> [!key] Four equivalent statements about a vector field $\mathbf{E}$ (on all of space, or any region without holes)
> 1. $\nabla\times\mathbf{E} = 0$ everywhere.
> 2. $\oint_C\mathbf{E}\cdot d\mathbf{l} = 0$ around every closed loop.
> 3. $\int_A^B\mathbf{E}\cdot d\mathbf{l}$ depends only on the endpoints $A$, $B$, not on the path.
> 4. There is a scalar function $V$ with $\mathbf{E} = -\nabla V$.
>
> (1)⟹(2) is Stokes' theorem; (2)⟺(3) is "two paths from $A$ to $B$ make a loop"; (3)⟹(4) is the definition $V(\mathbf{r}) = -\int_{\text{ref}}^{\mathbf{r}}\mathbf{E}\cdot d\mathbf{l}$; (4)⟹(1) is the identity $\nabla\times\nabla V = 0$.

**Why static E is conservative.** The Coulomb field of a point charge is radial, and every radial field $g(r)\mathbf{r}$ has zero curl; curl is linear, so any superposition of Coulomb fields — every electrostatic field — is curl-free ([[1-electrostatics/04-divergence-and-curl#6-static-electric-fields-are-curl-free|Lecture 4 §6]]). Time-varying $\mathbf{B}$ breaks this (Faraday, Lecture 14): the induced part of $\mathbf{E}$ has curl, and "voltage" then depends on the path.

**The potential** ([[1-electrostatics/05-curl-free-fields-and-the-electrostatic-potential|Lecture 5]]; full page: [[concepts/electrostatic-potential]]). With the minus-sign convention,
$$
V(A)-V(B) = \int_A^B\mathbf{E}\cdot d\mathbf{l},\qquad \mathbf{E} = -\nabla V,\qquad dV = -\mathbf{E}\cdot d\mathbf{l},
$$
$\mathbf{E}$ points from high $V$ to low $V$, and $V$ is defined up to a constant fixed by a reference point ($V=0$ at infinity when the charge is finite; at a chosen point otherwise).

> [!recipe] Finding V from a given E (the exam pattern)
> 1. Check $\nabla\times\mathbf{E} = 0$ (otherwise no $V$ exists).
> 2. Either spot the exact differential — $\mathbf{E}\cdot d\mathbf{l} = E_0(y^2\,dx+2xy\,dy+3z^2dz) = E_0\,d(xy^2+z^3)$, so $V = -E_0(xy^2+z^3)+\text{const}$ — or integrate along an axis-parallel path from the reference point: $(0,0,0)\to(x,0,0)\to(x,y,0)\to(x,y,z)$, using the *current* values of the coordinates already traversed.
> 3. Fix the constant from the reference condition.
> 4. Check: $-\nabla V$ must return $\mathbf{E}$.

> [!trap] Do not antidifferentiate each component separately and add
> $\int y^2\,dx+\int 2xy\,dy+\int 3z^2\,dz$ counts the $xy^2$ term twice. The line integral is one integral along one path.

**Where it appears.** [[1-electrostatics/04-divergence-and-curl#5-conservative-fields|Lecture 4 §5]], [[1-electrostatics/05-curl-free-fields-and-the-electrostatic-potential|Lecture 5]] (potential), [[1-electrostatics/06-circulation-and-boundary-conditions|Lecture 6]] (Kirchhoff's voltage law is statement 2), [[problems/curl-potential-and-charge-from-a-field]].

Related: [[concepts/curl]] · [[concepts/stokes-theorem]] · [[concepts/electrostatic-potential]].
