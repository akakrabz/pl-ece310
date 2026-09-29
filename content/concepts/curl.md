---
title: "Curl"
description: "Circulation per unit area at a point: a vector along the axis a tiny paddlewheel would spin about, nonzero where a field varies across its own direction."
tags: [concept, electrostatics, exam-1]
aliases: ["∇×", "rotation of a field"]
---

> [!key] Definition and Cartesian formula
> $$
> (\nabla\times\mathbf{A})\cdot\hat{n} \equiv \lim_{\Delta S\to0}\frac{\oint_C\mathbf{A}\cdot d\mathbf{l}}{\Delta S}
> \quad(\text{$C$ bounds $\Delta S$, right-hand rule})
> $$
> $$
> \nabla\times\mathbf{A} = \begin{vmatrix}\hat{x}&\hat{y}&\hat{z}\\[2pt]\partial_x&\partial_y&\partial_z\\[2pt]A_x&A_y&A_z\end{vmatrix}
> = \hat{x}(\partial_yA_z-\partial_zA_y)+\hat{y}(\partial_zA_x-\partial_xA_z)+\hat{z}(\partial_xA_y-\partial_yA_x)
> $$
> A **vector** field. Cylindrical and spherical forms on the [[0-toolkit/02-vector-calculus-cheatsheet|cheat sheet]].

**Picture.** A paddlewheel spins where the flow is faster on one side than the other: for a unidirectional field, curl is nonzero where the field varies **across its own direction** (arrows pointing $+y$ that grow with $x$). A uniform stream has none; a shear flow does, even though nothing "goes around". Curved flow adds a geometric part — rigid rotation has curl; the free vortex $\hat{\phi}/r$ (a wire's magnetic field) has none off its axis. The spin axis (right-hand rule) is the direction of $\nabla\times\mathbf{A}$. Equivalently: the line integral between two points depends on the path (1→2→3 ≠ 1→4→3 around a small rectangle).

**Physics.** $\nabla\times\mathbf{E} = -\partial\mathbf{B}/\partial t$ (Faraday) and $\nabla\times\mathbf{H} = \mathbf{J}+\partial\mathbf{D}/\partial t$ (Ampère–Maxwell): time-varying magnetic fields and currents make fields circulate. In **statics** $\nabla\times\mathbf{E} = 0$: every Coulomb field is curl-free, so by superposition every electrostatic field is → a potential exists ([[concepts/conservative-field]]).

**Examples.**
- $\mathbf{v} = -3y\hat{x}+2x\hat{y}$: $\nabla\times\mathbf{v} = \hat{z}(2-(-3)) = 5\hat{z}$, a uniform whirlpool.
- $\mathbf{E} = E_0(y^2\hat{x}+2xy\hat{y}+3z^2\hat{z})$: $x$- and $y$-components vanish by inspection (no $z$-dependence in $E_x,E_y$; $E_z$ depends only on $z$); $\hat{z}$ component $\partial_x(2xy)-\partial_y(y^2) = 2y-2y = 0$. Curl-free ⇒ a potential $V = -E_0(xy^2+z^3)$ exists.
- Any radial field $g(r)\,\mathbf{r}$ is curl-free ([[1-electrostatics/04-divergence-and-curl#6-static-electric-fields-are-curl-free|Lecture 4 §6]]).

> [!trap] Expanding the determinant
> The $\hat{y}$ term carries a minus sign when you expand along the top row, and the $\hat{z}$ term is $\partial_xA_y-\partial_yA_x$ in that order. Reversing that bracket is a common slip — it turns up even in official solution keys — and it goes unnoticed whenever the two partials happen to be equal. Expand mechanically, then sanity-check one component with the "across the flow" picture.

**Theorem.** [[concepts/stokes-theorem]]: $\oint_C\mathbf{A}\cdot d\mathbf{l} = \int_S(\nabla\times\mathbf{A})\cdot d\mathbf{S}$.

**Identities.** $\nabla\times\nabla f = 0$; $\nabla\cdot(\nabla\times\mathbf{A}) = 0$; $\nabla\times\nabla\times\mathbf{A} = \nabla(\nabla\cdot\mathbf{A})-\nabla^2\mathbf{A}$ (with $\nabla^2\mathbf{A}$ taken component-wise in Cartesian coordinates).

**Where it appears.** [[1-electrostatics/04-divergence-and-curl#3-curl-circulation-per-unit-area|Lecture 4 §3]], [[problems/curl-potential-and-charge-from-a-field]], [[1-electrostatics/06-circulation-and-boundary-conditions|Lecture 6]] (circulation and boundary conditions), Lectures 12–14.
