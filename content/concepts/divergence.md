---
title: "Divergence"
description: "Net outward flux per unit volume at a point: a scalar that is positive at sources, negative at sinks, zero where whatever flows in flows out."
tags: [concept, electrostatics, exam-1]
aliases: ["div", "∇·"]
---

> [!key] Definition and Cartesian formula
> $$
> \nabla\cdot\mathbf{A} \equiv \lim_{\Delta V\to0}\frac{\oint_{\partial\Delta V}\mathbf{A}\cdot d\mathbf{S}}{\Delta V}
> = \frac{\partial A_x}{\partial x}+\frac{\partial A_y}{\partial y}+\frac{\partial A_z}{\partial z}
> $$
> A **scalar** field. Cylindrical and spherical forms on the [[0-toolkit/02-vector-calculus-cheatsheet|cheat sheet]].

**Picture.** A water faucet: more flows out of a tiny box around the faucet than flows in. For a field whose arrows all point one way, divergence is nonzero where the field varies **along its own direction** — arrows pointing $+x$ that grow with $x$ — and zero if it only changes *across* its direction (that is [[concepts/curl|curl]]). When field lines spread or converge, geometry contributes as well: $Q\hat{r}/r^2$ weakens along its own direction yet has zero divergence off the charge. When in doubt, compute.

**Physics.** Gauss's law at a point: $\nabla\cdot\mathbf{D} = \rho$ — the divergence of $\mathbf{D}$ *is* the local charge density. $\nabla\cdot\mathbf{B} = 0$ — no magnetic sources. Continuity: $\nabla\cdot\mathbf{J} = -\partial\rho/\partial t$ — current flowing out of a point drains its charge.

**Examples.**
- $\mathbf{D} = 5x\hat{x}+12\hat{y}$: $\nabla\cdot\mathbf{D} = 5$ C/m³ everywhere ("a faucet of strength 5"). Constant components contribute nothing.
- $\mathbf{E} = E_0(y^2\hat{x}+2xy\hat{y}+3z^2\hat{z})$: $\nabla\cdot\mathbf{E} = E_0(0+2x+6z)$, so $\rho = \epsilon_0E_0(2x+6z)$ — only the "diagonal" derivatives $\partial_xE_x,\partial_yE_y,\partial_zE_z$ enter; the $y^2$ in $E_x$ is irrelevant.
- Point-charge field $\hat{r}Q/(4\pi\epsilon_0r^2)$: divergence zero for $r>0$ (all of it is concentrated in a δ at the origin).

> [!trap]
> - Differentiate each component with respect to **its own** coordinate only.
> - Don't forget $\epsilon_0$ (or $\epsilon$) when converting $\nabla\cdot\mathbf{E}$ to $\rho$.
> - "Curl-free" does **not** imply "divergence-free": a field can be conservative and still have charge as its source (the [[problems/curl-potential-and-charge-from-a-field|worked problem]] has both).

**Theorem.** [[concepts/divergence-theorem]]: $\int_V\nabla\cdot\mathbf{A}\,dV = \oint_S\mathbf{A}\cdot d\mathbf{S}$.

**Where it appears.** [[1-electrostatics/04-divergence-and-curl#2-divergence-flux-per-unit-volume|Lecture 4 §2]], [[problems/curl-potential-and-charge-from-a-field]], [[1-electrostatics/07-poisson-and-laplace|Lecture 7]] (Poisson: $\nabla\cdot\nabla V = \nabla^2V = -\rho/\epsilon$).
