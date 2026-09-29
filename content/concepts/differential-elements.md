---
title: "Differential elements dl, dS, dV"
description: "The tiny vectors every line and surface integral is built from; tangent for dl, normal (outward on closed surfaces) for dS, and the cross-product rule that sets orientation."
tags: [concept, toolkit]
aliases: ["dl", "dS", "differential length", "differential surface"]
---

> [!key] Definitions
> - $d\mathbf{l}$: an infinitesimal step along a curve; **tangent** to it. Cartesian: $d\mathbf{l} = \hat{x}dx+\hat{y}dy+\hat{z}dz$.
> - $d\mathbf{S}$: an infinitesimal patch of surface; **normal** to it, magnitude = area. Built from two tangent steps, $d\mathbf{S} = d\mathbf{l}_1\times d\mathbf{l}_2$; the order fixes the orientation. On a **closed** surface, $d\mathbf{S}$ points **outward**.
> - $dV$: an infinitesimal volume, $dx\,dy\,dz$ in Cartesian.

**What the integrals mean.**
- $\int_C\mathbf{F}\cdot d\mathbf{l}$ adds up the *tangential* component of $\mathbf{F}$ along $C$: work, EMF, voltage drop. A closed-loop version $\oint_C$ is a **circulation**.
- $\int_S\mathbf{F}\cdot d\mathbf{S}$ adds up the *normal* component through $S$: **flux**. Closed version $\oint_S$ = flux out.
- $\int_V f\,dV$ adds up a density: total charge from $\rho$.

**Orientation rules.**
- Closed surface: outward. (Wrong sign here flips a flux — see [[problems/flux-through-a-plane-from-two-charges]].)
- Open surface with boundary curve $C$: right-hand rule — fingers along $d\mathbf{l}$ of $C$, thumb along $d\mathbf{S}$. This links Stokes' theorem's two sides.
- Path integrals: $d\mathbf{l}$ follows the direction of travel from the start point to the end point; reversing the path negates the integral.

**In other coordinate systems** $d\mathbf{l}$ picks up metric factors ($r\,d\phi$, $r\sin\theta\,d\phi$) and $d\mathbf{S}$ on a sphere is $\hat{r}\,r^2\sin\theta\,d\theta\,d\phi$ — tabulated in [[0-toolkit/01-coordinates-and-differential-elements]].

**Where it appears.** [[1-electrostatics/01-fields-forces-and-the-maxwell-roadmap#5-the-two-new-tools-dl-and-ds|Lecture 1 §5]] introduces them; [[concepts/flux]], [[concepts/gauss-law]], [[concepts/stokes-theorem]] use them.
