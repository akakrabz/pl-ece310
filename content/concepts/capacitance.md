---
title: "Capacitance"
description: "C = Q/V for two conductors carrying ±Q, a property of geometry and dielectric alone. Computed by a chain (given V: Laplace → E → D → ρs → Q; given Q: Gauss → D → E → V). Plates εA/d, coax 2πε/ln(b/a) per length, spheres 4πε ab/(b − a); layers in series or parallel."
tags: [concept, electrostatics, exam-1]
aliases: ["capacitor", "per-unit-length capacitance", "series capacitors", "parallel capacitors"]
---

> [!key] Definition and the three geometries
> $$
> Q = CV,\qquad C = \frac{Q}{V}\ [\text{F}],\quad V = V_+ - V_- > 0 .
> $$
> | geometry | $C$ | condition |
> |---|---|---|
> | parallel plates, area $A$, gap $d$ | $\dfrac{\epsilon A}{d}$ | $d\ll\sqrt A$ |
> | coax, radii $a<b$, length $L$ | $\dfrac{2\pi\epsilon L}{\ln(b/a)}$; per unit length $\mathcal{C} = \dfrac{2\pi\epsilon}{\ln(b/a)}$ [F/m] | $L\gg b$ |
> | concentric spheres $a<b$ | $4\pi\epsilon\dfrac{ab}{b-a}$; isolated sphere $4\pi\epsilon a$ | — |

**How to compute it.** Never from memory alone — from a chain. *Given $V$:* Laplace's equation with the conductors' potentials → $\mathbf{E} = -\nabla V$ → $\mathbf{D} = \epsilon\mathbf{E}$ → $\rho_s = \hat{n}\cdot\mathbf{D}$ on one conductor → $Q = \int\rho_s\,dA$ → $C = Q/V$. *Given $Q$:* Gauss's law → $\mathbf{D}$ → $\mathbf{E} = \mathbf{D}/\epsilon$ → $V = -\int_-^+\mathbf{E}\cdot d\mathbf{l}$ (from the negative to the positive conductor, a positive number) → $C = Q/V$. Both routes pass through the same field pattern.

**Layers.** When the field crosses layers one after another (interfaces normal to $\mathbf{E}$), $\mathbf{D}$ is common and the voltages add: **series**, $1/C = \sum1/C_i$; e.g. plates with $d_1,\epsilon_1$ and $d_2,\epsilon_2$: $1/C = d_1/\epsilon_1A + d_2/\epsilon_2A$; a coax with $\epsilon_1$ for $a<r<c$ and $\epsilon_2$ for $c<r<b$: $1/\mathcal{C} = \ln(c/a)/2\pi\epsilon_1 + \ln(b/c)/2\pi\epsilon_2$. When layers sit side by side along the field, $\mathbf{E}$ is common and the charges add: **parallel**, $C = \sum C_i$; e.g. two dielectrics each filling half the angle of a coax: $\mathcal{C} = \pi(\epsilon_1+\epsilon_2)/\ln(b/a)$.

**What depends on what.** Inserting a dielectric of $\epsilon_r$ multiplies $C$ by $\epsilon_r$. With $Q$ fixed, $\mathbf{D}$ stays and $\mathbf{E}, V$ fall by $\epsilon_r$; with $V$ fixed, $\mathbf{E}$ stays and $\mathbf{D}, \rho_s, Q$ rise by $\epsilon_r$. Stored energy $W = \tfrac12CV^2 = Q^2/2C$ ([[concepts/electrostatic-energy]]); a lossy filling adds $G = (\sigma/\epsilon)C$ in parallel ([[concepts/conductance]]). A pn junction has $Q\propto\sqrt V$ and a small-signal $C = dQ/dV = \epsilon_0A/(W_1+W_2)$ that falls as $V^{-1/2}$.

**Two conditions.** No fringing ($d\ll\sqrt A$: the field is that of infinite plates) and quasi-static ($\sqrt A\ll c/f$: otherwise it is a transmission line).

> [!trap]
> - A negative $C$ means you used $V(-)-V(+)$; $V$ is the drop from the positive conductor to the negative one.
> - F versus F/m: a coax "capacitance" without a length is per unit length, $\mathcal{C}$.
> - Layers crossed in succession are in series, not parallel — the voltages add, the charge is common.
> - The outer conductor of a coax carries $-\lambda$ per unit length ($\lambda$ is the exam's symbol for $\rho_l$) on its inner surface; "conductors are neutral, so zero" is wrong.

**Where it appears.** [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]], [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]] (layered plates), [[problems/two-layer-coaxial-capacitor]] (FA26 4(e)); Unit 4 (per-unit-length $\mathcal{C}$ and $\mathcal{L}$ of a transmission line).

Related: [[concepts/electrostatic-potential]] · [[concepts/gauss-law]] · [[concepts/conductors]] · [[concepts/permittivity]] · [[concepts/conductance]] · [[concepts/electrostatic-energy]].
