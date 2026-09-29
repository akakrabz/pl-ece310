---
title: "Conductance and the leaky capacitor"
description: "A slightly conducting filling lets the same field that stores charge drive a current I = GV, with G = (σ/ε)C for any single-medium geometry. The element is C in parallel with R = 1/G; charge self-discharges with τ = RC = ε/σ."
tags: [concept, electrostatics, exam-1]
aliases: ["resistance", "leakage", "RC time constant", "current density", "conduction current"]
---

> [!key] Definition and the key identity
> $$
> I = GV,\qquad G = \frac{I}{V}\ [\text{S}],\qquad R = \frac1G,\qquad\qquad G = \frac{\sigma}{\epsilon}\,C,
> $$
> where $I = \displaystyle\int\mathbf{J}\cdot d\mathbf{S}$ is the current through any surface separating the two conductors and $\mathbf{J} = \sigma\mathbf{E}$. $G$ and $C$ share one geometric factor ($A/d$, $2\pi L/\ln(b/a)$, $4\pi ab/(b-a)$): $C = \epsilon\times$factor, $G = \sigma\times$factor. The identity needs a **single homogeneous** medium (uniform $\epsilon$ and $\sigma$) filling the field region.

**Recipe.** $V\to\mathbf{E} = -\nabla V\to\mathbf{J} = \sigma\mathbf{E}\to I = \int\mathbf{J}\cdot d\mathbf{S}\to G = I/V$; per unit length $\mathcal{G} = G/L$ [S/m]. Plates: $G = \sigma A/d$, $R = d/\sigma A$. Coax: $J_r = \sigma V_0/(r\ln(b/a))$, $I = J_r\cdot2\pi rL$ is independent of $r$ ($\nabla\cdot\mathbf{J} = 0$), $\mathcal{G} = 2\pi\sigma/\ln(b/a)$. Spheres: shells in series, $dR = dr/4\pi\sigma r^2$, $G = 4\pi\sigma ab/(b-a)$.

**The circuit.** Charge balance on the positive conductor, $C\,dV/dt = I - GV$, gives $I = C\,dV/dt + GV$: a capacitor in parallel with a resistor. Left alone, $Q(t) = Q(0)e^{-t/\tau}$ with
$$
\tau = RC = \frac{C}{G} = \frac{\epsilon}{\sigma},
$$
independent of geometry — the relaxation time of [[concepts/conductors]] again. Driven at $\omega$: resistor-like for $\omega\ll\sigma/\epsilon$, capacitor-like for $\omega\gg\sigma/\epsilon$.

**Example.** Concentric spheres, $a=1$ m, $b\to\infty$, $\epsilon = 2\epsilon_0$, $\sigma = 2\times10^{-6}$ S/m: $C = 8\pi\epsilon_0\approx222$ pF, $G = (\sigma/\epsilon)C = 8\pi\ \mu\text{S}\approx25\ \mu$S, $\tau = 8.85\ \mu$s; an initial 1 C decays as $e^{-t/\tau}$ through the medium itself, with no external wire.

> [!trap]
> - $G = (\sigma/\epsilon)C$ fails for layered fillings unless every layer has the same $\sigma/\epsilon$; then compute $G$ by the recipe.
> - $\sigma$ is conductivity here; surface charge is $\rho_s$.
> - $R = 1/G$ is the *leakage* resistance across the dielectric, not the resistance of the conductors.

**Where it appears.** [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]], [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]] (Ohm's law, relaxation), HW4 #4; Unit 4 (the per-unit-length $\mathcal{G}$ of a lossy transmission line).

Related: [[concepts/capacitance]] · [[concepts/conductors]] · [[concepts/permittivity]].
