---
title: "Polarization P and bound charge"
description: "P is the dipole moment per unit volume of a dielectric. Its ends leave bound surface charge ρsb = P·n̂ and its divergence leaves bound volume charge ρb = −∇·P. D = ε₀E + P absorbs the bound charge so that Gauss's law counts free charge only."
tags: [concept, electrostatics, exam-1]
aliases: ["bound charge", "polarization charge", "dielectric", "dipole moment per unit volume", "D = ε₀E + P"]
---

> [!key] Definition
> $$
> \mathbf{P} = N\mathbf{p}\ \ (\text{dipole moment per unit volume})\quad[\text{C/m}^2],\qquad
> \rho_{sb} = \mathbf{P}\cdot\hat{n},\qquad \rho_b = -\nabla\cdot\mathbf{P},
> \qquad
> \mathbf{D} \equiv \epsilon_0\mathbf{E}+\mathbf{P}.
> $$
> $\hat{n}$ is the outward normal of the dielectric's surface. Gauss's law counts all charge, $\nabla\cdot(\epsilon_0\mathbf{E}) = \rho_f+\rho_b$; moving $\rho_b = -\nabla\cdot\mathbf{P}$ to the left gives $\nabla\cdot\mathbf{D} = \rho_f$. The total bound charge of a neutral dielectric is always zero.

**Physics.** An applied field stretches each neutral atom into a small dipole $\mathbf{p} = q\mathbf{d}$ aligned with $\mathbf{E}$. Inside the material the head of one dipole cancels the tail of the next; only where $\mathbf{P}$ *ends* (a surface) or *changes* (a gradient) does uncompensated charge remain. The bound charge's own field opposes the applied field, so the field inside a dielectric is **reduced but not zero** — unlike a conductor, whose free charge keeps moving until the field vanishes. For an infinite slab polarized normal to its faces the bound field is $\mathbf{E}_p = -\mathbf{P}/\epsilon_0$ and $\mathbf{D}$ is the same inside and out; for a sphere it is $-\mathbf{P}/3\epsilon_0$, and for a slab polarized along its faces it is zero (then $\mathbf{E}$ is unchanged and $\mathbf{D}$ jumps), so "$\mathbf{D}$ unchanged" is a fact about one geometry, not a law. The general facts are $\nabla\cdot\mathbf{D} = \rho_f$ and the continuity of $D_n$ across free-charge-free interfaces.

**Linear media.** $\mathbf{P} = \epsilon_0\chi_e\mathbf{E}$, with $\mathbf{E}$ the *total* field, gives $\mathbf{D} = \epsilon\mathbf{E}$ and $\mathbf{P} = (\epsilon_r-1)\epsilon_0\mathbf{E} = \mathbf{D}-\epsilon_0\mathbf{E}$ ([[concepts/permittivity]]). In practice: get $\mathbf{D}$ from the free charge, $\mathbf{E} = \mathbf{D}/\epsilon$, and $\mathbf{P}$ last, if asked.

**Examples.** Slab in $18\hat{x}$ V/m with $3\hat{x}$ inside: $\mathbf{D} = 18\epsilon_0\hat{x}$ everywhere, $\mathbf{P} = 15\epsilon_0\hat{x}$ inside, faces $\pm15\epsilon_0$ C/m². Dielectric interface with $P_{1z} = -10\epsilon_0$ above and $P_{2z} = -16\epsilon_0$ below: $\rho_{sb} = P_{2n}-P_{1n} = -6\epsilon_0$ C/m². Graded $\epsilon(z)$ between plates: $P_z = \epsilon_0z/2$ ⇒ uniform $\rho_b = -\epsilon_0/2$ C/m³. Metal sphere $Q$ in a dielectric shell: $\rho_b = 0$, bound totals $\mp(1-1/\epsilon_r)Q$ on the inner and outer surfaces.

> [!trap]
> - $\mathbf{P} = \mathbf{D} - \epsilon\mathbf{E}$ is identically zero; it is $\mathbf{D} - \epsilon_0\mathbf{E}$.
> - $\mathbf{P} = \chi_e\mathbf{E}$ is missing an $\epsilon_0$ (units: C/m² ≠ V/m).
> - $\mathbf{P}$ responds to the total field, not the applied field alone.
> - Interface between two dielectrics: $\hat{n}\cdot(\mathbf{P}_1-\mathbf{P}_2) = -\rho_{sb}$ (minus sign; compare $\hat{n}\cdot(\mathbf{D}_1-\mathbf{D}_2) = +\rho_s$).

**Where it appears.** [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]] (derived), [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]] (used), [[problems/fields-across-a-dielectric-interface]] (FA26 3(c)), [[problems/two-layer-coaxial-capacitor]] (bound line charges as a bonus).

Related: [[concepts/electric-flux-density]] · [[concepts/permittivity]] · [[concepts/conductors]] · [[concepts/charge-density]].
