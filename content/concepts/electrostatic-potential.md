---
title: "Electrostatic potential V"
description: "Work per unit charge: V(b) − V(a) = −∫E·dl, path-independent because static E is curl-free; E = −∇V; the point-charge potential and scalar superposition; equipotentials."
tags: [concept, electrostatics, exam-1]
aliases: ["potential", "voltage", "electric potential", "equipotential", "gradient"]
---

> [!key] Definition and the two-way map
> $$
> V(b) - V(a) = -\int_a^b\mathbf{E}\cdot d\mathbf{l}\quad[\text{V} = \text{J/C}],
> \qquad\qquad
> \mathbf{E} = -\nabla V\quad[\text{V/m}].
> $$
> $V(b) - V(a)$ is the work per unit charge needed to carry a test charge from $a$ to $b$ against the field; $\int_a^b\mathbf{E}\cdot d\mathbf{l} = V(a) - V(b)$ is the **voltage drop** from $a$ to $b$. The integral is independent of path because $\nabla\times\mathbf{E} = 0$ ([[concepts/conservative-field]]), which is also why a single-valued $V$ exists. $V$ is defined up to a constant, fixed by choosing a **reference** ($V=0$ at infinity for finite charge distributions; at a stated point otherwise).

**Picture.** $V$ is altitude; $\nabla V$ points straight uphill, so $\mathbf{E} = -\nabla V$ points straight *downhill*, from high $V$ to low $V$, always perpendicular to the **equipotential** surfaces. Crowded equipotentials mean a strong field: $\lvert\mathbf{E}\rvert$ is the potential drop per unit distance. A positive charge released from rest accelerates toward lower $V$.

**From charge to V.** Point charge: $V = \dfrac{Q}{4\pi\epsilon_0r}$ (reference at infinity), falling as $1/r$ while $E\propto1/r^2$. Superposition is a **scalar** sum, $V(\mathbf{r}) = \dfrac{1}{4\pi\epsilon_0}\sum_i\dfrac{Q_i}{\lvert\mathbf{r}-\mathbf{r}_i\rvert}$ or $\dfrac{1}{4\pi\epsilon_0}\displaystyle\int\frac{\rho(\mathbf{r}')\,d^3\mathbf{r}'}{\lvert\mathbf{r}-\mathbf{r}'\rvert}$ — no components to track. Potential energy of a charge $q$ at $\mathbf{r}$: $U = qV(\mathbf{r})$; one electron-volt is $e\times1$ V $= 1.602\times10^{-19}$ J.

**From E to V.** Check the curl; then either recognize $\mathbf{E}\cdot d\mathbf{l}$ as an exact differential $-dV$ (e.g. $y\,dx + x\,dy = d(xy)$), or integrate along an axis-parallel staircase from the reference point, using the current values of the coordinates already traversed; fix the constant; verify $-\nabla V = \mathbf{E}$. In curvilinear coordinates, $\nabla V = \hat{r}\partial_rV + \hat{\phi}\frac1r\partial_\phi V + \hat{z}\partial_zV$ (cylindrical) and $\hat{r}\partial_rV + \hat{\theta}\frac1r\partial_\theta V + \hat{\phi}\frac{1}{r\sin\theta}\partial_\phi V$ (spherical).

**Examples.** $V = x^2-6y$ ⇒ $\mathbf{E} = -2x\hat{x}+6\hat{y}$. $\mathbf{E} = 2x\hat{x}+3z\hat{y}+3(y+1)\hat{z}$ with $V(0)=0$ ⇒ $V = -x^2-3(y+1)z$. Parallel plates $d$ apart with the top at $V_0$: $V = V_0z/d$, $\mathbf{E} = -(V_0/d)\hat{z}$. Coax with $V(a) = V_0$, $V(b) = 0$: $V = V_0\ln(b/r)/\ln(b/a)$.

> [!trap]
> - The minus sign. Away from a positive charge $V$ must *fall*; if yours rises, a sign is missing.
> - Do not reverse the limits *and* flip $d\mathbf{l}$; write $d\mathbf{l} = \hat{r}\,dr$ and let $\int_\infty^r$ carry the direction.
> - Do not antidifferentiate the three components separately and add: cross terms get counted twice. One integral along one path.
> - Reference at infinity is unavailable for infinite lines, sheets and slabs (the integral diverges); the problem must supply "$V=0$ at …".
> - $\nabla\cdot\nabla V = \nabla^2V$ is the Laplacian, not zero; the identity is $\nabla\times\nabla V = 0$.

**Where it appears.** [[1-electrostatics/05-curl-free-fields-and-the-electrostatic-potential|Lecture 5]] (everything above), [[1-electrostatics/07-poisson-and-laplace|Lecture 7]] (Poisson's equation for $V$), [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]] (conductors are equipotentials), [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]] ($C = Q/V$), [[problems/curl-potential-and-charge-from-a-field]], [[problems/charged-slab-with-a-power-law-profile]], [[problems/two-layer-coaxial-capacitor]].

Related: [[concepts/conservative-field]] · [[concepts/electric-field]] · [[concepts/poissons-equation]].
