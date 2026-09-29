---
title: "Worked problem — Curl, potential, and charge from a given field"
description: "Given a polynomial electric field: check it is conservative, find its potential, find the charge density that produces it, and evaluate a line integral by path independence. Modelled on FA26 Exam 1, problem 1, with different numbers."
tags: [problem, electrostatics, exam-1]
---

*Problem · style of FA26 Exam 1 #1 (25 pts) · uses [[1-electrostatics/04-divergence-and-curl|Lecture 4]] and Lecture 5 · concepts: [[concepts/curl]], [[concepts/conservative-field]], [[concepts/divergence]], [[concepts/gauss-law]]*

> [!question] Problem
> An electric field in free space is
> $$
> \mathbf{E} = E_0\left(y^2\,\hat{x} + 2xy\,\hat{y} + 3z^2\,\hat{z}\right)\quad[\text{V/m}].
> $$
> (a) Find $\nabla\times\mathbf{E}$. (b) Find the electrostatic potential $V(x,y,z)$, taking $V(0,0,0)=0$. (c) Find the volume charge density $\rho(x,y,z)$ that produces this field. (d) Evaluate $\displaystyle\int_A^B\mathbf{E}\cdot d\mathbf{l}$ from $A=(0,0,0)$ to $B=(1,1,0)$.

## Before computing anything

Look at the structure. $E_x$ depends only on $y$, $E_y$ on $x$ and $y$, $E_z$ only on $z$. So in the curl, every derivative involving $z$ vanishes and only the $\hat{z}$-component can survive; in the divergence, only $\partial_yE_y$ and $\partial_zE_z$ can be nonzero. That is half the work done by inspection.

## (a) Curl — is the field electrostatic at all?

Expand the determinant along the top row (middle term with a minus sign):

$$
\nabla\times\mathbf{E} = \hat{x}\Big(\frac{\partial E_z}{\partial y}-\frac{\partial E_y}{\partial z}\Big)
+\hat{y}\Big(\frac{\partial E_x}{\partial z}-\frac{\partial E_z}{\partial x}\Big)
+\hat{z}\Big(\frac{\partial E_y}{\partial x}-\frac{\partial E_x}{\partial y}\Big)
= E_0\big[\hat{x}(0-0)+\hat{y}(0-0)+\hat{z}(2y-2y)\big] = \mathbf{0}.
$$

> [!key] Answer (a)
> $\nabla\times\mathbf{E} = 0$. The field is curl-free, hence [[concepts/conservative-field|conservative]]: a potential exists and line integrals are path-independent. Parts (b) and (d) are now legitimate.

> [!trap] The order in the $\hat{z}$ bracket
> It is $\partial_xE_y - \partial_yE_x$. Here both are $2yE_0$, so a reversed bracket would still give zero — which is exactly how a sign error hides on an exam and then bites on the next problem. Expand mechanically every time.

## (b) The potential

Two ways.

**Exact differential** (the quickest route). $dV = -\mathbf{E}\cdot d\mathbf{l} = -E_0\,(y^2\,dx + 2xy\,dy + 3z^2\,dz)$. Notice $y^2\,dx + 2xy\,dy = d(xy^2)$ and $3z^2dz = d(z^3)$, so

$$
V = -E_0\left(xy^2 + z^3\right) + C,\qquad V(0,0,0) = 0 \Rightarrow C = 0 .
$$

**Explicit path** (if you don't spot the differential): go $(0,0,0)\to(x,0,0)\to(x,y,0)\to(x,y,z)$ along the axes.

$$
V(x,y,z) = -\int_{\text{path}}\mathbf{E}\cdot d\mathbf{l}
= -E_0\Big[\underbrace{\int_0^x (0)^2\,dx'}_{y=0\ \text{on leg 1}} + \int_0^y 2x\,y'\,dy' + \int_0^z 3z'^2\,dz'\Big]
= -E_0\left(xy^2+z^3\right).
$$

On leg 2, $x$ is held at its final value and $y'$ runs from 0 to $y$; on leg 1 the integrand is zero because $y=0$ there.

**Check:** $-\nabla V = E_0\left(y^2\hat{x} + 2xy\,\hat{y} + 3z^2\hat{z}\right) = \mathbf{E}$ ✓.

> [!key] Answer (b)
> $V(x,y,z) = -E_0\left(xy^2 + z^3\right)$ [V].

> [!trap] The classic wrong answer
> Antidifferentiating each component separately and adding — $\int y^2dx + \int 2xy\,dy + \int 3z^2dz = xy^2 + xy^2 + z^3$ — counts the $xy^2$ term twice. There is one line integral along one path, not three independent integrals.

## (c) The charge density

Gauss's law in differential form, $\rho = \nabla\cdot\mathbf{D} = \epsilon_0\nabla\cdot\mathbf{E}$:

$$
\nabla\cdot\mathbf{E} = \frac{\partial}{\partial x}(E_0y^2) + \frac{\partial}{\partial y}(2E_0xy) + \frac{\partial}{\partial z}(3E_0z^2) = E_0\,(0 + 2x + 6z).
$$

> [!key] Answer (c)
> $\rho(x,y,z) = \epsilon_0E_0\,(2x + 6z)$ [C/m³].

Units: $\epsilon_0$ [F/m] × $E_0$ [V/m³] × length [m] = C/m³ ✓ (for $\mathbf{E}$ to be in V/m, $E_0$ must be V/m³).

> [!trap]
> - Forgetting $\epsilon_0$ leaves you with V/m², not a charge density.
> - "The curl is zero, so there's no charge" — no. Curl-free tells you about *circulation*; charge is *divergence*. This field is conservative *and* has charge everywhere $2x+6z\neq0$.

## (d) The line integral, without integrating

The field is conservative, so $\int_A^B\mathbf{E}\cdot d\mathbf{l} = V(A) - V(B)$ — note the order: $V_A - V_B = \int_A^B\mathbf{E}\cdot d\mathbf{l}$ (potential *drops* along $\mathbf{E}$).

$$
\int_{(0,0,0)}^{(1,1,0)}\mathbf{E}\cdot d\mathbf{l} = V(0,0,0) - V(1,1,0) = 0 - \big[-E_0(1\cdot1^2 + 0)\big] = +E_0\ [\text{V}].
$$

Cross-check by two explicit paths (both must agree, since the field is conservative):

- *Via $(1,0,0)$*: leg 1 runs along $x$ at $y=0$, where $E_x = E_0y^2 = 0$ → contributes 0. Leg 2 runs along $y$ at $x=1$: $\int_0^1 E_y\,dy' = E_0\int_0^1 2y'\,dy' = E_0$. Total $E_0$.
- *Via $(0,1,0)$*: leg 1 runs along $y$ at $x=0$, where $E_y = 2E_0xy = 0$ → 0. Leg 2 runs along $x$ at $y=1$: $\int_0^1 E_x\,dx' = E_0\int_0^1 1^2\,dx' = E_0$. Total $E_0$. ✓

On each leg, write down *which* coordinate varies and hold the others at their current values — mislabelling a leg is the most common way this cross-check "fails" on an exam.

> [!key] Answer (d)
> $\displaystyle\int_A^B\mathbf{E}\cdot d\mathbf{l} = +E_0$ [V] — numerically $E_0\times(1\ \text{m})^3$.

> [!trap] Sign
> $V(B) - V(A) = -E_0$ is the negative of the line integral. Remember: moving *along* $\mathbf{E}$, the potential *falls*, and $\int\mathbf{E}\cdot d\mathbf{l}$ is positive.

## What this problem is really testing

| part | idea | lecture |
|---|---|---|
| (a) | curl in Cartesian coordinates; conservative ⟺ curl-free | [[1-electrostatics/04-divergence-and-curl#3-curl-circulation-per-unit-area|L4 §3]], [[1-electrostatics/04-divergence-and-curl#5-conservative-fields|L4 §5]] |
| (b) | $V$ from $\mathbf{E}$: exact differential or explicit path; reference point | Lecture 5 |
| (c) | $\rho = \nabla\cdot\mathbf{D}$; only diagonal derivatives | [[1-electrostatics/04-divergence-and-curl#2-divergence-flux-per-unit-volume|L4 §2]] |
| (d) | path independence, $V_A - V_B = \int_A^B\mathbf{E}\cdot d\mathbf{l}$ | Lecture 5 |

**Make your own variants.** Start from any potential $V = -E_0(\alpha x^my^n + \beta z^p)$ and differentiate: the field is conservative *by construction*, $\rho = \epsilon_0E_0[\alpha m(m-1)x^{m-2}y^n + \alpha n(n-1)x^my^{n-2} + \beta p(p-1)z^{p-2}]$, and $\int_A^B\mathbf{E}\cdot d\mathbf{l} = V(A)-V(B)$. A field written down at random, such as $E_0(a\,y\,\hat{x} + b\,x^2\hat{y})$, is almost never conservative — its curl is $E_0(2bx - a)\hat{z}$ — and then part (b) has no answer and part (d) depends on the path. The PrairieLearn question set in this repository generates fresh variants of this problem from the potential.
