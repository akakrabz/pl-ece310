---
title: "Vector calculus cheat sheet"
description: "Gradient, divergence, curl and Laplacian in the three coordinate systems; the identities and the integral theorems the course uses; and the vector-algebra facts assumed from MATH 241."
tags: [toolkit, reference]
---

*Toolkit · reference page · mirrors the "Vector Calculus" block of the Midterm 1 formula sheet, plus the cylindrical and spherical forms it omits*

## Vector algebra (assumed)

- $\mathbf{A}\cdot\mathbf{B} = A_xB_x+A_yB_y+A_zB_z = |\mathbf{A}||\mathbf{B}|\cos\theta$ — projection; zero for perpendicular vectors.
- $\mathbf{A}\times\mathbf{B} = \begin{vmatrix}\hat{x}&\hat{y}&\hat{z}\\A_x&A_y&A_z\\B_x&B_y&B_z\end{vmatrix}$, $|\mathbf{A}\times\mathbf{B}| = |\mathbf{A}||\mathbf{B}|\sin\theta$, direction by the right-hand rule; $\mathbf{A}\times\mathbf{B} = -\mathbf{B}\times\mathbf{A}$; zero for parallel vectors.
- Right-handed cycle: $\hat{x}\times\hat{y}=\hat{z}$, $\hat{y}\times\hat{z}=\hat{x}$, $\hat{z}\times\hat{x}=\hat{y}$ (and $\hat{r}\times\hat{\phi}=\hat{z}$ cylindrical; $\hat{r}\times\hat{\theta}=\hat{\phi}$ spherical).
- Triple products: $\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C})$ = signed volume of the parallelepiped (cyclic; its absolute value is the volume); $\mathbf{A}\times(\mathbf{B}\times\mathbf{C}) = \mathbf{B}(\mathbf{A}\cdot\mathbf{C})-\mathbf{C}(\mathbf{A}\cdot\mathbf{B})$ ("BAC–CAB").
- Separation vector from point 1 to point 2: $\mathbf{R}_{12} = \mathbf{r}_2-\mathbf{r}_1$; unit vector $\hat{R}_{12} = \mathbf{R}_{12}/|\mathbf{R}_{12}|$.

## The del operator in Cartesian coordinates

$$
\nabla = \hat{x}\frac{\partial}{\partial x}+\hat{y}\frac{\partial}{\partial y}+\hat{z}\frac{\partial}{\partial z}
$$

| operation | formula | type | meaning |
|---|---|---|---|
| gradient | $\nabla f = \hat{x}\dfrac{\partial f}{\partial x}+\hat{y}\dfrac{\partial f}{\partial y}+\hat{z}\dfrac{\partial f}{\partial z}$ | scalar → vector | direction and rate of steepest increase; $df = \nabla f\cdot d\mathbf{l}$ |
| divergence | $\nabla\cdot\mathbf{A} = \dfrac{\partial A_x}{\partial x}+\dfrac{\partial A_y}{\partial y}+\dfrac{\partial A_z}{\partial z}$ | vector → scalar | flux per unit volume; sources and sinks |
| curl | $\nabla\times\mathbf{A} = \hat{x}\Big(\dfrac{\partial A_z}{\partial y}-\dfrac{\partial A_y}{\partial z}\Big)+\hat{y}\Big(\dfrac{\partial A_x}{\partial z}-\dfrac{\partial A_z}{\partial x}\Big)+\hat{z}\Big(\dfrac{\partial A_y}{\partial x}-\dfrac{\partial A_x}{\partial y}\Big)$ | vector → vector | circulation per unit area; whirlpools |
| Laplacian | $\nabla^2 f = \nabla\cdot\nabla f = \dfrac{\partial^2 f}{\partial x^2}+\dfrac{\partial^2 f}{\partial y^2}+\dfrac{\partial^2 f}{\partial z^2}$ | scalar → scalar | how $f$ differs from its local average |

The curl as a determinant (expand along the top row, middle term with a minus sign):

$$
\nabla\times\mathbf{A} = \begin{vmatrix}\hat{x}&\hat{y}&\hat{z}\\[2pt]\dfrac{\partial}{\partial x}&\dfrac{\partial}{\partial y}&\dfrac{\partial}{\partial z}\\[6pt]A_x&A_y&A_z\end{vmatrix}
$$

## Cylindrical $(r,\phi,z)$

$$
\nabla f = \hat{r}\frac{\partial f}{\partial r}+\hat{\phi}\frac{1}{r}\frac{\partial f}{\partial\phi}+\hat{z}\frac{\partial f}{\partial z},\qquad
\nabla\cdot\mathbf{A} = \frac{1}{r}\frac{\partial (rA_r)}{\partial r}+\frac{1}{r}\frac{\partial A_\phi}{\partial\phi}+\frac{\partial A_z}{\partial z}
$$

$$
\nabla\times\mathbf{A} = \hat{r}\Big(\frac{1}{r}\frac{\partial A_z}{\partial\phi}-\frac{\partial A_\phi}{\partial z}\Big)
+\hat{\phi}\Big(\frac{\partial A_r}{\partial z}-\frac{\partial A_z}{\partial r}\Big)
+\hat{z}\,\frac{1}{r}\Big(\frac{\partial (rA_\phi)}{\partial r}-\frac{\partial A_r}{\partial\phi}\Big),\qquad
\nabla^2 f = \frac{1}{r}\frac{\partial}{\partial r}\Big(r\frac{\partial f}{\partial r}\Big)+\frac{1}{r^2}\frac{\partial^2 f}{\partial\phi^2}+\frac{\partial^2 f}{\partial z^2}
$$

Sanity check on the line-charge field $\mathbf{E} = \hat{r}\,\rho_l/(2\pi\epsilon_0 r)$: $\nabla\cdot\mathbf{E} = \frac{1}{r}\frac{\partial}{\partial r}\big(r\cdot\frac{\rho_l}{2\pi\epsilon_0 r}\big) = 0$ for $r>0$ (no charge off the axis), and every curl component vanishes. ✓

## Spherical $(r,\theta,\phi)$

$$
\nabla f = \hat{r}\frac{\partial f}{\partial r}+\hat{\theta}\frac{1}{r}\frac{\partial f}{\partial\theta}+\hat{\phi}\frac{1}{r\sin\theta}\frac{\partial f}{\partial\phi},\qquad
\nabla\cdot\mathbf{A} = \frac{1}{r^2}\frac{\partial (r^2A_r)}{\partial r}+\frac{1}{r\sin\theta}\frac{\partial(\sin\theta\,A_\theta)}{\partial\theta}+\frac{1}{r\sin\theta}\frac{\partial A_\phi}{\partial\phi}
$$

$$
\nabla\times\mathbf{A} = \hat{r}\,\frac{1}{r\sin\theta}\Big[\frac{\partial(\sin\theta\,A_\phi)}{\partial\theta}-\frac{\partial A_\theta}{\partial\phi}\Big]
+\hat{\theta}\,\frac{1}{r}\Big[\frac{1}{\sin\theta}\frac{\partial A_r}{\partial\phi}-\frac{\partial(rA_\phi)}{\partial r}\Big]
+\hat{\phi}\,\frac{1}{r}\Big[\frac{\partial(rA_\theta)}{\partial r}-\frac{\partial A_r}{\partial\theta}\Big]
$$

$$
\nabla^2 f = \frac{1}{r^2}\frac{\partial}{\partial r}\Big(r^2\frac{\partial f}{\partial r}\Big)+\frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\Big(\sin\theta\frac{\partial f}{\partial\theta}\Big)+\frac{1}{r^2\sin^2\theta}\frac{\partial^2 f}{\partial\phi^2}
$$

Sanity check on the point-charge field $\mathbf{E} = \hat{r}\,Q/(4\pi\epsilon_0 r^2)$: $\nabla\cdot\mathbf{E} = \frac{1}{r^2}\frac{\partial}{\partial r}\big(r^2\cdot\frac{Q}{4\pi\epsilon_0 r^2}\big) = 0$ for $r>0$. All the divergence is concentrated at the origin — that is the δ-function of [[concepts/charge-density]].

## Identities (valid for any smooth fields)

$$
\nabla\times(\nabla f) = 0,\qquad
\nabla\cdot(\nabla\times\mathbf{A}) = 0,\qquad
\nabla\times(\nabla\times\mathbf{A}) = \nabla(\nabla\cdot\mathbf{A})-\nabla^2\mathbf{A}
$$

The last one defines the vector Laplacian, $\nabla^2\mathbf{A} \equiv \hat{x}\nabla^2A_x+\hat{y}\nabla^2A_y+\hat{z}\nabla^2A_z$ — **Cartesian components only**. Also useful:

$$
\nabla(fg) = f\nabla g+g\nabla f,\qquad
\nabla\cdot(f\mathbf{A}) = f\,\nabla\cdot\mathbf{A}+\mathbf{A}\cdot\nabla f,\qquad
\nabla\cdot(\mathbf{A}\times\mathbf{B}) = \mathbf{B}\cdot(\nabla\times\mathbf{A})-\mathbf{A}\cdot(\nabla\times\mathbf{B}).
$$

## Integral theorems

| theorem | statement | used for |
|---|---|---|
| gradient theorem | $\displaystyle\int_a^b\nabla f\cdot d\mathbf{l} = f(b)-f(a)$ | $V_A - V_B = \int_A^B\mathbf{E}\cdot d\mathbf{l}$ |
| Stokes | $\displaystyle\oint_C\mathbf{A}\cdot d\mathbf{l} = \int_S(\nabla\times\mathbf{A})\cdot d\mathbf{S}$ | Faraday, Ampère ↔ differential form |
| divergence (Gauss) | $\displaystyle\oint_S\mathbf{A}\cdot d\mathbf{S} = \int_V\nabla\cdot\mathbf{A}\,dV$ | Gauss's law ↔ differential form |

Orientation conventions: $d\mathbf{S}$ outward on a closed surface; on an open surface bounded by $C$, $d\mathbf{l}$ and $d\mathbf{S}$ are related by the right-hand rule (curl the fingers along $C$, the thumb gives $d\mathbf{S}$).

Related: [[1-electrostatics/04-divergence-and-curl]] · [[concepts/divergence]] · [[concepts/curl]] · [[concepts/divergence-theorem]] · [[concepts/stokes-theorem]].
