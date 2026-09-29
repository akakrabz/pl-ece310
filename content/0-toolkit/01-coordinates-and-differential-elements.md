---
title: "Coordinates and differential elements"
description: "Cartesian, cylindrical and spherical coordinates: unit vectors, dl, dS and dV, with the one warning that matters — only Cartesian unit vectors are constant."
tags: [toolkit, reference]
---

*Toolkit · reference page · used by every lecture from [[1-electrostatics/01-fields-forces-and-the-maxwell-roadmap|Lecture 1]] onward*

Every integral in the course is built from a differential element: a step $d\mathbf{l}$ along a curve, a patch $d\mathbf{S}$ on a surface, a cell $dV$ in a volume. Pick the coordinate system that matches the symmetry of the problem, read the elements off this page, and the integral usually becomes one-dimensional.

## Cartesian $(x,y,z)$

Unit vectors $\hat{x},\hat{y},\hat{z}$ are **constant** — the same at every point — which is why they can be pulled out of integrals.

$$
d\mathbf{l} = \hat{x}\,dx+\hat{y}\,dy+\hat{z}\,dz,\qquad
d\mathbf{S} = \pm\hat{x}\,dy\,dz,\ \pm\hat{y}\,dz\,dx,\ \pm\hat{z}\,dx\,dy,\qquad
dV = dx\,dy\,dz .
$$

## Cylindrical $(r,\phi,z)$

$r$ = distance from the $z$-axis, $\phi$ = azimuth measured from $+x$ toward $+y$. Conversions $x = r\cos\phi$, $y = r\sin\phi$.

$$
\hat{r} = \hat{x}\cos\phi+\hat{y}\sin\phi,\qquad \hat{\phi} = -\hat{x}\sin\phi+\hat{y}\cos\phi,\qquad \hat{r}\times\hat{\phi} = \hat{z}.
$$

$$
d\mathbf{l} = \hat{r}\,dr + \hat{\phi}\,r\,d\phi + \hat{z}\,dz,\qquad dV = r\,dr\,d\phi\,dz .
$$

| surface | $d\mathbf{S}$ (outward for the closed cylinder) |
|---|---|
| curved side, radius $r$ | $\hat{r}\,r\,d\phi\,dz$ → total area $2\pi rL$ |
| top / bottom cap | $\pm\hat{z}\,r\,dr\,d\phi$ → area $\pi r^2$ |
| half-plane $\phi=$ const | $\pm\hat{\phi}\,dr\,dz$ |

## Spherical $(r,\theta,\phi)$

$r$ = distance from the origin, $\theta$ = polar angle from $+z$, $\phi$ = azimuth. Conversions $x = r\sin\theta\cos\phi$, $y = r\sin\theta\sin\phi$, $z = r\cos\theta$.

$$
d\mathbf{l} = \hat{r}\,dr + \hat{\theta}\,r\,d\theta + \hat{\phi}\,r\sin\theta\,d\phi,\qquad
dV = r^2\sin\theta\,dr\,d\theta\,d\phi .
$$

| surface | $d\mathbf{S}$ |
|---|---|
| sphere of radius $r$ | $\hat{r}\,r^2\sin\theta\,d\theta\,d\phi$ → total area $4\pi r^2$ |
| cone $\theta=$ const | $\pm\hat{\theta}\,r\sin\theta\,dr\,d\phi$ |
| half-plane $\phi=$ const | $\pm\hat{\phi}\,r\,dr\,d\theta$ |

The solid angle of a patch is $d\Omega = \sin\theta\,d\theta\,d\phi$; a full sphere has $4\pi$ steradians, a hemisphere $2\pi$ — which is why a point charge sends exactly half its flux through any plane that does not contain it ([[1-electrostatics/03-gauss-law-at-work#2-gausss-law-restated-for-use|Lecture 3]]).

> [!trap] $\hat{r}$, $\hat{\phi}$, $\hat{\theta}$ change from point to point
> Only Cartesian unit vectors are constant. In a superposition integral, the direction $\hat{R}$ from each source element to the field point is different for each element, so you may **not** write "$\hat{r}\int(\dots)$" unless you have first shown (by symmetry) that the surviving component is along a direction fixed *at the field point*. When in doubt, resolve into $\hat{x},\hat{y},\hat{z}$ before integrating. See the finite-line-charge example in [[1-electrostatics/02-coulombs-law-superposition-and-gauss#4-continuous-distributions--the-5-step-program|Lecture 2]].

> [!tip] Which system?
> Match the source: point/sphere → spherical; line/cylinder/coax → cylindrical; sheet/slab/parallel plates → Cartesian. Gauss's law only pays off when $|\mathbf{D}|$ is constant on a coordinate surface of the system you chose ([[concepts/gauss-law]]).

## The three products of the course

| object | formula | meaning |
|---|---|---|
| line integral | $\displaystyle\int_C\mathbf{F}\cdot d\mathbf{l}$ | adds the tangential component; work, EMF, voltage |
| surface integral | $\displaystyle\int_S\mathbf{F}\cdot d\mathbf{S}$ | adds the normal component; flux |
| volume integral | $\displaystyle\int_V f\,dV$ | adds a density; total charge |

Related: [[concepts/differential-elements]] · [[0-toolkit/02-vector-calculus-cheatsheet]].
