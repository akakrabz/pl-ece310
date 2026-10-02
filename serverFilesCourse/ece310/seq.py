"""Finite sequences with an explicit start index (values[i] is x[start + i])."""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from .fmt import Q


def trim_seq(values: Sequence, start: int) -> tuple[list, int]:
    """Drop zeros at both ends, keeping track of the start index. [] -> ([0], start)."""
    vals = [Q(v) for v in values]
    while len(vals) > 1 and vals[0] == 0:
        vals.pop(0)
        start += 1
    while len(vals) > 1 and vals[-1] == 0:
        vals.pop()
    return (vals if vals else [Fraction(0)]), start


def conv(x: Sequence, nx: int, h: Sequence, nh: int) -> tuple[list, int]:
    """Exact convolution: (y, ny) with ny = nx + nh and len(y) = len(x) + len(h) - 1."""
    x = [Q(v) for v in x]
    h = [Q(v) for v in h]
    y = [Fraction(0)] * (len(x) + len(h) - 1)
    for i, a in enumerate(x):
        for j, b in enumerate(h):
            y[i + j] += a * b
    return y, nx + nh


def value(values: Sequence, start: int, n: int):
    i = n - start
    return Q(values[i]) if 0 <= i < len(values) else Fraction(0)
