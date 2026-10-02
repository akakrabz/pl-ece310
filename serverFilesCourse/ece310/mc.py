"""Helpers for pl-multiple-choice / pl-checkbox options generated in server.py.

Use in question.html (PrairieLearn renders the template with Mustache first):

    <pl-multiple-choice answers-name="roc">
      {{#params.roc_choices}}<pl-answer correct="{{correct}}">{{text}}</pl-answer>{{/params.roc_choices}}
    </pl-multiple-choice>

``correct`` is a JSON boolean, which Mustache renders as "true"/"false".
"""

from __future__ import annotations

import random
import re
from typing import Iterable


def _norm(s: str) -> str:
    return re.sub(r"\s+", "", s)


def choices(correct: str, distractors: Iterable[str], n: int = 4, shuffle: bool = True,
            feedback: dict | None = None) -> list[dict]:
    """One correct option plus up to n-1 distinct distractors (duplicates of the correct text or of
    each other are dropped). feedback maps option text -> feedback string (optional)."""
    seen = {_norm(correct)}
    out = [{"text": correct, "correct": True}]
    for d in distractors:
        if len(out) >= n:
            break
        k = _norm(d)
        if k in seen:
            continue
        seen.add(k)
        out.append({"text": d, "correct": False})
    if feedback:
        for o in out:
            o["feedback"] = feedback.get(o["text"], "")
    if shuffle:
        random.shuffle(out)
    return out


def yes_no(answer: bool, yes: str = "Yes", no: str = "No") -> list[dict]:
    """Fixed-order two-option list (use order="fixed" on the element)."""
    return [{"text": yes, "correct": bool(answer)}, {"text": no, "correct": not answer}]


def true_false(answer: bool) -> list[dict]:
    return yes_no(answer, "True", "False")
