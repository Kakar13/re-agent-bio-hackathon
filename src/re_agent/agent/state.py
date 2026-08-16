"""Typed LangGraph state shared by runtime and UI."""

from __future__ import annotations

import operator
from typing import Annotated, Any, Literal, TypedDict

from langgraph.graph import MessagesState


class DirectScreenRequest(TypedDict):
    """A UI-triggered screen that must bypass model tool selection."""

    sequence: str
    candidate_id: str


class ScientificAgentState(MessagesState):
    """Append-only scientific artifacts and reviewer findings."""

    artifacts: Annotated[list[dict[str, Any]], operator.add]
    reviews: Annotated[list[dict[str, Any]], operator.add]
    screening_profile: Literal[
        "mhc_ii_standard",
        "mhc_ii_plus_chao1",
    ]
    direct_screen_request: DirectScreenRequest | None
