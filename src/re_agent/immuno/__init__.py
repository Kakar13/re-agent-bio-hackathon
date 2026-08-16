"""Semi-supervised MHC class II immunogenicity model for de novo binders.

A small attention-pooling head is trained on labeled natural peptides (IEDB) over
frozen ESM-2 embeddings, then adapted to the de novo distribution with a Mean
Teacher consistency objective over unlabeled designed sequences. No pseudo-labels
and no distillation from an external predictor: the teacher is an EMA of the student.
"""

from re_agent.immuno.config import PATHS, WINDOW, TrainConfig

__all__ = ["PATHS", "WINDOW", "TrainConfig"]
