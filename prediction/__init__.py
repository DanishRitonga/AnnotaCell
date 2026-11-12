"""Cell type prediction for AnnotaCell.

Stage one scores every leaf cell type in a tissue from /query expression data;
stage two annotates the shortlist with /markers evidence.
"""

from prediction.predict import CellTypePredictor
from prediction.schemas import (
    GeneEvidence,
    LineagePrediction,
    Prediction,
    PredictionResult,
    Verification,
)
from prediction.verify import MarkerVerifier

__all__ = [
    "CellTypePredictor",
    "GeneEvidence",
    "LineagePrediction",
    "MarkerVerifier",
    "Prediction",
    "PredictionResult",
    "Verification",
]
