"""Result types for the prediction pipeline."""

from pydantic import BaseModel, ConfigDict, Field


class GeneEvidence(BaseModel):
    model_config = ConfigDict(extra="ignore")

    gene_id: str
    gene_symbol: str | None = None
    in_markers: bool | None = None
    marker_score: float | None = None
    specificity: float | None = None
    rank: int | None = None


class Verification(BaseModel):
    """Stage-two evidence for one predicted cell type.

    ``coverage`` is None when the API returned no marker list at all: absence of
    a marker list is not evidence against the prediction.
    """

    model_config = ConfigDict(extra="ignore")

    coverage: float | None = None
    matched: int = 0
    total: int = 0
    verified: bool = False
    evidence: list[GeneEvidence] = Field(default_factory=list)


class Prediction(BaseModel):
    model_config = ConfigDict(extra="ignore")

    cell_id: str
    cell_name: str | None = None
    score: float
    pc_wilson: float
    me_norm: float
    vote_frac: float
    n_present: int
    eligible: bool
    rank: int
    verification: Verification | None = None


class Confidence(BaseModel):
    """Label-free signals for how trustworthy the top answer is.

    ``margin`` is the top1-top2 log-odds gap; a small margin means the runner-up
    (a different lineage, since the winner is a leaf) scores almost as well,
    which is the signature of a mixed cluster. ``n_present_ratio`` is how much
    of the input the winner actually expressed; ``vote_frac`` the winner's share
    of confidently-expressed genes.
    """

    model_config = ConfigDict(extra="ignore")

    margin: float
    n_present_ratio: float
    vote_frac: float


class AlternativeLineage(BaseModel):
    """A top-ranked prediction outside the winner's subtree.

    Populated in rank order; the presence of a high-scoring alternative is the
    signature of a mixed cluster, and ``lca_depth`` says how deep the two
    lineages diverge (0 = they only share the root).
    """

    model_config = ConfigDict(extra="ignore")

    cell_id: str
    cell_name: str | None = None
    score: float
    rank: int
    lca_id: str | None = None
    lca_name: str | None = None
    lca_depth: int | None = None


class PredictionResult(BaseModel):
    model_config = ConfigDict(extra="ignore")

    tissue_id: str
    tissue_name: str | None = None
    input_symbols: list[str] = Field(default_factory=list)
    input_genes: list[str] = Field(default_factory=list)
    missing_genes: list[str] = Field(default_factory=list)
    snapshot_id: str | None = None
    min_present: int = 0
    n_candidates: int = 0
    n_eligible: int = 0
    gate_relaxed: bool = False
    lineage: list[LineagePrediction] = Field(default_factory=list)
    alternatives: list[AlternativeLineage] = Field(default_factory=list)
    confidence: Confidence | None = None
    predictions: list[Prediction] = Field(default_factory=list)


class LineagePrediction(BaseModel):
    """One rung of the predicted ontology chain, root first.

    ``support`` is the summed leaf score beneath the term; ``rival`` is the
    strongest sibling at the same level; ``margin`` is the difference.
    """

    model_config = ConfigDict(extra="ignore")

    cell_id: str
    cell_name: str | None = None
    depth: int
    support: float
    rival: float
    margin: float
    is_leaf: bool = False
