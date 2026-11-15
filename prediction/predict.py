"""Stage one: a cluster's genes plus a known tissue -> ranked leaf cell types.

    predictor = CellTypePredictor(CellxGeneClient())
    result = predictor.predict(["KRT14", "KRT5", ...], "skin of body")
"""

import json
import logging
from pathlib import Path

import numpy as np
from sklearn.naive_bayes import GaussianNB

from dataIngestion.client import CellxGeneClient
from prediction.features import (
    DEFAULT_MIN_PRESENT_FRAC,
    build_lineage,
    decision_scores,
    minimum_present,
    order_candidates,
    per_gene_max_from_summary,
    score_cells,
    sigmoid,
)
from prediction.schemas import LineagePrediction, Prediction, PredictionResult

log = logging.getLogger(__name__)

BATCH_SIZE = 6
ENSEMBL_PREFIX = "ENSG"


def load_model(path):
    """Rebuild the fitted estimator from the plain-JSON model artifact."""
    payload = json.loads(Path(path).read_text())
    model = GaussianNB()
    model.classes_ = np.asarray(payload["classes"], dtype=int)
    model.theta_ = np.asarray(payload["theta"], dtype=float)
    model.var_ = np.asarray(payload["var"], dtype=float)
    model.class_prior_ = np.asarray(payload["prior"], dtype=float)
    model.n_features_in_ = model.theta_.shape[1]
    return payload, model


class CellTypePredictor:
    """Predict which cell types a cluster's marker genes describe, within one tissue."""

    def __init__(self, client: CellxGeneClient, model_path=None):
        path = Path(model_path) if model_path else Path(__file__).with_name("model.json")
        if not path.exists():
            raise FileNotFoundError(f"model artifact not found: {path}")
        self.client = client
        self.payload, self.model = load_model(path)
        self.min_present_frac = float(
            self.payload.get("min_present_frac", DEFAULT_MIN_PRESENT_FRAC)
        )
        self._positive = list(self.model.classes_).index(1)
        self._known = None

    @property
    def known_cells(self):
        if self._known is None:
            self._known = {row[0] for row in self.client._ids("cell")}
        return self._known

    def predict(self, genes, tissue, top_n=10) -> PredictionResult:
        input_symbols = [g.strip() for g in genes if g and g.strip()]
        if not input_symbols:
            raise ValueError("no genes supplied")

        tissue_id, tissue_name = self._resolve_tissue(tissue)
        resolved, missing = self._resolve_genes(input_symbols)
        if not resolved:
            raise ValueError(f"none of the supplied genes are known: {', '.join(missing)}")

        current = self.client.snapshot_id()
        trained = self.payload.get("snapshot_id")
        if trained and current and trained != current:
            log.warning(
                "model was trained on snapshot %s but the cache holds %s", trained, current
            )

        expression = {}
        for start in range(0, len(resolved), BATCH_SIZE):
            summary = self.client.get_expression(resolved[start : start + BATCH_SIZE])
            for gene, tissues in summary.expression_summary.items():
                expression.setdefault(gene, {}).update(tissues)

        present = [
            gene
            for gene in resolved
            if any(cell != "tissue_stats" for cells in expression.get(gene, {}).values() for cell in cells)
        ]
        if not present:
            raise ValueError("the supplied genes have no expression data in this snapshot")

        min_present = minimum_present(len(present), self.min_present_frac)
        gmax = per_gene_max_from_summary(expression, present)
        candidates = self.client.tissue_leaves(tissue_id)
        features = score_cells(expression, present, tissue_id, gmax, allowed=candidates)
        if not features:
            raise ValueError(
                f"no cell types in {tissue_name!r} express any of the supplied genes"
            )

        matrix = np.array(
            [[f.pc_wilson, f.me_norm, f.vote_frac] for f in features.values()], dtype=float
        )
        scores = sigmoid(decision_scores(self.model, matrix))
        score_map = dict(zip(features, scores))
        ordered = order_candidates(features, score_map, min_present)

        predictions = []
        for rank, (cell_id, eligible) in enumerate(ordered[:top_n], start=1):
            feature = features[cell_id]
            predictions.append(
                Prediction(
                    cell_id=cell_id,
                    cell_name=self.client.name_for("cell", cell_id),
                    score=float(score_map[cell_id]),
                    pc_wilson=feature.pc_wilson,
                    me_norm=feature.me_norm,
                    vote_frac=feature.vote_frac,
                    n_present=feature.n_present,
                    eligible=eligible,
                    rank=rank,
                )
            )

        n_eligible = sum(1 for _cell, eligible in ordered if eligible)
        lineage = []
        if ordered:
            winner = ordered[0][0]
            lineage = [
                LineagePrediction(
                    cell_id=node.cell_id,
                    cell_name=self.client.label(node.cell_id),
                    depth=node.depth,
                    support=node.support,
                    rival=node.rival,
                    margin=node.margin,
                    is_leaf=node.is_leaf,
                )
                for node in build_lineage(
                    score_map,
                    winner,
                    self.client.ancestors,
                    self.client.lineage_path,
                    self.client.children,
                    self.known_cells,
                    self.client.is_leaf,
                )
            ]

        return PredictionResult(
            tissue_id=tissue_id,
            tissue_name=tissue_name,
            input_symbols=input_symbols,
            input_genes=resolved,
            missing_genes=missing,
            snapshot_id=current,
            min_present=min_present,
            n_candidates=len(features),
            n_eligible=n_eligible,
            gate_relaxed=n_eligible == 0,
            lineage=lineage,
            predictions=predictions,
        )

    def _resolve_genes(self, symbols):
        resolved, missing, seen = [], [], set()
        for symbol in symbols:
            if symbol.upper().startswith(ENSEMBL_PREFIX):
                ensembl = symbol.upper() if self.client.name_for("gene", symbol.upper()) else None
            else:
                ensembl = self.client.gene_to_ensembl(symbol)
            if ensembl is None:
                missing.append(symbol)
                continue
            if ensembl not in seen:
                seen.add(ensembl)
                resolved.append(ensembl)
        return resolved, missing

    def _resolve_tissue(self, tissue):
        text = (tissue or "").strip()
        if not text:
            raise ValueError("tissue is required")

        if text.upper().startswith("UBERON:"):
            tissue_id = text.upper()
            name = self.client.name_for("tissue", tissue_id)
            if name is None:
                raise ValueError(f"unknown tissue id: {tissue_id}")
            return tissue_id, name

        tissue_id = self.client.id_for("tissue", text)
        if tissue_id:
            return tissue_id, self.client.name_for("tissue", tissue_id)

        matches = self.client.find("tissue", text)
        if not matches:
            raise ValueError(f"unknown tissue: {text!r}")
        if len(matches) > 1:
            options = ", ".join(f"{name} ({term_id})" for term_id, name in matches)
            raise ValueError(f"ambiguous tissue {text!r}; candidates: {options}")
        tissue_id, name = matches[0]
        return tissue_id, name
