"""Stage two: check a prediction against /markers evidence.

The marker list is used only as a verifier, never as a scoring feature, and it
never reorders the predictions. A cell type with no marker list gets
``coverage=None``: absent evidence is not negative evidence.
"""

import logging

from dataIngestion.client import FULL_MARKERS, CellxGeneClient
from prediction.schemas import GeneEvidence, PredictionResult, Verification

log = logging.getLogger(__name__)


class MarkerVerifier:
    """Annotate stage-one predictions with how many input genes are markers of the cell type."""

    def __init__(self, client: CellxGeneClient, n_markers=FULL_MARKERS):
        self.client = client
        self.n_markers = n_markers

    def verify(self, result: PredictionResult, top_n=None) -> PredictionResult:
        predictions = result.predictions if top_n is None else result.predictions[:top_n]
        for prediction in predictions:
            prediction.verification = self._verify_one(result, prediction)
        return result

    def _verify_one(self, result, prediction) -> Verification:
        markers = self.client.get_markers(
            result.tissue_id, prediction.cell_id, self.n_markers
        )
        if not markers:
            log.info(
                "no marker list for %s in %s", prediction.cell_id, result.tissue_id
            )
            return Verification(verified=False)

        by_gene = {
            marker.gene_ontology_term_id: (position + 1, marker)
            for position, marker in enumerate(markers)
        }

        evidence = []
        matched = 0
        for gene in result.input_genes:
            symbol = self.client.ensembl_to_gene(gene)
            hit = by_gene.get(gene)
            if hit is None:
                evidence.append(
                    GeneEvidence(gene_id=gene, gene_symbol=symbol, in_markers=False)
                )
                continue
            rank, marker = hit
            matched += 1
            evidence.append(
                GeneEvidence(
                    gene_id=gene,
                    gene_symbol=symbol,
                    in_markers=True,
                    marker_score=marker.marker_score,
                    specificity=marker.specificity,
                    rank=rank,
                )
            )

        total = len(result.input_genes)
        return Verification(
            coverage=matched / total if total else None,
            matched=matched,
            total=total,
            verified=matched > 0,
            evidence=evidence,
        )
