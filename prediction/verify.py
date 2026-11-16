"""Stage two: check a prediction against /markers evidence.

The marker list is used only as a verifier, never as a scoring feature, and it
never reorders the predictions. A cell type with no marker list gets
``coverage=None``: absent evidence is not negative evidence.
"""

import logging
from contextlib import closing
from math import comb

from dataIngestion.client import FULL_MARKERS, CellxGeneClient
from prediction.schemas import GeneEvidence, PredictionResult, Verification

log = logging.getLogger(__name__)


def hyper_tail(population, successes, draws, observed) -> float:
    """P(X >= observed) for Hypergeometric(successes of population, draws)."""
    if not 0 <= observed <= draws:
        return 1.0
    total = comb(population, draws)
    if total == 0:
        return 1.0
    lo, hi = max(observed, 0), min(draws, successes)
    tail = sum(
        comb(successes, i) * comb(population - successes, draws - i)
        for i in range(lo, hi + 1)
    )
    return min(1.0, tail / total)


class MarkerVerifier:
    """Annotate stage-one predictions with how many input genes are markers of the cell type.

    ``coverage`` (fraction of input genes that are markers) is size-dependent: more
    input genes mean more chances to overlap the fixed marker list. ``enrichment``
    divides the observed matches by the matches expected when the same number of
    genes is drawn at random from all genes expressed in the tissue, and
    ``p_value`` is the exact hypergeometric tail of that null. Both are
    comparable across input sizes; neither reorders the predictions.
    """

    def __init__(self, client: CellxGeneClient, n_markers=FULL_MARKERS):
        self.client = client
        self.n_markers = n_markers
        self._tissue_universe: dict[str, set[str]] = {}

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
        universe = self._universe(result.tissue_id)
        draws = sum(1 for gene in result.input_genes if gene in universe)
        enrichment = p_value = None
        if draws and universe:
            expected = draws * len(markers) / len(universe)
            if expected > 0:
                enrichment = matched / expected
            p_value = hyper_tail(len(universe), len(markers), draws, matched)
        return Verification(
            coverage=matched / total if total else None,
            matched=matched,
            total=total,
            verified=matched > 0,
            enrichment=enrichment,
            p_value=p_value,
            evidence=evidence,
        )

    def _universe(self, tissue_id) -> set[str]:
        """Every gene with expression data in the tissue: the random-null population."""
        cached = self._tissue_universe.get(tissue_id)
        if cached is None:
            with closing(self.client._connect()) as conn:
                cached = {
                    row[0]
                    for row in conn.execute(
                        "SELECT DISTINCT gene_id FROM expression WHERE tissue_id = ?",
                        (tissue_id,),
                    )
                }
            self._tissue_universe[tissue_id] = cached
        return cached
