from typing import Any

from pydantic import BaseModel, ConfigDict, Field

TISSUE_STATS = "tissue_stats"


class Marker(BaseModel):
    model_config = ConfigDict(extra="ignore")

    gene_ontology_term_id: str
    marker_score: float | None = None
    specificity: float | None = None


class MarkerList(BaseModel):
    model_config = ConfigDict(extra="ignore")

    marker_genes: list[Marker] = Field(default_factory=list)
    snapshot_id: str | None = None


class ExpressionMetric(BaseModel):
    model_config = ConfigDict(extra="ignore")

    me: float | None = None
    n: int | None = None
    pc: float | None = None
    tpc: float | None = None


class AggregatedCell(BaseModel):
    model_config = ConfigDict(extra="ignore")

    aggregated: ExpressionMetric = Field(default_factory=ExpressionMetric)


class ExpressionRow(BaseModel):
    model_config = ConfigDict(extra="ignore")

    gene: str
    tissue: str
    cell: str
    expression: float | None
    cell_count: int | None
    cell_percentage: float | None
    tissue_composition: float | None


class CellLabelAggregate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    cell_type_ontology_term_id: str | None = None
    name: str | None = None


class CellLabel(BaseModel):
    model_config = ConfigDict(extra="ignore")

    aggregated: CellLabelAggregate = Field(default_factory=CellLabelAggregate)


TERM_ID_LABELS = "term_id_labels"


class ExpressionSummary(BaseModel):
    model_config = ConfigDict(extra="ignore")

    expression_summary: dict[str, dict[str, dict[str, AggregatedCell]]] = Field(default_factory=dict)
    term_id_labels: dict[str, Any] = Field(default_factory=dict)
    snapshot_id: str | None = None

    def cell_labels(self) -> dict[str, str]:
        """Map every observed cell type id to its name, skipping tissue-level rollup rows."""
        labels: dict[str, str] = {}
        for cells in (self.term_id_labels.get("cell_types") or {}).values():
            for cell_id, raw in cells.items():
                if cell_id == TISSUE_STATS:
                    continue
                label = CellLabel.model_validate(raw)
                if label.aggregated.name:
                    labels[cell_id] = label.aggregated.name
        return labels

    def to_rows(self) -> list[ExpressionRow]:
        rows: list[ExpressionRow] = []
        for gene, tissues in self.expression_summary.items():
            for tissue, cells in tissues.items():
                for cell, data in cells.items():
                    m = data.aggregated
                    rows.append(
                        ExpressionRow(
                            gene=gene,
                            tissue=tissue,
                            cell=cell,
                            expression=m.me,
                            cell_count=m.n,
                            cell_percentage=m.pc,
                            tissue_composition=m.tpc,
                        )
                    )
        return rows
