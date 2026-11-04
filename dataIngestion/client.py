import logging
import sqlite3
from contextlib import closing

import requests

try:
    from .models import (
        TISSUE_STATS,
        AggregatedCell,
        ExpressionMetric,
        ExpressionSummary,
        Marker,
        MarkerList,
    )
except ImportError:
    from models import (
        TISSUE_STATS,
        AggregatedCell,
        ExpressionMetric,
        ExpressionSummary,
        Marker,
        MarkerList,
    )

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

FILTERS_URL = "https://api.cellxgene.cziscience.com/wmg/v2/filters"
PRIMARY_FILTER_URL = (
    "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions"
)
QUERY_URL = "https://api.cellxgene.cziscience.com/wmg/v2/query"
MARKERS_URL = "https://api.cellxgene.cziscience.com/wmg/v2/markers"

ORGANISM_HUMAN = "NCBITaxon:9606"

SCHEMA = """
CREATE TABLE IF NOT EXISTS tissue (
    id   TEXT PRIMARY KEY,
    name TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS cell (
    id   TEXT PRIMARY KEY,
    name TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS gene (
    id   TEXT PRIMARY KEY,
    name TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS marker (
    tissue_id   TEXT NOT NULL REFERENCES tissue(id) ON DELETE CASCADE,
    cell_id     TEXT NOT NULL REFERENCES cell(id)   ON DELETE CASCADE,
    gene_id     TEXT NOT NULL REFERENCES gene(id)   ON DELETE CASCADE,
    marker_score REAL,
    specificity  REAL,
    snapshot_id TEXT NOT NULL,
    PRIMARY KEY (tissue_id, cell_id, gene_id)
);
CREATE TABLE IF NOT EXISTS expression (
    gene_id     TEXT NOT NULL REFERENCES gene(id)   ON DELETE CASCADE,
    tissue_id   TEXT NOT NULL REFERENCES tissue(id) ON DELETE CASCADE,
    cell_id     TEXT NOT NULL REFERENCES cell(id)   ON DELETE CASCADE,
    me  REAL,
    n   INTEGER,
    pc  REAL,
    tpc REAL,
    snapshot_id TEXT NOT NULL,
    PRIMARY KEY (gene_id, tissue_id, cell_id)
);
CREATE INDEX IF NOT EXISTS idx_gene_name   ON gene(name);
CREATE INDEX IF NOT EXISTS idx_cell_name   ON cell(name);
CREATE INDEX IF NOT EXISTS idx_tissue_name ON tissue(name);
CREATE INDEX IF NOT EXISTS idx_marker_snapshot     ON marker(snapshot_id);
CREATE INDEX IF NOT EXISTS idx_expression_snapshot ON expression(snapshot_id);
CREATE INDEX IF NOT EXISTS idx_expression_gene     ON expression(gene_id, snapshot_id);
"""


class CellxGeneClient:
    def __init__(
        self,
        db_path="ontology.sqlite",
        timeout=120,
        session=None,
        organism=ORGANISM_HUMAN,
    ):
        self.db_path = db_path
        self.timeout = timeout
        self.session = session or requests.Session()
        self.organism = organism
        self._init_db()

    def _init_db(self):
        with closing(self._connect()) as conn:
            self._migrate(conn)
            conn.executescript(SCHEMA)

    @staticmethod
    def _migrate(conn):
        """Drop derived caches created before snapshot tracking so they can be rebuilt."""
        for table in ("marker", "expression"):
            cols = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
            if cols and "snapshot_id" not in cols:
                conn.execute(f"DROP TABLE {table}")
        conn.commit()

    def _connect(self):
        conn = sqlite3.connect(self.db_path)
        conn.execute("PRAGMA foreign_keys = ON")
        return conn

    def _post(self, url, payload):
        response = self.session.post(url, json=payload, timeout=self.timeout)
        response.raise_for_status()
        return response.json()

    def _get(self, url):
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        return response.json()

    def get_ontology(self, refresh=False):
        """Fetch tissue/cell/gene vocabularies from CELLxGENE and persist to SQLite."""
        with closing(self._connect()) as conn:
            stored = self._get_meta(conn, "organism")
            if not refresh:
                if stored and stored != self.organism:
                    raise ValueError(
                        f"Database {self.db_path} holds ontology for organism {stored}, "
                        f"but client was configured for {self.organism}. "
                        "Use a different db_path or call get_ontology(refresh=True)."
                    )
                counts = self._counts(conn)
                if all(counts.values()) and stored == self.organism:
                    logging.info(
                        "Ontology already populated %(tissues)d tissues, %(cells)d cell types, %(genes)d genes; skipping fetch.",
                        counts,
                    )
                    return counts

            logging.info("Fetching ontology from CELLxGENE...")
            filters = self._post(
                FILTERS_URL, {"filter": {"organism_ontology_term_id": self.organism}}
            )
            dims = filters.get("filter_dims", {})
            tissues = [
                (k, v) for d in dims.get("tissue_terms", []) for k, v in d.items()
            ]
            cells = [
                (k, v) for d in dims.get("cell_type_terms", []) for k, v in d.items()
            ]

            primary = self._get(PRIMARY_FILTER_URL)
            genes = [
                (k, v)
                for d in primary.get("gene_terms", {}).get(self.organism, [])
                for k, v in d.items()
            ]

            conn.execute("DELETE FROM tissue")
            conn.execute("DELETE FROM cell")
            conn.execute("DELETE FROM gene")
            conn.execute("DELETE FROM marker")
            conn.execute("DELETE FROM expression")
            conn.executemany(
                "INSERT OR REPLACE INTO tissue(id, name) VALUES (?, ?)", tissues
            )

            conn.executemany(
                "INSERT OR REPLACE INTO cell(id, name) VALUES (?, ?)", cells
            )
            conn.executemany(
                "INSERT OR REPLACE INTO gene(id, name) VALUES (?, ?)", genes
            )
            conn.execute(
                "INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
                ("organism", self.organism),
            )
            if primary.get("snapshot_id") is not None:
                conn.execute(
                    "INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
                    ("snapshot_id", str(primary["snapshot_id"])),
                )
            conn.commit()

            counts = self._counts(conn)
            logging.info(
                "Stored %(tissues)d tissues, %(cells)d cell types, %(genes)d genes.",
                counts,
            )
            return counts

    @staticmethod
    def _counts(conn):
        def n(table):
            return conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]

        return {"tissues": n("tissue"), "cells": n("cell"), "genes": n("gene")}

    @staticmethod
    def _get_meta(conn, key):
        row = conn.execute("SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
        return row[0] if row else None

    def snapshot_id(self, conn=None):
        """Return the CELLxGENE snapshot id the cached ontology was built from."""
        if conn is None:
            with closing(self._connect()) as own:
                return self._get_meta(own, "snapshot_id")
        return self._get_meta(conn, "snapshot_id")

    def id_for(self, table, name):
        with closing(self._connect()) as conn:
            row = conn.execute(
                f"SELECT id FROM {table} WHERE name = ? COLLATE NOCASE", (name,)
            ).fetchone()
        return row[0] if row else None

    def name_for(self, table, term_id):
        with closing(self._connect()) as conn:
            row = conn.execute(
                f"SELECT name FROM {table} WHERE id = ?", (term_id,)
            ).fetchone()
        return row[0] if row else None

    def gene_to_ensembl(self, symbol):
        return self.id_for("gene", symbol)

    def ensembl_to_gene(self, ensembl_id):
        return self.name_for("gene", ensembl_id)

    def marker_to_gene(self, ensembl_id):
        return self.name_for("gene", ensembl_id)

    def get_markers(
        self, tissue_id, cell_id, n_markers=50, test="ttest"
    ) -> list[Marker]:
        """Fetch marker genes for a (tissue, cell type) pair, validate and cache them in SQLite."""
        cached = self._cached_markers(tissue_id, cell_id)
        if cached is not None:
            return cached

        payload = {
            "tissue": tissue_id,
            "celltype": cell_id,
            "n_markers": n_markers,
            "organism": self.organism,
            "test": test,
        }
        parsed = MarkerList.model_validate(self._post(MARKERS_URL, payload))
        markers = parsed.marker_genes
        with closing(self._connect()) as conn:
            current = self.snapshot_id(conn)
            snapshot = parsed.snapshot_id or current
            conn.execute(
                "DELETE FROM marker WHERE tissue_id = ? AND cell_id = ?",
                (tissue_id, cell_id),
            )
            conn.executemany(
                "INSERT OR REPLACE INTO marker(tissue_id, cell_id, gene_id, marker_score, specificity, snapshot_id) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                [
                    (
                        tissue_id,
                        cell_id,
                        m.gene_ontology_term_id,
                        m.marker_score,
                        m.specificity,
                        snapshot,
                    )
                    for m in markers
                ],
            )
            conn.commit()
        return markers

    def _cached_markers(self, tissue_id, cell_id) -> list[Marker] | None:
        with closing(self._connect()) as conn:
            current = self.snapshot_id(conn)
            if current is None:
                return None
            rows = conn.execute(
                "SELECT gene_id, marker_score, specificity, snapshot_id FROM marker "
                "WHERE tissue_id = ? AND cell_id = ? ORDER BY marker_score DESC",
                (tissue_id, cell_id),
            ).fetchall()
        if not rows or any(r[3] != current for r in rows):
            return None
        return [
            Marker(gene_ontology_term_id=r[0], marker_score=r[1], specificity=r[2])
            for r in rows
        ]

    def get_expression(self, gene_ids, is_rollup=True) -> ExpressionSummary:
        """Return the expression summary for a batch of Ensembl gene ids.

        Served from the SQLite cache when available; otherwise fetched, sanitized to the
        stored shape, and persisted. Both paths return identical structures.
        """
        gene_ids = list(gene_ids)
        cached = self._cached_expression(gene_ids, is_rollup)
        if cached is not None:
            return cached

        payload = {
            "filter": {
                "dataset_ids": [],
                "development_stage_ontology_term_ids": [],
                "disease_ontology_term_ids": [],
                "gene_ontology_term_ids": gene_ids,
                "organism_ontology_term_id": self.organism,
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": [],
                "publication_citations": [],
            },
            "is_rollup": is_rollup,
        }
        fetched = ExpressionSummary.model_validate(self._post(QUERY_URL, payload))
        self._observe_cell_types(fetched)

        with closing(self._connect()) as conn:
            current = self.snapshot_id(conn)
            snapshot = fetched.snapshot_id or current
            if current is not None and snapshot is not None and snapshot != current:
                logging.warning(
                    "Expression snapshot %s differs from cached ontology snapshot %s; "
                    "returning sanitized rows without caching.",
                    snapshot,
                    current,
                )
                return self._sanitize(fetched)

            rows = self._expression_rows(fetched, snapshot)
            if rows:
                conn.executemany(
                    "INSERT OR REPLACE INTO expression(gene_id, tissue_id, cell_id, me, n, pc, tpc, snapshot_id) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    rows,
                )
                conn.commit()

        return self._cached_expression(gene_ids, is_rollup) or self._sanitize(fetched)

    def _sanitize(self, summary: ExpressionSummary) -> ExpressionSummary:
        """Reduce an API summary to the rows the cache can represent."""
        known = {
            "gene": {r[0] for r in self._ids("gene")},
            "tissue": {r[0] for r in self._ids("tissue")},
            "cell": {r[0] for r in self._ids("cell")},
        }
        kept: dict[str, dict[str, dict[str, AggregatedCell]]] = {}
        for row in summary.to_rows():
            if (
                row.cell == TISSUE_STATS
                or row.gene not in known["gene"]
                or row.tissue not in known["tissue"]
                or row.cell not in known["cell"]
            ):
                continue
            kept.setdefault(row.gene, {}).setdefault(row.tissue, {})[row.cell] = (
                AggregatedCell(
                    aggregated=ExpressionMetric(
                        me=row.expression,
                        n=row.cell_count,
                        pc=row.cell_percentage,
                        tpc=row.tissue_composition,
                    )
                )
            )
        return ExpressionSummary(
            expression_summary=kept,
            term_id_labels=summary.term_id_labels,
            snapshot_id=summary.snapshot_id,
        )

    def _cached_expression(self, gene_ids, is_rollup) -> ExpressionSummary | None:
        """Rebuild an ExpressionSummary from the cache when every gene is present at this snapshot."""
        if not is_rollup or not gene_ids:
            return None
        with closing(self._connect()) as conn:
            current = self.snapshot_id(conn)
            if current is None:
                return None
            placeholders = ",".join("?" * len(gene_ids))
            rows = conn.execute(
                "SELECT gene_id, tissue_id, cell_id, me, n, pc, tpc, snapshot_id FROM expression "
                f"WHERE gene_id IN ({placeholders})",
                gene_ids,
            ).fetchall()

        if not rows or any(r[7] != current for r in rows):
            return None
        if {r[0] for r in rows} != set(gene_ids):
            return None

        summary: dict[str, dict[str, dict[str, AggregatedCell]]] = {}
        for gene, tissue, cell, me, n, pc, tpc, _ in rows:
            summary.setdefault(gene, {}).setdefault(tissue, {})[cell] = AggregatedCell(
                aggregated=ExpressionMetric(me=me, n=n, pc=pc, tpc=tpc)
            )
        return ExpressionSummary(expression_summary=summary, snapshot_id=current)

    def _observe_cell_types(self, summary: ExpressionSummary) -> None:
        """Persist cell types the /query rollup observed but /filters never listed."""
        observed = summary.cell_labels()
        if not observed:
            return
        with closing(self._connect()) as conn:
            conn.executemany(
                "INSERT OR IGNORE INTO cell(id, name) VALUES (?, ?)",
                list(observed.items()),
            )
            conn.commit()

    def _expression_rows(self, summary: ExpressionSummary, snapshot) -> list[tuple]:
        """Flatten the summary into FK-valid expression rows.

        tissue_stats rows are dropped because they are tissue-level rollups, not cell types.
        """
        known = {
            "gene": {r[0] for r in self._ids("gene")},
            "tissue": {r[0] for r in self._ids("tissue")},
            "cell": {r[0] for r in self._ids("cell")},
        }
        rows = []
        dropped = 0
        for row in summary.to_rows():
            if row.cell == TISSUE_STATS:
                dropped += 1
                continue
            if (
                row.gene not in known["gene"]
                or row.tissue not in known["tissue"]
                or row.cell not in known["cell"]
            ):
                dropped += 1
                continue
            rows.append(
                (
                    row.gene,
                    row.tissue,
                    row.cell,
                    row.expression,
                    row.cell_count,
                    row.cell_percentage,
                    row.tissue_composition,
                    snapshot,
                )
            )
        if dropped:
            logging.info(
                "Dropped %d expression rows (tissue_stats rollups or terms absent from the ontology).",
                dropped,
            )
        return rows

    def _ids(self, table):
        with closing(self._connect()) as conn:
            return conn.execute(f"SELECT id FROM {table}").fetchall()
