"""Offline trainer: /markers -> features -> Gaussian Naive Bayes -> model.json.

The trainer never talks to the network. It reads the marker and expression
tables already cached in SQLite, derives a training set from the markers of leaf
cell types, and evaluates the model on a group-level hold-out split.

    uv run python -m prediction.train --db cellxgene.sqlite --out prediction/model.json
"""

import argparse
import json
import random
import sqlite3
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.naive_bayes import GaussianNB

from dataIngestion.client import CellxGeneClient
from prediction.features import (
    DEFAULT_MIN_PRESENT_FRAC,
    DEFAULT_Z,
    CellFeature,
    decision_scores,
    minimum_present,
    order_candidates,
    sigmoid,
    score_rows,
)

TOP_K_MARKERS = 5
TRAIN_FRACTION = 0.7
MODEL_VERSION = 1

NOISE_MODES = {
    "clean": (0, None),
    "swap1_ambient": (1, "ambient"),
    "swap2_ambient": (2, "ambient"),
    "swap1_crosstype": (1, "crosstype"),
    "swap2_crosstype": (2, "crosstype"),
}

ZERO = CellFeature("", 0.0, 0.0, 0)


def load_markers(con, tissue_leaves, expr_genes, k=TOP_K_MARKERS):
    """Top-K marker genes per (tissue, leaf cell type), restricted to genes with expression data."""
    by_pair = defaultdict(list)
    for tissue, cell, gene in con.execute(
        "SELECT tissue_id, cell_id, gene_id FROM marker "
        "ORDER BY tissue_id, cell_id, marker_score DESC"
    ):
        by_pair[(tissue, cell)].append(gene)

    inputs = {}
    for (tissue, cell), genes in by_pair.items():
        if cell not in tissue_leaves.get(tissue, ()):
            continue
        usable = [g for g in genes if g in expr_genes][:k]
        if len(usable) == k:
            inputs[(tissue, cell)] = usable
    return inputs


def load_expression(con, tissues, tissue_leaves):
    """Leaf-cell expression rows grouped by (tissue, gene), plus per-gene peak expression."""
    store = defaultdict(list)
    gmax = {}
    placeholders = ",".join("?" * len(tissues))
    for tissue, gene, cell, me, n, pc in con.execute(
        "SELECT tissue_id, gene_id, cell_id, me, n, pc FROM expression "
        f"WHERE tissue_id IN ({placeholders})",
        tissues,
    ):
        if cell not in tissue_leaves.get(tissue, ()):
            continue
        store[(tissue, gene)].append((cell, me, n, pc))
        if me is not None and me > gmax.get(gene, float("-inf")):
            gmax[gene] = me
    return store, gmax


def load_present(con, tissue_leaves):
    """Leaf cell types with at least one expression row, per tissue."""
    present = defaultdict(set)
    for tissue, cell in con.execute("SELECT DISTINCT tissue_id, cell_id FROM expression"):
        if cell in tissue_leaves.get(tissue, ()):
            present[tissue].add(cell)
    return present


def load_tissue_genes(con, tissues):
    genes = defaultdict(list)
    placeholders = ",".join("?" * len(tissues))
    for tissue, gene in con.execute(
        f"SELECT DISTINCT tissue_id, gene_id FROM expression WHERE tissue_id IN ({placeholders})",
        tissues,
    ):
        genes[tissue].append(gene)
    return genes


def make_inputs(mode, tissue, truth, genes, rng, by_pair, expr_genes, tissue_genes, present):
    """Apply one noise mode to a group's marker genes."""
    swaps, kind = NOISE_MODES[mode]
    if not swaps:
        return list(genes)

    genes = list(genes)
    used = set(genes)
    for _ in range(swaps):
        replacement = None
        if kind == "ambient":
            own = set(by_pair.get((tissue, truth), ()))
            pool = [g for g in tissue_genes if g not in used and g not in own]
            if pool:
                replacement = rng.choice(pool)
        else:
            donors = [
                c for c in present
                if c != truth and (tissue, c) in by_pair
            ]
            if donors:
                pool = [g for g in by_pair[(tissue, rng.choice(donors))] if g not in used]
                if pool:
                    replacement = rng.choice(pool)
        if replacement is None:
            continue
        position = rng.randrange(len(genes))
        used.discard(genes[position])
        genes[position] = replacement
        used.add(replacement)
    return genes


def build_samples(groups, inputs, store, gmax, present, tissue_genes, by_pair, expr_genes, rng):
    """One sample per (group, noise mode): the input genes and every candidate's features."""
    samples = []
    for tissue, truth in groups:
        candidates = present.get(tissue, set()) | {truth}
        for mode in NOISE_MODES:
            genes = make_inputs(
                mode, tissue, truth, inputs[(tissue, truth)], rng,
                by_pair, expr_genes, tissue_genes.get(tissue, []), present.get(tissue, set()),
            )
            rows = []
            for gene in genes:
                for cell, me, n, pc in store.get((tissue, gene), ()):
                    rows.append((gene, cell, me, n, pc))
            scored = score_rows(rows, genes, gmax)
            features = {}
            for cell in candidates:
                features[cell] = scored.get(cell) or CellFeature(cell, 0.0, 0.0, 0)
            samples.append(((tissue, truth), mode, genes, features))
    return samples


def evaluate(model, samples, test_keys):
    """Rank held-out groups under every noise mode and aggregate the metrics."""
    results = {}
    for mode in NOISE_MODES:
        ranks = []
        for key, sample_mode, genes, features in samples:
            if sample_mode != mode or key not in test_keys:
                continue
            truth = key[1]
            matrix = np.array(
                [[f.pc_wilson, f.me_norm, f.vote_frac] for f in features.values()],
                dtype=float,
            )
            scores = sigmoid(decision_scores(model, matrix))
            score_map = dict(zip(features, scores))
            ordered = order_candidates(features, score_map, minimum_present(len(genes)))
            ids = [cell for cell, _ in ordered]
            rank = ids.index(truth) + 1
            ranks.append(rank)
        if not ranks:
            continue
        results[mode] = {
            "groups": len(ranks),
            "top1": sum(r == 1 for r in ranks) / len(ranks),
            "r5": sum(r <= 5 for r in ranks) / len(ranks),
            "r10": sum(r <= 10 for r in ranks) / len(ranks),
            "mrr": sum(1.0 / r for r in ranks) / len(ranks),
            "medrank": statistics.median(ranks),
        }
    return results


def main():
    parser = argparse.ArgumentParser(description="Train the cell type predictor.")
    parser.add_argument("--db", default="cellxgene.sqlite")
    parser.add_argument("--out", default=str(Path(__file__).with_name("model.json")))
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--markers-per-group", type=int, default=TOP_K_MARKERS)
    args = parser.parse_args()

    client = CellxGeneClient(db_path=args.db)
    con = sqlite3.connect(args.db)

    marker_tissues = [
        row[0] for row in con.execute("SELECT DISTINCT tissue_id FROM marker")
    ]
    tissue_leaves = {t: set(client.tissue_leaves(t)) for t in marker_tissues}

    expr_genes = {row[0] for row in con.execute("SELECT DISTINCT gene_id FROM expression")}
    k = args.markers_per_group
    inputs = load_markers(con, tissue_leaves, expr_genes, k)
    if not inputs:
        raise SystemExit(
            f"no training groups: no leaf cell type has {k} markers with expression data"
        )

    tissues = sorted({tissue for tissue, _cell in inputs})
    store, gmax = load_expression(con, tissues, tissue_leaves)
    present = load_present(con, tissue_leaves)
    tissue_genes = load_tissue_genes(con, tissues)

    groups = sorted(inputs)
    shuffled = list(groups)
    random.Random(args.seed).shuffle(shuffled)
    cut = int(len(shuffled) * TRAIN_FRACTION)
    train_keys = set(shuffled[:cut])
    test_keys = set(shuffled[cut:])

    samples = build_samples(
        groups, inputs, store, gmax, present, tissue_genes, inputs, expr_genes,
        random.Random(args.seed),
    )

    X, y = [], []
    for key, _mode, _genes, features in samples:
        if key not in train_keys:
            continue
        for feature in features.values():
            X.append([feature.pc_wilson, feature.me_norm, feature.vote_frac])
            y.append(1 if feature.cell_id == key[1] else 0)

    model = GaussianNB()
    model.fit(np.asarray(X, dtype=float), np.asarray(y, dtype=int))

    payload = {
        "version": MODEL_VERSION,
        "snapshot_id": client.snapshot_id(),
        "features": ["pc_wilson", "me_norm", "vote_frac"],
        "z": DEFAULT_Z,
        "min_present_frac": DEFAULT_MIN_PRESENT_FRAC,
        "markers_per_group": k,
        "classes": model.classes_.tolist(),
        "theta": model.theta_.tolist(),
        "var": model.var_.tolist(),
        "prior": model.class_prior_.tolist(),
        "trained_on": {
            "groups": len(train_keys),
            "rows": len(X),
            "positives": int(sum(y)),
            "tissues": len(tissues),
        },
        "evaluation": evaluate(model, samples, test_keys),
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n")

    print(f"wrote {args.out}")
    print(f"groups {len(groups)} train {len(train_keys)} test {len(test_keys)} tissues {len(tissues)}")
    print(f"rows {len(X)} positives {sum(y)}")
    print(f"theta {model.theta_.tolist()} prior {model.class_prior_.tolist()}")
    for mode, metrics in payload["evaluation"].items():
        print(
            f"  {mode:<16} top1 {metrics['top1']:.3f} r@5 {metrics['r5']:.3f} "
            f"r@10 {metrics['r10']:.3f} mrr {metrics['mrr']:.3f} medrank {metrics['medrank']}"
        )


if __name__ == "__main__":
    main()
