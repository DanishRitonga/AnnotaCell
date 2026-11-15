"""Feature construction for cell type prediction.

No I/O happens here: both the offline trainer and the online predictor feed the
same row tuples through :func:`score_rows`, so a model can never be trained on
features computed differently from the ones it is served.
"""

from dataclasses import dataclass
from math import ceil, exp, log, pi, sqrt

from dataIngestion.models import TISSUE_STATS

DEFAULT_Z = 1.959963984540054
DEFAULT_MIN_PRESENT_FRAC = 0.6
VOTE_THRESHOLD = 0.5


def decision_scores(model, x):
    """Stable log-likelihood ratio P(class=1|x)/P(class=0|x) for a fitted GaussianNB.

    ``predict_proba`` underflows to exactly 0.0 for unlikely candidates, which makes their
    ordering arbitrary. Ranking on this ratio is mathematically identical but never saturates.
    """
    classes = list(model.classes_)
    pos = classes.index(1)
    other = 0 if pos == 1 else 1
    theta, var, prior = model.theta_, model.var_, model.class_prior_
    half_log_2pi = log(2.0 * pi) / 2.0
    out = []
    for row in x:
        joint = []
        for c in (other, pos):
            total = log(prior[c])
            for value, mean, variance in zip(row, theta[c], var[c]):
                total -= half_log_2pi + log(variance) / 2.0 + (value - mean) ** 2 / (2.0 * variance)
            joint.append(total)
        out.append(joint[1] - joint[0])
    return out


def sigmoid(x):
    """Numerically stable logistic, used to turn log-likelihood ratios into probabilities."""

    def one(value):
        if value >= 0:
            return 1.0 / (1.0 + exp(-value))
        e = exp(value)
        return e / (1.0 + e)

    if isinstance(x, (list, tuple)):
        return [one(value) for value in x]
    return one(x)


def wilson(count: float, total: float, z: float = DEFAULT_Z) -> float:
    """Lower bound of the Wilson score interval for ``count / total``.

    Used instead of the raw proportion because ``pc`` saturates at 1.0 for cell
    types with a handful of cells: one expressing cell out of one is not the
    same evidence as five hundred out of five hundred.
    """
    if not total or total <= 0:
        return 0.0
    p = count / total
    denom = 1.0 + z * z / total
    centre = p + z * z / (2.0 * total)
    margin = z * sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total))
    return max(0.0, (centre - margin) / denom)


def minimum_present(k: int, frac: float = DEFAULT_MIN_PRESENT_FRAC) -> int:
    """How many of ``k`` input genes must be expressed by a candidate to be eligible."""
    if k <= 0:
        return 0
    return min(k, ceil(frac * k))


def per_gene_max(rows) -> dict[str, float]:
    """Largest mean expression observed for each gene, over every tissue and cell type."""
    gmax: dict[str, float] = {}
    for gene, _cell, me, _n, _pc in rows:
        if me is None:
            continue
        if me > gmax.get(gene, float("-inf")):
            gmax[gene] = me
    return gmax


def per_gene_max_from_summary(expression_summary, genes) -> dict[str, float]:
    """Same as :func:`per_gene_max` but read from an ``expression_summary`` mapping."""
    gmax: dict[str, float] = {}
    for gene in genes:
        best = float("-inf")
        for cells in expression_summary.get(gene, {}).values():
            for cell_id, node in cells.items():
                if cell_id == TISSUE_STATS:
                    continue
                me = node.aggregated.me
                if me is not None and me > best:
                    best = me
        if best > float("-inf"):
            gmax[gene] = best
    return gmax


@dataclass(frozen=True)
class CellFeature:
    """Aggregated evidence for one candidate cell type."""

    cell_id: str
    pc_wilson: float
    me_norm: float
    n_present: int
    vote_frac: float = 0.0


def score_rows(rows, genes, gmax) -> dict[str, CellFeature]:
    """Aggregate ``(gene, cell, me, n, pc)`` rows into one feature vector per cell type.

    ``pc_wilson`` and ``me_norm`` are means over the input genes the cell type
    actually expresses; ``n_present`` is how many those were and ``vote_frac``
    the share of those genes clearing :data:`VOTE_THRESHOLD`. The vote makes the
    score robust to a few contaminant genes that pull the mean elsewhere.
    """
    wanted = set(genes)
    sums: dict[str, list[float]] = {}
    for gene, cell, me, n, pc in rows:
        if gene not in wanted or not pc or pc <= 0 or not n or n <= 0:
            continue
        acc = sums.get(cell)
        if acc is None:
            acc = sums[cell] = [0.0, 0.0, 0.0, 0.0]
        bound = wilson(n, n / pc)
        acc[0] += bound
        peak = gmax.get(gene)
        acc[1] += (me or 0.0) / peak if peak and peak > 0 else 0.0
        acc[2] += 1.0 if bound >= VOTE_THRESHOLD else 0.0
        acc[3] += 1.0

    features = {}
    for cell, (pw, mn, votes, count) in sums.items():
        features[cell] = CellFeature(
            cell_id=cell,
            pc_wilson=pw / count,
            me_norm=mn / count,
            n_present=int(count),
            vote_frac=votes / count,
        )
    return features


def score_cells(expression_summary, genes, tissue, gmax, allowed=None) -> dict[str, CellFeature]:
    """Collect one tissue's rows out of an ``expression_summary`` mapping and score them."""
    rows = []
    for gene in genes:
        cells = expression_summary.get(gene, {}).get(tissue, {})
        for cell_id, node in cells.items():
            if cell_id == TISSUE_STATS:
                continue
            if allowed is not None and cell_id not in allowed:
                continue
            metric = node.aggregated
            rows.append((gene, cell_id, metric.me, metric.n, metric.pc))
    return score_rows(rows, genes, gmax)


def order_candidates(
    features: dict[str, CellFeature], scores: dict[str, float], min_present: int
):
    """Order candidates by evidence, keeping gated ones below eligible ones.

    Returns a list of ``(cell_id, eligible)`` tuples. Gating is a ranking tier
    rather than a filter so a request can never come back empty.
    """
    ordered = []
    for cell_id, feature in features.items():
        eligible = feature.n_present >= min_present
        ordered.append((not eligible, -scores[cell_id], -feature.n_present, cell_id))
    ordered.sort()
    return [(cell_id, not blocked) for blocked, _s, _n, cell_id in ordered]


@dataclass(frozen=True)
class LineageNode:
    """One level of a prediction's ontology chain, root first."""

    cell_id: str
    depth: int
    support: float
    rival: float
    margin: float
    is_leaf: bool


def build_lineage(
    scores: dict[str, float],
    winner: str,
    ancestors_of,
    path_of,
    children_of,
    known: set[str],
    is_leaf,
) -> list[LineageNode]:
    """Aggregate leaf scores up the winner's ancestor chain.

    Leaves are scored once; a term's ``support`` is the summed score of the
    scored leaves beneath it, so terms at the same level are compared on the
    same scale instead of against subtree-inflated ancestor rows. ``rival`` is
    the strongest sibling at that level, and ``margin`` the difference: how
    decisively this branch beats its alternatives.
    """

    total = sum(scores.values()) or 1.0

    def support(term: str) -> float:
        inside = sum(
            value
            for cell, value in scores.items()
            if cell == term or term in ancestors_of(cell)
        )
        return inside / total

    chain = list(path_of(winner))
    nodes = []
    for depth, term in enumerate(chain):
        path = path_of(term)
        rival = 0.0
        if len(path) > 1:
            for sibling in children_of(path[-2]):
                if sibling == term or sibling not in known:
                    continue
                rival = max(rival, support(sibling))
        value = support(term)
        nodes.append(
            LineageNode(
                cell_id=term,
                depth=depth,
                support=value,
                rival=rival,
                margin=value - rival,
                is_leaf=is_leaf(term),
            )
        )
    return nodes
