"""Fetch expression for every gene that appears as a marker somewhere.

Training uses marker genes as ground truth, so a model trained on a subset of
them is biased toward whichever tissues that subset happened to cover. This
fills the cache with all of them, resuming from whatever is already stored.

    uv run python -m dataIngestion.fetch_all_genes --db cellxgene.sqlite
"""

import argparse
import logging
import time
from contextlib import closing

from dataIngestion.client import CellxGeneClient

log = logging.getLogger("fetch_all_genes")

BATCH_SIZE = 6
RETRY_SLEEP = 5.0


def missing_genes(client, limit=None):
    with closing(client._connect()) as conn:
        rows = conn.execute(
            """
            SELECT DISTINCT gene_id FROM marker
            EXCEPT
            SELECT DISTINCT gene_id FROM expression
            ORDER BY gene_id
            """
        ).fetchall()
    genes = [row[0] for row in rows]
    return genes[:limit] if limit else genes


def run(client, batch_size=BATCH_SIZE, limit=None, sleep=0.0):
    genes = missing_genes(client, limit)
    total = len(genes)
    log.info("fetching %d genes in batches of %d", total, batch_size)
    started = time.time()
    fetched = 0
    failures = 0

    for start in range(0, total, batch_size):
        batch = genes[start : start + batch_size]
        try:
            client.get_expression(batch)
            fetched += len(batch)
        except Exception as exc:
            failures += len(batch)
            log.warning("batch %d failed (%s): %s", start // batch_size, type(exc).__name__, exc)
            time.sleep(RETRY_SLEEP)

        done = start + len(batch)
        if done % (batch_size * 10) < batch_size:
            elapsed = time.time() - started
            rate = done / elapsed if elapsed else 0.0
            remaining = (total - done) / rate if rate else 0.0
            log.info(
                "%d/%d genes | %.1f min elapsed | %.2f genes/s | ~%.0f min left",
                done,
                total,
                elapsed / 60,
                rate,
                remaining / 60,
            )
        if sleep:
            time.sleep(sleep)

    elapsed = time.time() - started
    log.info(
        "done: %d fetched, %d failed, %.1f min",
        fetched,
        failures,
        elapsed / 60,
    )
    return fetched, failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="cellxgene.sqlite")
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--sleep", type=float, default=0.0)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    run(
        CellxGeneClient(db_path=args.db),
        batch_size=args.batch_size,
        limit=args.limit,
        sleep=args.sleep,
    )


if __name__ == "__main__":
    main()
