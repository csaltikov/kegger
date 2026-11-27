"""Batch downloads from KEGG with ``kegger.batch_tools``.

``basic_usage.py`` fetches one record per call. This example uses the async
batch layer, which fetches many records concurrently behind a persistent
cache -- useful when you need whole genomes or dozens of pathways.

Run it twice: the second run is served from ``kegg_async_cache.sqlite`` and
returns almost instantly.

    python examples/batch_usage.py
"""

import logging
import time

import pandas as pd

import kegger as kg
from kegger.batch_tools import get_paths, org_annotation

# Failed downloads are logged and skipped, never raised. Turn the log on so
# you can see which identifiers were dropped and why.
logging.basicConfig(level=logging.WARNING, format="[%(levelname)s] %(message)s")

# Show every column of the wide annotation tables.
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 160)
pd.set_option("display.max_colwidth", 45)

ORGS = ["eco", "shn", "bsu"]  # E. coli K-12, Shewanella sp. ANA-3, B. subtilis

# The last one is deliberately bogus, to show that a bad identifier
# does not take down the rest of the batch.
PATHWAYS = ["shn03020", "shn02024", "shn00910", "shn99999"]


def annotate_genomes() -> pd.DataFrame:
    """Download gene annotations for several organisms at once."""
    print(f"--- Annotating {len(ORGS)} genomes: {', '.join(ORGS)} ---")

    start = time.perf_counter()
    genomes = org_annotation(ORGS)
    elapsed = time.perf_counter() - start

    print(f"{len(genomes):,} genes in {elapsed:.1f}s")
    print("\nGenes per organism:")
    print(genomes["org"].value_counts().to_string())
    print("\nSame columns as kg.get_org(), plus 'org':")
    print(genomes.head())

    # Cached responses cost nothing, so a repeat call is effectively free.
    start = time.perf_counter()
    org_annotation(ORGS)
    print(f"\nSecond call (from cache): {time.perf_counter() - start:.1f}s")

    return genomes


def summarize_pathways() -> list[dict]:
    """Download pathway records and parse them into structured entries."""
    print(f"\n--- Downloading {len(PATHWAYS)} pathway records ---")

    records = get_paths(PATHWAYS)
    skipped = len(PATHWAYS) - len(records)
    print(f"{len(records)} of {len(PATHWAYS)} downloaded, {skipped} skipped "
          f"(bad identifiers are reported on the warning log)")

    # Each record holds the raw flat file under 'raw'; kegg_parser turns it
    # into a dict of KEGG fields (ENTRY, NAME, GENE, ORTHOLOG, ...).
    entries = [kg.kegg_parser(record["raw"]) for record in records]

    print(f"\n{'Pathway':<12} {'Genes':>6}  Name")
    for entry in entries:
        pathid = entry["ENTRY"][0]
        print(f"{pathid:<12} {len(entry['GENE']):>6}  {entry['NAME']}")

    return entries


def genes_in_pathway(genomes: pd.DataFrame, entry: dict) -> pd.DataFrame:
    """Join one parsed pathway back onto the batch annotation table."""
    pathid, name = entry["ENTRY"][0], entry["NAME"]
    org = pathid[:3]

    print(f"\n--- Annotations for the {len(entry['GENE'])} genes in {pathid} ---")
    print(f"{name}\n")

    # 'gene' in the annotation table and 'GENE' in the parsed record both use
    # locus tags, so they join directly.
    pathway_genes = genomes[
        (genomes["org"] == org) & (genomes["gene"].isin(entry["GENE"]))
    ]
    print(pathway_genes[["gene", "position", "annotation"]].head(10).to_string(index=False))

    return pathway_genes


def main() -> None:
    genomes = annotate_genomes()
    entries = summarize_pathways()
    genes_in_pathway(genomes, entries[0])

    # In a notebook an event loop is already running, so call the async
    # versions instead and await them:
    #
    #     from kegger.batch_tools import org_annotation_async
    #     genomes = await org_annotation_async(ORGS)


if __name__ == "__main__":
    main()
