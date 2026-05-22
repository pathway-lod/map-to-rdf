#!/usr/bin/env python3
"""Export BGC↔pathway gene links to a TSV for supplementary tables.

Loads the BGC RDF (output_ttl/) and the pathway RDF, finds all cases where
a BGC cluster member gene is also present in a metabolic pathway, and writes
a TSV with full metadata per link.

Usage
-----
    python scripts/export_bgc_pathway_links.py \
        --pathways input_ttl/all-plantcyc17.0.0-gpml2021.ttl \
        --output summaries/bgc_pathway_links.tsv
"""

from __future__ import annotations

import argparse
import csv
import warnings
from pathlib import Path

from rdflib import Graph, Namespace

warnings.filterwarnings("ignore")  # suppress rdflib URI warnings for old TTL

# ============================================================
# Namespaces
# ============================================================

WP        = Namespace("http://vocabularies.wikipathways.org/wp#")
DCT       = Namespace("http://purl.org/dc/terms/")
RDFS      = Namespace("http://www.w3.org/2000/01/rdf-schema#")
PMW       = Namespace("http://rdf-plantmetwiki.bioinformatics.nl/vocab/")
OBO       = Namespace("http://purl.obolibrary.org/obo/")

DC = Namespace("http://purl.org/dc/elements/1.1/")

INIT_NS = {
    "wp":      WP,
    "dcterms": DCT,
    "dc":      DC,
    "rdfs":    RDFS,
    "pmw":     PMW,
    "obo":     OBO,
}

# ============================================================
# Queries
# ============================================================

# Direct match: BGC member gene IRI == pathway gene IRI (tair.locus)
# Works when both pathway RDF and BGC RDF use the same identifiers.org IRI.
# Two sub-patterns tried because different WikiPathways RDF versions use
# different predicates to link DataNode → pathway:
#   (a) wp:isAbout + dcterms:isPartOf  (GPMLRDF / newer WP model)
#   (b) direct: gene IRI itself has dcterms:isPartOf (older WP model)

QUERY_DIRECT_VIA_DATANODE = """
SELECT DISTINCT
    ?gene ?pathway ?pathway_title ?bgc ?bgc_source ?bgc_species ?bgc_taxon
WHERE {
    # Gene in pathway (via DataNode → wp:isAbout)
    ?dataNode a wp:GeneProduct ;
              wp:isAbout ?gene ;
              dcterms:isPartOf ?pathway .

    # Pathway title
    OPTIONAL { ?pathway rdfs:label ?pathway_title . }

    # Gene in BGC
    ?bgc obo:RO_0000051 ?gene .
    ?bgc dcterms:source ?bgc_source .
    OPTIONAL { ?bgc wp:organismName ?bgc_species . }
    OPTIONAL { ?bgc wp:organism    ?bgc_taxon   . }

    # Restrict to tair.locus genes (Arabidopsis direct matches)
    FILTER(STRSTARTS(STR(?gene), "https://identifiers.org/tair.name/"))
}
ORDER BY ?bgc ?gene
"""

QUERY_DIRECT_GENE_ISPARTOF = """
SELECT DISTINCT
    ?gene ?pathway ?pathway_title ?bgc ?bgc_source ?bgc_species ?bgc_taxon
WHERE {
    # Gene in pathway (gene IRI directly has dcterms:isPartOf)
    ?gene a wp:GeneProduct ;
          dcterms:isPartOf ?pathway .

    # Only keep pathway-level IRIs, not interaction-level sub-IRIs
    FILTER(STRSTARTS(STR(?pathway), "http://rdf-plantmetwiki.bioinformatics.nl/pathways/"))

    # Pathway title (try multiple predicates used in WikiPathways RDF)
    OPTIONAL { ?pathway dc:title     ?pathway_title . }
    OPTIONAL { ?pathway rdfs:label   ?pathway_title . }
    OPTIONAL { ?pathway dcterms:title ?pathway_title . }

    # Gene in BGC
    ?bgc obo:RO_0000051 ?gene .
    ?bgc dcterms:source ?bgc_source .
    OPTIONAL { ?bgc wp:organismName ?bgc_species . }
    OPTIONAL { ?bgc wp:organism    ?bgc_taxon   . }

    FILTER(STRSTARTS(STR(?gene), "https://identifiers.org/tair.name/"))
}
ORDER BY ?bgc ?gene
"""

# ============================================================
# Helpers
# ============================================================

def short(iri: str) -> str:
    """Return last segment of an IRI (after / or #)."""
    return iri.rstrip("/").rsplit("/", 1)[-1].rsplit("#", 1)[-1]


def tair_id(gene_iri: str) -> str:
    """AT1G23500 from https://identifiers.org/tair.name/AT1G23500."""
    return gene_iri.rsplit("/", 1)[-1]


def ncbi_id(taxon_iri: str) -> str:
    """3702 from http://purl.obolibrary.org/obo/NCBITaxon_3702."""
    return taxon_iri.rsplit("_", 1)[-1] if "_" in taxon_iri else taxon_iri


def run_query(g: Graph, sparql: str, label: str) -> list[dict]:
    rows = []
    for r in g.query(sparql, initNs=INIT_NS):
        src = str(r.bgc_source) if r.bgc_source else ""
        rows.append({
            "gene_iri":       str(r.gene),
            "gene_id":        tair_id(str(r.gene)),
            "pathway_iri":    str(r.pathway),
            "pathway_id":     short(str(r.pathway)),
            "pathway_title":  str(r.pathway_title) if r.pathway_title else "",
            "bgc_iri":        str(r.bgc),
            "bgc_id":         short(str(r.bgc)),
            "bgc_source":     src,
            "bgc_url":        bgc_url(str(r.bgc), src),
            "bgc_species":    str(r.bgc_species) if r.bgc_species else "",
            "bgc_ncbi_taxon": ncbi_id(str(r.bgc_taxon)) if r.bgc_taxon else "",
            "link_type":      label,
        })
    return rows


# ============================================================
# Main
# ============================================================

FIELDNAMES = [
    "bgc_id", "bgc_source", "bgc_url", "bgc_species", "bgc_ncbi_taxon",
    "gene_id", "gene_iri",
    "pathway_id", "pathway_title", "pathway_iri",
    "link_type",
]

MIBIG_URL_BASE       = "https://mibig.secondarymetabolites.org/repository/"
PLANTISMASH_URL_BASE = "https://plantismash.bioinformatics.nl/precalc/v2/"


def bgc_url(bgc_iri: str, bgc_source: str) -> str:
    """Construct a human-clickable URL for a BGC."""
    if bgc_source == "MIBIG":
        # bioregistry IRI: https://bioregistry.io/mibig:BGC0000670
        # extract BGC ID from end of IRI
        bgc_id = bgc_iri.rsplit(":", 1)[-1].rsplit("/", 1)[-1]
        return f"{MIBIG_URL_BASE}{bgc_id}"
    # plantiSMASH: the IRI IS the URL
    return bgc_iri


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pathways", type=Path,
                   default=Path("input_ttl/all-plantcyc17.0.0-gpml2021.ttl"),
                   help="Pathway RDF TTL (default: input_ttl/all-plantcyc17.0.0-gpml2021.ttl)")
    p.add_argument("--bgc-plantismash", type=Path,
                   default=Path("output_ttl/plantismash.ttl"))
    p.add_argument("--bgc-mibig", type=Path,
                   default=Path("output_ttl/mibig.ttl"))
    p.add_argument("--output", type=Path,
                   default=Path("summaries/bgc_pathway_links.tsv"))
    return p.parse_args()


def main() -> int:
    args = parse_args()

    g = Graph()

    for ttl in (args.bgc_plantismash, args.bgc_mibig, args.pathways):
        if not ttl.exists():
            print(f"[SKIP] {ttl} not found")
            continue
        print(f"Loading {ttl} ...", end=" ", flush=True)
        g.parse(str(ttl), format="turtle")
        print(f"{len(g):,} triples")

    print(f"\nTotal triples: {len(g):,}")
    print("Running queries ...\n")

    all_rows: list[dict] = []

    # Try both predicate patterns
    for label, sparql in [
        ("direct_via_datanode", QUERY_DIRECT_VIA_DATANODE),
        ("direct_gene_ispartof", QUERY_DIRECT_GENE_ISPARTOF),
    ]:
        rows = run_query(g, sparql, label)
        print(f"  [{label}]: {len(rows)} links")
        all_rows.extend(rows)

    # Deduplicate (same gene/pathway/bgc might appear from both patterns)
    seen = set()
    deduped = []
    for r in all_rows:
        key = (r["gene_iri"], r["pathway_iri"], r["bgc_iri"])
        if key not in seen:
            seen.add(key)
            deduped.append(r)

    print(f"\n  Total unique links: {len(deduped)}")

    if not deduped:
        print("\n  No links found. Check that the pathway and BGC IRIs share the same gene namespace.")
        return 0

    # Summary statistics
    from collections import Counter
    by_source  = Counter(r["bgc_source"]  for r in deduped)
    by_species = Counter(r["bgc_species"] for r in deduped)
    unique_bgcs     = len({r["bgc_id"]      for r in deduped})
    unique_genes    = len({r["gene_id"]     for r in deduped})
    unique_pathways = len({r["pathway_id"]  for r in deduped})

    print(f"\n  Unique BGCs linked:     {unique_bgcs}")
    print(f"  Unique genes linked:    {unique_genes}")
    print(f"  Unique pathways linked: {unique_pathways}")
    print(f"  By source:  {dict(by_source)}")
    print(f"  By species: {dict(by_species)}")

    # Write TSV
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sorted(deduped, key=lambda r: (r["bgc_source"], r["bgc_id"], r["gene_id"])))

    print(f"\n✔ Wrote {len(deduped)} rows → {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
