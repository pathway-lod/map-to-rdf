#!/usr/bin/env python3
"""Local SPARQL query tester and BridgeDb-assisted BGC↔pathway linker.

Loads the generated BGC RDF (output_ttl/) and optionally the pathway RDF
(input_ttl/all_pathways.ttl) into an in-memory rdflib graph, then runs:

  A) Basic sanity queries against the BGC graph
  B) BridgeDb-assisted pathway→BGC matching
       For each TAIR locus in the pathway graph, fetch BridgeDb xrefs
       (NCBI Gene IDs), then check if those xrefs appear as BGC member IRIs.

Usage
-----
    # BGC queries only (no pathway graph required)
    python scripts/test_queries.py

    # Include pathway graph + BridgeDb-assisted linking
    python scripts/test_queries.py --pathways input_ttl/all_pathways.ttl

    # Limit BridgeDb calls (useful for quick checks)
    python scripts/test_queries.py --pathways input_ttl/all_pathways.ttl --max-genes 20
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from urllib.parse import quote

import requests
from rdflib import Graph, Namespace

# ============================================================
# Namespaces — must match convert_bgc_to_rdf.py
# ============================================================

PMW       = Namespace("http://rdf-plantmetwiki.bioinformatics.nl/vocab/")
WP        = Namespace("http://vocabularies.wikipathways.org/wp#")
DCT       = Namespace("http://purl.org/dc/terms/")
FOAF      = Namespace("http://xmlns.com/foaf/0.1/")
NCBITAXON = Namespace("http://purl.obolibrary.org/obo/NCBITaxon_")

TAIR_LOCUS = Namespace("https://identifiers.org/tair.locus/")
NCBIGENE   = Namespace("https://identifiers.org/ncbigene:")

BRIDGEDB_BASE = "https://webservice.bridgedb.org"

INIT_NS = {
    "pmw":     PMW,
    "wp":      WP,
    "dcterms": DCT,
    "obo":     Namespace("http://purl.obolibrary.org/obo/"),
    "foaf":    FOAF,
    "ncbi":    NCBITAXON,
}

# ============================================================
# Helpers
# ============================================================

def sep(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def run(g: Graph, sparql: str, title: str, limit: int = 50) -> list:
    sep(title)
    results = list(g.query(sparql, initNs=INIT_NS))
    for row in results[:limit]:
        print(" ", row)
    if len(results) > limit:
        print(f"  ... ({len(results)} total, showing {limit})")
    return results


# ============================================================
# BridgeDb direction: TAIR locus → NCBI Gene IDs → BGC members
# ============================================================

def bridgedb_ncbi_xrefs(tair_id: str, timeout: int = 8) -> list[str]:
    """Return NCBI Gene IDs for *tair_id* via BridgeDb.

    Calls: GET /Arabidopsis thaliana/xrefs/A/{tair_id}
    Filters for system code 'Ec' (Entrez Gene / NCBI Gene).
    Returns a list of NCBI gene ID strings (e.g. ['831234']).
    """
    organism = "Arabidopsis%20thaliana"
    url = f"{BRIDGEDB_BASE}/{organism}/xrefs/A/{quote(tair_id, safe='')}"
    try:
        r = requests.get(url, timeout=timeout)
        if not r.ok:
            return []
        ids = []
        for line in r.text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1] == "Ec":
                ids.append(parts[0].strip())
        return ids
    except Exception:
        return []


def find_bgc_links_via_bridgedb(
    g: Graph,
    max_genes: int | None = None,
) -> list[dict]:
    """For every TAIR locus gene in the pathway graph, use BridgeDb to find
    matching BGC cluster members (via NCBI Gene IDs).

    Returns a list of match dicts: {pathway, gene_tair, ncbi_id, cluster, source}.
    """
    # 1. Get all TAIR locus gene IRIs from pathway graph
    pathway_genes_q = """
    PREFIX tair: <https://identifiers.org/tair.locus/>
    SELECT DISTINCT ?gene ?pathway
    WHERE {
        ?pathway a wp:Pathway .
        ?pathway wp:containsElement ?gene .
        FILTER(STRSTARTS(STR(?gene), "https://identifiers.org/tair.locus/"))
    }
    """
    # Fallback: also catch genes linked via dcterms/wp predicates
    pathway_genes_q2 = """
    SELECT DISTINCT ?gene
    WHERE {
        ?gene a wp:GeneProduct .
        FILTER(STRSTARTS(STR(?gene), "https://identifiers.org/tair.locus/"))
    }
    """

    tair_genes: set[str] = set()
    for row in g.query(pathway_genes_q, initNs=INIT_NS):
        tair_genes.add(str(row.gene))
    for row in g.query(pathway_genes_q2, initNs=INIT_NS):
        tair_genes.add(str(row.gene))

    if not tair_genes:
        print("  No TAIR locus gene IRIs found in the graph.")
        return []

    print(f"  Found {len(tair_genes)} unique TAIR locus genes in the pathway graph.")
    if max_genes:
        tair_genes_list = sorted(tair_genes)[:max_genes]
        print(f"  Limiting BridgeDb calls to {max_genes} genes.")
    else:
        tair_genes_list = sorted(tair_genes)

    # 2. Collect all BGC member IRIs that are NCBI Gene URIs
    bgc_members_q = """
    SELECT ?cluster ?member ?source
    WHERE {
        ?cluster obo:RO_0000051 ?member .
        ?cluster dcterms:source ?source .
        FILTER(STRSTARTS(STR(?member), "https://identifiers.org/ncbigene:"))
    }
    """
    ncbigene_in_bgc: dict[str, list[tuple[str, str]]] = {}
    for row in g.query(bgc_members_q, initNs=INIT_NS):
        ncbi_id = str(row.member).replace("https://identifiers.org/ncbigene:", "")
        ncbigene_in_bgc.setdefault(ncbi_id, []).append(
            (str(row.cluster), str(row.source))
        )

    if not ncbigene_in_bgc:
        print("  No NCBI Gene IRI members found in BGC graph — "
              "BridgeDb matching requires plantiSMASH data with LOC IDs.")
        return []

    print(f"  BGC graph contains {len(ncbigene_in_bgc)} unique NCBI Gene member IRIs.")

    # 3. For each TAIR gene, get NCBI Gene xrefs and check for BGC membership
    matches = []
    for i, gene_iri in enumerate(tair_genes_list, 1):
        tair_id = gene_iri.replace("https://identifiers.org/tair.locus/", "")
        ncbi_ids = bridgedb_ncbi_xrefs(tair_id)
        time.sleep(0.05)  # be polite

        for ncbi_id in ncbi_ids:
            if ncbi_id in ncbigene_in_bgc:
                for cluster_uri, source in ncbigene_in_bgc[ncbi_id]:
                    matches.append({
                        "gene_tair":  tair_id,
                        "ncbi_id":    ncbi_id,
                        "cluster":    cluster_uri,
                        "source":     source,
                    })

        if i % 10 == 0:
            print(f"  [{i}/{len(tair_genes_list)}] checked, {len(matches)} matches so far ...")

    return matches


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--pathways",
        type=Path,
        default=None,
        help="Path to pathway TTL (e.g. input_ttl/all_pathways.ttl). "
             "Required for BridgeDb-assisted linking.",
    )
    p.add_argument(
        "--max-genes",
        type=int,
        default=None,
        help="Limit BridgeDb calls to this many TAIR genes (useful for quick tests).",
    )
    p.add_argument(
        "--test",
        action="store_true",
        help="Load test BGC TTLs from output_ttl/ (generated by --test run of convert script).",
    )
    return p.parse_args()


# ============================================================
# Main
# ============================================================

def main() -> int:
    args = parse_args()

    g = Graph()

    # Load BGC graphs
    for ttl in (Path("output_ttl/plantismash.ttl"), Path("output_ttl/mibig.ttl")):
        if ttl.exists():
            g.parse(str(ttl), format="turtle")
            print(f"Loaded {ttl} ({len(g)} triples so far)")
        else:
            print(f"[SKIP] {ttl} not found — run convert_bgc_to_rdf.py first")

    # Load pathway graph if provided
    has_pathways = False
    if args.pathways:
        if args.pathways.exists():
            g.parse(str(args.pathways), format="turtle")
            print(f"Loaded {args.pathways} ({len(g)} triples total)")
            has_pathways = True
        else:
            print(f"[WARN] Pathway file not found: {args.pathways}")

    print(f"\nTotal triples in graph: {len(g)}\n")

    # ---- A) Basic BGC sanity queries ----

    run(g, """
        SELECT ?cluster ?gene
        WHERE { ?cluster obo:RO_0000051 ?gene . }
        ORDER BY ?cluster ?gene
        LIMIT 20
    """, "A1) Sample gene→BGC links")

    run(g, """
        SELECT (COUNT(*) AS ?nLinks)
        WHERE { ?cluster obo:RO_0000051 ?gene . }
    """, "A2) Total gene→BGC links")

    run(g, """
        SELECT ?source (COUNT(*) AS ?nLinks)
        WHERE {
            ?cluster obo:RO_0000051 ?gene .
            ?cluster dcterms:source ?source .
        }
        GROUP BY ?source
        ORDER BY DESC(?nLinks)
    """, "A3) Links by source (MIBiG vs plantiSMASH)")

    run(g, """
        SELECT ?species (COUNT(DISTINCT ?cluster) AS ?nClusters)
        WHERE {
            ?cluster a pmw:BiosyntheticGeneCluster ;
                     wp:organismName ?species .
        }
        GROUP BY ?species
        ORDER BY DESC(?nClusters)
    """, "A4) BGC clusters by species")

    run(g, """
        SELECT (COUNT(DISTINCT ?cluster) AS ?nUntyped)
        WHERE {
            ?cluster a pmw:BiosyntheticGeneCluster .
            FILTER NOT EXISTS { ?cluster wp:organism ?taxon . }
        }
    """, "A5) Untyped BGC clusters (no species assigned)")

    # ---- B1) Direct Arabidopsis join (no BridgeDb needed) ----

    if has_pathways:
        run(g, """
            SELECT (COUNT(DISTINCT ?pathway) AS ?nPathways)
                   (COUNT(DISTINCT ?gene)    AS ?nGenes)
                   (COUNT(DISTINCT ?cluster) AS ?nClusters)
            WHERE {
                ?gene dcterms:isPartOf ?pathway .
                FILTER(STRSTARTS(STR(?gene), "https://identifiers.org/tair.locus/"))
                ?cluster obo:RO_0000051 ?gene .
                ?cluster dcterms:source ?source .
            }
        """, "B1a) Count: distinct pathways / genes / BGC clusters linked")

        run(g, """
            SELECT DISTINCT ?pathway ?gene ?cluster ?source
            WHERE {
                ?gene dcterms:isPartOf ?pathway .
                FILTER(STRSTARTS(STR(?gene), "https://identifiers.org/tair.locus/"))
                ?cluster obo:RO_0000051 ?gene .
                ?cluster dcterms:source ?source .
            }
            ORDER BY ?cluster ?gene
        """, "B1b) Direct Arabidopsis pathway↔BGC links (gene in pathway AND in BGC)")

    # ---- B2) BridgeDb-assisted pathway→BGC linking (tomato) ----

    if has_pathways:
        sep("B2) BridgeDb-assisted pathway→BGC linking")
        print("  Strategy: TAIR locus → BridgeDb xrefs (NCBI Gene IDs) → BGC members")

        matches = find_bgc_links_via_bridgedb(g, max_genes=args.max_genes)

        if matches:
            print(f"\n  ✔ Found {len(matches)} pathway gene → BGC link(s):\n")
            for m in matches:
                print(f"  {m['gene_tair']} (NCBI:{m['ncbi_id']}) "
                      f"→ {m['cluster'].split('/')[-1]} [{m['source']}]")
        else:
            print("\n  No cross-source links found via BridgeDb.")
            print("  This is expected if the pathway and BGC graphs cover different genes.")
    else:
        sep("B) BridgeDb-assisted pathway→BGC linking")
        print("  Skipped — provide --pathways <file.ttl> to enable.")
        print("  Example: python scripts/test_queries.py "
              "--pathways input_ttl/all_pathways.ttl --max-genes 50")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
