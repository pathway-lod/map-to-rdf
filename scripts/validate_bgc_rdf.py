#!/usr/bin/env python3
"""Validate the generated BGC RDF files and supporting artefacts.

Checks performed
----------------
1. TTL syntax     — rdflib parses plantismash.ttl and mibig.ttl without errors
2. Triple count   — both files have at least MIN_TRIPLES triples
3. Required types — every cluster is typed as pmw:BiosyntheticGeneCluster
4. Provenance     — every cluster has dcterms:source
5. Gene links     — RO:0000051 (has_part) triples are present
6. Species        — at least one cluster carries wp:organism
7. VoID           — void-bgc.ttl parses and declares both void:Dataset URIs
8. Link table     — bgc_pathway_links.tsv exists and has the required columns

Exit codes
----------
0  all checks passed
1  one or more checks failed
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF, DCTERMS

# ── Namespaces ────────────────────────────────────────────────────────────
PMW   = Namespace("http://rdf-plantmetwiki.bioinformatics.nl/vocab/")
WP    = Namespace("http://vocabularies.wikipathways.org/wp#")
VOID  = Namespace("http://rdfs.org/ns/void#")
RO_HAS_PART = URIRef("http://purl.obolibrary.org/obo/RO_0000051")

# ── Thresholds ────────────────────────────────────────────────────────────
MIN_TRIPLES = {
    "plantismash.ttl": 1_000,
    "mibig.ttl":         100,
}

REQUIRED_TSV_COLUMNS = {
    "bgc_id", "bgc_source", "bgc_url", "bgc_species",
    "gene_id", "pathway_id", "pathway_title", "link_type",
}

REQUIRED_VOID_DATASETS = {
    "http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/plantismash-v2",
    "http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/mibig-4.0",
}

# ── Helpers ───────────────────────────────────────────────────────────────

PASS = "\033[32m✔\033[0m"
FAIL = "\033[31m✗\033[0m"


class Results:
    def __init__(self):
        self._failures: list[str] = []

    def ok(self, label: str) -> None:
        print(f"  {PASS} {label}")

    def fail(self, label: str, detail: str = "") -> None:
        msg = f"{label}" + (f": {detail}" if detail else "")
        print(f"  {FAIL} {msg}")
        self._failures.append(msg)

    def check(self, condition: bool, ok_label: str,
               fail_label: str, detail: str = "") -> None:
        if condition:
            self.ok(ok_label)
        else:
            self.fail(fail_label, detail)

    @property
    def failed(self) -> bool:
        return bool(self._failures)

    def summary(self) -> None:
        print()
        if self._failures:
            print(f"FAILED — {len(self._failures)} check(s) failed:")
            for f in self._failures:
                print(f"  • {f}")
        else:
            print("All checks passed.")


# ── Individual checks ─────────────────────────────────────────────────────

def check_ttl(path: Path, r: Results) -> Graph | None:
    """Parse TTL and return graph, or record failure."""
    label = path.name
    if not path.exists():
        r.fail(f"{label} exists", "file not found")
        return None
    try:
        g = Graph()
        g.parse(str(path), format="turtle")
        r.ok(f"{label} parses without errors")
        return g
    except Exception as exc:
        r.fail(f"{label} parses without errors", str(exc))
        return None


def check_triples(g: Graph, path: Path, r: Results) -> None:
    n = len(g)
    minimum = MIN_TRIPLES.get(path.name, 10)
    r.check(
        n >= minimum,
        f"{path.name}: {n:,} triples (≥ {minimum:,})",
        f"{path.name} triple count too low",
        f"{n} < {minimum}",
    )


def check_cluster_types(g: Graph, label: str, r: Results) -> None:
    clusters = set(g.subjects(RDF.type, PMW.BiosyntheticGeneCluster))
    r.check(
        len(clusters) > 0,
        f"{label}: {len(clusters)} pmw:BiosyntheticGeneCluster resources",
        f"{label}: no pmw:BiosyntheticGeneCluster found",
    )


def check_dcterms_source(g: Graph, label: str, r: Results) -> None:
    clusters  = set(g.subjects(RDF.type, PMW.BiosyntheticGeneCluster))
    with_src  = {s for s in clusters if list(g.objects(s, DCTERMS.source))}
    missing   = clusters - with_src
    r.check(
        not missing,
        f"{label}: all clusters have dcterms:source",
        f"{label}: {len(missing)} cluster(s) missing dcterms:source",
    )


def check_gene_links(g: Graph, label: str, r: Results) -> None:
    links = list(g.triples((None, RO_HAS_PART, None)))
    r.check(
        len(links) > 0,
        f"{label}: {len(links):,} RO:0000051 (has_part) gene links",
        f"{label}: no RO:0000051 gene links found",
    )


def check_species(g: Graph, label: str, r: Results) -> None:
    with_sp = list(g.triples((None, WP.organism, None)))
    r.check(
        len(with_sp) > 0,
        f"{label}: {len(with_sp)} cluster(s) with wp:organism",
        f"{label}: no wp:organism annotations found",
    )


def check_void(path: Path, r: Results) -> None:
    if not path.exists():
        r.fail(f"{path.name} exists", "file not found")
        return
    try:
        g = Graph()
        g.parse(str(path), format="turtle")
    except Exception as exc:
        r.fail(f"{path.name} parses", str(exc))
        return
    r.ok(f"{path.name} parses without errors")
    datasets = {str(s) for s in g.subjects(RDF.type, VOID.Dataset)}
    for uri in REQUIRED_VOID_DATASETS:
        r.check(
            uri in datasets,
            f"void:Dataset <{uri.split('/')[-1]}>",
            f"missing void:Dataset",
            uri,
        )


def check_link_table(path: Path, r: Results) -> None:
    if not path.exists():
        r.fail(f"{path.name} exists", "file not found")
        return
    r.ok(f"{path.name} exists")
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = set(reader.fieldnames or [])
        missing = REQUIRED_TSV_COLUMNS - cols
        r.check(
            not missing,
            f"{path.name}: all required columns present",
            f"{path.name}: missing columns",
            ", ".join(sorted(missing)),
        )
        rows = list(reader)
    r.check(
        len(rows) > 0,
        f"{path.name}: {len(rows)} link rows",
        f"{path.name}: empty (no rows)",
    )


# ── Main ──────────────────────────────────────────────────────────────────

def main() -> int:
    # Resolve paths relative to repo root (script may be called from anywhere)
    root = Path(__file__).parent.parent

    r = Results()

    for fname in ("plantismash.ttl", "mibig.ttl"):
        path = root / "output_ttl" / fname
        print(f"\n── {fname} ──")
        g = check_ttl(path, r)
        if g is not None:
            check_triples(g, path, r)
            check_cluster_types(g, fname, r)
            check_dcterms_source(g, fname, r)
            check_gene_links(g, fname, r)
            check_species(g, fname, r)

    print("\n── void-bgc.ttl ──")
    check_void(root / "output_ttl" / "void-bgc.ttl", r)

    print("\n── bgc_pathway_links.tsv ──")
    check_link_table(root / "summaries" / "bgc_pathway_links.tsv", r)

    r.summary()
    return 1 if r.failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
