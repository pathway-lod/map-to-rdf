#!/usr/bin/env python3
"""Convert biosynthetic gene cluster (BGC) data from MIBiG and plantiSMASH to RDF.

Usage
-----
    # Full run (reads from input/)
    python scripts/convert_bgc_to_rdf.py

    # Quick test with root-level test JSON files
    python scripts/convert_bgc_to_rdf.py --test

    # Only one source
    python scripts/convert_bgc_to_rdf.py --source plantismash

    # Skip BridgeDb lookups (offline / CI)
    python scripts/convert_bgc_to_rdf.py --no-bridgedb
"""

from __future__ import annotations

import argparse
import json
import re
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

import requests
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, DCTERMS

# ============================================================
# Repo layout
# ============================================================

INPUT_DIR  = Path("input")
OUTPUT_TTL_DIR = Path("output_ttl")
SUMMARY_DIR    = Path("summaries")

# Full-run inputs
INPUT_MIBIG_GENES_JSON  = INPUT_DIR / "mibig_genes.json"
INPUT_PLANTISMASH_JSON  = INPUT_DIR / "plantismash_v2_clusters_minimal.json"

# Test-run inputs (root of repo, small samples)
TEST_PLANTISMASH_JSON = Path("testPlantIsMash.json")
TEST_MIBIG_JSON       = Path("testMibig.json")

OUT_PLANTISMASH_TTL = OUTPUT_TTL_DIR / "plantismash.ttl"
OUT_MIBIG_TTL       = OUTPUT_TTL_DIR / "mibig.ttl"
OUT_SUMMARY         = SUMMARY_DIR    / "bgc_conversion_summary.json"

# ============================================================
# Namespaces / vocab
# ============================================================

# Aligned with gpml-to-rdf: http://rdf-plantmetwiki.bioinformatics.nl/vocab/
PMW       = Namespace("http://rdf-plantmetwiki.bioinformatics.nl/vocab/")
WP        = Namespace("http://vocabularies.wikipathways.org/wp#")
FOAF      = Namespace("http://xmlns.com/foaf/0.1/")
NCBITAXON = Namespace("http://purl.obolibrary.org/obo/NCBITaxon_")

# Canonical RO has_part term
RO_HAS_PART = URIRef("http://purl.obolibrary.org/obo/RO_0000051")

# MIBiG stable identifier via Bioregistry
MIBIG_BIOREGISTRY_PREFIX = "https://bioregistry.io/mibig:"

# plantiSMASH cluster namespace
PLANTISMASH = Namespace("https://plantismash.bioinformatics.nl/precalc/v2/")

# Internal URI for unresolvable member identifiers
PMW_MEMBER_BASE = "http://rdf-plantmetwiki.bioinformatics.nl/id/member/"

# ============================================================
# Species detection and gene IRI patterns
# ============================================================

RE_AT_CLUSTER = re.compile(r"Arabidopsis_thaliana")
RE_CROSS      = re.compile(r"_x_")

# Standard TAIR locus ID: AT{1-5,M,C}G{5 digits}, optional .N isoform
RE_AT_GENE = re.compile(r"^AT[1-5MC]G\d{5}(\.\d+)?$", re.IGNORECASE)

# plantiSMASH isoform notation: AT1G2351_1 → AT1G23510
# (drops trailing 0 from locus number and appends underscore+isoform)
RE_AT_PLANTISMASH_ISOFORM = re.compile(
    r"^(?P<prefix>AT[1-5MC]G)(?P<digits>\d{4})_\d+$", re.IGNORECASE
)

# MIBiG tomato pattern (Solyc gene model)
RE_TOMATO_SOLYC = re.compile(r"^Solyc\d{2}g\d{6}(\.\d+)+$", re.IGNORECASE)

# Identifier namespaces
# tair.name is used by the WikiPathways/PlantMetWiki pathway RDF (from GPML Xref
# dataSource="TAIR gene name"), so we must match it here for direct SPARQL joins.
# tair.locus is kept as foaf:page for discoverability.
TAIR_LOCUS    = Namespace("https://identifiers.org/tair.locus/")
TAIR_NAME     = Namespace("https://identifiers.org/tair.name/")
ENSEMBL_PLANT = Namespace("http://identifiers.org/ensembl.plant:")

TAIR_SEARCH          = "https://www.arabidopsis.org/results?mainType=general&category=genes&searchText="
ENSEMBL_TOMATO_SEARCH = "https://plants.ensembl.org/Solanum_lycopersicum/Search/Results?species=Solanum_lycopersicum;idx=;q="

# ============================================================
# BridgeDb integration
# ============================================================

BRIDGEDB_BASE = "https://webservice.bridgedb.org"
_BRIDGEDB_CACHE: dict[str, str | None] = {}

# System codes to try when resolving an unknown identifier to TAIR locus.
# "A" = TAIR locus (target); others are potential source codes.
_SOURCE_CODES = ("L", "X", "En", "Uc", "Ag")  # locus tag, RefSeq, Ensembl, UniGene, Affy


def bridgedb_to_tair(gene_id: str, organism: str = "Arabidopsis thaliana",
                     timeout: int = 8) -> str | None:
    """Try to resolve *gene_id* to a TAIR locus ID via the BridgeDb REST API.

    Returns the first TAIR locus hit (e.g. 'AT1G23500'), or None.
    Results are cached in-process to avoid redundant HTTP calls.
    """
    cache_key = f"{organism}|{gene_id}"
    if cache_key in _BRIDGEDB_CACHE:
        return _BRIDGEDB_CACHE[cache_key]

    organism_enc = quote(organism)
    gene_enc     = quote(gene_id, safe="")

    for src_code in _SOURCE_CODES:
        url = f"{BRIDGEDB_BASE}/{organism_enc}/xrefs/{src_code}/{gene_enc}"
        try:
            r = requests.get(url, timeout=timeout)
            if not r.ok:
                continue
            for line in r.text.strip().splitlines():
                parts = line.split("\t")
                if len(parts) >= 2 and parts[1] == "A":
                    tair_id = parts[0].strip()
                    _BRIDGEDB_CACHE[cache_key] = tair_id
                    return tair_id
        except Exception:
            continue
        time.sleep(0.05)  # be polite to the BridgeDb API

    _BRIDGEDB_CACHE[cache_key] = None
    return None


# ============================================================
# ID normalisation helpers
# ============================================================

def normalise_plantismash_at_id(gene_id: str) -> str | None:
    """Fix plantiSMASH's non-standard Arabidopsis isoform notation.

    plantiSMASH writes AT1G2351_1 meaning isoform 1 of AT1G23510
    (4-digit locus number + underscore suffix instead of the canonical
    5-digit locus + dot suffix). Returns the canonical locus string
    (e.g. 'AT1G23510') or None if the pattern does not match.
    """
    m = RE_AT_PLANTISMASH_ISOFORM.fullmatch(gene_id)
    if m:
        return f"{m.group('prefix')}{m.group('digits')}0"
    return None


def species_from_plantismash_cluster_id(cluster_id: str) -> tuple[str, str] | None:
    """Return (species_label, ncbitaxon_id) for supported species, else None.

    Currently restricted to Arabidopsis thaliana: plantiSMASH gene models for
    Solanum lycopersicum use a different genome annotation (NCBI RefSeq) than
    the PlantCyc pathway graph (Phytozome/ITAG), so they cannot be joined to
    pathway genes without a BridgeDb-assisted crosswalk (see README). That
    work is preserved on the `plantismash-solanum-lycopersicum` branch.
    """
    if RE_CROSS.search(cluster_id):
        return None
    if RE_AT_CLUSTER.search(cluster_id):
        return ("Arabidopsis thaliana", "3702")
    return None


def mint_arabidopsis_gene_iri(gene_id: str) -> URIRef:
    """Mint a tair.name IRI to match the WikiPathways/PlantMetWiki pathway RDF.

    The pathway RDF uses identifiers.org/tair.name/ (from GPML Xref dataSource
    'TAIR gene name'). Using the same namespace enables direct SPARQL joins.
    tair.locus is kept as foaf:page for cross-reference discoverability.
    """
    core = gene_id.split(".")[0].upper()
    return URIRef(TAIR_NAME[core])


def add_gene_pages_arabidopsis(g: Graph, gene: URIRef, gene_id: str) -> None:
    core = gene_id.split(".")[0].upper()
    g.add((gene, FOAF.page, URIRef(TAIR_SEARCH + core)))
    g.add((gene, FOAF.page, URIRef(TAIR_LOCUS[core])))


def mint_member_identifier_iri(member_id: str) -> URIRef:
    return URIRef(PMW_MEMBER_BASE + quote(member_id, safe=""))


# ============================================================
# Shared RDF helpers
# ============================================================

def add_common_cluster_metadata(
    g: Graph,
    cluster: URIRef,
    source: str,
    species_label: str | None = None,
    taxon_id: str | None = None,
) -> None:
    g.add((cluster, RDF.type, PMW.BiosyntheticGeneCluster))
    g.add((cluster, DCTERMS.source, Literal(source)))
    if species_label and taxon_id:
        g.add((cluster, WP.organismName, Literal(species_label)))
        g.add((cluster, WP.organism, URIRef(NCBITAXON[taxon_id])))
    label = str(cluster).split("/")[-1]
    g.add((cluster, RDFS.label, Literal(label)))


def add_common_gene_metadata(
    g: Graph, gene: URIRef, species_label: str, taxon_id: str
) -> None:
    g.add((gene, RDF.type, WP.GeneProduct))
    g.add((gene, WP.organismName, Literal(species_label)))
    g.add((gene, WP.organism, URIRef(NCBITAXON[taxon_id])))


# ============================================================
# plantiSMASH conversion
# ============================================================

def convert_plantismash(
    plantismash_json: dict, use_bridgedb: bool = True
) -> tuple[Graph, dict]:
    g = Graph()
    g.bind("pmw",     PMW)
    g.bind("wp",      WP)
    g.bind("dcterms", DCTERMS)
    g.bind("foaf",    FOAF)
    g.bind("ncbi",    NCBITAXON)
    g.bind("obo",     Namespace("http://purl.obolibrary.org/obo/"))
    g.bind("planti",  PLANTISMASH)

    summary: dict = {
        "source": "plantiSMASH",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "species": defaultdict(lambda: {"bgcs": 0, "gene_links": 0, "unique_genes": set()}),
        "bridgedb_resolved": 0,
    }

    for cluster_id, gene_list in plantismash_json.items():
        sp = species_from_plantismash_cluster_id(cluster_id)
        if sp is None:
            continue
        species_label, taxon_id = sp

        cluster = URIRef(PLANTISMASH[cluster_id])
        add_common_cluster_metadata(g, cluster, "plantiSMASH", species_label, taxon_id)
        summary["species"][species_label]["bgcs"] += 1

        for raw_gene_id in gene_list:
            gene_id = raw_gene_id.strip()
            if not gene_id:
                continue

            # ---- Arabidopsis ----
            resolved_id = gene_id

            # 1. Already a valid TAIR locus ID
            if not RE_AT_GENE.fullmatch(gene_id):
                # 2. plantiSMASH isoform notation (AT1G2351_1 → AT1G23510)
                fixed = normalise_plantismash_at_id(gene_id)
                if fixed:
                    resolved_id = fixed
                elif use_bridgedb:
                    # 3. Try BridgeDb as last resort
                    tair = bridgedb_to_tair(gene_id)
                    if tair:
                        resolved_id = tair
                        summary["bridgedb_resolved"] += 1

            if RE_AT_GENE.fullmatch(resolved_id):
                gene = mint_arabidopsis_gene_iri(resolved_id)
                add_common_gene_metadata(g, gene, species_label, taxon_id)
                add_gene_pages_arabidopsis(g, gene, resolved_id)
                g.add((cluster, RO_HAS_PART, gene))
                summary["species"][species_label]["gene_links"] += 1
                summary["species"][species_label]["unique_genes"].add(str(gene))
            else:
                member = mint_member_identifier_iri(gene_id)
                g.add((member, RDF.type, PMW.MemberIdentifier))
                g.add((member, RDFS.label, Literal(gene_id)))
                g.add((cluster, RO_HAS_PART, member))

    for sp_label in list(summary["species"].keys()):
        summary["species"][sp_label]["unique_genes"] = len(
            summary["species"][sp_label]["unique_genes"]
        )

    return g, summary


# ============================================================
# MIBiG conversion
# ============================================================

def convert_mibig(
    mibig_json: dict, use_bridgedb: bool = True
) -> tuple[Graph, dict]:
    """Convert MIBiG gene cluster membership to RDF.

    Species is only assigned when identifiers match reliable patterns.
    XP_/NP_/KNA_ etc. are kept as pmw:MemberIdentifier without species.
    """
    g = Graph()
    g.bind("pmw",     PMW)
    g.bind("wp",      WP)
    g.bind("dcterms", DCTERMS)
    g.bind("foaf",    FOAF)
    g.bind("ncbi",    NCBITAXON)
    g.bind("obo",     Namespace("http://purl.obolibrary.org/obo/"))
    g.bind("mibig",   Namespace(MIBIG_BIOREGISTRY_PREFIX))

    summary: dict = {
        "source": "MIBiG",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "species": defaultdict(lambda: {"bgcs": set(), "gene_links": 0, "unique_genes": set()}),
        "untyped": {"bgcs": set(), "member_links": 0, "unique_members": set()},
        "bridgedb_resolved": 0,
    }

    for bgc_id, member_list in mibig_json.items():
        cluster = URIRef(f"{MIBIG_BIOREGISTRY_PREFIX}{bgc_id}")
        add_common_cluster_metadata(g, cluster, "MIBIG", None, None)
        typed_this_cluster = False

        for raw_member_id in member_list:
            member_id = str(raw_member_id).strip()
            if not member_id:
                continue

            # ---- Arabidopsis (safe: TAIR locus regex) ----
            if RE_AT_GENE.fullmatch(member_id):
                species_label, taxon_id = "Arabidopsis thaliana", "3702"
                gene = mint_arabidopsis_gene_iri(member_id)
                add_common_gene_metadata(g, gene, species_label, taxon_id)
                add_gene_pages_arabidopsis(g, gene, member_id)
                add_common_cluster_metadata(g, cluster, "MIBIG", species_label, taxon_id)
                g.add((cluster, RO_HAS_PART, gene))
                summary["species"][species_label]["bgcs"].add(str(cluster))
                summary["species"][species_label]["gene_links"] += 1
                summary["species"][species_label]["unique_genes"].add(str(gene))
                typed_this_cluster = True
                continue

            # ---- Tomato (safe: Solyc gene model) ----
            if RE_TOMATO_SOLYC.fullmatch(member_id):
                species_label, taxon_id = "Solanum lycopersicum", "4081"
                gene = URIRef(ENSEMBL_PLANT[member_id])
                add_common_gene_metadata(g, gene, species_label, taxon_id)
                g.add((gene, FOAF.page, URIRef(ENSEMBL_TOMATO_SEARCH + member_id)))
                add_common_cluster_metadata(g, cluster, "MIBIG", species_label, taxon_id)
                g.add((cluster, RO_HAS_PART, gene))
                summary["species"][species_label]["bgcs"].add(str(cluster))
                summary["species"][species_label]["gene_links"] += 1
                summary["species"][species_label]["unique_genes"].add(str(gene))
                typed_this_cluster = True
                continue

            # ---- Unknown: preserve membership, skip species annotation ----
            member = mint_member_identifier_iri(member_id)
            g.add((member, RDF.type, PMW.MemberIdentifier))
            g.add((member, RDFS.label, Literal(member_id)))
            g.add((cluster, RO_HAS_PART, member))
            summary["untyped"]["bgcs"].add(str(cluster))
            summary["untyped"]["member_links"] += 1
            summary["untyped"]["unique_members"].add(str(member))

        if not typed_this_cluster:
            summary["untyped"]["bgcs"].add(str(cluster))

    # Finalise counts
    final_species = {
        sp: {
            "bgcs":         len(s["bgcs"]),
            "gene_links":   s["gene_links"],
            "unique_genes": len(s["unique_genes"]),
        }
        for sp, s in summary["species"].items()
    }
    final_summary = {
        "source":           summary["source"],
        "generated_at":     summary["generated_at"],
        "species":          final_species,
        "bridgedb_resolved": summary["bridgedb_resolved"],
        "untyped": {
            "bgcs":           len(summary["untyped"]["bgcs"]),
            "member_links":   summary["untyped"]["member_links"],
            "unique_members": len(summary["untyped"]["unique_members"]),
        },
        "notes": [
            "Species only assigned when identifiers match reliable patterns.",
            "XP_/NP_/KNA_ etc. kept as pmw:MemberIdentifier (no species claim).",
        ],
    }
    return g, final_summary


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--test",
        action="store_true",
        help="Use root-level test JSON files (testPlantIsMash.json, testMibig.json)",
    )
    p.add_argument(
        "--source",
        choices=["both", "plantismash", "mibig"],
        default="both",
        help="Which source to convert (default: both)",
    )
    p.add_argument(
        "--no-bridgedb",
        action="store_true",
        help="Skip BridgeDb lookups (offline / CI mode)",
    )
    return p.parse_args()


# ============================================================
# Main
# ============================================================

def main() -> int:
    args = parse_args()
    use_bridgedb = not args.no_bridgedb

    OUTPUT_TTL_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    # Resolve input paths
    if args.test:
        plantismash_path = TEST_PLANTISMASH_JSON
        mibig_path       = TEST_MIBIG_JSON
        print("[TEST MODE] using root-level test JSON files")
    else:
        plantismash_path = INPUT_PLANTISMASH_JSON
        mibig_path       = INPUT_MIBIG_GENES_JSON

    if use_bridgedb:
        print("[INFO] BridgeDb lookups enabled (use --no-bridgedb to skip)")
    else:
        print("[INFO] BridgeDb lookups disabled")

    summaries: dict = {"generated_at": datetime.now().isoformat(timespec="seconds")}

    # ---- plantiSMASH ----
    if args.source in ("both", "plantismash"):
        if not plantismash_path.exists():
            print(f"[STOP] Missing: {plantismash_path}")
            return 1
        data = json.loads(plantismash_path.read_text(encoding="utf-8"))
        g, s = convert_plantismash(data, use_bridgedb=use_bridgedb)
        g.serialize(str(OUT_PLANTISMASH_TTL), format="turtle")
        print(f"✔ Wrote {OUT_PLANTISMASH_TTL} ({len(g)} triples)")
        summaries["plantiSMASH"] = s

    # ---- MIBiG ----
    if args.source in ("both", "mibig"):
        if not mibig_path.exists():
            print(f"[STOP] Missing: {mibig_path}")
            return 1
        data = json.loads(mibig_path.read_text(encoding="utf-8"))
        g, s = convert_mibig(data, use_bridgedb=use_bridgedb)
        g.serialize(str(OUT_MIBIG_TTL), format="turtle")
        print(f"✔ Wrote {OUT_MIBIG_TTL} ({len(g)} triples)")
        summaries["MIBIG"] = s

    OUT_SUMMARY.write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    print(f"✔ Wrote {OUT_SUMMARY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
