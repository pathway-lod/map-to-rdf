#!/usr/bin/env python3

import json
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, DCTERMS

# ============================================================
# Repo layout
# ============================================================

INPUT_DIR = Path("input")
OUTPUT_TTL_DIR = Path("output_ttl")
SUMMARY_DIR = Path("summaries")

INPUT_MIBIG_GENES_JSON = INPUT_DIR / "mibig_genes.json"
INPUT_PLANTISMASH_JSON = INPUT_DIR / "plantismash_v2_clusters_minimal.json"

OUT_PLANTISMASH_TTL = OUTPUT_TTL_DIR / "plantismash.ttl"
OUT_MIBIG_TTL = OUTPUT_TTL_DIR / "mibig.ttl"
OUT_SUMMARY = SUMMARY_DIR / "bgc_conversion_summary.json"

# ============================================================
# Namespaces / vocab
# ============================================================

PMW = Namespace("https://plantmetwiki.bioinformatics.nl/vocab/")
WP = Namespace("http://vocabularies.wikipathways.org/wp#")
FOAF = Namespace("http://xmlns.com/foaf/0.1/")
NCBITAXON = Namespace("http://purl.obolibrary.org/obo/NCBITaxon_")

# Use canonical RO term (has_part)
RO_HAS_PART = URIRef("http://purl.obolibrary.org/obo/RO_0000051")

# MiBIG stable identifier via Bioregistry (as requested)
MIBIG_BIOREGISTRY_PREFIX = "https://bioregistry.io/mibig:"

# plantiSMASH cluster namespace
PLANTISMASH = Namespace("https://plantismash.bioinformatics.nl/precalc/v2/")

# Our internal URI space for “raw member identifiers” that we cannot confidently type as genes
PMW_MEMBER_BASE = "https://plantmetwiki.bioinformatics.nl/id/member/"

# ============================================================
# Species detection and gene IRI minting
# ============================================================

# cluster_id strings contain species like "Arabidopsis_thaliana" or "Solanum_lycopersicum"
RE_AT_CLUSTER = re.compile(r"Arabidopsis_thaliana")
RE_SL_CLUSTER = re.compile(r"Solanum_lycopersicum")
RE_CROSS = re.compile(r"_x_")

# Tomato gene ID patterns in plantiSMASH
RE_SL_LOC = re.compile(r"^LOC(?P<geneid>\d+)$")
RE_SL_LOC_UNMATCHING = re.compile(r"^LOC\d+_\d+$")

# Arabidopsis gene IDs are usually like AT1G01010, AT1G01010.1 etc.
RE_AT_GENE = re.compile(r"^AT[1-5MC]G\d{5}(\.\d+)?$", re.IGNORECASE)

# Tomato gene model IDs (only treat as tomato if we see these; NOT XP_/NP_)
RE_TOMATO_SOLYC = re.compile(r"^Solyc\d{2}g\d{6}(\.\d+)+$", re.IGNORECASE)

# Prefer identifiers.org style for gene nodes (better for joining later)
TAIR_LOCUS = Namespace("https://identifiers.org/tair.locus/")
TAIR_NAME = Namespace("https://identifiers.org/tair.name/")

NCBIGENE = Namespace("https://identifiers.org/ncbigene:")
ENSEMBL_PLANT = Namespace("http://identifiers.org/ensembl.plant:")

# Useful “page” URLs (optional)
TAIR_SEARCH = "https://www.arabidopsis.org/results?mainType=general&category=genes&searchText="
ENSEMBL_TOMATO_SEARCH = "https://plants.ensembl.org/Solanum_lycopersicum/Search/Results?species=Solanum_lycopersicum;idx=;q="


def species_from_plantismash_cluster_id(cluster_id: str) -> tuple[str, str] | None:
    """Return (species_label, ncbitaxon_id) for supported species, else None."""
    if RE_CROSS.search(cluster_id):
        return None
    if RE_AT_CLUSTER.search(cluster_id):
        return ("Arabidopsis thaliana", "3702")
    if RE_SL_CLUSTER.search(cluster_id):
        return ("Solanum lycopersicum", "4081")
    return None


def mint_arabidopsis_gene_iri(gene_id: str) -> URIRef:
    """
    Mint a stable gene IRI for Arabidopsis.
    Prefer tair.locus because it commonly appears in pathway RDF exports.
    """
    core = gene_id.split(".")[0]
    return URIRef(TAIR_LOCUS[core])


def add_gene_pages_arabidopsis(g: Graph, gene: URIRef, gene_id: str) -> None:
    """Add linkouts (optional, but helpful)."""
    core = gene_id.split(".")[0]
    g.add((gene, FOAF.page, URIRef(TAIR_SEARCH + core)))
    g.add((gene, FOAF.page, URIRef(TAIR_NAME[core])))


def mint_tomato_gene_iri_from_plantismash(gene_id: str) -> tuple[URIRef | None, URIRef | None]:
    """
    Tomato: if LOC123 -> map to NCBIGene:123.
    Else -> treat as Ensembl Plant identifier string (best effort).
    Returns (gene_iri, optional_page_iri) or (None, None) if we skip.
    """
    m = RE_SL_LOC.fullmatch(gene_id)
    if m:
        return (URIRef(NCBIGENE[m.group("geneid")]), None)

    if RE_SL_LOC_UNMATCHING.search(gene_id):
        # skip junk LOC formats
        return (None, None)

    # "gene symbol / ensembl-ish identifier" – this is still a best-effort choice
    gene_iri = URIRef(ENSEMBL_PLANT[gene_id])
    page = URIRef(ENSEMBL_TOMATO_SEARCH + gene_id)
    return (gene_iri, page)


def mint_member_identifier_iri(member_id: str) -> URIRef:
    """Mint a PMW-local URI for raw member identifiers (safe for XP_/NP_/etc.)."""
    return URIRef(PMW_MEMBER_BASE + quote(member_id, safe=""))


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


def add_common_gene_metadata(g: Graph, gene: URIRef, species_label: str, taxon_id: str) -> None:
    g.add((gene, RDF.type, WP.GeneProduct))
    g.add((gene, WP.organismName, Literal(species_label)))
    g.add((gene, WP.organism, URIRef(NCBITAXON[taxon_id])))


# ============================================================
# Conversion functions
# ============================================================

def convert_plantismash(plantismash_json: dict) -> tuple[Graph, dict]:
    g = Graph()
    g.bind("pmw", PMW)
    g.bind("wp", WP)
    g.bind("dcterms", DCTERMS)
    g.bind("foaf", FOAF)
    g.bind("ncbi", NCBITAXON)
    g.bind("planti", PLANTISMASH)

    summary = {
        "source": "plantiSMASH",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "species": defaultdict(lambda: {"bgcs": 0, "gene_links": 0, "unique_genes": set()}),
    }

    for cluster_id, gene_list in plantismash_json.items():
        sp = species_from_plantismash_cluster_id(cluster_id)
        if sp is None:
            continue
        species_label, taxon_id = sp

        cluster = URIRef(PLANTISMASH[cluster_id])
        add_common_cluster_metadata(g, cluster, "plantiSMASH", species_label, taxon_id)

        summary["species"][species_label]["bgcs"] += 1

        for gene_id in gene_list:
            gene_id = gene_id.strip()
            if not gene_id:
                continue

            if species_label == "Arabidopsis thaliana":
                # Only mint as TAIR locus if it matches TAIR-ish pattern; otherwise keep as raw member
                if RE_AT_GENE.fullmatch(gene_id):
                    gene = mint_arabidopsis_gene_iri(gene_id)
                    add_common_gene_metadata(g, gene, species_label, taxon_id)
                    add_gene_pages_arabidopsis(g, gene, gene_id)
                    g.add((cluster, RO_HAS_PART, gene))

                    summary["species"][species_label]["gene_links"] += 1
                    summary["species"][species_label]["unique_genes"].add(str(gene))
                else:
                    member = mint_member_identifier_iri(gene_id)
                    g.add((member, RDF.type, PMW.MemberIdentifier))
                    g.add((member, RDFS.label, Literal(gene_id)))
                    g.add((cluster, RO_HAS_PART, member))
                continue

            if species_label == "Solanum lycopersicum":
                gene, page = mint_tomato_gene_iri_from_plantismash(gene_id)
                if gene is None:
                    # still preserve as raw member node
                    member = mint_member_identifier_iri(gene_id)
                    g.add((member, RDF.type, PMW.MemberIdentifier))
                    g.add((member, RDFS.label, Literal(gene_id)))
                    g.add((cluster, RO_HAS_PART, member))
                    continue

                add_common_gene_metadata(g, gene, species_label, taxon_id)
                if page is not None:
                    g.add((gene, FOAF.page, page))
                g.add((cluster, RO_HAS_PART, gene))

                summary["species"][species_label]["gene_links"] += 1
                summary["species"][species_label]["unique_genes"].add(str(gene))
                continue

    # Convert sets to counts for JSON
    for sp_label in list(summary["species"].keys()):
        summary["species"][sp_label]["unique_genes"] = len(summary["species"][sp_label]["unique_genes"])

    return g, summary


def convert_mibig(mibig_json: dict) -> tuple[Graph, dict]:
    """
    Convert mibig_genes.json (derived from knownclusters.txt) into RDF.

    Key behavior change vs Denise’s original approach:
    - DO NOT assume XP_/NP_/etc. implies tomato (or any species).
    - Only assign species when we have a reliable ID rule.
    - Always keep membership info by minting PMW.MemberIdentifier nodes for unknown IDs.
    """
    g = Graph()
    g.bind("pmw", PMW)
    g.bind("wp", WP)
    g.bind("dcterms", DCTERMS)
    g.bind("foaf", FOAF)
    g.bind("ncbi", NCBITAXON)
    g.bind("mibig", Namespace(MIBIG_BIOREGISTRY_PREFIX))

    summary = {
        "source": "MIBiG",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "species": defaultdict(lambda: {"bgcs": set(), "gene_links": 0, "unique_genes": set()}),
        "untyped": {"bgcs": set(), "member_links": 0, "unique_members": set()},
    }

    for bgc_id, member_list in mibig_json.items():
        cluster = URIRef(f"{MIBIG_BIOREGISTRY_PREFIX}{bgc_id}")

        # Always add minimal cluster metadata (no species assumptions here)
        add_common_cluster_metadata(g, cluster, "MIBIG", None, None)

        typed_this_cluster = False

        for member_id in member_list:
            member_id = str(member_id).strip()
            if not member_id:
                continue

            # ---- Arabidopsis (safe) ----
            if RE_AT_GENE.fullmatch(member_id):
                species_label, taxon_id = ("Arabidopsis thaliana", "3702")
                gene = mint_arabidopsis_gene_iri(member_id)

                # add gene + species metadata (safe)
                add_common_gene_metadata(g, gene, species_label, taxon_id)
                add_gene_pages_arabidopsis(g, gene, member_id)

                # optionally annotate cluster with species if we have at least one reliable member
                add_common_cluster_metadata(g, cluster, "MIBIG", species_label, taxon_id)

                g.add((cluster, RO_HAS_PART, gene))

                summary["species"][species_label]["bgcs"].add(str(cluster))
                summary["species"][species_label]["gene_links"] += 1
                summary["species"][species_label]["unique_genes"].add(str(gene))
                typed_this_cluster = True
                continue

            # ---- Tomato (only if really Solyc gene model; safe-ish) ----
            if RE_TOMATO_SOLYC.fullmatch(member_id):
                species_label, taxon_id = ("Solanum lycopersicum", "4081")
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

            # ---- Unknown identifier type (XP_/NP_/BAF_/etc.): keep, but don't type species ----
            member = mint_member_identifier_iri(member_id)
            g.add((member, RDF.type, PMW.MemberIdentifier))
            g.add((member, RDFS.label, Literal(member_id)))
            g.add((cluster, RO_HAS_PART, member))

            summary["untyped"]["bgcs"].add(str(cluster))
            summary["untyped"]["member_links"] += 1
            summary["untyped"]["unique_members"].add(str(member))

        # If nothing matched a reliable species rule, cluster stays untyped (but still present)
        if not typed_this_cluster:
            summary["untyped"]["bgcs"].add(str(cluster))

    # Finalize summary counts
    final_species = {}
    for sp_label, stats in summary["species"].items():
        final_species[sp_label] = {
            "bgcs": len(stats["bgcs"]),
            "gene_links": stats["gene_links"],
            "unique_genes": len(stats["unique_genes"]),
        }

    final_untyped = {
        "bgcs": len(summary["untyped"]["bgcs"]),
        "member_links": summary["untyped"]["member_links"],
        "unique_members": len(summary["untyped"]["unique_members"]),
    }

    final_summary = {
        "source": summary["source"],
        "generated_at": summary["generated_at"],
        "species": final_species,
        "untyped": final_untyped,
        "notes": [
            "MIBiG membership is preserved for all identifiers.",
            "Species annotations are only added when identifiers match reliable patterns (e.g., TAIR locus IDs).",
            "XP_/NP_/etc. are NOT assumed to be tomato (or any species); they are kept as pmw:MemberIdentifier nodes.",
        ],
    }

    return g, final_summary


# ============================================================
# Main
# ============================================================

def main() -> int:
    OUTPUT_TTL_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_PLANTISMASH_JSON.exists():
        print(f"[STOP] Missing input: {INPUT_PLANTISMASH_JSON}")
        return 1
    if not INPUT_MIBIG_GENES_JSON.exists():
        print(f"[STOP] Missing input: {INPUT_MIBIG_GENES_JSON}")
        return 1

    plantismash_json = json.loads(INPUT_PLANTISMASH_JSON.read_text(encoding="utf-8"))
    mibig_json = json.loads(INPUT_MIBIG_GENES_JSON.read_text(encoding="utf-8"))

    plant_g, plant_summary = convert_plantismash(plantismash_json)
    plant_g.serialize(str(OUT_PLANTISMASH_TTL), format="turtle")
    print(f"✔ Wrote {OUT_PLANTISMASH_TTL} ({len(plant_g)} triples)")

    mibig_g, mibig_summary = convert_mibig(mibig_json)
    mibig_g.serialize(str(OUT_MIBIG_TTL), format="turtle")
    print(f"✔ Wrote {OUT_MIBIG_TTL} ({len(mibig_g)} triples)")

    full_summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "outputs": {
            "plantiSMASH": str(OUT_PLANTISMASH_TTL),
            "MIBIG": str(OUT_MIBIG_TTL),
        },
        "plantiSMASH": plant_summary,
        "MIBIG": mibig_summary,
    }

    OUT_SUMMARY.write_text(json.dumps(full_summary, indent=2), encoding="utf-8")
    print(f"✔ Wrote {OUT_SUMMARY}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())