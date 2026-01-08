#!/usr/bin/env python3

import json
import re
import sys
import requests
from datetime import datetime
from pathlib import Path
from typing import Optional

from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, DCTERMS

# ============================================================
# Namespaces
# ============================================================

PMW = Namespace("https://plantmetwiki.bioinformatics.nl/vocab/")
RO = Namespace("http://purl.obolibrary.org/obo/RO_")

RO_HAS_PART = RO["0000051"]   # has_part
RO_PART_OF  = RO["0000050"]   # part_of (optional but nice)

# ============================================================
# External services
# ============================================================

BRIDGEDB_BASE = "https://webservice.bridgedb.org"
ORGANISM = "Arabidopsis thaliana"

# ============================================================
# Identifier conventions
# ============================================================

TAIR_URI_PREFIX = "https://identifiers.org/tair.name/"
WP_DATANODE = URIRef("http://vocabularies.wikipathways.org/wp#DataNode")

# MIBiG via Bioregistry (stable!)
MIBIG_BIOREGISTRY_PREFIX = "https://bioregistry.io/mibig:"

# ============================================================
# Repo layout
# ============================================================

INPUT_DIR = Path("input")
INPUT_TTL_DIR = Path("input_ttl")
OUTPUT_TTL_DIR = Path("output_ttl")
SUMMARY_DIR = Path("summaries")

INPUT_REACTIONS_TTL = INPUT_TTL_DIR / "all_pathways.ttl"
INPUT_MIBIG_GENES_JSON = INPUT_DIR / "mibig_genes.json"
INPUT_PLANTISMASH_JSON = INPUT_DIR / "plantismash_v2_clusters_minimal.json"

OUT_BGC_LINKS_TTL = OUTPUT_TTL_DIR / "bgc_links.ttl"
OUT_MERGED_TTL = OUTPUT_TTL_DIR / "reactions_with_bgc_links.ttl"
OUT_SUMMARY_LATEST = SUMMARY_DIR / "linking_summary_latest.json"
OUT_SUMMARY_LOG = SUMMARY_DIR / "linking_summary.log"

# ============================================================
# HTTP session reuse
# ============================================================

SESSION = requests.Session()


# ============================================================
# BridgeDb helpers
# ============================================================

def http_get_text(url: str, timeout: int = 30) -> str:
    r = SESSION.get(url, timeout=timeout, headers={"Accept": "text/plain"})
    if r.status_code >= 400:
        raise RuntimeError(f"HTTP {r.status_code} for {url}\n{r.text[:300]}")
    return r.text


def bridgedb_source_data_sources() -> list[str]:
    url = f"{BRIDGEDB_BASE}/{requests.utils.quote(ORGANISM)}/sourceDataSources"
    txt = http_get_text(url)
    return sorted({line.split("\t")[0] for line in txt.splitlines() if line.strip()})


def bridgedb_xrefs(system_code: str, identifier: str) -> list[tuple[str, str]]:
    url = f"{BRIDGEDB_BASE}/{requests.utils.quote(ORGANISM)}/xrefs/{system_code}/{identifier}"
    txt = http_get_text(url)
    out = []
    for line in txt.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            out.append((parts[1].strip(), parts[0].strip()))
    return out


def pick_working_source_code(tair_id: str, candidates: list[str]) -> Optional[str]:
    for code in candidates:
        try:
            if bridgedb_xrefs(code, tair_id):
                print(f"[OK] BridgeDb source '{code}' works for {tair_id}")
                return code
        except Exception as e:
            print(f"[FAIL] {code}: {e}")
    return None


# ============================================================
# RDF helpers
# ============================================================

def add_cluster_metadata(g: Graph, cluster_uri: str, source: str) -> None:
    c = URIRef(cluster_uri)
    g.add((c, RDF.type, PMW.BiosyntheticGeneCluster))

    # nicer label
    if source == "MIBIG":
        label = cluster_uri.split("mibig:", 1)[1]
    else:
        label = cluster_uri.rstrip("/").split("/")[-1]

    g.add((c, RDFS.label, Literal(label)))
    g.add((c, DCTERMS.source, Literal(source)))


# ============================================================
# Main
# ============================================================

def main() -> int:
    OUTPUT_TTL_DIR.mkdir(exist_ok=True)
    SUMMARY_DIR.mkdir(exist_ok=True)

    # ---- load inputs -------------------------------------------------------
    mibig = json.loads(INPUT_MIBIG_GENES_JSON.read_text())
    plantismash = json.loads(INPUT_PLANTISMASH_JSON.read_text())

    # MIBiG: gene → clusters (via Bioregistry!)
    mibig_by_gene: dict[str, set[str]] = {}
    for bgc_id, genes in mibig.items():
        cluster_uri = f"{MIBIG_BIOREGISTRY_PREFIX}{bgc_id}"
        for gene in genes:
            mibig_by_gene.setdefault(gene, set()).add(cluster_uri)

    # plantiSMASH: gene → clusters
    plant_by_gene: dict[str, set[str]] = {}
    for cid, genes in plantismash.items():
        cluster_uri = "https://plantismash.bioinformatics.nl/precalc/v2/" + cid
        for gene in genes:
            plant_by_gene.setdefault(gene, set()).add(cluster_uri)

    # ---- parse pathway RDF -------------------------------------------------
    base = Graph()
    base.parse(INPUT_REACTIONS_TTL, format="turtle")

    tair_genes = {
        str(s)
        for s in base.subjects(RDF.type, WP_DATANODE)
        if str(s).startswith(TAIR_URI_PREFIX)
    }

    print(f"TAIR genes in pathways: {len(tair_genes)}")

    # ---- BridgeDb setup ----------------------------------------------------
    supported = bridgedb_source_data_sources()
    candidates = [c for c in ["A", "L", "T", "En"] if c in supported or c == "A"]

    working_source = pick_working_source_code("AT5G48000", candidates)
    if not working_source:
        print("[STOP] No working BridgeDb source code")
        return 1

    # ---- build output graph -----------------------------------------------
    out = Graph()
    out.bind("pmw", PMW)
    out.bind("ro", RO)

    tair_re = re.compile(r"tair\.name/([^/]+)$")

    seen_clusters = set()
    has_part_triples = 0

    direct_mibig_hits = 0
    mapped_mibig_hits = 0
    mapped_plant_hits = 0
    mapped_total = 0

    for gene_uri in sorted(tair_genes):
        m = tair_re.search(gene_uri)
        if not m:
            continue

        tair_id = m.group(1).split(".")[0]
        gene_ref = URIRef(gene_uri)

        # ---- direct MIBiG gene match --------------------------------------
        if tair_id in mibig_by_gene:
            direct_mibig_hits += 1
            for cluster_uri in mibig_by_gene[tair_id]:
                c = URIRef(cluster_uri)
                out.add((c, RO_HAS_PART, gene_ref))
                out.add((gene_ref, RO_PART_OF, c))
                has_part_triples += 1
                if cluster_uri not in seen_clusters:
                    add_cluster_metadata(out, cluster_uri, "MIBIG")
                    seen_clusters.add(cluster_uri)

        # ---- BridgeDb mappings --------------------------------------------
        for _, mapped_id in bridgedb_xrefs(working_source, tair_id):
            mapped_total += 1

            if mapped_id in mibig_by_gene:
                mapped_mibig_hits += 1
                for cluster_uri in mibig_by_gene[mapped_id]:
                    c = URIRef(cluster_uri)
                    out.add((c, RO_HAS_PART, gene_ref))
                    out.add((gene_ref, RO_PART_OF, c))
                    has_part_triples += 1
                    if cluster_uri not in seen_clusters:
                        add_cluster_metadata(out, cluster_uri, "MIBIG")
                        seen_clusters.add(cluster_uri)

            if mapped_id in plant_by_gene:
                mapped_plant_hits += 1
                for cluster_uri in plant_by_gene[mapped_id]:
                    c = URIRef(cluster_uri)
                    out.add((c, RO_HAS_PART, gene_ref))
                    out.add((gene_ref, RO_PART_OF, c))
                    has_part_triples += 1
                    if cluster_uri not in seen_clusters:
                        add_cluster_metadata(out, cluster_uri, "plantiSMASH")
                        seen_clusters.add(cluster_uri)

    # ---- write outputs ----------------------------------------------------
    out.serialize(OUT_BGC_LINKS_TTL, format="turtle")

    merged = Graph()
    merged.parse(INPUT_REACTIONS_TTL, format="turtle")
    merged.parse(OUT_BGC_LINKS_TTL, format="turtle")
    merged.serialize(OUT_MERGED_TTL, format="turtle")

    summary = {
        "date": datetime.now().isoformat(timespec="seconds"),
        "organism": ORGANISM,
        "bridgedb_source_code": working_source,
        "tair_genes_in_ttl": len(tair_genes),
        "direct_mibig_gene_hits": direct_mibig_hits,
        "mapped_ids_hit_mibig": mapped_mibig_hits,
        "mapped_ids_hit_plantismash": mapped_plant_hits,
        "mapped_xref_pairs_total": mapped_total,
        "ro_has_part_triples": has_part_triples,
        "output_ttl": str(OUT_BGC_LINKS_TTL),
    }

    OUT_SUMMARY_LATEST.write_text(json.dumps(summary, indent=2))
    OUT_SUMMARY_LOG.open("a").write(json.dumps(summary) + "\n")

    print("✔ Crosslinks written")
    print(json.dumps(summary, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())