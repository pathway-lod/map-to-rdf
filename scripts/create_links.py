import json
import re
import sys
import requests
from datetime import datetime
from pathlib import Path
from typing import Optional

from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, DCTERMS

# ---- Constants / Paths (repo layout) ----
PMW = Namespace("https://plantmetwiki.bioinformatics.nl/vocab/")

BRIDGEDB_BASE = "https://webservice.bridgedb.org"
ORGANISM = "Arabidopsis thaliana"

TAIR_URI_PREFIX = "https://identifiers.org/tair.locus/"
WP_DATANODE = URIRef("http://vocabularies.wikipathways.org/wp#DataNode")

INPUT_DIR = Path("input")
INPUT_TTL_DIR = Path("input_ttl")
OUTPUT_TTL_DIR = Path("output_ttl")
SUMMARY_DIR = Path("summaries")

INPUT_REACTIONS_TTL = INPUT_TTL_DIR / "reactions.ttl"
INPUT_MIBIG_GENES_JSON = INPUT_DIR / "mibig_genes.json"
INPUT_PLANTISMASH_JSON = INPUT_DIR / "plantismash_v2_clusters_minimal.json"

OUT_BGC_LINKS_TTL = OUTPUT_TTL_DIR / "bgc_links.ttl"
OUT_MERGED_TTL = OUTPUT_TTL_DIR / "reactions_with_bgc_links.ttl"
OUT_SUMMARY_LATEST = SUMMARY_DIR / "linking_summary_latest.json"
OUT_SUMMARY_LOG = SUMMARY_DIR / "linking_summary.log"

# Reuse HTTP connections (faster)
SESSION = requests.Session()


def http_get_text(url: str, timeout: int = 30) -> str:
    r = SESSION.get(url, timeout=timeout, headers={"Accept": "text/plain"})
    if r.status_code >= 400:
        raise RuntimeError(f"HTTP {r.status_code} for {url}\nBody:\n{r.text[:500]}")
    return r.text


def bridgedb_source_data_sources() -> list[str]:
    """Return list of source system codes supported for this organism by this webservice."""
    url = f"{BRIDGEDB_BASE}/{requests.utils.quote(ORGANISM)}/sourceDataSources"
    txt = http_get_text(url)
    codes: list[str] = []
    for line in txt.strip().splitlines():
        if not line.strip():
            continue
        codes.append(line.split("\t")[0].strip())
    return sorted(set(codes))


def bridgedb_xrefs(system_code: str, identifier: str) -> list[tuple[str, str]]:
    """Returns list of (mapped_code, mapped_id) from BridgeDb webservice."""
    url = f"{BRIDGEDB_BASE}/{requests.utils.quote(ORGANISM)}/xrefs/{system_code}/{identifier}"
    txt = http_get_text(url)
    out: list[tuple[str, str]] = []
    for line in txt.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            mapped_id = parts[0].strip()
            mapped_code = parts[1].strip()
            out.append((mapped_code, mapped_id))
    return out


def pick_working_source_code(tair_id: str, preferred: list[str]) -> Optional[str]:
    """Try plausible source codes until we get non-empty xrefs."""
    for code in preferred:
        try:
            xrefs = bridgedb_xrefs(code, tair_id)
            if xrefs:
                print(f"[OK] Source code '{code}' returns {len(xrefs)} xrefs for {tair_id}")
                return code
            else:
                print(f"[NO XREFS] Source code '{code}' returns 0 xrefs for {tair_id}")
        except Exception as e:
            print(f"[FAIL] Source code '{code}' failed for {tair_id}: {e}")
    return None


def add_bgc_metadata(outg: Graph, bgc_uri: str, source: str) -> None:
    """Add minimal metadata for BGC/cluster nodes."""
    bgc_ref = URIRef(bgc_uri)
    outg.add((bgc_ref, RDF.type, PMW.BiosyntheticGeneCluster))

    label = bgc_uri.rstrip("/").split("/")[-1]
    outg.add((bgc_ref, RDFS.label, Literal(label)))

    outg.add((bgc_ref, DCTERMS.source, Literal(source)))


def write_summary_files(
    working_source: str,
    tair_genes_in_ttl: int,
    mapped_xref_pairs_total: int,
    direct_mibig_gene_hits: int,
    mapped_ids_hit_mibig: int,
    mapped_ids_hit_plantismash: int,
    triples_written: int,
    output_ttl: Path,
) -> None:
    summary = {
        "date": datetime.now().isoformat(timespec="seconds"),
        "organism": ORGANISM,
        "bridgedb_source_code": working_source,
        "tair_genes_in_ttl": tair_genes_in_ttl,
        "mapped_xref_pairs_total": mapped_xref_pairs_total,
        "direct_mibig_gene_hits": direct_mibig_gene_hits,
        "mapped_ids_hit_mibig": mapped_ids_hit_mibig,
        "mapped_ids_hit_plantismash": mapped_ids_hit_plantismash,
        "triples_written": triples_written,
        "output_ttl": str(output_ttl),
    }

    OUT_SUMMARY_LATEST.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    with OUT_SUMMARY_LOG.open("a", encoding="utf-8") as fh:
        fh.write(
            (
                f"=== Summary ({summary['date']}) ===\n"
                f"Organism: {summary['organism']}\n"
                f"Working BridgeDb source code used: {summary['bridgedb_source_code']}\n"
                f"TAIR genes found in TTL: {summary['tair_genes_in_ttl']}\n"
                f"Mapped xref pairs seen total: {summary['mapped_xref_pairs_total']}\n"
                f"Direct MIBIG gene hits (no BridgeDb): {summary['direct_mibig_gene_hits']}\n"
                f"Mapped IDs that matched MIBIG members: {summary['mapped_ids_hit_mibig']}\n"
                f"Mapped IDs that matched plantiSMASH genes: {summary['mapped_ids_hit_plantismash']}\n"
                f"Triples written (gene → belongsToBGC → bgc): {summary['triples_written']}\n"
                f"Wrote {summary['output_ttl']}\n\n"
            )
        )


def main() -> int:
    # Ensure directories exist
    OUTPUT_TTL_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    # Sanity check input files
    for p in [INPUT_REACTIONS_TTL, INPUT_MIBIG_GENES_JSON, INPUT_PLANTISMASH_JSON]:
        if not p.exists():
            print(f"[STOP] Missing required input file: {p}")
            return 1

    # 1) Load maps
    with INPUT_MIBIG_GENES_JSON.open() as f:
        mibig = json.load(f)  # {BGC: [genes...]}

    with INPUT_PLANTISMASH_JSON.open() as f:
        plantismash = json.load(f)  # {cluster: [genes...]}

    # invert: member_id -> set(bgc_uri)
    mibig_inv: dict[str, set[str]] = {}
    for bgc_id, members in mibig.items():
        bgc_uri = f"https://mibig.secondarymetabolites.org/repository/{bgc_id}"
        for m in members:
            mibig_inv.setdefault(m, set()).add(bgc_uri)

    plant_inv: dict[str, set[str]] = {}
    for cluster_id, genes in plantismash.items():
        cluster_uri = "https://plantismash.bioinformatics.nl/precalc/v2/" + cluster_id
        for g in genes:
            plant_inv.setdefault(g, set()).add(cluster_uri)

    print(f"MIBIG unique gene IDs: {len(mibig_inv)}")
    print(f"plantiSMASH unique gene IDs: {len(plant_inv)}")

    # 2) Parse TTL and collect TAIR URIs
    g = Graph()
    g.parse(str(INPUT_REACTIONS_TTL), format="turtle")

    tair_gene_uris = set()
    for s in g.subjects(RDF.type, WP_DATANODE):
        if str(s).startswith(TAIR_URI_PREFIX):
            tair_gene_uris.add(str(s))

    print(f"TAIR genes found in TTL: {len(tair_gene_uris)}")
    print("Example TAIR genes:", list(sorted(tair_gene_uris))[:5])

    if not tair_gene_uris:
        print(
            "\n[STOP] No TAIR genes were found. "
            "Likely the RDF type is different in this TTL, or TAIR URIs appear as objects not subjects.\n"
        )
        return 2

    # 3) Discover supported sources (BridgeDb)
    try:
        supported_sources = bridgedb_source_data_sources()
        print(f"BridgeDb supported source system codes for {ORGANISM}: {len(supported_sources)}")
    except Exception as e:
        print(f"[STOP] Could not query BridgeDb sourceDataSources: {e}")
        return 3

    # 4) Choose a working source code
    candidates = ["A", "L", "T", "En", "At"]
    candidates = [c for c in candidates if (c == "A" or c in supported_sources)]

    test_id = "AT5G48000"
    working_source = pick_working_source_code(test_id, candidates)
    if working_source is None:
        print(
            "\n[STOP] None of the tested source codes returned xrefs.\n"
            "This usually means the BridgeDb instance you're calling doesn't have a suitable mapping DB loaded.\n"
        )
        return 4

    # 5) Build overlay triples
    outg = Graph()
    outg.bind("pmw", PMW)

    tair_id_re = re.compile(r"^https://identifiers\.org/tair\.locus/([^/]+)$")

    mapped_total = 0
    direct_mibig_gene_hits = 0
    mapped_ids_hit_mibig = 0
    mapped_ids_hit_plantismash = 0

    seen_bgc_nodes: set[str] = set()

    for gene_uri in sorted(tair_gene_uris):
        m = tair_id_re.match(gene_uri)
        if not m:
            continue

        tair_id = m.group(1)
        tair_id_core = tair_id.split(".")[0]  # strip isoform suffix

        # ---- Direct MIBIG gene match (fast, recommended) ----
        if tair_id_core in mibig_inv:
            direct_mibig_gene_hits += 1
            for bgc_uri in mibig_inv[tair_id_core]:
                outg.add((URIRef(gene_uri), PMW.belongsToBGC, URIRef(bgc_uri)))
                if bgc_uri not in seen_bgc_nodes:
                    add_bgc_metadata(outg, bgc_uri, "MIBIG")
                    seen_bgc_nodes.add(bgc_uri)

        # ---- BridgeDb mappings (needed for plantiSMASH / other id types) ----
        try:
            xrefs = bridgedb_xrefs(working_source, tair_id_core)
        except Exception as e:
            print(f"[BridgeDb ERROR] {tair_id_core}: {e}")
            continue

        for mapped_code, mapped_id in xrefs:
            mapped_total += 1

            # MIBIG via mapped IDs (rare if using gene-based JSON, but keep for robustness)
            if mapped_id in mibig_inv:
                mapped_ids_hit_mibig += 1
                for bgc_uri in mibig_inv[mapped_id]:
                    outg.add((URIRef(gene_uri), PMW.belongsToBGC, URIRef(bgc_uri)))
                    if bgc_uri not in seen_bgc_nodes:
                        add_bgc_metadata(outg, bgc_uri, "MIBIG")
                        seen_bgc_nodes.add(bgc_uri)

            # plantiSMASH matches
            if mapped_id in plant_inv:
                mapped_ids_hit_plantismash += 1
                for cluster_uri in plant_inv[mapped_id]:
                    outg.add((URIRef(gene_uri), PMW.belongsToBGC, URIRef(cluster_uri)))
                    if cluster_uri not in seen_bgc_nodes:
                        add_bgc_metadata(outg, cluster_uri, "plantiSMASH")
                        seen_bgc_nodes.add(cluster_uri)

    print("\n=== Summary ===")
    print(f"Working BridgeDb source code used: {working_source}")
    print(f"Mapped xref pairs seen total: {mapped_total}")
    print(f"Direct MIBIG gene hits (no BridgeDb): {direct_mibig_gene_hits}")
    print(f"Mapped IDs that matched MIBIG members: {mapped_ids_hit_mibig}")
    print(f"Mapped IDs that matched plantiSMASH genes: {mapped_ids_hit_plantismash}")
    print(f"Triples written (gene → belongsToBGC → bgc): {len(outg)}")

    outg.serialize(str(OUT_BGC_LINKS_TTL), format="turtle")
    print(f"Wrote {OUT_BGC_LINKS_TTL}")

    merged = Graph()
    merged.parse(str(INPUT_REACTIONS_TTL), format="turtle")
    merged.parse(str(OUT_BGC_LINKS_TTL), format="turtle")
    merged.serialize(str(OUT_MERGED_TTL), format="turtle")
    print(f"Wrote {OUT_MERGED_TTL}")

    write_summary_files(
        working_source=working_source,
        tair_genes_in_ttl=len(tair_gene_uris),
        mapped_xref_pairs_total=mapped_total,
        direct_mibig_gene_hits=direct_mibig_gene_hits,
        mapped_ids_hit_mibig=mapped_ids_hit_mibig,
        mapped_ids_hit_plantismash=mapped_ids_hit_plantismash,
        triples_written=len(outg),
        output_ttl=OUT_BGC_LINKS_TTL,
    )

    if len(outg) == 0:
        print(
            "\n[NOTE] bgc_links.ttl is empty.\n"
            "No overlaps were found between pathway genes and MIBIG/plantiSMASH mappings.\n"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())