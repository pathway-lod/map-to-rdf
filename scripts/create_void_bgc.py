#!/usr/bin/env python3
"""Create a VoID metadata file for the BGC RDF datasets (MIBiG + plantiSMASH).

Follows the same structure as the VoID produced by the gpml-to-rdf pipeline,
so both VoID files can be loaded into the same named graph in Virtuoso.

Usage
-----
    python scripts/create_void_bgc.py --output output_ttl/void-bgc.ttl
"""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path


PREFIXES = """\
@prefix void:    <http://rdfs.org/ns/void#> .
@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix pav:     <http://purl.org/pav/> .
@prefix dcat:    <http://www.w3.org/ns/dcat#> .
@prefix foaf:    <http://xmlns.com/foaf/0.1/> .
@prefix xsd:     <http://www.w3.org/2001/XMLSchema#> .
@prefix cc:      <http://creativecommons.org/ns#> .
"""

# -----------------------------------------------------------------------
# Publisher organization (shared across all PlantMetWiki VoID files)
# -----------------------------------------------------------------------

PUBLISHER = """\
# ── Publisher ────────────────────────────────────────────────────────────
<http://rdf-plantmetwiki.bioinformatics.nl/organization/wur-plant-sciences> a foaf:Organization ;
    foaf:name "Wageningen University & Research, Department of Plant Sciences"@en ;
    foaf:homepage <https://www.wur.nl/> .
"""

# -----------------------------------------------------------------------
# Licence documents
# -----------------------------------------------------------------------

LICENSES = """\
# ── MIBiG licence ──────────────────────────────────────────────────────
<http://rdf-plantmetwiki.bioinformatics.nl/license/cc-by-4.0> a dcterms:LicenseDocument ;
    dcterms:title "Creative Commons Attribution 4.0 International (CC BY 4.0)"@en ;
    dcterms:description "You are free to share and adapt this material for any purpose, provided you give appropriate credit."@en ;
    foaf:page <https://creativecommons.org/licenses/by/4.0/> .

# ── plantiSMASH licence ─────────────────────────────────────────────────
<http://rdf-plantmetwiki.bioinformatics.nl/license/gpl-3.0> a dcterms:LicenseDocument ;
    dcterms:title "GNU General Public License v3.0 (GPL-3.0)"@en ;
    dcterms:description "Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications."@en ;
    foaf:page <https://www.gnu.org/licenses/gpl-3.0.en.html> .
"""

# -----------------------------------------------------------------------
# Dataset descriptions
# -----------------------------------------------------------------------

DATASETS = """\
# ── MIBiG dataset ───────────────────────────────────────────────────────
<http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/mibig-4.0> a void:Dataset ;
    dcterms:title "PlantMetWiki BGC crosslinks — MIBiG 4.0"@en ;
    dcterms:description "RDF representation of curated biosynthetic gene cluster (BGC) membership data from MIBiG 4.0. Clusters are typed as pmw:BiosyntheticGeneCluster and linked to their constituent genes via RO:0000051 (has_part). Gene identifiers are normalised to identifiers.org/tair.name/ for Arabidopsis thaliana. Unknown identifiers are preserved as pmw:MemberIdentifier nodes."@en ;
    dcterms:hasVersion "mibig-4.0" ;
    pav:version "mibig-4.0" ;
    pav:createdOn "{today}"^^xsd:date ;
    dcterms:modified "{today}"^^xsd:date ;
    dcterms:publisher <http://rdf-plantmetwiki.bioinformatics.nl/organization/wur-plant-sciences> ;
    pav:createdWith <https://github.com/pathway-lod/map-to-rdf> ;
    pav:derivedFrom <https://github.com/mite-standard/mite_data> ;
    dcterms:source <https://github.com/plantismash/plantismash/blob/master/antismash/generic_modules/knownclusterblast/knownclusters.txt> ;
    dcterms:license <http://rdf-plantmetwiki.bioinformatics.nl/license/cc-by-4.0> ;
    dcterms:rights "MIBiG data is available under the Creative Commons Attribution 4.0 International (CC BY 4.0) license. Please cite: Terlouw et al. 2023 (doi:10.1093/nar/gkac1049)."@en ;
    void:sparqlEndpoint <https://plantmetwiki.bioinformatics.nl/sparql> ;
    void:vocabulary <http://rdf-plantmetwiki.bioinformatics.nl/vocab/> ,
                    <http://vocabularies.wikipathways.org/wp#> ,
                    <http://purl.obolibrary.org/obo/RO_> ,
                    <http://purl.org/dc/terms/> ;
    foaf:homepage <https://plantmetwiki.bioinformatics.nl/> ;
    foaf:page <https://mibig.secondarymetabolites.org/> ;
    foaf:page <https://doi.org/10.5281/zenodo.20345133> ;
    void:triples {mibig_triples} ;
    dcat:distribution <http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/mibig-4.0/distribution/mibig.ttl> .

<http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/mibig-4.0/distribution/mibig.ttl> a dcat:Distribution ;
    dcterms:title "mibig.ttl" ;
    dcat:byteSize {mibig_bytes} .

# ── plantiSMASH dataset ─────────────────────────────────────────────────
<http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/plantismash-v2> a void:Dataset ;
    dcterms:title "PlantMetWiki BGC crosslinks — plantiSMASH v2"@en ;
    dcterms:description "RDF representation of predicted biosynthetic gene cluster (BGC) membership data from the plantiSMASH pre-calculated database (v2), currently scoped to Arabidopsis thaliana clusters. Clusters are typed as pmw:BiosyntheticGeneCluster and linked to their constituent genes via RO:0000051 (has_part). Gene identifiers use identifiers.org/tair.name/ IRIs enabling direct SPARQL joins with the PlantMetWiki pathway graph."@en ;
    dcterms:hasVersion "plantismash-v2" ;
    pav:version "plantismash-v2" ;
    pav:createdOn "{today}"^^xsd:date ;
    dcterms:modified "{today}"^^xsd:date ;
    dcterms:publisher <http://rdf-plantmetwiki.bioinformatics.nl/organization/wur-plant-sciences> ;
    pav:createdWith <https://github.com/pathway-lod/map-to-rdf> ;
    pav:derivedFrom <https://github.com/plantismash/plantismash-database> ;
    dcterms:source <https://raw.githubusercontent.com/plantismash/plantismash-database/main/data/plantismash_v2_clusters_minimal.json> ;
    dcterms:license <http://rdf-plantmetwiki.bioinformatics.nl/license/gpl-3.0> ;
    dcterms:rights "plantiSMASH data is distributed under the GNU General Public License v3.0. Please cite: Kautsar et al. 2017 (doi:10.1093/nar/gkx338)."@en ;
    void:sparqlEndpoint <https://plantmetwiki.bioinformatics.nl/sparql> ;
    void:vocabulary <http://rdf-plantmetwiki.bioinformatics.nl/vocab/> ,
                    <http://vocabularies.wikipathways.org/wp#> ,
                    <http://purl.obolibrary.org/obo/RO_> ,
                    <http://purl.org/dc/terms/> ;
    foaf:homepage <https://plantmetwiki.bioinformatics.nl/> ;
    foaf:page <https://plantismash.bioinformatics.nl/> ;
    foaf:page <https://doi.org/10.5281/zenodo.20345133> ;
    void:triples {plantismash_triples} ;
    dcat:distribution <http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/plantismash-v2/distribution/plantismash.ttl> .

<http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/plantismash-v2/distribution/plantismash.ttl> a dcat:Distribution ;
    dcterms:title "plantismash.ttl" ;
    dcat:byteSize {plantismash_bytes} .

# ── Linkset: BGC genes ↔ PlantMetWiki pathway genes ────────────────────
<http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/linkset/bgc-pathway-genes> a void:Linkset ;
    dcterms:title "BGC–pathway gene crosslinks (Arabidopsis thaliana)"@en ;
    dcterms:description "Links between biosynthetic gene cluster member genes (plantiSMASH / MIBiG) and PlantMetWiki pathway gene nodes, established via shared identifiers.org/tair.name/ IRIs for Arabidopsis thaliana. Derived by SPARQL join across the BGC and pathway graphs."@en ;
    void:subjectsTarget <http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/plantismash-v2> ;
    void:objectsTarget  <http://rdf-plantmetwiki.bioinformatics.nl/dataset/plantcyc17.0.0-gpml2021/core> ;
    void:linkPredicate  <http://purl.obolibrary.org/obo/RO_0000051> ;
    dcterms:isPartOf <http://rdf-plantmetwiki.bioinformatics.nl/dataset/bgc/plantismash-v2> ;
    void:triples {link_triples} ;
    pav:createdOn "{today}"^^xsd:date .
"""


def count_triples(ttl_path: Path) -> int:
    """Count triples by counting '.' at end of lines (fast, avoids full parse)."""
    try:
        from rdflib import Graph
        g = Graph()
        g.parse(str(ttl_path), format="turtle")
        return len(g)
    except Exception:
        # Fallback: count lines ending in ' .' or '\t.'
        count = 0
        with ttl_path.open(encoding="utf-8", errors="replace") as f:
            for line in f:
                s = line.strip()
                if s.endswith(" .") or s.endswith("\t.") or s == ".":
                    count += 1
        return count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plantismash", type=Path, default=Path("output_ttl/plantismash.ttl"))
    parser.add_argument("--mibig",       type=Path, default=Path("output_ttl/mibig.ttl"))
    parser.add_argument("--links",       type=Path, default=Path("summaries/bgc_pathway_links.tsv"),
                        help="BGC–pathway links TSV (to count link triples)")
    parser.add_argument("--output",      type=Path, default=Path("output_ttl/void-bgc.ttl"))
    args = parser.parse_args()

    today = date.today().isoformat()

    # File sizes
    mibig_bytes       = args.mibig.stat().st_size       if args.mibig.exists()       else 0
    plantismash_bytes = args.plantismash.stat().st_size if args.plantismash.exists() else 0

    # Triple counts
    print(f"Counting triples in {args.mibig} ...")
    mibig_triples = count_triples(args.mibig) if args.mibig.exists() else 0
    print(f"  {mibig_triples:,} triples")

    print(f"Counting triples in {args.plantismash} ...")
    plantismash_triples = count_triples(args.plantismash) if args.plantismash.exists() else 0
    print(f"  {plantismash_triples:,} triples")

    # Link count (rows in TSV = unique gene-pathway-BGC links)
    link_triples = 0
    if args.links.exists():
        import csv
        with args.links.open(encoding="utf-8") as f:
            link_triples = sum(1 for _ in csv.DictReader(f, delimiter="\t"))
    print(f"Link triples (from TSV): {link_triples}")

    # Render
    content = (
        PREFIXES + "\n"
        + PUBLISHER + "\n"
        + LICENSES + "\n"
        + DATASETS.format(
            today=today,
            mibig_bytes=mibig_bytes,
            mibig_triples=mibig_triples,
            plantismash_bytes=plantismash_bytes,
            plantismash_triples=plantismash_triples,
            link_triples=link_triples,
        )
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(content, encoding="utf-8")
    print(f"\n✔ Wrote {args.output}")
    print(f"\nLoad into Virtuoso alongside the pathway VoID:")
    print(f"  LOAD <file://{args.output.resolve()}> INTO GRAPH <http://rdf-plantmetwiki.bioinformatics.nl/void>")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
