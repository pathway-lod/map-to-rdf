# map-to-rdf

Mapping biosynthetic gene clusters (BGCs) to PlantMetWiki / WikiPathways RDF

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20345133-blue)](https://doi.org/10.5281/zenodo.20345133)
[![GitHub Release](https://img.shields.io/github/v/release/pathway-lod/map-to-rdf?include_prereleases&sort=semver&display_name=tag&logo=github)](https://github.com/pathway-lod/map-to-rdf/releases)
[![PlantMetWiki](https://img.shields.io/static/v1?label=webserver&message=PlantMetWiki&color=green)](https://plantmetwiki.bioinformatics.nl/)

Visit the PlantMetWiki webserver at: https://plantmetwiki.bioinformatics.nl/

This repository converts biosynthetic gene cluster (BGC) data from MIBiG and plantiSMASH into interoperable RDF, suitable for integration with WikiPathways / PlantMetWiki and querying via SPARQL.

⚠️ Important change in scope
This repository now only converts BGC sources to RDF.
Linking BGCs to pathways is intentionally done at the SPARQL/query layer, not during RDF generation.

## What this repository produces

The pipeline generates:
	•	output_ttl/plantismash.ttl
→ RDF representation of plantiSMASH-predicted BGCs
	•	output_ttl/mibig.ttl
→ RDF representation of MIBiG curated BGCs
	•	summaries/bgc_conversion_summary.json
→ per-source, per-species summary of generated content

These RDF files can then be loaded into Virtuoso together with pathway RDF and queried jointly.

## Overview of the approach 

1. BGC-first modeling (not pathway-first)

Biosynthetic gene clusters are treated as first-class RDF resources, independent of pathways.

``` 
BGC  ── ro:0000051 (has_part) ──▶ gene / member
``` 

Pathways are not modified at conversion time.

2. Conservative species assignment (no guessing)

We only assign species when identifiers are reliable like in A. thaliana.

3. Known limitation: tomato annotation incompatibility

Direct BGC↔pathway linking currently works only for *Arabidopsis thaliana*.
For *Solanum lycopersicum* (tomato), the two data sources use gene models from
**different genome annotations**:

| Source | Genome annotation | Gene ID format |
|---|---|---|
| plantiSMASH | NCBI RefSeq (`GCF_036512215.1`) | Gene symbols (`ABCG1`, `loxC`) and LOC IDs |
| PlantCyc / pathway graph | Phytozome / ITAG | Solyc IDs (`Solyc01g105890.2.1`) |

NCBI RefSeq and the ITAG/Phytozome annotation generate different gene models
and identifiers for the same genome. Gene symbols (e.g. `ABCG1`) are shared
in principle but not systematically assigned, and the mapping is not 1:1.
BridgeDb maintains a partial crosswalk between ITAG (Solyc) and NCBI Gene
(LOC) IDs, but coverage is incomplete due to differences in gene model
boundaries between annotation versions.

As a result, tomato BGC clusters from plantiSMASH (1,025 gene links across
63 clusters in the current dataset) cannot be directly joined to PlantCyc
pathway genes without a BridgeDb-assisted identifier resolution step.
This is a known data limitation, not a pipeline bug.

4. Untyped MiBIG clusters are preserved (by design)

Many MIBiG clusters contain only protein accessions or non-locus identifiers.

Instead of discarding them:
	•	The cluster is kept
	•	Gene membership is preserved
	•	No wp:organism is added

These are reported as “untyped” BGCs in the summary.

This avoids:
	•	incorrect biological claims
	•	loss of curated knowledge


## Ontologies and predicates used

### Core model


```
BGC
 ├─ ro:0000051 → gene / member
 ├─ rdf:type → pmw:BiosyntheticGeneCluster
 └─ dcterms:source → "MIBIG" | "plantiSMASH"
```
Ontologies used

| Purpose  | Ontology / Term  | 
|---|---|
|  Gene–cluster relation |  RO:0000051 (has_part) |  
| Cluster type |  pmw:BiosyntheticGeneCluster |   
| Gene type | wp:GeneProduct | 
| Provenance |  dcterms:source | 
| Species | wp:organism, wp:organismName | 


### Stable identifiers for MIBiG (via Bioregistry) 

MIBIG mapping files to enzymes are available on MITE at https://github.com/mite-standard/mite_data/blob/dev/mite_data/mibig/mibig_proteins.json

The list of MIBiG 4.0 BGCs is available from the plantiSMASH gitHub repository: https://github.com/plantismash/plantismash/blob/fea399bdd040f49f7cfd9018e8cbec51f9ba6684/antismash/generic_modules/knownclusterblast/knownclusters.txt#L4 

They can be obtained with: 

```
mkdir input 
cd input
wget https://raw.githubusercontent.com/mite-standard/mite_data/main/mite_data/mibig/mibig_proteins.json
# or from plantiSMASH a txt file with all the gene identifiers 
wget https://raw.githubusercontent.com/plantismash/plantismash/refs/heads/master/antismash/generic_modules/knownclusterblast/knownclusters.txt
```

This pipeline uses gene-based MIBiG cluster membership, derived from plantiSMASH knownclusters list (MIBiG 4.x)
https://github.com/plantismash/plantismash/blob/master/antismash/generic_modules/knownclusterblast/knownclusters.txt

This file is parsed into a JSON mapping:
```
 {
  "BGC0000670": ["AT5G48000", "AT5G48010", "..."],
  ...
}
```

The identifiers should be linked to the prefix `https://mibig.secondarymetabolites.org/repository/` 

Instead of hardcoding MIBiG website URLs, this pipeline uses Bioregistry IRIs:
`https://bioregistry.io/mibig:BGC0002906`

Benefits:

	•	Stable, canonical identifier
	•	Automatically redirects to the current MIBiG landing page
	•	Future-proof if MIBiG changes its URL structure again

## plantiSMASH URI mapping to rdf 

The mapping files for plantiSMASH are available here: 

in the input directory 
```
wget https://raw.githubusercontent.com/plantismash/plantismash-database/main/data/plantismash_v2_clusters_minimal.json
```

The identifiers should be linked to the prefix: 
`https://plantismash.bioinformatics.nl/precalc/v2/<cluster-id-from-json>` 

```
{
  "<cluster-id>": ["gene1", "gene2", ...]
}
```

# Workflow 

### 1. Direct MIBiG gene matching 

If a TAIR gene ID from the pathway appears directly in the MIBiG gene list:
`AT5G48000 ∈ MIBiG(BGC0000670)`

Then the following RDF is created:
```
<https://bioregistry.io/mibig:BGC0000670>
  ro:0000051 <https://identifiers.org/tair.locus/AT5G48000> .
```

This is fast, robust, and does not require BridgeDb

### 2. BridgeDb mapping (required for plantiSMASH)

plantiSMASH clusters often use:

	•	locus tags
	•	alternative gene identifiers
	•	non-TAIR identifiers

To resolve these, the pipeline uses the BridgeDb webservice:

`TAIR locus → BridgeDb → mapped identifiers → cluster gene lists`

Example BridgeDb call:

`https://webservice.bridgedb.org/Arabidopsis%20thaliana/xrefs/A/AT5G48000`

BridgeDb allows standardized identifier mapping without maintaining custom crosswalk tables.

BridgeDb has a REST webservice at https://webservice.bridgedb.org/ with OpenAPI docs.

Tutorial examples show the pattern:
/{Organism}/xrefs/{systemCode}/{identifier} 

So for Arabidopsis you’ll do something like (illustrative):
	•	https://webservice.bridgedb.org/Arabidopsis%20thaliana/xrefs/A/AT5G48000

Where:
	•	Arabidopsis thaliana is the organism name the service recognizes
	•	A is the TAIR locus system code (you already have this in your datasource row)

To find the target system codes you need (e.g., RefSeq protein, GenBank protein, UniProt, etc.), check the BridgeDb sdatasource registry: https://github.com/bridgedb/datasources/blob/main/datasources.tsv


### 3. Run python script to make the links 

#### 3.1 Installation

This project uses Python and requires a small set of RDF-related libraries.

Create a Conda environment We recommend using Conda (or Mamba) to ensure reproducibility.

```bash
conda env create -f environment.yml
conda activate map-to-rdf

``` 

#### 3.2 (Optional) Download the pathway RDF for SPARQL testing

The PlantMetWiki pathway RDF bundle is permanently archived on Zenodo:
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17967619-blue)](https://doi.org/10.5281/zenodo.17967619)

```bash
python scripts/download_pathways.py
# Downloads all_pathways.ttl, reactions.ttl, all.ttl into input_ttl/
```

#### 3.3 Run the script

```bash
# Full run (reads from input/)
python scripts/convert_bgc_to_rdf.py

# Quick test with bundled sample data
python scripts/convert_bgc_to_rdf.py --test

# One source only
python scripts/convert_bgc_to_rdf.py --source plantismash

# Offline / CI (skip BridgeDb HTTP calls)
python scripts/convert_bgc_to_rdf.py --no-bridgedb
```

#### 3.4 Validate the generated RDF (run locally before committing)

Always run the validation script before committing changes to the TTL files:

```bash
python scripts/validate_bgc_rdf.py
```

This checks 8 conditions across all output artefacts:

| Check | What it verifies |
|---|---|
| TTL syntax | `plantismash.ttl` and `mibig.ttl` parse without errors |
| Triple count | ≥ 1,000 triples (plantiSMASH) and ≥ 100 (MIBiG) |
| Cluster types | All clusters typed as `pmw:BiosyntheticGeneCluster` |
| Provenance | All clusters have `dcterms:source` |
| Gene links | `RO:0000051` (has_part) triples are present |
| Species | At least one `wp:organism` annotation |
| VoID | `void-bgc.ttl` parses and declares both `void:Dataset` URIs |
| Link table | `bgc_pathway_links.tsv` exists with all required columns |

The same validation runs automatically on GitHub via the **Validate BGC RDF** Actions workflow:
- On every push to `main` / `dev_edp` (path-filtered)
- On every pull request to `main`
- On every release tag (`bgc-*`, `v*`) — also diffs the VoID against a fresh generation

#### 3.5 Check the summaries generated

[/summaries](./summaries/)


#### 3.5 Optional: run SPARQL queries and BridgeDb-assisted pathway→BGC linking

Run basic BGC queries, or include the pathway graph for BridgeDb-assisted linking:

```bash
# BGC queries only
python scripts/test_queries.py

# BGC + pathway graph + BridgeDb pathway→BGC matching
python scripts/test_queries.py --pathways input_ttl/all_pathways.ttl --max-genes 50
```

Cheat sheet using the OBO Relations Ontology model: 
```
PREFIX ro: <http://purl.obolibrary.org/obo/RO_>

# cluster → genes
?cluster ro:0000051 ?gene .

# gene → cluster
?gene ro:0000050 ?cluster .
``` 


## LICENSE 

This repository modifies data from Plant Metabolic Network (PMN), PlantCyc and AraCyc repository. 
It includes modified data derived from the PlantCyc / PMN databases. Original data © Carnegie Institution for Science. Licensed under the PMN Open Database License:  https://plantcyc.org/?webform=license-agreement

The License for the code contained in this repository is available at [LICENSE](./LICENSE). 


## Execute the script with Rstudio (optional for developers preferring this)

We've tested the Python script in Rstudio (version 2025.09.2+418 "Cucumberleaf Sunflower" Release (12f6d5e22720bd78dbd926bb344efe12d0dce83d, 2025-10-20) for windows).
If you want to run this code locally, you can use your own favorite GUI, or follow these steps in Rstudio:

Open a Terminal (can be done inside Rstudio:, next to the Console button in the bottom left window)
```bash
python --version
python.exe -m pip install --upgrade pip
python -m pip install rdflib

pip install requests
```