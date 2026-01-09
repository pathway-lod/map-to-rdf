# map-to-rdf

Mapping biosynthetic gene clusters (BGCs) to PlantMetWiki / WikiPathways RDF

Visit the PlantMetWiki webserver at: https://plantmetwiki.bioinformatics.nl/

This repository links biosynthetic gene clusters (BGCs) from MIBiG and plantiSMASH to plant metabolic pathway genes represented in WikiPathways RDF, producing interoperable RDF crosslinks that can be queried via SPARQL.

The resulting RDF enables questions such as:

	•	Which genes in this pathway are part of a known BGC?
	•	Which MIBiG or plantiSMASH clusters overlap with plant metabolic pathways?
	•	Which pathways contain genes from a specific BGC (e.g. thalianol / arabidiol)?

## Overview of the approach 

Overview of the approach

The pipeline:

	1.	Reads pathway RDF (WikiPathways GPML → RDF)
	2.	Loads BGC gene membership data from:
        •	MIBiG (gene-based, curated)
        •	plantiSMASH (predicted clusters)
	3.	Matches pathway genes to BGC genes
        •	Directly (when identifiers match)
        •	Via BridgeDb identifier mapping
	4.	Creates RDF triples where BGCs are first-class nodes
	5.	Uses standard ontologies and stable identifiers
	6.	Merges the new crosslinks back into the pathway RDF

## Identifier and ontology design choices

## Cluster-centric modeling
BGCs are represented as explicit RDF nodes, and genes are linked as parts of clusters:

```
BGC  ── RO:has_part (RO:0000051) ──▶ Gene
Gene ── RO:part_of (RO:0000050) ──▶ BGC
```
Ontologies used

| Purpose  | Ontology / Term  | 
|---|---|
|  Gene–cluster relation |  RO:0000051 (has_part) |  
| Reverse relation  | RO:0000050 (part_of)  |   
| Cluster type |  pmw:BiosyntheticGeneCluster |   
| Provenance |  dcterms:source | 

The model is 
```
BGC
 ├─ has_part → gene
 ├─ rdf:type → pmw:BiosyntheticGeneCluster
 └─ dcterms:source → "MIBIG" / "plantiSMASH"
```

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

#### 3.2 Run the script

```
python scripts/create_links.py
``` 

#### 3.3 Check the summaries generated 

[/summaries](./summaries/)

and the output rdf including the crosslinks 

[./output_ttl/reactions_with_bgc_links.ttl](./output_ttl/reactions_with_bgc_links.ttl)

#### 3.4 Optional: check sample SPARQL queries 

Adapt the script by pasting your test sample queries and run: 

```
python scripts/test_queries.py
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