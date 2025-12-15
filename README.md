# map-to-rdf

Mapping biosynthetic gene cluster (BGC) repositories to RDF files for Plant Wikipathways using BridgeDB. 

## MIBIG mapping to rdf 

MIBIG mapping files to enzymes are available on MITE at https://github.com/mite-standard/mite_data/blob/dev/mite_data/mibig/mibig_proteins.json


They can be obtained with: 

```
mkdir input 
cd input
wget https://raw.githubusercontent.com/mite-standard/mite_data/main/mite_data/mibig/mibig_proteins.json
# or from plantiSMASH a txt file with all the gene identifiers 
wget https://raw.githubusercontent.com/plantismash/plantismash/refs/heads/master/antismash/generic_modules/knownclusterblast/knownclusters.txt
```

The identifiers should be linked to the prefix `https://mibig.secondarymetabolites.org/repository/` 

## plantiSMASH mapping to rdf 

The mapping files for plantiSMASH are available here: 

in the input directory 
```
wget https://raw.githubusercontent.com/plantismash/plantismash-database/main/data/plantismash_v2_clusters_minimal.json
```

The identifiers should be linked to the prefix: 
`https://plantismash.bioinformatics.nl/precalc/v2/` 

# Workflow 

### 1. Create the BGC URIs 

Create BGC URIs using your desired prefixes:
	•	MIBIG BGC: https://mibig.secondarymetabolites.org/repository/BGC0000001
	•	plantiSMASH cluster: https://plantismash.bioinformatics.nl/precalc/v2/Abeliophyllum_distichum_GCA_043235775.1/#cluster-1

Then add one triple per membership, e.g.:
```
@prefix pmw: <https://plantmetwiki.bioinformatics.nl/vocab/> .
@prefix dcterms: <http://purl.org/dc/terms/> .

<https://identifiers.org/tair.locus/AT5G48000.1>
  pmw:belongsToBGC <https://plantismash.bioinformatics.nl/precalc/v2/Arabidopsis_thaliana.../#cluster-thalianol> . 

```

Where pmw:belongsToBGC is your own predicate (recommended), because dcterms:isPartOf is already heavily used in WP RDF for graph membership, and you don’t want to overload it.

### 2. Using BridgeDb 

Your JSON files are:
	•	MIBIG: BGC → [protein accessions / gene names]
	•	plantiSMASH: cluster → [locus tags / gene IDs]

Your pathway genes are (example): TAIR locus IDs.

So the pipeline is:

TAIR locus → (BridgeDb maps) → protein accession / locus tag → (match against JSON) → BGC/cluster

BridgeDb is ideal here because you avoid writing tons of ad-hoc crosswalk logic yourself. BridgeDb provides standardized identifier mapping via databases + system codes

### 3. Use BridgeDb webservice to map TAIR → other identifier types

BridgeDb has a REST webservice at https://webservice.bridgedb.org/ with OpenAPI docs.

Tutorial examples show the pattern:
/{Organism}/xrefs/{systemCode}/{identifier} 

So for Arabidopsis you’ll do something like (illustrative):
	•	https://webservice.bridgedb.org/Arabidopsis%20thaliana/xrefs/A/AT5G48000

Where:
	•	Arabidopsis thaliana is the organism name the service recognizes
	•	A is the TAIR locus system code (you already have this in your datasource row)

To find the target system codes you need (e.g., RefSeq protein, GenBank protein, UniProt, etc.), check the BridgeDb sdatasource registry: https://github.com/bridgedb/datasources/blob/main/datasources.tsv

### 4. Run python script to make the links 

## Installation

This project uses Python and requires a small set of RDF-related libraries.

### Create a Conda environment

We recommend using Conda (or Mamba) to ensure reproducibility.

```bash
conda env create -f environment.yml
conda activate map-to-rdf
# run the script
python ./scripts/create_links.py
``` 


## LICENSE 

This repository modifies data from Plant Metabolic Network (PMN), PlantCyc and AraCyc repository. 
It includes modified data derived from the PlantCyc / PMN databases. Original data © Carnegie Institution for Science. Licensed under the PMN Open Database License:  https://plantcyc.org/?webform=license-agreement