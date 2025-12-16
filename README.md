# map-to-rdf

Mapping biosynthetic gene cluster (BGC) repositories to RDF files for Plant Wikipathways. 

## MIBIG mapping to rdf 

MIBIG mapping files to enzymes are available on MITE at https://github.com/mite-standard/mite_data/blob/dev/mite_data/mibig/mibig_proteins.json


They can be obtained with: 

```
wget https://raw.githubusercontent.com/mite-standard/mite_data/dev/mite_data/mibig/mibig_proteins.json
```

## plantiSMASH mapping to rdf 

The mapping files for plantiSMASH are available here: 

https://github.com/plantismash/plantismash-database/blob/main/data/plantismash_v2_clusters_minimal.json

## Execute the script:
We've tested the Python script in Rstudio (version 2025.09.2+418 "Cucumberleaf Sunflower" Release (12f6d5e22720bd78dbd926bb344efe12d0dce83d, 2025-10-20) for windows).
If you want to run this code locally, you can use your own favorite GUI, or follow these steps in Rstudio:

Open a Terminal (can be done inside Rstudio:, next to the Console button in the bottom left window)
```bash
python --version
python.exe -m pip install --upgrade pip
python -m pip install rdflib
```