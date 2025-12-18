import json
import re
from rdflib import Graph, Namespace, URIRef, Literal
#import subprocess
import requests

url = "https://raw.githubusercontent.com/plantismash/plantismash-database/main/data/plantismash_v2_clusters_minimal.json"
out = "plantismash.json"

r = requests.get(url)
r.raise_for_status()

data = r.json()

#Simple test dataset to check the script is working as expected:
#with open("test.json") as f:
#    data = json.load(f)

#Namespace MIBIG: https://mibig.secondarymetabolites.org/repository/
#Namespace PlantisMash: https://plantismash.bioinformatics.nl/precalc/v2/Abeliophyllum_distichum_GCA_043235775.1/#cluster-1

#Example for TAIR (Arabidopsis) that should be in a PW: AT5G48010

g = Graph()

##Linking genes to cluster using 'has_part' from relation ontology:
RO = Namespace("http://purl.obolibrary.org/obo/RO_")
g.bind("obo", RO)

##PlantisMash data:

plantismash = Namespace("https://plantismash.bioinformatics.nl/precalc/v2/")
g.bind("ex", plantismash)

##Load data for arabidopsis only
pattern = re.compile(r"Arabidopsis_thaliana")
patternCrossbreeding = re.compile(r"_x_")

arabidopsis = Namespace("https://www.arabidopsis.org/results?mainType=general&category=genes&searchText=")
g.bind("at", arabidopsis)

foafPage = Namespace("http://xmlns.com/foaf/0.1/")
g.bind("foaf", foafPage)

arabidopsisIDs = Namespace("https://identifiers.org/tair.name/")
g.bind("atid", arabidopsisIDs)

# Iterate over clusters and genes
for cluster_id, gene_list in data.items():
  
  # Only keep clusters matching the regex for arabidopsis
    if not pattern.search(cluster_id):
        continue
    if patternCrossbreeding.search(cluster_id):
        continue
      
    cluster_uri = URIRef(plantismash[cluster_id])

    for gene_id in gene_list:
        gene_uri = URIRef(arabidopsis[gene_id])
        gene_uri_page = URIRef(arabidopsisIDs[gene_id])

        # cluster has_Part gene
        g.add((cluster_uri, RO.has_part, gene_uri))
        
        g.add((gene_uri, foafPage.page, gene_uri_page))

# Serialize
g.serialize("output.ttl", format="turtle")


##TODO:
  ## Other species are a combination of resources. GPML example has ASAT1 as a label (for tomato, Solanum Lycopersicum), the xRef is G18C3-27. ASAT1 is in PlantisMash (but not an official ID), while the PlantCyc URI is https://pmn.plantcyc.org/gene?orgid=PLANT&id=G18C3-27 . This will cause a mismatch between data from PlantisMash and the PW data we have. MIBIG does not have this ID to begin with.
  ## Define what a BCG actually is (so we can query how many there are); same for genes in clusters
 
