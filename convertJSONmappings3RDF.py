import json
import re
from rdflib import Graph, Namespace, URIRef, Literal

with open("test.json") as f:
    data = json.load(f)

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

arabidopsis = Namespace("https://www.arabidopsis.org/results?mainType=general&category=genes&searchText=")
g.bind("at", arabidopsis)

# Iterate over clusters and genes
for cluster_id, gene_list in data.items():
  
  # Only keep clusters matching the regex for arabidopsis
    if not pattern.search(cluster_id):
        continue
      
    cluster_uri = URIRef(plantismash[cluster_id])

    for gene_id in gene_list:
        gene_uri = URIRef(arabidopsis[gene_id])

        # cluster has_Part gene
        g.add((cluster_uri, RO.has_part, gene_uri))

# Serialize
g.serialize("output.ttl", format="turtle")


##TODO:
  ##Other species are a combination of resources. GPML example has ASAT1 as a label (for tomato, Solanum Lycopersicum), the xRef is G18C3-27. ASAT1 is in PlantisMash (but not an official ID), while the PlantCyc URI is https://pmn.plantcyc.org/gene?orgid=PLANT&id=G18C3-27 . This will cause a mismatch between data from PlantisMash and the PW data we have. MIBIG does not have this ID to begin with.
 ## CUrrent IRIs for dc_identifiers of TAIR not working, use foaf:page to link to working URI
