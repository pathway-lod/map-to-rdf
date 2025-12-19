import json
import re
from rdflib import Graph, Namespace, URIRef, Literal
import requests


###################### Import PlantIsMash Data and Convert to RDF ######################

##Declare if you want to use the test dataset or the full dataset
input = "full" #Options: 'test' or 'full'

if input == "test":
  #Simple test dataset to check the script is working as expected:
  with open("testPlantIsMash.json") as f:
    data = json.load(f)
elif input == "full":
  #Download the full json from GitHub
  url = "https://raw.githubusercontent.com/plantismash/plantismash-database/main/data/plantismash_v2_clusters_minimal.json"
  out = "plantismash.json"

  r = requests.get(url)
  r.raise_for_status()

  data = r.json()
else:
  print("Select 'test' or 'full' for variable 'input' to run this script on the PlantIsMash dataset")
  
#Create the RDF graph
g = Graph()

##Describing data types in this dataset:
##RDF:type (or a) predicate
rdfType = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
g.bind("rdf", rdfType)
##PMW vocabulary (to define the BioSynthetic Gene Clusters)
pmwType = Namespace("https://plantmetwiki.bioinformatics.nl/vocab/")
g.bind("pmwvocab", pmwType)

##WP vocabulary (to define the genes)
wpType = Namespace("http://vocabularies.wikipathways.org/wp#")
g.bind("wp", wpType)

##NCBI prefix for taxon IRIs
ncbiTaxon = Namespace("http://purl.obolibrary.org/obo/NCBITaxon_")
g.bind("ncbi", ncbiTaxon)

##Linking genes to cluster using 'has_part' from relation ontology:
RO = Namespace("http://purl.obolibrary.org/obo/RO_")
g.bind("obo", RO)

#Define namespace for database
plantismash = Namespace("https://plantismash.bioinformatics.nl/precalc/v2/")
g.bind("ex", plantismash)

##Load data for arabidopsis (filter out other species with regex)
pattern_At_speciesName = re.compile(r"Arabidopsis_thaliana")
##Filter out data of crossbreeds to Arabidopsis
patternCrossbreeding = re.compile(r"_x_")

##Define the IRIs for the TAIR database
arabidopsisIDs = Namespace("https://identifiers.org/tair.name/")
g.bind("atid", arabidopsisIDs)

##Define the URIs for finding data in TAIR website
arabidopsis = Namespace("https://www.arabidopsis.org/results?mainType=general&category=genes&searchText=")
g.bind("at", arabidopsis)

##Load data for arabidopsis (filter out other species with regex)
pattern_Sl_speciesName = re.compile(r"Solanum_lycopersicum")

###Other species have a mixture of IDs that need to be matched:

##Load data for solanum lycopersicum (tomato) if LOC (NCBI LOC gene IDs); assuming that all others are "regular" gene names
pattern_Sl_LOC = re.compile(r"^LOC(?P<geneid>\d+)$") 
#### Data for solanum lycopersicum (tomato) if gene name can match to several entries in EnsemblPlants --> 
#### Besides identifiers.org, also link to general webpage (will not resolve for all, working example is ABCG1)
#### https://plants.ensembl.org/Solanum_lycopersicum/Search/Results?species=Solanum_lycopersicum;idx=;q=ABCG1

##Some of the LOCs are not following the correct structure it seems, filtere these out for now.
pattern_Sl_LOC_unmatching = re.compile(r"^LOC\d+_\d+$")

##Define the IRIs for the LOC terms from NCBI
ncbiGene = Namespace("https://identifiers.org/ncbigene:")
g.bind("slidNCBI", ncbiGene)

##Define the IRIs for the EnsemblPlants database
solanumIDsEns = Namespace("http://identifiers.org/ensembl.plant:")
g.bind("slidEns", solanumIDsEns)

##Deine the URIs for finding data in EnsemblPlants website
solanumEns = Namespace("https://plants.ensembl.org/Solanum_lycopersicum/Search/Results?species=Solanum_lycopersicum;idx=;q=")
g.bind("slEns", solanumEns)

##Connect to data in WPRDF
foafPage = Namespace("http://xmlns.com/foaf/0.1/")
g.bind("foaf", foafPage)

##Iterate over clusters and genes in json data for arabidopsis
for cluster_id, gene_list in data.items():
  
    ##Only keep clusters matching the regex for arabidopsis
    if not pattern_At_speciesName.search(cluster_id):
        continue
    ##Ignore crossbreeding clusters
    if patternCrossbreeding.search(cluster_id):
        continue
    
    ##Define which part of the data contains the cluster IDs  
    cluster_uri = URIRef(plantismash[cluster_id])
    
    ##Loop over the gene IDs in each cluster and add these to the graph
    for gene_id in gene_list:
      
        ##Define which part of the data contains the gene IDs  
        gene_uri = URIRef(arabidopsis[gene_id])
        ##Define which part of the data contains the gene URL links  
        gene_uri_page = URIRef(arabidopsisIDs[gene_id])

        ##Add information on data types in this RDF model
        g.add((cluster_uri, rdfType.type, pmwType.BiosyntheticGeneCluster))
        ##Add organism name arabidopsis for cluster to specify better and avoid mismatches
        g.add((cluster_uri, wpType.organismName, Literal("Arabidopsis thaliana")))
        ##Add taxonomy ID for species arabidopsis:
        g.add((cluster_uri, wpType.organism, ncbiTaxon['3702']))#3702
        ##Add to RDF: cluster has_Part geneID
        g.add((cluster_uri, RO.has_part, gene_uri))
        
        ##Add to RDF: geneID foaf:page geneURL
        g.add((gene_uri, foafPage.page, gene_uri_page))
        ##Add type to RDF:
        g.add((gene_uri, rdfType.type, wpType.GeneProduct))
        ##Add organism name arabidopsis for cluster to specify better and avoid mismatches
        g.add((gene_uri, wpType.organismName, Literal("Arabidopsis thaliana")))
        ##Add taxonomy ID for species arabidopsis:
        g.add((gene_uri, wpType.organism, ncbiTaxon['3702']))#3702


##Iterate over clusters and genes in json data for tomato
for cluster_id, gene_list in data.items():
  
    ##Only keep clusters matching the regex for tomato
    if not pattern_Sl_speciesName.search(cluster_id):
        continue
    ##Ignore crossbreeding clusters
    if patternCrossbreeding.search(cluster_id):
        continue
    
    ##Define which part of the data contains the cluster IDs  
    cluster_uri = URIRef(plantismash[cluster_id])
    
    ##Loop over the gene IDs in each cluster and add these to the graph
    for gene_id in gene_list:
      m = pattern_Sl_LOC.fullmatch(gene_id)
      if m:
        gene_numeric_id = m.group("geneid")   # trims 'LOC'
        ##Define which part of the data contains the gene IDs  
        gene_uri = ncbiGene[gene_numeric_id]
          
        ##Define which part of the data contains the gene URL links  
        #gene_uri_page = URIRef(arabidopsisIDs[gene_id])

        ##Add information on data types in this RDF model
        g.add((cluster_uri, rdfType.type, pmwType.BiosyntheticGeneCluster))
        g.add((cluster_uri, wpType.organismName, Literal("Solanum lycopersicum")))
        ##Add taxonomy ID for species tomato:
        g.add((cluster_uri, wpType.organism, ncbiTaxon['4081']))#4081
        ##Add to RDF: cluster has_Part geneID
        g.add((cluster_uri, RO.has_part, gene_uri))
        
        ##Add to RDF: geneID foaf:page geneURL
        #g.add((gene_uri, foafPage.page, gene_uri_page))
        ##Add type to RDF:
        g.add((gene_uri, rdfType.type, wpType.GeneProduct))
        ##Add organism name tomato for cluster to specify better and avoid mismatches
        g.add((gene_uri, wpType.organismName, Literal("Solanum lycopersicum")))
        ##Add taxonomy ID for species tomato:
        g.add((gene_uri, wpType.organism, ncbiTaxon['4081']))#4081
      elif pattern_Sl_LOC_unmatching.search(gene_id):
        continue
      else: ##When there are most likely gene names
        ##Define which part of the data contains the gene IDs  
        gene_uri = URIRef(solanumIDsEns[gene_id])
        ##Define which part of the data contains the gene URL links  
        gene_uri_page = URIRef(solanumEns[gene_id])
          
        ##Add information on data types in this RDF model
        g.add((cluster_uri, rdfType.type, pmwType.BiosyntheticGeneCluster))
        g.add((cluster_uri, wpType.organismName, Literal("Solanum lycopersicum")))
        ##Add taxonomy ID for species tomato:
        g.add((cluster_uri, wpType.organism, ncbiTaxon['4081']))#4081
        ##Add to RDF: cluster has_Part geneID
        g.add((cluster_uri, RO.has_part, gene_uri))
        
        ##Add to RDF: geneID foaf:page geneURL
        g.add((gene_uri, foafPage.page, gene_uri_page))
        ##Add type to RDF:
        g.add((gene_uri, rdfType.type, wpType.GeneProduct))
        ##Add organism name tomato for cluster to specify better and avoid mismatches
        g.add((gene_uri, wpType.organismName, Literal("Solanum lycopersicum")))
        ##Add taxonomy ID for species tomato:
        g.add((gene_uri, wpType.organism, ncbiTaxon['4081']))#4081


# Serialize data into output file
if input == "test":
  #Simple test dataset to check the script is working as expected:
  g.serialize("plantIsMashTestOutput.ttl", format="turtle")
elif input == "full":
  #full json from GitHub
  g.serialize("plantIsMash.ttl", format="turtle")
else:
  print("Select 'test' or 'full' for variable 'input' to run this script on the PlantIsMash dataset")



###################### Import MIBIG Data and Convert to RDF ######################


##Declare if you want to use the test dataset or the full dataset
#input = "test" #Options: 'test' or 'full'

if input == "test":
  #Simple test dataset to check the script is working as expected:
  with open("testMibig.json") as f:
    data = json.load(f)
elif input == "full":
  #Download the full json from GitHub
  url = "https://raw.githubusercontent.com/mite-standard/mite_data/main/mite_data/mibig/mibig_proteins.json"
  out = "mibig.json"

  r = requests.get(url)
  r.raise_for_status()

  data = r.json()
else:
  print("Select 'test' or 'full' for variable 'input' to run this script on the MIBIG dataset")

#Create the RDF graph
g = Graph()

##Describing data types in this dataset:
##RDF:type (or a) predicate
rdfType = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
g.bind("rdf", rdfType)

##PMW vocabulary (to define the BioSynthetic Gene Clusters)
pmwType = Namespace("https://plantmetwiki.bioinformatics.nl/vocab/")
g.bind("pmwvocab", pmwType)

##WP vocabulary (to define the genes)
wpType = Namespace("http://vocabularies.wikipathways.org/wp#")
g.bind("wp", wpType)

##NCBI prefix for taxon IRIs
ncbiTaxon = Namespace("http://purl.obolibrary.org/obo/NCBITaxon_")
g.bind("ncbi", ncbiTaxon)

####Data model relevant predicates:

##Linking genes to cluster using 'has_part' from relation ontology:
RO = Namespace("http://purl.obolibrary.org/obo/RO_")
g.bind("obo", RO)

#Define namespace for database
mibig = Namespace("https://mibig.secondarymetabolites.org/repository/")
g.bind("ex", mibig)

##Load data for arabidopsis only (doesn't match to anything in MIBIG it seems...)
pattern_At = re.compile(r"^AT.G[0-9]{5}([._][0-9]+)?$") #Based on identifiers.org (but rewritten to be more restricted and including isoforms:) ^AT.G[0-9]{5}$

##Load data for solanum lycopersicum (tomato) ^XP_\d+\.\d+$
pattern_Sl = re.compile(r"^XP_\d+\.\d+$") 
#### These XP_IDs or LOC IDs cannot be found in EnsemblPlants, so will not add linkouts via BridgeDb mapping file. 
#### Example: 	XP_004243624.1 (in json) is also registered in MIBIG with LOC101250233
#### The sr.iu.a.u-tokyo database is now linked to the XP_ entries, which provide links to UniProt (A0A3Q7HCW9)
#### UniProt has a SPARQL endpoint, and does contain findable IDs for EnsemblPlants (Solyc07g043680.3.1)

##Define the IRIs for the TAIR database
arabidopsisIDs = Namespace("https://identifiers.org/tair.name/")
g.bind("atid", arabidopsisIDs)

##Deine the URIs for finding data in TAIR website
arabidopsis = Namespace("https://www.arabidopsis.org/results?mainType=general&category=genes&searchText=")
g.bind("at", arabidopsis)

##Define the IRIs for a Tomato hosted at (since there doesn't seem to be anything better at the moment)
tomatoIDs = Namespace("https://sr.iu.a.u-tokyo.ac.jp/db/protein.pl?mibig_accession=BGC0002405&protein_id=")
g.bind("sl", tomatoIDs)

##Connect to data in WPRDF
foafPage = Namespace("http://xmlns.com/foaf/0.1/")
g.bind("foaf", foafPage)

##Iterate over clusters and genes in json data for arabidopsis
for cluster_id, gene_list in data.items():
  
    ##Define which part of the data contains the cluster IDs  
    cluster_uri = URIRef(mibig[cluster_id])
    
    ##Loop over the gene IDs in each cluster and add these to the graph
    for gene_id in gene_list:
      ##Only keep clusters matching the regex for arabidopsis
      if not pattern_At.search(gene_id):
        continue
      
      ##Define which part of the data contains the gene IDs  
      gene_uri = URIRef(arabidopsis[gene_id])
      ##Define which part of the data contains the gene URL links  
      gene_uri_page = URIRef(arabidopsisIDs[gene_id])

      ##Add information on data types in this RDF model
      g.add((cluster_uri, rdfType.type, pmwType.BiosyntheticGeneCluster))
      ##Add organism name arabidopsis for cluster to specify better and avoid mismatches
      g.add((cluster_uri, wpType.organismName, Literal("Arabidopsis thaliana")))
      ##Add taxonomy ID for species arabidopsis:
      g.add((cluster_uri, wpType.organism, ncbiTaxon['3702']))#3702
      ##Add to RDF: cluster has_Part geneID
      g.add((cluster_uri, RO.has_part, gene_uri))
      
      ##Add to RDF: geneID foaf:page geneURL
      g.add((gene_uri, foafPage.page, gene_uri_page))
      ##Add type to RDF:
      g.add((gene_uri, rdfType.type, wpType.GeneProduct))
      ##Add organism name arabidopsis for cluster to specify better and avoid mismatches
      g.add((gene_uri, wpType.organismName, Literal("Arabidopsis thaliana")))
      ##Add taxonomy ID for species arabidopsis:
      g.add((gene_uri, wpType.organism, ncbiTaxon['3702']))#3702
       
##Iterate over clusters and genes in json data for tomato
for cluster_id, gene_list in data.items():
  
        ##Define which part of the data contains the cluster IDs  
    cluster_uri = URIRef(mibig[cluster_id])
    
    ##Loop over the gene IDs in each cluster and add these to the graph
    for gene_id in gene_list:
      ##Only keep genes matching the regex for tomato
      if not pattern_Sl.search(gene_id):
          continue
      
      
      ##Define which part of the data contains the gene IDs  
      gene_uri = URIRef(tomatoIDs[gene_id])
      ##Define which part of the data contains the gene URL links (could be switches to website from Tokyo later)  
      #gene_uri_page = URIRef(tomatoIDs[gene_id])

      ##Add information on data types in this RDF model
      g.add((cluster_uri, rdfType.type, pmwType.BiosyntheticGeneCluster))
      ##Add organism name tomato for cluster to specify better and avoid mismatches
      g.add((cluster_uri, wpType.organismName, Literal("Solanum lycopersicum")))
      ##Add taxonomy ID for species tomato:
      g.add((cluster_uri, wpType.organism, ncbiTaxon['4081']))#4081
      ##Add to RDF: cluster has_Part geneID
      g.add((cluster_uri, RO.has_part, gene_uri))
      
      ##Add to RDF: geneID foaf:page geneURL
      #gr.add((gene_uri, foafPage.page, gene_uri_page))
      ##Add type of datato RDF:
      g.add((gene_uri, rdfType.type, wpType.GeneProduct))
      ##Add organism name tomato for cluster to specify better and avoid mismatches
      g.add((gene_uri, wpType.organismName, Literal("Solanum lycopersicum")))
      ##Add taxonomy ID for species tomato:
      g.add((gene_uri, wpType.organism, ncbiTaxon['4081']))#4081

# Serialize data into output file
if input == "test":
  #Simple test dataset to check the script is working as expected:
  g.serialize("mibigTestOutput.ttl", format="turtle")
elif input == "full":
  #full json from GitHub
  g.serialize("mibig.ttl", format="turtle")
else:
  print("Select 'test' or 'full' for variable 'input' to run this script on the MIBIG dataset")



##TODO:
  ## Other species are a combination of resources. GPML example has ASAT1 as a label (for tomato, Solanum Lycopersicum), the xRef is G18C3-27. ASAT1 is in PlantisMash (but not an official ID), while the PlantCyc URI is https://pmn.plantcyc.org/gene?orgid=PLANT&id=G18C3-27 . This will cause a mismatch between data from PlantisMash and the PW data we have. MIBIG does not have this ID to begin with.
  ## add metadata as separate void header file.
 
