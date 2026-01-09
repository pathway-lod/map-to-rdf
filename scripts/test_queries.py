from rdflib import Graph, Namespace

PMW = Namespace("https://plantmetwiki.bioinformatics.nl/vocab/")
WP  = Namespace("http://vocabularies.wikipathways.org/wp#")
DCT = Namespace("http://purl.org/dc/terms/")
RO  = Namespace("http://purl.obolibrary.org/obo/RO_")

g = Graph()
g.parse("input_ttl/all_pathways.ttl", format="turtle")
g.parse("output_ttl/bgc_links.ttl", format="turtle")

print("Triples in graph:", len(g))

def run(q, title):
    print("\n" + "="*80)
    print(title)
    print("="*80)
    res = g.query(q, initNs={"pmw": PMW, "wp": WP, "dcterms": DCT, "ro": RO})
    for row in res:
        print(row)

# --- Put queries below ---

# list all gene → BGC links
run("""
PREFIX ro: <http://purl.obolibrary.org/obo/RO_>
SELECT ?cluster ?gene
WHERE { ?cluster ro:0000051 ?gene . }
ORDER BY ?cluster ?gene
LIMIT 50
""", "A) Gene → BGC links")

# count total links 
run("""
SELECT (COUNT(*) AS ?nLinks)
WHERE { ?cluster ro:0000051 ?gene . }
""", "B) Count Total Links")

# count links by source 
run("""
PREFIX ro: <http://purl.obolibrary.org/obo/RO_>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT ?source (COUNT(*) AS ?nLinks)
WHERE {
  ?cluster ro:0000051 ?gene .
  ?cluster dcterms:source ?source .
}
GROUP BY ?source
ORDER BY DESC(?nLinks)
""", "C) Link count by source")