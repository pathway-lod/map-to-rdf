from rdflib import Graph, Namespace

PMW = Namespace("https://plantmetwiki.bioinformatics.nl/vocab/")
WP  = Namespace("http://vocabularies.wikipathways.org/wp#")
DCT = Namespace("http://purl.org/dc/terms/")

g = Graph()
g.parse("reactions.ttl", format="turtle")
g.parse("bgc_links.ttl", format="turtle")  # or reactions_with_bgc_links.ttl

print("Triples in graph:", len(g))

def run(q, title):
    print("\n" + "="*80)
    print(title)
    print("="*80)
    res = g.query(q, initNs={"pmw": PMW, "wp": WP, "dcterms": DCT})
    for row in res:
        print(row)

# --- Put queries below ---

# list all gene → BGC links
run("""
SELECT ?gene ?bgc
WHERE { ?gene pmw:belongsToBGC ?bgc . }
ORDER BY ?gene ?bgc
LIMIT 50
""", "A) Gene → BGC links")

# count total links 
run("""
SELECT (COUNT(*) AS ?nLinks)
WHERE { ?gene pmw:belongsToBGC ?bgc . }
""", "B) Total link count")

# count links by source 
run("""
SELECT ?source (COUNT(*) AS ?nLinks)
WHERE {
  ?gene pmw:belongsToBGC ?bgc .
  ?bgc dcterms:source ?source .
}
GROUP BY ?source
ORDER BY DESC(?nLinks)
""", "C) Link count by source")