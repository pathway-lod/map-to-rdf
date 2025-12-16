import json
from pathlib import Path

INPUT_TXT = "input/knownclusters.txt"
OUTPUT_JSON = "input/mibig_genes.json"

mibig = {}

with open(INPUT_TXT, encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")

        # We need at least 5 columns to get gene IDs
        if len(parts) < 5:
            print(f"[SKIP] malformed line: {line}")
            continue

        bgc_full = parts[0]        # e.g. BGC0000670.4
        gene_field = parts[4]      # ← THIS is the gene column

        # Strip version: BGC0000670.4 → BGC0000670
        bgc_id = bgc_full.split(".", 1)[0]

        genes = [
            g.strip()
            for g in gene_field.split(";")
            if g.strip()
        ]

        # Merge safely if BGC appears more than once
        mibig.setdefault(bgc_id, set()).update(genes)

# Convert sets to sorted lists for JSON
mibig = {
    bgc_id: sorted(genes)
    for bgc_id, genes in mibig.items()
}

print(f"Parsed {len(mibig)} MIBIG BGCs")

Path(OUTPUT_JSON).write_text(
    json.dumps(mibig, indent=2),
    encoding="utf-8"
)

print(f"Wrote {OUTPUT_JSON}")