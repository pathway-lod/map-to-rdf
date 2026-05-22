#!/usr/bin/env python3
"""Download the PlantMetWiki pathway RDF bundle from Zenodo into input_ttl/.

The permanent Zenodo record is: https://doi.org/10.5281/zenodo.19928985
This resolves to the latest deposited version automatically.

Usage
-----
    python scripts/download_pathways.py

Files downloaded to input_ttl/:
    all_pathways.ttl   — all pathway RDF (use with test_queries.py --pathways)
    reactions.ttl      — individual reaction RDF
    all.ttl            — pathways + reactions combined
"""

from __future__ import annotations

import gzip
import shutil
import sys
from pathlib import Path

import requests

CONCEPT_DOI = "10.5281/zenodo.17967619"
ZENODO_API  = "https://zenodo.org/api/records"
OUTPUT_DIR  = Path("input_ttl")

# Filename patterns to download (matched by substring).
# Versioned names like all-plantcyc17.0.0-gpml2021.ttl.gz are matched by prefix.
WANTED_PREFIXES = ("all-", "all_gpml_taxonomy_extra-", "all_gpml_properties_extra-")


def resolve_concept_doi(concept_doi: str) -> dict:
    """Resolve a concept DOI to the latest Zenodo record metadata."""
    # Strip DOI prefix if present
    record_id = concept_doi.replace("10.5281/zenodo.", "")
    url = f"{ZENODO_API}/{record_id}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def download_file(url: str, dest: Path, chunk_size: int = 1 << 20) -> None:
    """Stream-download *url* to *dest*, showing progress."""
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with dest.open("wb") as fh:
            for chunk in r.iter_content(chunk_size=chunk_size):
                fh.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    print(f"  {dest.name}: {pct:.0f}%", end="\r")
    print(f"  ✔ {dest.name} ({downloaded / 1e6:.1f} MB)")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Resolving Zenodo record: https://doi.org/{CONCEPT_DOI}")
    try:
        record = resolve_concept_doi(CONCEPT_DOI)
    except Exception as e:
        print(f"[ERROR] Could not fetch Zenodo metadata: {e}")
        return 1

    version = record.get("metadata", {}).get("version", "unknown")
    title   = record.get("metadata", {}).get("title", "")
    print(f"Record: {title} (version {version})")
    print(f"Files available:")

    files = record.get("files", [])
    for f in files:
        print(f"  {f['key']}  ({f['size'] / 1e6:.1f} MB)")

    to_download = [f for f in files
                   if any(f["key"].startswith(p) for p in WANTED_PREFIXES)]
    if not to_download:
        print(f"[WARN] No files matching prefixes {WANTED_PREFIXES} found.")
        print("Available files:", [f["key"] for f in files])
        return 1

    print(f"\nDownloading {len(to_download)} file(s) to {OUTPUT_DIR}/")
    downloaded_ttls = []
    for f in to_download:
        key  = f["key"]
        dest = OUTPUT_DIR / key

        if dest.exists() or (dest.with_suffix("") if key.endswith(".gz") else dest).exists():
            print(f"  [SKIP] {key} already exists")
            downloaded_ttls.append(dest.with_suffix("") if key.endswith(".gz") else dest)
            continue

        download_file(f["links"]["self"], dest)

        # Decompress .gz files in-place
        if key.endswith(".gz"):
            ttl_dest = dest.with_suffix("")  # strip .gz
            print(f"  Decompressing → {ttl_dest.name} ...", end=" ")
            with gzip.open(dest, "rb") as gz_in, ttl_dest.open("wb") as ttl_out:
                shutil.copyfileobj(gz_in, ttl_out)
            dest.unlink()  # remove .gz
            print("done")
            downloaded_ttls.append(ttl_dest)
        else:
            downloaded_ttls.append(dest)

    print(f"\nFiles in {OUTPUT_DIR}/:")
    for p in sorted(OUTPUT_DIR.glob("*.ttl")):
        print(f"  {p.name}  ({p.stat().st_size / 1e6:.1f} MB)")

    # Find the main combined TTL for the usage hint
    main_ttl = next((p for p in downloaded_ttls if p.name.startswith("all-")), None)
    if main_ttl:
        print(f"\nRun queries with:")
        print(f"  python scripts/test_queries.py --pathways {main_ttl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
