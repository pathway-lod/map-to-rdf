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

import sys
from pathlib import Path

import requests

CONCEPT_DOI   = "10.5281/zenodo.19928985"
ZENODO_API    = "https://zenodo.org/api/records"
OUTPUT_DIR    = Path("input_ttl")

# Files to download from the record (subset — add more if needed)
WANTED_FILES = {"all_pathways.ttl", "reactions.ttl", "all.ttl"}


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

    to_download = [f for f in files if f["key"] in WANTED_FILES]
    if not to_download:
        print(f"[WARN] None of {WANTED_FILES} found in this record.")
        return 1

    print(f"\nDownloading {len(to_download)} file(s) to {OUTPUT_DIR}/")
    for f in to_download:
        dest = OUTPUT_DIR / f["key"]
        if dest.exists():
            print(f"  [SKIP] {f['key']} already exists")
            continue
        download_file(f["links"]["self"], dest)

    print(f"\nDone. Run queries with:")
    print(f"  python scripts/test_queries.py --pathways {OUTPUT_DIR}/all_pathways.ttl")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
