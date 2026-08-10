#!/usr/bin/env python
"""
Fetch and cache external source data for the genera project.

Phase 1 needs only WALS `languages.csv` (Dryer genera, used as calibration and
validation ground truth). ASJP is fetched in phase 2.

Downloads are cached under ``data/sources/`` and skipped if already present
(use ``--force`` to refetch). Network access goes through the environment's
proxy if one is configured; TLS uses a CA bundle from ``SSL_CERT_FILE`` /
``REQUESTS_CA_BUNDLE`` / the CCR bundle when present.
"""

from __future__ import annotations

import argparse
import os
import ssl
import sys
import urllib.request
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
SOURCES = BASE / "data" / "sources"

# Pinned to a specific WALS CLDF commit for reproducibility. Bump deliberately.
WALS_REF = "v2020.1"
SOURCES_SPEC = {
    "wals_languages.csv": (
        f"https://raw.githubusercontent.com/cldf-datasets/wals/{WALS_REF}"
        "/cldf/languages.csv"
    ),
}


def _ssl_context() -> ssl.SSLContext:
    """Build an SSL context, honoring an explicit CA bundle if one is set."""
    for var in ("SSL_CERT_FILE", "REQUESTS_CA_BUNDLE"):
        path = os.environ.get(var)
        if path and Path(path).exists():
            return ssl.create_default_context(cafile=path)
    ccr = Path("/root/.ccr/ca-bundle.crt")
    if ccr.exists():
        return ssl.create_default_context(cafile=str(ccr))
    return ssl.create_default_context()


def fetch(url: str, dest: Path, force: bool = False) -> None:
    if dest.exists() and not force:
        print(f"  cached  {dest.name} ({dest.stat().st_size:,} bytes)")
        return
    print(f"  fetch   {dest.name}  <-  {url}")
    # urllib honors http(s)_proxy env vars via the default opener.
    req = urllib.request.Request(url, headers={"User-Agent": "genera-fetch/0.1"})
    with urllib.request.urlopen(req, context=_ssl_context(), timeout=60) as resp:
        data = resp.read()
    dest.write_bytes(data)
    print(f"          wrote {len(data):,} bytes")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="refetch even if cached")
    args = ap.parse_args()

    SOURCES.mkdir(parents=True, exist_ok=True)
    print(f"Caching sources into {SOURCES}")
    for name, url in SOURCES_SPEC.items():
        try:
            fetch(url, SOURCES / name, force=args.force)
        except Exception as exc:  # noqa: BLE001 - report and continue
            print(f"  ERROR   {name}: {exc}", file=sys.stderr)
            return 1
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
