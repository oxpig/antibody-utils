"""Fetch NCBI's BLOSUM62 matrix into the package data directory.

The matrix is not tracked in git.  Run this script once before building or
testing:

    uv run python scripts/fetch_blosum62.py

The download is verified against a pinned SHA-256 checksum, so a build can only
ever ship the exact file published by NCBI.
"""

from __future__ import annotations

import hashlib
import logging
import sys
import urllib.request
from pathlib import Path

URL = "https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62"
SHA256 = "85510d3846ee6d5f4778e425cf8daf6e0dbb889b306f2d13434e1254780efb40"
DESTINATION = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "antibody_utils"
    / "data"
    / "BLOSUM62"
)

logger = logging.getLogger(__name__)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    """Download and verify the matrix, skipping the download if already present.

    Returns:
        The process exit status.
    """
    if DESTINATION.exists() and _sha256(DESTINATION.read_bytes()) == SHA256:
        logger.info("%s is present and verified.", DESTINATION)
        return 0

    with urllib.request.urlopen(URL, timeout=60) as response:
        data = response.read()

    if (digest := _sha256(data)) != SHA256:
        logger.error(
            "Checksum mismatch for %s: expected %s, got %s.", URL, SHA256, digest
        )
        return 1

    DESTINATION.parent.mkdir(parents=True, exist_ok=True)
    DESTINATION.write_bytes(data)
    logger.info("Fetched %s to %s.", URL, DESTINATION)
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    sys.exit(main())
