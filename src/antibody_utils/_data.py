"""Load the data files shipped in `antibody_utils/data`."""

from __future__ import annotations

import json
import tomllib
from functools import cache
from importlib.resources import files
from typing import Any


def _read_bytes(name: str) -> bytes:
    return files("antibody_utils").joinpath("data", name).read_bytes()


@cache
def load_toml(name: str) -> dict[str, Any]:
    """Load a TOML data file.

    Args:
        name: The file name, relative to `antibody_utils/data`.

    Returns:
        The parsed document.  It is cached, so callers must not modify it.
    """
    return tomllib.loads(_read_bytes(name).decode())


@cache
def load_json(name: str) -> Any:
    """Load a JSON data file.

    Args:
        name: The file name, relative to `antibody_utils/data`.

    Returns:
        The parsed document.  It is cached, so callers must not modify it.
    """
    return json.loads(_read_bytes(name))


@cache
def load_blosum62() -> dict[tuple[str, str], int]:
    """Load NCBI's BLOSUM62 substitution matrix.

    Returns:
        The score for each ordered pair of residue codes, including the
        ambiguity codes `B`, `Z` and `X` and the stop `*`.

    Raises:
        FileNotFoundError: If the matrix has not been fetched, which can only
            happen in a source checkout.
    """
    try:
        text = _read_bytes("BLOSUM62").decode()
    except FileNotFoundError:
        raise FileNotFoundError(
            "The BLOSUM62 matrix is missing.  In a source checkout, fetch it "
            "with: uv run --no-project python scripts/fetch_blosum62.py"
        ) from None
    rows = [line.split() for line in text.splitlines() if not line.startswith("#")]
    header, *body = filter(None, rows)
    return {
        (row[0], column): int(score)
        for row in body
        for column, score in zip(header, row[1:], strict=True)
    }
