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
