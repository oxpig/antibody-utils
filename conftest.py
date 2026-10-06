"""Collect the interactive examples in the Markdown docs and README as tests.

Docstring examples are collected by pytest's own `--doctest-modules`, but its
plain-text doctest parser can't tell where a Markdown code fence ends, so
Sybil handles the Markdown files instead.

Only ```` ```pycon ```` blocks (interactive sessions with `>>>` prompts) are
run; ```` ```python ```` blocks are illustrative and are not executed.

Each Markdown document runs in a temporary directory holding the example
structure files, and can skip examples that need ANARCII with
`% skip: start if(not HAS_ANARCII, reason="...")` and `% skip: end`.
"""

import os
import tempfile
from collections.abc import Iterable
from doctest import NORMALIZE_WHITESPACE
from importlib.util import find_spec
from pathlib import Path
from typing import Any

from sybil import Document, Region, Sybil
from sybil.evaluators.doctest import DocTestEvaluator
from sybil.parsers.abstract.doctest import DocTestStringParser
from sybil.parsers.myst import CodeBlockParser, SkipParser


class PyconParser:
    """Run each ```` ```pycon ```` code block as a doctest."""

    def __init__(self, doctest_optionflags: int = 0) -> None:
        """Initialise the parser.

        Args:
            doctest_optionflags: `doctest` option flags for every example.
        """
        self.code_blocks = CodeBlockParser(language="pycon")
        self.doctests = DocTestStringParser(DocTestEvaluator(doctest_optionflags))

    def __call__(self, document: Document) -> Iterable[Region]:
        """Find the doctest examples in a document's `pycon` blocks.

        Args:
            document: The Markdown document to parse.

        Yields:
            One region per doctest example, positioned within the document.
        """
        for block in self.code_blocks(document):
            source = block.parsed
            for region in self.doctests(source, document.path):
                region.adjust(block, source)
                yield region


# The ANARCII wrapper's doctests need the `numbering` extra.
collect_ignore = (
    [] if find_spec("anarcii") else ["src/antibody_utils/numbering/anarcii.py"]
)

# Structure files that the documentation's examples read, by the name they use.
EXAMPLE_FILES = {
    "12e8.pdb.gz": "tests/data/pdb/12e8.pdb.gz",
    "pdb_000012e8_H_L_ab.cif.gz": "tests/data/sabdab/pdb_000012e8_H_L_ab.cif.gz",
    "pdb_000012e8_sabdab.cif.gz": "tests/data/sabdab/pdb_000012e8_sabdab.cif.gz",
}


def _setup(namespace: dict[str, Any]) -> None:
    """Run a document's examples in a directory holding the example files."""
    root = Path(__file__).parent
    directory = tempfile.TemporaryDirectory()
    for name, source in EXAMPLE_FILES.items():
        (Path(directory.name) / name).symlink_to(root / source)
    namespace["_directory"] = directory
    namespace["_previous_directory"] = Path.cwd()
    namespace["HAS_ANARCII"] = find_spec("anarcii") is not None
    os.chdir(directory.name)


def _teardown(namespace: dict[str, Any]) -> None:
    os.chdir(namespace["_previous_directory"])
    namespace["_directory"].cleanup()


pytest_collect_file = Sybil(
    parsers=[PyconParser(doctest_optionflags=NORMALIZE_WHITESPACE), SkipParser()],
    patterns=["*.md"],
    excludes=["docs/_build/*"],
    setup=_setup,
    teardown=_teardown,
).pytest()
