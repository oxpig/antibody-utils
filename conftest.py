"""Collect the interactive examples in the Markdown docs and README as tests.

Docstring examples are collected by pytest's own `--doctest-modules`, but its
plain-text doctest parser can't tell where a Markdown code fence ends, so
Sybil handles the Markdown files instead.

Only ```` ```pycon ```` blocks (interactive sessions with `>>>` prompts) are
run; ```` ```python ```` blocks are illustrative and are not executed.
"""

from collections.abc import Iterable
from doctest import NORMALIZE_WHITESPACE

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


pytest_collect_file = Sybil(
    parsers=[PyconParser(doctest_optionflags=NORMALIZE_WHITESPACE), SkipParser()],
    patterns=["*.md"],
    excludes=["docs/_build/*"],
).pytest()
