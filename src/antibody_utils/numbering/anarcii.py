"""Number antibody sequences with ANARCII.

This module needs the `numbering` extra: `pip install "antibody-utils[numbering]"`.
ANARCII (and with it PyTorch) is imported only when sequences are first
numbered, so importing `antibody_utils` stays quick without it.

For large inputs, use `iter_number`, which numbers sequences in batches and
yields the results as it goes, so that memory use depends on the batch size
rather than on the number of sequences.
"""

from __future__ import annotations

import contextlib
import io
import logging
import re
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from itertools import islice
from typing import TYPE_CHECKING, Any

from antibody_utils.numbering.positions import Chain, Position, Scheme

if TYPE_CHECKING:
    from anarcii import Anarcii

    from antibody_utils.regions import Definition, Region

__all__ = [
    "NumberedSequence",
    "NumberingError",
    "iter_number",
    "number",
    "number_sequence",
]

logger = logging.getLogger(__name__)

# Only letters: ANARCII would treat `-`, `/` and `\` as separators between
# paired chains, and split the sequence.
_SEQUENCE = re.compile(r"[A-Za-z]+")

# Sequences per call to ANARCII.  Each call reloads the model weights (about a
# second), and ANARCII holds about 35 KiB per sequence until we convert its
# results, so this keeps both overheads small.  ANARCII writes its results to a
# file instead of returning them above 102,400 sequences per call.
DEFAULT_BATCH_SIZE = 10_000
_MAX_BATCH_SIZE = 100_000


class NumberingError(ValueError):
    """A sequence could not be numbered as an antibody variable domain."""


@dataclass(frozen=True, slots=True)
class NumberedSequence:
    """An antibody variable domain, numbered in one scheme.

    Attributes:
        name: The name of the input sequence.
        chain_type: The chain type that ANARCII assigned: `"H"` (heavy), `"K"`
            (kappa) or `"L"` (lambda).
        scheme: The numbering scheme.
        positions: The position of each residue of the numbered domain, in
            sequence order.  Positions with no residue are omitted.
        sequence: The residues of the numbered domain.
        score: ANARCII's confidence score for the numbering.
        start: The index of the domain's first residue in the input sequence.
        end: The index just past the domain's last residue in the input
            sequence, so that `input[start:end]` is the numbered domain.

    Raises:
        ValueError: If `positions` and `sequence` differ in length.
    """

    name: str
    chain_type: str
    scheme: Scheme
    positions: tuple[Position, ...]
    sequence: str
    score: float
    start: int
    end: int

    def __post_init__(self) -> None:
        """Check that there is one residue per position."""
        if len(self.positions) != len(self.sequence):
            raise ValueError(
                f"{len(self.positions)} positions but {len(self.sequence)} residues"
            )

    @property
    def chain(self) -> Chain:
        """The chain: heavy, or light for kappa and lambda chains."""
        return Chain.HEAVY if self.chain_type == "H" else Chain.LIGHT

    @property
    def numbering(self) -> list[tuple[Position, str]]:
        """The `(position, residue)` pairs, in sequence order."""
        return list(self)

    def __iter__(self) -> Iterator[tuple[Position, str]]:
        """Iterate over the `(position, residue)` pairs."""
        return zip(self.positions, self.sequence, strict=True)

    def __len__(self) -> int:
        """The number of numbered residues."""
        return len(self.sequence)

    def annotate_regions(
        self, definition: Definition | str
    ) -> list[tuple[Position, str, Region | None]]:
        """Assign each residue to its region.

        Args:
            definition: The region definition to apply.

        Returns:
            A `(position, residue, region)` triple for each residue; see
            `antibody_utils.regions.annotate_regions`.
        """
        # Imported here: `antibody_utils.regions` imports this package.
        from antibody_utils.regions import annotate_regions

        return annotate_regions(
            self, self.chain, scheme=self.scheme, definition=definition
        )


def _import_anarcii() -> type[Anarcii]:
    try:
        from anarcii import Anarcii
    except ImportError as error:
        raise ImportError(
            "Numbering sequences needs ANARCII, which is not installed.  "
            'Install the numbering extra: pip install "antibody-utils[numbering]"'
        ) from error
    return Anarcii


@cache
def _model(cpu: bool) -> Anarcii:
    """Load ANARCII's antibody model once per device choice."""
    Anarcii = _import_anarcii()
    with contextlib.redirect_stdout(io.StringIO()):
        return Anarcii(seq_type="antibody", mode="accuracy", cpu=cpu, verbose=False)


def _run(model: Anarcii, sequences: dict[str, str], scheme: Scheme) -> dict[str, Any]:
    # ANARCII prints progress messages, even when not verbose.
    with contextlib.redirect_stdout(io.StringIO()):
        results = model.number(sequences)
        if scheme is not Scheme.IMGT:
            results = model.to_scheme(scheme.value)
    return results


def _position(number: int | str, insertion: str) -> Position:
    """Share one `Position` object between all sequences with that position.

    The key is normalised first, so that, for example, `(111, " ")` and
    `("111", "")` give the same object.

    This is safe only because `Position` and `NumberedSequence` are immutable,
    and `NumberedSequence.positions` is a tuple; see
    `test_shared_positions_cannot_be_modified`.
    """
    return _interned_position(int(number), insertion.strip())


@cache
def _interned_position(number: int, insertion: str) -> Position:
    return Position(number, insertion)


def _convert(name: str, result: dict[str, Any], scheme: Scheme) -> NumberedSequence:
    if result["numbering"] is None:
        raise NumberingError(f"Could not number {name!r}: {result['error']}")
    residues = [pair for pair in result["numbering"] if pair[1] != "-"]
    return NumberedSequence(
        name=name,
        chain_type=result["chain_type"],
        scheme=scheme,
        positions=tuple(_position(*position) for position, _ in residues),
        sequence="".join(residue for _, residue in residues),
        score=float(result["score"]),
        start=result["query_start"],
        end=result["query_end"] + 1,
    )


def _validated(
    sequences: Mapping[str, str] | Iterable[str],
) -> Iterator[tuple[str, str]]:
    if isinstance(sequences, str):
        raise TypeError(
            "Expected a collection of sequences, not a single sequence string; "
            "use number_sequence to number one sequence"
        )
    items: Iterable[tuple[str, str]] = (
        sequences.items()
        if isinstance(sequences, Mapping)
        else ((str(index), sequence) for index, sequence in enumerate(sequences))
    )
    for name, sequence in items:
        if not _SEQUENCE.fullmatch(sequence):
            raise ValueError(
                f"Sequence {name!r} must contain only one-letter amino-acid codes"
            )
        yield name, sequence


def _iter_number(
    sequences: Mapping[str, str] | Iterable[str],
    scheme: Scheme | str,
    cpu: bool,
    batch_size: int,
) -> Iterator[tuple[str, NumberedSequence | NumberingError]]:
    scheme = Scheme(scheme)
    if not 1 <= batch_size <= _MAX_BATCH_SIZE:
        raise ValueError(f"batch_size must be between 1 and {_MAX_BATCH_SIZE:,}")
    items = _validated(sequences)
    model = None
    while batch := dict(islice(items, batch_size)):
        model = model or _model(cpu)
        results = _run(model, batch, scheme)
        for name in batch:
            try:
                yield name, _convert(name, results[name], scheme)
            except NumberingError as error:
                yield name, error


def iter_number(
    sequences: Mapping[str, str] | Iterable[str],
    scheme: Scheme | str = Scheme.IMGT,
    *,
    cpu: bool = False,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> Iterator[tuple[str, NumberedSequence | None]]:
    """Number antibody variable-domain sequences with ANARCII, lazily.

    Sequences are read from `sequences`, and numbered `batch_size` at a time,
    only as the results are consumed.  This suits inputs too large to hold in
    memory, such as a generator over a FASTA file.

    Each sequence should contain one variable domain; any residues before or
    after it (such as a constant domain) are left unnumbered.

    Args:
        sequences: Sequences of one-letter amino-acid codes, as a mapping from
            names to sequences, or as an iterable of sequences, which are
            named by their index (`"0"`, `"1"`, …).
        scheme: The numbering scheme.
        cpu: Run on the CPU even if a GPU is available.
        batch_size: How many sequences to pass to ANARCII at a time.

    Yields:
        A `(name, numbered sequence)` pair for each input sequence, in input
        order.  A sequence that ANARCII could not number gives `None`, and a
        warning is logged.

    Raises:
        TypeError: If `sequences` is a single sequence string.
        ValueError: If a sequence contains anything other than letters.  (In
            particular, ANARCII would treat `-` and `/` as separators between
            paired chains.)  The error is raised when that sequence's batch is
            reached, after the results of earlier batches have been yielded.
        ImportError: If ANARCII is not installed.
    """
    for name, result in _iter_number(sequences, scheme, cpu, batch_size):
        if isinstance(result, NumberingError):
            logger.warning("%s", result)
            yield name, None
        else:
            yield name, result


def number(
    sequences: Mapping[str, str] | Sequence[str],
    scheme: Scheme | str = Scheme.IMGT,
    *,
    cpu: bool = False,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> dict[str, NumberedSequence | None]:
    """Number a collection of antibody variable-domain sequences with ANARCII.

    This returns every result at once, at about 1.2 KiB per sequence, so it
    takes only collections that are already in memory: a mapping or a
    sequence such as a list.  To number a stream of sequences, such as a
    generator over a large FASTA file, use `iter_number`, which yields each
    result as it goes.

    Args:
        sequences: Sequences of one-letter amino-acid codes, as a mapping from
            names to sequences, or as a sequence (such as a list) of them,
            which are named by their index (`"0"`, `"1"`, …).
        scheme: The numbering scheme.
        cpu: Run on the CPU even if a GPU is available.
        batch_size: How many sequences to pass to ANARCII at a time.

    Returns:
        The numbered sequences, by name and in input order.  A sequence that
        ANARCII could not number maps to `None`, and a warning is logged.

    Raises:
        TypeError: If `sequences` is not a mapping or a sequence, or is a
            single sequence string.
        ValueError: If a sequence contains anything other than letters.  (In
            particular, ANARCII would treat `-` and `/` as separators between
            paired chains.)
        ImportError: If ANARCII is not installed.
    """
    if not isinstance(sequences, Mapping | Sequence):
        raise TypeError(
            f"number() takes a mapping or a sequence (such as a list), not "
            f"{type(sequences).__name__}; use iter_number to number a stream of "
            f"sequences without holding every result in memory"
        )
    return dict(iter_number(sequences, scheme, cpu=cpu, batch_size=batch_size))


def number_sequence(
    sequence: str, scheme: Scheme | str = Scheme.IMGT, *, cpu: bool = False
) -> NumberedSequence:
    """Number a single antibody variable-domain sequence with ANARCII.

    Args:
        sequence: The sequence, as one-letter amino-acid codes.
        scheme: The numbering scheme.
        cpu: Run on the CPU even if a GPU is available.

    Returns:
        The numbered sequence, named `"sequence"`.

    Raises:
        NumberingError: If ANARCII could not number the sequence.
        ValueError: If the sequence contains anything other than letters.
        ImportError: If ANARCII is not installed.

    >>> numbered = number_sequence(
    ...     "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSV"
    ...     "KGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
    ...     cpu=True,
    ... )
    >>> numbered.chain_type, len(numbered)
    ('H', 120)
    >>> "".join(
    ...     residue
    ...     for _, residue, region in numbered.annotate_regions("imgt")
    ...     if region == "cdrh3"
    ... )
    'SRWGGDGFYAMDY'
    """
    [(_, result)] = _iter_number({"sequence": sequence}, scheme, cpu, 1)
    if isinstance(result, NumberingError):
        raise result
    return result
