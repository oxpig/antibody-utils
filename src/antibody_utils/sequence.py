"""Extract regions from numbered antibody sequences, and compare them.

Comparisons work position by position on numbered sequences, so they need no
alignment: residues at the same position in the same numbering scheme are
compared.  Only like chains are compared (heavy with heavy, light with light).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from antibody_utils._data import load_blosum62
from antibody_utils.numbering.anarcii import NumberedSequence
from antibody_utils.numbering.positions import Chain, Position, Scheme
from antibody_utils.regions import Definition, Region, RegionSelector

if TYPE_CHECKING:
    from collections.abc import Iterable

__all__ = [
    "cdr_lengths",
    "cdr_sequences",
    "identity",
    "region_sequences",
    "similarity",
]


def region_sequences(
    numbered: NumberedSequence, definition: Definition | str
) -> dict[Region, str]:
    """Split a numbered domain into its framework and CDR sequences.

    Args:
        numbered: The numbered domain.
        definition: The region definition to apply.

    Returns:
        The sequence of each region present, in sequence order.  Residues
        beyond the numbered variable domain are left out.

    >>> from antibody_utils.numbering import NumberedSequence, Position
    >>> numbered = NumberedSequence(
    ...     name="example",
    ...     chain_type="H",
    ...     scheme="imgt",
    ...     positions=tuple(map(Position.parse, ["1", "2", "27", "28", "39"])),
    ...     sequence="EVGFW",
    ...     score=30.0,
    ...     start=0,
    ...     end=5,
    ... )
    >>> {
    ...     str(region): residues
    ...     for region, residues in region_sequences(numbered, "imgt").items()
    ... }
    {'fwh1': 'EV', 'cdrh1': 'GF', 'fwh2': 'W'}
    """
    sequences: dict[Region, str] = {}
    for _, residue, region in numbered.annotate_regions(definition):
        if region is not None:
            sequences[region] = sequences.get(region, "") + residue
    return sequences


def cdr_sequences(
    numbered: NumberedSequence, definition: Definition | str
) -> dict[Region, str]:
    """Extract the three CDR sequences of a numbered domain.

    Args:
        numbered: The numbered domain.
        definition: The CDR definition to apply.

    Returns:
        CDRs 1, 2 and 3 of the domain's chain, in order.  A CDR with no
        residues (for example, in a truncated sequence) maps to `""`.
    """
    regions = region_sequences(numbered, definition)
    c = numbered.chain.lower()
    return {
        cdr: regions.get(cdr, "") for cdr in (Region(f"cdr{c}{n}") for n in (1, 2, 3))
    }


def cdr_lengths(
    numbered: NumberedSequence, definition: Definition | str
) -> dict[Region, int]:
    """Measure the three CDRs of a numbered domain.

    Args:
        numbered: The numbered domain.
        definition: The CDR definition to apply.

    Returns:
        The number of residues in each of CDRs 1, 2 and 3, in order.
    """
    return {
        cdr: len(sequence)
        for cdr, sequence in cdr_sequences(numbered, definition).items()
    }


def _by_chain(
    domains: NumberedSequence | Iterable[NumberedSequence],
) -> dict[Chain, NumberedSequence]:
    # A single NumberedSequence is itself iterable (over its residues).
    if isinstance(domains, NumberedSequence):
        domains = [domains]
    by_chain: dict[Chain, NumberedSequence] = {}
    for domain in domains:
        if domain.chain in by_chain:
            raise ValueError(f"More than one {domain.chain.name.lower()} chain given")
        by_chain[domain.chain] = domain
    return by_chain


def _paired_residues(
    a: NumberedSequence | Iterable[NumberedSequence],
    b: NumberedSequence | Iterable[NumberedSequence],
    definition: Definition | str,
    regions: Iterable[Region | str] | None,
) -> list[tuple[str, str]]:
    """The residue pairs at the selected positions present in both sequences."""
    a_chains, b_chains = _by_chain(a), _by_chain(b)
    common = [chain for chain in Chain if chain in a_chains and chain in b_chains]
    schemes = {
        Scheme(domain.scheme) for domain in (*a_chains.values(), *b_chains.values())
    }
    if len(schemes) > 1:
        raise ValueError("Sequences numbered in different schemes cannot be compared")
    if not common:
        return []
    selector = RegionSelector(
        ["fv"] if regions is None else regions,
        scheme=schemes.pop(),
        definition=definition,
    )
    pairs = []
    for chain in common:
        b_residues: dict[Position, str] = dict(b_chains[chain])
        for position, residue in a_chains[chain]:
            if position in b_residues and selector.accepts(position, chain):
                pairs.append((residue, b_residues[position]))
    return pairs


def identity(
    a: NumberedSequence | Iterable[NumberedSequence],
    b: NumberedSequence | Iterable[NumberedSequence],
    *,
    definition: Definition | str,
    regions: Iterable[Region | str] | None = None,
) -> float:
    """Calculate the sequence identity between numbered sequences.

    Each of `a` and `b` is one domain, or several (such as a VH and its VL).
    Like chains are compared, and the identity is pooled over all the positions
    compared.  A position counts only if both sequences have a residue there.

    Region membership is decided position by position (as by
    `antibody_utils.regions.get_region`), so that both sequences are compared
    over the same positions.

    Args:
        a: The first sequence, or sequences.
        b: The second sequence, or sequences.
        definition: The region definition used to select positions.
        regions: The regions, or region groups such as `"cdrs"`, to compare
            over.  By default, the whole variable domain (`"fv"`).

    Returns:
        The fraction of compared positions with the same residue, or 0.0 if no
        positions are compared.

    Raises:
        ValueError: If the sequences are numbered in different schemes, or if
            `a` or `b` contains two domains of the same chain.

    >>> from antibody_utils.numbering import NumberedSequence, Position
    >>> def heavy(sequence):
    ...     return NumberedSequence(
    ...         name=sequence,
    ...         chain_type="H",
    ...         scheme="imgt",
    ...         positions=tuple(Position(n) for n in (1, 2, 27, 28)),
    ...         sequence=sequence,
    ...         score=30.0,
    ...         start=0,
    ...         end=4,
    ...     )
    >>> identity(heavy("EVGF"), heavy("QVGF"), definition="imgt")
    0.75
    >>> identity(heavy("EVGF"), heavy("QVGF"), definition="imgt", regions=["cdrs"])
    1.0
    >>> identity(heavy("EVGF"), heavy("QVGF"), definition="imgt", regions=["fwh1"])
    0.5
    """
    pairs = _paired_residues(a, b, definition, regions)
    if not pairs:
        return 0.0
    return sum(x == y for x, y in pairs) / len(pairs)


def similarity(
    a: NumberedSequence | Iterable[NumberedSequence],
    b: NumberedSequence | Iterable[NumberedSequence],
    *,
    definition: Definition | str,
    regions: Iterable[Region | str] | None = None,
    normalise: bool = False,
) -> float:
    """Calculate the BLOSUM62 similarity between numbered sequences.

    Positions are selected as for `identity`.  Residues with no BLOSUM62 score
    (anything other than the 20 amino acids, `B`, `Z`, `X` and `*`) are
    skipped.

    Args:
        a: The first sequence, or sequences.
        b: The second sequence, or sequences.
        definition: The region definition used to select positions.
        regions: The regions, or region groups such as `"cdrs"`, to compare
            over.  By default, the whole variable domain (`"fv"`).
        normalise: Divide the total score by the number of positions compared.

    Returns:
        The total BLOSUM62 score over the compared positions, or its mean if
        `normalise` is set; 0.0 if no positions are compared.

    Raises:
        ValueError: If the sequences are numbered in different schemes, or if
            `a` or `b` contains two domains of the same chain.
    """
    blosum62 = load_blosum62()
    scores = [
        blosum62[pair]
        for pair in _paired_residues(a, b, definition, regions)
        if pair in blosum62
    ]
    if not scores:
        return 0.0
    return sum(scores) / len(scores) if normalise else float(sum(scores))
