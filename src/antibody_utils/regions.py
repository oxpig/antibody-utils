"""Framework and CDR regions of antibody variable domains.

A *numbering scheme* (IMGT, Kabat, Chothia or Martin) assigns a number to each
residue of a variable domain.  A *region definition* (IMGT, Kabat, Chothia,
Contact or North) says which positions belong to each complementarity-
determining region (CDR) and framework region (FW).  The two are independent:
any definition can be applied to sequences numbered in any scheme, so long as
positions can be translated between schemes.

Each definition is stored as boundaries on IMGT positions, and positions in the
other schemes are mapped onto IMGT positions.  That mapping is exact except
where the schemes place insertions differently around CDR-H1 and CDR-L2, where
`annotate_regions` uses heuristics to reproduce the CDR that the definition
gives in its own scheme.  `get_region` looks up a single position without those
heuristics, so prefer `annotate_regions` when the definition and the scheme
differ.
"""

from __future__ import annotations

from collections.abc import Iterable
from enum import StrEnum
from functools import cache
from itertools import islice
from string import ascii_uppercase

from antibody_utils._data import load_json, load_toml
from antibody_utils.numbering.positions import (
    Chain,
    Position,
    Scheme,
    _CaseInsensitiveStrEnum,
)

__all__ = [
    "REGION_GROUPS",
    "Definition",
    "Region",
    "RegionSelector",
    "annotate_regions",
    "get_region",
]


class Definition(_CaseInsensitiveStrEnum):
    """A definition of the CDR and framework regions.

    The North definition is that of North B, Lehmann A and Dunbrack RL Jr,
    `"A new clustering of antibody CDR loop conformations"
    <https://doi.org/10.1016/j.jmb.2010.10.030>`__, *J Mol Biol* 406:228–256
    (2011).

    >>> Definition("North")
    <Definition.NORTH: 'north'>
    """

    IMGT = "imgt"
    KABAT = "kabat"
    CHOTHIA = "chothia"
    CONTACT = "contact"
    NORTH = "north"


class Region(StrEnum):
    """A framework or CDR region of a heavy or light variable domain.

    >>> Region.CDRH3 == "cdrh3"
    True
    >>> Region.CDRH3.chain, Region.CDRH3.is_cdr
    (<Chain.HEAVY: 'H'>, True)
    """

    FWH1 = "fwh1"
    CDRH1 = "cdrh1"
    FWH2 = "fwh2"
    CDRH2 = "cdrh2"
    FWH3 = "fwh3"
    CDRH3 = "cdrh3"
    FWH4 = "fwh4"
    FWL1 = "fwl1"
    CDRL1 = "cdrl1"
    FWL2 = "fwl2"
    CDRL2 = "cdrl2"
    FWL3 = "fwl3"
    CDRL3 = "cdrl3"
    FWL4 = "fwl4"

    @property
    def chain(self) -> Chain:
        """The chain that the region belongs to."""
        return Chain(self.value[-2])

    @property
    def is_cdr(self) -> bool:
        """Whether the region is a CDR, rather than a framework region."""
        return self.value.startswith("cdr")


def _group(*names: str) -> frozenset[Region]:
    return frozenset(Region(name) for name in names)


#: Named groups of regions, accepted wherever region names are.
REGION_GROUPS: dict[str, frozenset[Region]] = {
    "hframework": _group("fwh1", "fwh2", "fwh3", "fwh4"),
    "hcdrs": _group("cdrh1", "cdrh2", "cdrh3"),
    "lframework": _group("fwl1", "fwl2", "fwl3", "fwl4"),
    "lcdrs": _group("cdrl1", "cdrl2", "cdrl3"),
}
REGION_GROUPS |= {
    "framework": REGION_GROUPS["hframework"] | REGION_GROUPS["lframework"],
    "cdrs": REGION_GROUPS["hcdrs"] | REGION_GROUPS["lcdrs"],
    "vh": REGION_GROUPS["hframework"] | REGION_GROUPS["hcdrs"],
    "vl": REGION_GROUPS["lframework"] | REGION_GROUPS["lcdrs"],
}
REGION_GROUPS["fv"] = REGION_GROUPS["vh"] | REGION_GROUPS["vl"]


def _expand(regions: Iterable[Region | str]) -> set[Region]:
    expanded: set[Region] = set()
    for region in regions:
        name = region.lower()
        if name in REGION_GROUPS:
            expanded |= REGION_GROUPS[name]
        else:
            try:
                expanded.add(Region(name))
            except ValueError:
                message = f"Unknown region or region group: {region!r}"
                raise ValueError(message) from None
    return expanded


def _as_position(position: Position | str) -> Position:
    return position if isinstance(position, Position) else Position.parse(position)


@cache
def _imgt_positions(scheme: Scheme, chain: Chain) -> dict[int, int]:
    """Map residue numbers in `scheme` to IMGT positions."""
    table = load_json("scheme_to_imgt.json")[scheme][chain]
    return {int(number): imgt for number, imgt in table.items()}


@cache
def _regions_by_imgt_position(
    definition: Definition, chain: Chain
) -> dict[int, Region]:
    """Map each IMGT position (1-128) to its region under `definition`."""
    boundaries = load_toml("regions.toml")[definition][chain]
    return {
        imgt: Region(name)
        for name, (start, end) in boundaries.items()
        for imgt in range(start, end + 1)
    }


def get_region(
    position: Position | str,
    chain: Chain | str,
    *,
    scheme: Scheme | str,
    definition: Definition | str,
) -> Region | None:
    """Find the region that a single numbered position belongs to.

    This ignores the rest of the sequence, so it cannot account for insertions
    elsewhere in the domain.  When the definition was specified in a different
    scheme from `scheme`, it can therefore misassign positions near CDR-H1 and
    CDR-L2; use `annotate_regions` for whole sequences.

    Args:
        position: The position, as a `Position` or a string such as `"100A"`.
        chain: The chain type.
        scheme: The scheme that the position is numbered in.
        definition: The region definition to apply.

    Returns:
        The region, or `None` if the position lies outside the numbered
        variable domain.

    >>> get_region("100A", "H", scheme="chothia", definition="chothia")
    <Region.CDRH3: 'cdrh3'>
    >>> get_region("27", "L", scheme="imgt", definition="kabat")
    <Region.CDRL1: 'cdrl1'>
    """
    position = _as_position(position)
    chain, scheme, definition = Chain(chain), Scheme(scheme), Definition(definition)
    number, insertion = position.number, position.insertion

    # Kabat's CDR-H1 ends at H35B in Kabat numbering, so further insertions on
    # H35 fall in the framework.
    if (
        definition is Definition.KABAT
        and scheme is Scheme.KABAT
        and chain is Chain.HEAVY
        and 31 <= number <= 35
    ):
        if number == 35 and insertion not in ("", "A", "B"):
            return Region.FWH2
        return Region.CDRH1

    # Kabat numbering puts CDR-H1 insertions on H35, after the end of the
    # Chothia and IMGT definitions of CDR-H1.
    if scheme is Scheme.KABAT and chain is Chain.HEAVY:
        if definition is Definition.CHOTHIA and 33 <= number <= 35:
            return Region.FWH2
        if definition is Definition.IMGT and 34 <= number <= 35:
            return Region.FWH2

    imgt = _imgt_positions(scheme, chain).get(number)
    if imgt is None:
        return None
    return _regions_by_imgt_position(definition, chain)[imgt]


class RegionSelector:
    """Select positions by region, with optional extra and excluded positions.

    Args:
        regions: Region names (such as `"cdrh3"`) or group names (such as
            `"framework"`; see `REGION_GROUPS`) to select.
        scheme: The scheme that positions are numbered in.
        definition: The region definition to apply.

    Raises:
        ValueError: If a region or group name is not recognised.

    >>> selector = RegionSelector(["hcdrs"], scheme="imgt", definition="imgt")
    >>> selector.accepts("111A", "H")
    True
    >>> selector.exclude_positions(["111A"], "H")
    >>> selector.accepts("111A", "H")
    False
    """

    def __init__(
        self,
        regions: Iterable[Region | str] = (),
        *,
        scheme: Scheme | str,
        definition: Definition | str,
    ) -> None:
        """Create a selector for the given regions."""
        self.scheme = Scheme(scheme)
        self.definition = Definition(definition)
        self.regions: set[Region] = _expand(regions)
        self._extra: dict[Chain, set[Position]] = {chain: set() for chain in Chain}
        self._excluded: dict[Chain, set[Position]] = {chain: set() for chain in Chain}

    def add_regions(self, regions: Iterable[Region | str]) -> None:
        """Add regions or groups of regions to the selection.

        Args:
            regions: Region or group names.

        Raises:
            ValueError: If a region or group name is not recognised.
        """
        self.regions |= _expand(regions)

    def add_positions(
        self, positions: Iterable[Position | str], chain: Chain | str
    ) -> None:
        """Select individual positions, whatever their region.

        Args:
            positions: The positions to select.
            chain: The chain that the positions are on.
        """
        self._extra[Chain(chain)].update(map(_as_position, positions))

    def exclude_positions(
        self, positions: Iterable[Position | str], chain: Chain | str
    ) -> None:
        """Exclude individual positions, overriding any other selection.

        Args:
            positions: The positions to exclude.
            chain: The chain that the positions are on.
        """
        self._excluded[Chain(chain)].update(map(_as_position, positions))

    def accepts(self, position: Position | str, chain: Chain | str) -> bool:
        """Whether the selection includes a position.

        Args:
            position: The position.
            chain: The chain that the position is on.

        Returns:
            `True` if the position is selected and not excluded.
        """
        position, chain = _as_position(position), Chain(chain)
        if position in self._excluded[chain]:
            return False
        if position in self._extra[chain]:
            return True
        region = get_region(
            position, chain, scheme=self.scheme, definition=self.definition
        )
        return region in self.regions


def _insertions(number: int, insertions: str = ascii_uppercase) -> list[Position]:
    """Positions on `number` with each of the given insertion codes."""
    return [Position(number, insertion) for insertion in insertions]


def _with_insertions(number: int) -> list[Position]:
    """Position `number`, followed by its insertions in order."""
    return [Position(number), *_insertions(number)]


def _count_insertions(residues: dict[Position, str], number: int) -> int:
    """Count the consecutive insertions on `number` that hold a residue."""
    count = 0
    for insertion in ascii_uppercase:
        residue = residues.get(Position(number, insertion))
        if residue is None:
            break
        if residue != "-":
            count += 1
    return count


def annotate_regions(
    numbering: Iterable[tuple[Position | str, str]],
    chain: Chain | str,
    *,
    scheme: Scheme | str,
    definition: Definition | str,
) -> list[tuple[Position, str, Region | None]]:
    """Assign each residue of a numbered sequence to its region.

    Unlike `get_region`, this takes account of the insertions in the whole
    sequence, so that a definition applied in a different scheme from its own
    gives the same CDRs as in its own scheme.  (The Contact definition cannot
    be reproduced exactly in Kabat numbering.)

    Framework regions are assigned by their order relative to the CDRs: a
    residue before CDR 1 is in framework 1, one between CDRs 1 and 2 is in
    framework 2, and so on.  The numbering must therefore start at the
    beginning of the domain, or at least before CDR 1.

    Args:
        numbering: The `(position, residue)` pairs of the sequence, in sequence
            order.  A residue of `"-"` marks a gap.
        chain: The chain type.
        scheme: The scheme that the sequence is numbered in.
        definition: The region definition to apply.

    Returns:
        A `(position, residue, region)` triple for each residue, where the
        region is `None` for residues beyond the numbered variable domain.

    A sparse IMGT-numbered heavy chain, with one residue in each region:

    >>> numbering = [
    ...     ("1", "E"),
    ...     ("27", "G"),
    ...     ("39", "M"),
    ...     ("56", "I"),
    ...     ("66", "Y"),
    ...     ("111A", "D"),
    ...     ("118", "W"),
    ... ]
    >>> [
    ...     str(region)
    ...     for *_, region in annotate_regions(
    ...         numbering, "H", scheme="imgt", definition="imgt"
    ...     )
    ... ]
    ['fwh1', 'cdrh1', 'fwh2', 'cdrh2', 'fwh3', 'cdrh3', 'fwh4']
    """
    chain, scheme, definition = Chain(chain), Scheme(scheme), Definition(definition)
    numbering = [(_as_position(position), residue) for position, residue in numbering]
    residues = dict(numbering)
    c = chain.lower()

    def selector(n: int) -> RegionSelector:
        return RegionSelector([f"cdr{c}{n}"], scheme=scheme, definition=definition)

    cdr1, cdr2, cdr3 = selector(1), selector(2), selector(3)

    def first(positions: Iterable[Position], n: int) -> list[Position]:
        return list(islice(positions, n))

    if chain is Chain.HEAVY:
        if scheme is Scheme.IMGT and definition is Definition.KABAT:
            # Residues at IMGT H31-H34, and insertions on H33, become insertions
            # on H35 in Kabat numbering.  Extend CDR-H1 back by that many
            # positions; beyond two, the Kabat CDR-H1 (which ends at H35B) loses
            # positions at its C-terminal end instead.
            count = sum(
                residues.get(Position(number), "-") != "-" for number in range(31, 35)
            ) + _count_insertions(residues, 33)
            if count:
                extend = [Position(number) for number in (35, 34, 33, 32)]
                extend += _insertions(33)
                cdr1.add_positions(first(extend, count), chain)
                if count > 2:
                    trim = [Position(number) for number in range(40, 33, -1)]
                    trim += _insertions(33, ascii_uppercase[5:])
                    cdr1.exclude_positions(first(trim, count - 2), chain)
        elif scheme is Scheme.KABAT and definition in (
            Definition.CHOTHIA,
            Definition.IMGT,
        ):
            # Kabat numbering puts CDR-H1 insertions on H35, after the end of
            # the Chothia and IMGT CDR-H1s, so add positions from H33 (Chothia)
            # or H34 (IMGT) onwards, one for each insertion.
            count = _count_insertions(residues, 35)
            if count:
                first_added = 33 if definition is Definition.CHOTHIA else 34
                extend = [Position(number) for number in range(first_added, 35)]
                extend += _with_insertions(35)
                cdr1.add_positions(first(extend, count), chain)
        elif (
            scheme in (Scheme.CHOTHIA, Scheme.MARTIN) and definition is Definition.KABAT
        ):
            # Chothia numbering puts CDR-H1 insertions on H31.  Kabat's CDR-H1
            # has at most two insertions (H35A and H35B), so for each further
            # insertion, drop a position from its C-terminal end.
            count = _count_insertions(residues, 31)
            if count > 2:
                trim = [Position(number) for number in range(35, 31, -1)]
                trim += _insertions(31, ascii_uppercase[6:])
                cdr1.exclude_positions(first(trim, count - 2), chain)
    elif (
        scheme in (Scheme.KABAT, Scheme.CHOTHIA, Scheme.MARTIN)
        and definition is Definition.IMGT
    ):
        # Without insertions on L54, the IMGT CDR-L2 ends at L52 in Kabat-like
        # numbering; extend it by one position for each insertion.
        count = _count_insertions(residues, 54)
        if count:
            extend = [Position(53), *_with_insertions(54)]
            cdr2.add_positions(first(extend, count), chain)

    current = Region(f"fw{c}1")
    c_terminus = max(_imgt_positions(scheme, chain))
    annotated: list[tuple[Position, str, Region | None]] = []
    for position, residue in numbering:
        region: Region | None
        if cdr1.accepts(position, chain):
            region = Region(f"cdr{c}1")
            current = Region(f"fw{c}2")
        elif cdr2.accepts(position, chain):
            region = Region(f"cdr{c}2")
            current = Region(f"fw{c}3")
        elif cdr3.accepts(position, chain):
            region = Region(f"cdr{c}3")
            current = Region(f"fw{c}4")
        elif position.number <= c_terminus:
            region = current
        else:
            region = None
        annotated.append((position, residue, region))
    return annotated
