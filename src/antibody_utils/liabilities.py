"""Find sequence liabilities in antibody variable domains.

Liabilities are sequence motifs that may affect developability: chemical
modification such as oxidation, deamidation or isomerisation, glycosylation,
fragmentation, or unwanted binding.  The motifs, and the regions in which they
count, are listed in `antibody_utils/data/liabilities.toml`.

Matches are found in IMGT-numbered domains, and assigned to regions under the
North definition (`North, Lehmann and Dunbrack 2011
<https://doi.org/10.1016/j.jmb.2010.10.030>`__).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import cache
from typing import TYPE_CHECKING

from antibody_utils._data import load_toml
from antibody_utils.numbering.positions import Chain, Position, Scheme
from antibody_utils.regions import RegionSelector

if TYPE_CHECKING:
    from antibody_utils.numbering.anarcii import NumberedSequence

__all__ = ["Liability", "find_liabilities"]


@dataclass(frozen=True, slots=True)
class Liability:
    """A match to a sequence-liability motif.

    Attributes:
        name: The liability, such as `"Asn deamidation (NG NS NT)"`.
        chain: The chain that the match is on.
        residues: The matched residues.
        start: The index of the match's first residue in the chain's numbered
            domain sequence.
        end: The index just past the match's last residue.
        positions: The IMGT position of each matched residue.
    """

    name: str
    chain: Chain
    residues: str
    start: int
    end: int
    positions: tuple[Position, ...]


@dataclass(frozen=True, slots=True)
class _Motif:
    name: str
    pattern: re.Pattern[str]
    regions: tuple[str, ...]
    position_sets: tuple[str, ...]
    ignore_positions: frozenset[Position]
    requires_both_chains: bool


@cache
def _motifs() -> tuple[tuple[_Motif, ...], dict[str, dict[Chain, list[Position]]]]:
    data = load_toml("liabilities.toml")
    position_sets = {
        name: {
            Chain(chain): [Position.parse(p) for p in positions]
            for chain, positions in by_chain.items()
        }
        for name, by_chain in data["positions"].items()
    }
    motifs = []
    for entry in data["liability"]:
        regions = entry["regions"]
        motifs.append(
            _Motif(
                name=entry["name"],
                pattern=re.compile(entry["motif"]),
                regions=tuple(r for r in regions if r not in position_sets),
                position_sets=tuple(r for r in regions if r in position_sets),
                ignore_positions=frozenset(
                    Position.parse(p) for p in entry.get("ignore_positions", ())
                ),
                requires_both_chains=entry.get("requires_both_chains", False),
            )
        )
    return tuple(motifs), position_sets


def _domain(
    domain: NumberedSequence | str | None, chain: Chain, cpu: bool
) -> NumberedSequence | None:
    if domain is None:
        return None
    if isinstance(domain, str):
        from antibody_utils.numbering.anarcii import number_sequence

        domain = number_sequence(domain, Scheme.IMGT, cpu=cpu)
    if Scheme(domain.scheme) is not Scheme.IMGT:
        raise ValueError("Liabilities are found in IMGT-numbered domains")
    if domain.chain is not chain:
        raise ValueError(
            f"Expected a {chain.name.lower()} chain, but {domain.name!r} is "
            f"{domain.chain.name.lower()}"
        )
    return domain


def _matches(motif: _Motif, domain: NumberedSequence, position_sets) -> list[Liability]:
    chain = domain.chain
    selector = RegionSelector(motif.regions, scheme=Scheme.IMGT, definition="north")
    for name in motif.position_sets:
        selector.add_positions(position_sets[name][chain], chain)
    found = []
    for match in motif.pattern.finditer(domain.sequence):
        start, end = match.span()
        first = domain.positions[start]
        if first in motif.ignore_positions or not selector.accepts(first, chain):
            continue
        found.append(
            Liability(
                name=motif.name,
                chain=chain,
                residues=match.group(),
                start=start,
                end=end,
                positions=domain.positions[start:end],
            )
        )
    return found


def find_liabilities(
    heavy: NumberedSequence | str | None = None,
    light: NumberedSequence | str | None = None,
    *,
    cpu: bool = False,
) -> list[Liability]:
    """Find sequence liabilities in a heavy and a light variable domain.

    Either domain may be omitted, for example for a single-domain antibody,
    but liabilities that need both chains (such as N-terminal glutamate on
    both) are then never found.

    Args:
        heavy: The heavy-chain domain, numbered in IMGT, or its sequence to be
            numbered with ANARCII (which needs the `numbering` extra).
        light: The light-chain domain, likewise.
        cpu: When numbering sequences, run on the CPU even if a GPU is
            available.

    Returns:
        The liabilities found, in the order of `liabilities.toml`, then heavy
        before light, then by position.

    Raises:
        ValueError: If a domain is not IMGT-numbered, or is of the wrong chain.
        antibody_utils.numbering.NumberingError: If a sequence could not be
            numbered.
    """
    domains = [
        domain
        for domain in (
            _domain(heavy, Chain.HEAVY, cpu),
            _domain(light, Chain.LIGHT, cpu),
        )
        if domain is not None
    ]
    motifs, position_sets = _motifs()
    found: list[Liability] = []
    for motif in motifs:
        by_chain = [_matches(motif, domain, position_sets) for domain in domains]
        if motif.requires_both_chains and not (len(by_chain) == 2 and all(by_chain)):
            continue
        for matches in by_chain:
            found.extend(matches)
    return found
