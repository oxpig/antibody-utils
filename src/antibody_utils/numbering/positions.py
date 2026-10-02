"""Numbered positions in antibody variable domains."""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum
from typing import Self


class _CaseInsensitiveStrEnum(StrEnum):
    """A string enum whose lookup by value ignores case."""

    @classmethod
    def _missing_(cls, value: object) -> Self | None:
        if isinstance(value, str):
            folded = value.casefold()
            for member in cls:
                if member.value.casefold() == folded:
                    return member
        return None


class Chain(_CaseInsensitiveStrEnum):
    """An antibody variable-domain chain type.

    Members compare equal to their one-letter codes, and lookup is
    case-insensitive.

    >>> Chain("h")
    <Chain.HEAVY: 'H'>
    >>> Chain.LIGHT == "L"
    True
    """

    HEAVY = "H"
    LIGHT = "L"


class Scheme(_CaseInsensitiveStrEnum):
    """An antibody numbering scheme.

    >>> Scheme("Chothia")
    <Scheme.CHOTHIA: 'chothia'>
    """

    IMGT = "imgt"
    KABAT = "kabat"
    CHOTHIA = "chothia"
    MARTIN = "martin"


_POSITION = re.compile(r"\s*(\d+)\s*([A-Za-z]?)\s*")
_CHAIN_POSITION = re.compile(r"\s*([HhLl])\s*(\d+\s*[A-Za-z]?)\s*")


@dataclass(frozen=True, order=True, slots=True)
class Position:
    """A numbered position: a residue number and an optional insertion code.

    Positions order by number, then by insertion code, which is the order of
    residues in the Kabat, Chothia and Martin schemes.  (IMGT orders the
    insertions at position 112 in reverse, so sort IMGT positions by sequence
    order rather than by `Position`.)

    The `insertion` attribute is `""` when there is no insertion code.  Gemmi,
    Biopython, ANARCI and ANARCII write that as a single space instead, which
    `Position` accepts on input; use `to_tuple` to convert back.

    Args:
        number: The residue number.
        insertion: The insertion code: an upper-case letter, or `""` for none.
            Whitespace, such as the single space used by other tools, also
            means none.

    Raises:
        ValueError: If the insertion code is not a single letter.

    >>> Position(100, "A")
    Position(number=100, insertion='A')
    >>> Position(35) < Position(35, "A") < Position(36)
    True
    """

    number: int
    insertion: str = ""

    def __post_init__(self) -> None:
        """Normalise and validate the insertion code."""
        if not self.insertion.strip():
            object.__setattr__(self, "insertion", "")
        if self.insertion and not (
            len(self.insertion) == 1 and "A" <= self.insertion <= "Z"
        ):
            raise ValueError(
                f"Insertion code must be a single upper-case letter or empty, "
                f"not {self.insertion!r}"
            )

    @classmethod
    def parse(cls, text: str) -> Self:
        """Parse a position such as `"100A"`.

        Args:
            text: The residue number, optionally followed by an insertion code.

        Returns:
            The position.

        Raises:
            ValueError: If the text is not a valid position.

        >>> Position.parse("100a")
        Position(number=100, insertion='A')
        >>> Position.parse("35")
        Position(number=35, insertion='')
        """
        if not (match := _POSITION.fullmatch(text)):
            raise ValueError(f"Not a valid position: {text!r}")
        number, insertion = match.groups()
        return cls(int(number), insertion.upper())

    @classmethod
    def from_tuple(cls, position: tuple[int | str, str]) -> Self:
        """Convert a `(number, insertion)` tuple, as used by ANARCI and Biopython.

        Args:
            position: The residue number (an integer, or a string of digits as
                ANARCII sometimes gives) and the insertion code, where `" "`
                means no insertion.

        Returns:
            The position.

        >>> Position.from_tuple((52, " "))
        Position(number=52, insertion='')
        >>> Position.from_tuple(("111", "A"))
        Position(number=111, insertion='A')
        """
        number, insertion = position
        return cls(int(number), insertion)

    def to_tuple(self) -> tuple[int, str]:
        """Convert to a `(number, insertion)` tuple, with `" "` for no insertion.

        This is the convention of Gemmi, Biopython, ANARCI and ANARCII.

        Returns:
            The residue number and insertion code.

        >>> Position(52).to_tuple()
        (52, ' ')
        >>> Position(52, "A").to_tuple()
        (52, 'A')
        """
        return self.number, self.insertion or " "

    def __str__(self) -> str:
        """Format the position as its number followed by any insertion code.

        >>> str(Position(100, "A"))
        '100A'
        """
        return f"{self.number}{self.insertion}"


def parse_chain_position(text: str) -> tuple[Chain, Position]:
    """Parse a chain-qualified position such as `"H100A"`.

    Args:
        text: A chain letter (`H` or `L`) followed by a position.

    Returns:
        The chain and the position.

    Raises:
        ValueError: If the text is not a valid chain-qualified position.

    >>> parse_chain_position("H100A")
    (<Chain.HEAVY: 'H'>, Position(number=100, insertion='A'))
    """
    if not (match := _CHAIN_POSITION.fullmatch(text)):
        raise ValueError(f"Not a valid chain-qualified position: {text!r}")
    chain, position = match.groups()
    return Chain(chain), Position.parse(position)
