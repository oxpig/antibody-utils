"""Antibody chains in a structure, and the Fvs that pair them.

These classes wrap a `gemmi.Structure` rather than copying it.  Each
`AntibodyChain` maps the numbered positions of its variable domain to residues
of the structure, which keep their original numbering.  Treat the structure as
read-only while these objects are in use: adding or removing chains or
residues invalidates the mapping.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from antibody_utils.numbering.positions import Chain, Position, Scheme

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    import gemmi

    from antibody_utils.numbering.anarcii import NumberedSequence
    from antibody_utils.regions import Definition, Region

__all__ = ["AntibodyChain", "AntibodyStructure", "Fv"]


class AntibodyChain:
    """A chain of a structure with a numbered antibody variable domain.

    Only the residues of the variable domain are numbered.  Any others on the
    chain (such as a constant domain, a tag, ligands or water) are left out
    here, but remain in the structure.
    """

    __slots__ = ("_by_position", "_chain", "_indices", "_numbered", "_structure")

    def __init__(
        self, chain: gemmi.Chain, numbered: NumberedSequence, indices: Sequence[int]
    ) -> None:
        """Wrap a chain.  Use `read_structure` rather than calling this.

        Args:
            chain: The chain in the structure.
            numbered: The chain's numbered variable domain.
            indices: The index in `chain` of each residue of `numbered`.

        Raises:
            ValueError: If `indices` and `numbered` differ in length.
        """
        if len(indices) != len(numbered):
            raise ValueError(
                f"{len(indices)} residues but {len(numbered)} numbered positions"
            )
        self._chain = chain
        self._indices = tuple(indices)
        self._by_position = dict(zip(numbered.positions, self._indices, strict=True))
        self._structure: AntibodyStructure | None = None
        self._numbered = numbered

    def __repr__(self) -> str:
        """Describe the chain."""
        return (
            f"<AntibodyChain {self.name} {self.chain.name.lower()}, "
            f"{len(self)} residues numbered in the {self.scheme.value} scheme>"
        )

    @property
    def numbered(self) -> NumberedSequence:
        """The numbered variable domain."""
        return self._numbered

    @property
    def name(self) -> str:
        """The chain's name in the structure."""
        return self._chain.name

    @property
    def chain(self) -> Chain:
        """Whether this is a heavy or a light chain."""
        return self.numbered.chain

    @property
    def scheme(self) -> Scheme:
        """The numbering scheme."""
        return Scheme(self.numbered.scheme)

    @property
    def gemmi_chain(self) -> gemmi.Chain:
        """The whole chain in the structure, including unnumbered residues."""
        return self._chain

    @property
    def residues(self) -> list[gemmi.Residue]:
        """The residues of the numbered domain, in sequence order."""
        return [self._chain[i] for i in self._indices]

    def __len__(self) -> int:
        """The number of numbered residues."""
        return len(self._indices)

    def __iter__(self) -> Iterator[tuple[Position, gemmi.Residue]]:
        """Iterate over the `(position, residue)` pairs, in sequence order."""
        for position, i in zip(self.numbered.positions, self._indices, strict=True):
            yield position, self._chain[i]

    def __contains__(self, position: object) -> bool:
        """Whether a position, such as `"100A"`, holds a residue."""
        try:
            return _as_position(position) in self._by_position
        except (TypeError, ValueError):
            return False

    def __getitem__(self, position: Position | str) -> gemmi.Residue:
        """Get the residue at a position, such as `"100A"`.

        Raises:
            KeyError: If there is no residue at that position.
        """
        try:
            return self._chain[self._by_position[_as_position(position)]]
        except KeyError:
            raise KeyError(position) from None

    def region_residues(
        self, definition: Definition | str
    ) -> dict[Region, list[gemmi.Residue]]:
        """Group the residues of the domain by region.

        Regions are assigned with `antibody_utils.regions.annotate_regions`,
        as for `antibody_utils.sequence.region_sequences`.

        Args:
            definition: The region definition to apply.

        Returns:
            The residues of each region present, in sequence order.
        """
        residues = self.residues
        grouped: dict[Region, list[gemmi.Residue]] = {}
        for residue, (_, _, region) in zip(
            residues, self.numbered.annotate_regions(definition), strict=True
        ):
            if region is not None:
                grouped.setdefault(region, []).append(residue)
        return grouped


def _as_position(position: object) -> Position:
    if isinstance(position, Position):
        return position
    if isinstance(position, str):
        return Position.parse(position)
    raise TypeError(f"Expected a Position or a string, not {type(position).__name__}")


class Fv:
    """A variable fragment: a heavy and a light chain that pair.

    Chains are not paired automatically: choose the pair yourself, for example
    from the structure's metadata: `Fv(structure["H"], structure["L"])`.
    """

    __slots__ = ("_heavy", "_light")

    def __init__(self, heavy: AntibodyChain, light: AntibodyChain) -> None:
        """Pair a heavy and a light chain.

        Args:
            heavy: The heavy chain.
            light: The light chain.

        Raises:
            ValueError: If the chains are of the wrong types, are numbered in
                different schemes, or come from different structures.
        """
        if heavy.chain is not Chain.HEAVY:
            raise ValueError(f"Chain {heavy.name} is not a heavy chain")
        if light.chain is not Chain.LIGHT:
            raise ValueError(f"Chain {light.name} is not a light chain")
        if heavy.scheme is not light.scheme:
            raise ValueError("The chains are numbered in different schemes")
        if heavy._structure is not light._structure:
            raise ValueError("The chains come from different structures")
        self._heavy = heavy
        self._light = light

    def __repr__(self) -> str:
        """Describe the Fv."""
        return f"<Fv {self.name} heavy={self.heavy.name} light={self.light.name}>"

    @property
    def heavy(self) -> AntibodyChain:
        """The heavy chain."""
        return self._heavy

    @property
    def light(self) -> AntibodyChain:
        """The light chain."""
        return self._light

    @property
    def name(self) -> str:
        """The heavy and light chain names, joined, such as `"HL"`."""
        return self.heavy.name + self.light.name

    @property
    def scheme(self) -> Scheme:
        """The numbering scheme."""
        return self.heavy.scheme

    @property
    def numbered(self) -> tuple[NumberedSequence, NumberedSequence]:
        """The numbered heavy and light domains.

        These can be passed straight to `antibody_utils.sequence.identity` and
        `antibody_utils.sequence.similarity`.
        """
        return self.heavy.numbered, self.light.numbered

    def __iter__(self) -> Iterator[AntibodyChain]:
        """Iterate over the heavy chain, then the light chain."""
        yield self.heavy
        yield self.light

    def region_residues(
        self, definition: Definition | str
    ) -> dict[Region, list[gemmi.Residue]]:
        """Group the residues of both domains by region.

        Args:
            definition: The region definition to apply.

        Returns:
            The residues of each region present: the heavy chain's, then the
            light chain's.
        """
        return self.heavy.region_residues(definition) | self.light.region_residues(
            definition
        )


class AntibodyStructure:
    """A structure with its antibody chains identified and numbered.

    Create one with `antibody_utils.structure.read_structure` or
    `antibody_utils.structure.from_gemmi`.
    """

    __slots__ = ("_by_name", "_model_index", "chains", "scheme", "structure")

    def __init__(
        self,
        structure: gemmi.Structure,
        chains: Sequence[AntibodyChain],
        scheme: Scheme,
        model: int,
    ) -> None:
        """Collect the antibody chains of one model of a structure.

        Args:
            structure: The structure.
            chains: The antibody chains, from `structure[model]`.
            scheme: The numbering scheme of every chain.
            model: The index of the model that the chains belong to.

        Raises:
            ValueError: If two antibody chains have the same name, or a chain
                is not numbered in `scheme`.
        """
        #: The underlying GEMMI structure, with its original numbering.
        self.structure = structure
        #: The antibody chains, in the order they appear in the structure.
        self.chains = tuple(chains)
        #: The numbering scheme.
        self.scheme = scheme
        self._model_index = model
        self._by_name = {chain.name: chain for chain in self.chains}
        if len(self._by_name) != len(self.chains):
            raise ValueError("Two antibody chains have the same name")
        for chain in self.chains:
            if chain.scheme is not scheme:
                raise ValueError(
                    f"Chain {chain.name} is numbered in {chain.scheme.value}, "
                    f"not {scheme.value}"
                )
        for chain in self.chains:
            chain._structure = self

    def __repr__(self) -> str:
        """Describe the structure."""
        names = ", ".join(chain.name for chain in self.chains) or "none"
        return f"<AntibodyStructure {self.structure.name} antibody chains: {names}>"

    @property
    def model(self) -> gemmi.Model:
        """The model that the antibody chains belong to."""
        return self.structure[self._model_index]

    @property
    def heavy_chains(self) -> tuple[AntibodyChain, ...]:
        """The heavy chains."""
        return tuple(c for c in self.chains if c.chain is Chain.HEAVY)

    @property
    def light_chains(self) -> tuple[AntibodyChain, ...]:
        """The light chains."""
        return tuple(c for c in self.chains if c.chain is Chain.LIGHT)

    @property
    def other_chains(self) -> list[gemmi.Chain]:
        """The chains of the model that are not antibody chains."""
        return [chain for chain in self.model if chain.name not in self._by_name]

    def __len__(self) -> int:
        """The number of antibody chains."""
        return len(self.chains)

    def __iter__(self) -> Iterator[AntibodyChain]:
        """Iterate over the antibody chains."""
        return iter(self.chains)

    def __getitem__(self, name: str) -> AntibodyChain:
        """Get an antibody chain by name.

        Raises:
            KeyError: If there is no antibody chain of that name.
        """
        return self._by_name[name]
