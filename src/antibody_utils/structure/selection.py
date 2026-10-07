"""Copy selected regions and atoms of antibody chains into a new structure."""

from __future__ import annotations

from typing import TYPE_CHECKING

import gemmi

from antibody_utils.regions import Definition, Region, _expand
from antibody_utils.structure.models import AntibodyChain, Fv

if TYPE_CHECKING:
    from collections.abc import Iterable

__all__ = ["BACKBONE_ATOMS", "select"]

#: The names of the backbone atoms of an amino-acid residue.
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})


def select(
    source: AntibodyChain | Fv | Iterable[AntibodyChain],
    regions: Iterable[Region | str] = ("fv",),
    *,
    definition: Definition | str,
    exclude: Iterable[Region | str] = (),
    atoms: Iterable[str] | None = None,
) -> gemmi.Structure:
    """Copy regions of antibody chains into a new structure.

    The copy holds one model, with a chain for each antibody chain selected.
    Its residues are numbered in the chains' numbering scheme, so that, for
    example, `select(fv, definition="imgt").write_pdb("fv.pdb")` writes an
    IMGT-numbered Fv.  The original structure is left unchanged.

    Regions are assigned as by `AntibodyChain.region_residues`.

    Args:
        source: An antibody chain, an Fv, or several antibody chains from one
            structure.
        regions: The regions, or region groups such as `"cdrs"`, to copy.  By
            default, the whole variable domain (`"fv"`).
        definition: The region definition to apply.
        exclude: Regions, or region groups, to leave out, such as `["cdrh3"]`.
        atoms: The names of the atoms to copy, such as `BACKBONE_ATOMS`.  By
            default, every atom.

    Returns:
        The new structure, with the original's name, unit cell and space
        group.

    Raises:
        ValueError: If a region or group name is not recognised, or the chains
            come from different structures.
    """
    chains = [source] if isinstance(source, AntibodyChain) else list(source)
    wanted = _expand(regions) - _expand(exclude)
    atom_names = None if atoms is None else frozenset(atoms)

    structures = {id(chain._structure) for chain in chains}
    if len(structures) > 1:
        raise ValueError("The chains come from different structures")
    original = chains[0]._structure.structure if chains else gemmi.Structure()

    copy = gemmi.Structure()
    copy.name = original.name
    copy.cell = original.cell
    copy.spacegroup_hm = original.spacegroup_hm
    model = gemmi.Model(1)
    for antibody_chain in chains:
        new_chain = gemmi.Chain(antibody_chain.name)
        for (position, _, region), residue in zip(
            antibody_chain.numbered.annotate_regions(definition),
            antibody_chain.residues,
            strict=True,
        ):
            if region not in wanted:
                continue
            new_residue = residue.clone()
            new_residue.seqid = gemmi.SeqId(position.number, position.insertion or " ")
            if atom_names is not None:
                for i in reversed(range(len(new_residue))):
                    if new_residue[i].name not in atom_names:
                        del new_residue[i]
            new_chain.add_residue(new_residue)
        if len(new_chain):
            model.add_chain(new_chain)
    copy.add_model(model)
    copy.setup_entities()
    return copy
