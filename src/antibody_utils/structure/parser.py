"""Read structures and identify their numbered antibody chains.

There are two ways to number a structure's antibody chains:

- By default, every protein chain is numbered with ANARCII, which needs the `numbering`
  extra.  Chains that ANARCII cannot number are treated as non-antibody chains.
- If the file is already numbered, as SAbDab's files are, name its antibody chains with
  `chains` to use that numbering as it is, without renumberin with ANARCII.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import gemmi

from antibody_utils.numbering.anarcii import (
    DEFAULT_BATCH_SIZE,
    NumberedSequence,
    NumberingError,
    _iter_number,
)
from antibody_utils.numbering.positions import Chain, Position, Scheme
from antibody_utils.regions import get_region
from antibody_utils.structure.models import AntibodyChain, AntibodyStructure

if TYPE_CHECKING:
    from collections.abc import Mapping

__all__ = ["from_gemmi", "read_structure"]


def read_structure(
    path: str | os.PathLike[str],
    scheme: Scheme | str = Scheme.IMGT,
    *,
    chains: Mapping[str, Chain | str] | None = None,
    model: int = 0,
    cpu: bool = False,
) -> AntibodyStructure:
    """Read a structure file and number its antibody chains.

    Args:
        path: A PDB, mmCIF or other file that GEMMI can read, optionally
            gzipped.
        scheme: The numbering scheme.  With `chains`, this is the scheme that
            the file is already numbered in.
        chains: The antibody chains, as a mapping from chain names to their
            types (`"H"` or `"L"`), whose numbering in the file is used as it
            is.  By default, every chain is numbered with ANARCII instead.
        model: The index of the model to use.
        cpu: When numbering with ANARCII, run on the CPU even if a GPU is
            available.

    Returns:
        The structure, with its antibody chains.

    Raises:
        KeyError: If a chain named in `chains` is not in the structure.
        ValueError: If a chain named in `chains` has no residues numbered in
            the variable domain, or its numbering is out of order.
        ImportError: If `chains` is not given and ANARCII is not installed.
    """
    structure = gemmi.read_structure(os.fspath(path))
    return _antibody_structure(structure, Scheme(scheme), chains, model, cpu)


def from_gemmi(
    structure: gemmi.Structure,
    scheme: Scheme | str = Scheme.IMGT,
    *,
    chains: Mapping[str, Chain | str] | None = None,
    model: int = 0,
    cpu: bool = False,
) -> AntibodyStructure:
    """Number the antibody chains of a structure that is already loaded.

    The structure is copied, so later changes to `structure` do not affect the
    result.  The arguments are as for `read_structure`.

    Args:
        structure: The structure.
        scheme: The numbering scheme.
        chains: The antibody chains whose numbering to use as it is.
        model: The index of the model to use.
        cpu: When numbering with ANARCII, run on the CPU.

    Returns:
        The structure, with its antibody chains.
    """
    return _antibody_structure(structure.clone(), Scheme(scheme), chains, model, cpu)


def _antibody_structure(
    structure: gemmi.Structure,
    scheme: Scheme,
    chains: Mapping[str, Chain | str] | None,
    model: int,
    cpu: bool,
) -> AntibodyStructure:
    # Mark which residues belong to polymers, if the file didn't say.
    structure.setup_entities()
    gemmi_model = structure[model]
    if chains is None:
        antibody_chains = _number_with_anarcii(gemmi_model, scheme, cpu)
    else:
        antibody_chains = [
            _from_file_numbering(gemmi_model, name, Chain(chain_type), scheme)
            for name, chain_type in chains.items()
        ]
    order = {chain.name: i for i, chain in enumerate(gemmi_model)}
    antibody_chains.sort(key=lambda chain: order[chain.name])
    return AntibodyStructure(structure, antibody_chains, scheme, model)


def _amino_acids(chain: gemmi.Chain) -> tuple[list[int], str]:
    """Find the amino-acid residues of a chain's polymer.

    Returns:
        Each residue's index in the chain, and its one-letter code.  Residues
        that GEMMI doesn't know but that have a Cα atom are `X`.  Only the
        first of several alternative residues at one position is kept.
    """
    indices: list[int] = []
    codes: list[str] = []
    previous = None
    for i, residue in enumerate(chain):
        if residue.entity_type != gemmi.EntityType.Polymer:
            continue
        if residue.seqid == previous:
            continue
        info = gemmi.find_tabulated_residue(residue.name)
        if info is not None and info.found():
            if not info.is_amino_acid():
                continue
            code = info.one_letter_code.upper()
        elif residue.find_atom("CA", "*") is not None:
            code = "X"
        else:
            continue
        previous = residue.seqid
        indices.append(i)
        codes.append(code if code.isalpha() else "X")
    return indices, "".join(codes)


def _number_with_anarcii(
    model: gemmi.Model, scheme: Scheme, cpu: bool
) -> list[AntibodyChain]:
    residues: dict[str, tuple[gemmi.Chain, list[int]]] = {}
    sequences: dict[str, str] = {}
    for i, chain in enumerate(model):
        indices, sequence = _amino_acids(chain)
        if sequence:
            # Chain names need not be unique, so key by the chain's index.
            residues[str(i)] = chain, indices
            sequences[str(i)] = sequence

    antibody_chains = []
    for key, result in _iter_number(sequences, scheme, cpu, DEFAULT_BATCH_SIZE):
        if isinstance(result, NumberingError):
            continue
        chain, indices = residues[key]
        numbered = NumberedSequence(
            name=chain.name,
            chain_type=result.chain_type,
            scheme=result.scheme,
            positions=result.positions,
            sequence=result.sequence,
            score=result.score,
            start=result.start,
            end=result.end,
        )
        antibody_chains.append(
            AntibodyChain(chain, numbered, indices[result.start : result.end])
        )
    return antibody_chains


def _from_file_numbering(
    model: gemmi.Model, name: str, chain_type: Chain, scheme: Scheme
) -> AntibodyChain:
    chain = model.find_chain(name)
    if chain is None:
        raise KeyError(f"No chain {name!r} in the structure")
    indices, sequence = _amino_acids(chain)
    domain = []
    for n, i in enumerate(indices):
        seqid = chain[i].seqid
        position = Position(seqid.num, seqid.icode)
        region = get_region(position, chain_type, scheme=scheme, definition="imgt")
        if region is not None:
            domain.append((n, position))
    if not domain:
        raise ValueError(
            f"Chain {name!r} has no residues numbered in the {scheme.value} "
            f"{chain_type.name.lower()}-chain variable domain"
        )
    start, end = domain[0][0], domain[-1][0] + 1
    if end - start != len(domain):
        raise ValueError(
            f"Chain {name!r} has residues outside the variable domain among "
            f"those numbered in it"
        )
    positions = tuple(position for _, position in domain)
    if len(set(positions)) != len(positions):
        raise ValueError(f"Chain {name!r} has more than one residue per position")
    numbered = NumberedSequence(
        name=name,
        chain_type=chain_type.value,
        scheme=scheme,
        positions=positions,
        sequence=sequence[start:end],
        score=float("nan"),
        start=start,
        end=end,
    )
    return AntibodyChain(chain, numbered, indices[start:end])
