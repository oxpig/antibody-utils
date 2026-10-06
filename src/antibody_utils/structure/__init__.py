"""Antibody structure tools.

Read a structure with `read_structure`, which identifies and numbers its
antibody chains, then pair a heavy and a light chain into an `Fv`, and copy
regions of them with `select`.
"""

from antibody_utils.structure.models import AntibodyChain, AntibodyStructure, Fv
from antibody_utils.structure.parser import from_gemmi, read_structure
from antibody_utils.structure.selection import BACKBONE_ATOMS, select

__all__ = [
    "BACKBONE_ATOMS",
    "AntibodyChain",
    "AntibodyStructure",
    "Fv",
    "from_gemmi",
    "read_structure",
    "select",
]
