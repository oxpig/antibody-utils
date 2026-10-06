"""Measure the relative orientation of the VH and VL domains with ABangle.

ABangle describes the orientation of an Fv's VH and VL domains with five
angles and a distance.  It was introduced by Dunbar J, Fuchs A, Shi J and
Deane CM, `"ABangle: characterising the VH–VL orientation in antibodies"
<https://doi.org/10.1093/protein/gzt020>`__, *Protein Eng Des Sel* 26:611–620
(2013).  Each domain is given a coordinate frame by superposing a consensus
domain on to the CA atoms of a core set of framework positions.  A vector C
joins the two domains, and each domain has two vectors, 1 and 2, in a plane
fitted to it:

- HL: the torsion angle between H1 and L1 about C.
- HC1 and HC2: the angles between H1 and H2, and C.
- LC1 and LC2: the angles between L1 and L2, and C.
- dc: the length of C.

The core sets are defined in Chothia numbering, so the Fv must be numbered in
the Chothia or IMGT scheme.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from functools import cache

import numpy as np

from antibody_utils._data import load_toml
from antibody_utils.geometry.superposition import superpose
from antibody_utils.numbering.positions import Position, Scheme
from antibody_utils.structure.models import AntibodyChain, Fv

__all__ = ["OrientationAngles", "abangle"]

# The schemes in which the core sets can be located: the key of each core
# position in that scheme.
_SCHEMES = {Scheme.CHOTHIA: "chothia", Scheme.IMGT: "imgt"}


@dataclass(frozen=True, slots=True)
class OrientationAngles:
    """The ABangle measures of an Fv.

    Attributes:
        hl: The torsion angle HL, in degrees, from -180 to 180.
        hc1: The bend angle HC1, in degrees.
        lc1: The bend angle LC1, in degrees.
        hc2: The twist angle HC2, in degrees.
        lc2: The twist angle LC2, in degrees.
        dc: The distance dc between the domains, in angstroms.
    """

    hl: float
    hc1: float
    lc1: float
    hc2: float
    lc2: float
    dc: float

    def as_dict(self) -> dict[str, float]:
        """The measures, keyed by their names in the ABangle paper (`"HL"` …)."""
        return {
            name.upper() if name != "dc" else name: value
            for name, value in asdict(self).items()
        }


@dataclass(frozen=True, slots=True)
class _Domain:
    """ABangle's definition of one domain."""

    consensus: dict[str, dict[Position, np.ndarray]]  # by scheme key
    point: np.ndarray  # the domain's end of C, in the consensus frame
    plane: tuple[np.ndarray, np.ndarray]  # points along vectors 1 and 2


@cache
def _domain(chain: str) -> _Domain:
    data = load_toml("abangle.toml")[chain]
    vectors = [np.array(v) for v in data["plane_vectors"]]
    a, b = data["separation"]
    point = np.array(data["point"]) + a * vectors[0] + b * vectors[1]
    consensus: dict[str, dict[Position, np.ndarray]] = {"chothia": {}, "imgt": {}}
    for chothia, entry in data["core"].items():
        coordinates = np.array(entry["ca"])
        consensus["chothia"][Position(int(chothia))] = coordinates
        consensus["imgt"][Position(entry["imgt"])] = coordinates
    return _Domain(consensus, point, (point + vectors[0], point + vectors[1]))


def _frame(chain: AntibodyChain) -> np.ndarray:
    """Place a domain's C end and plane points on the structure.

    Returns:
        The C end, and the points along vectors 1 and 2, as rows.
    """
    domain = _domain(chain.chain.value)
    consensus = domain.consensus[_SCHEMES[chain.scheme]]
    mobile, target = [], []
    for position, residue in chain:
        if position in consensus:
            atom = residue.find_atom("CA", "*")
            if atom is not None:
                mobile.append(consensus[position])
                target.append([atom.pos.x, atom.pos.y, atom.pos.z])
    if len(target) < 3:
        raise ValueError(
            f"Chain {chain.name} has CA atoms at only {len(target)} of the "
            f"{len(consensus)} ABangle core positions"
        )
    superposition = superpose(mobile, target)
    return superposition.apply([domain.point, *domain.plane])


def _unit(vector: np.ndarray) -> np.ndarray:
    return vector / np.linalg.norm(vector)


def _degrees(cosine: float) -> float:
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def abangle(fv: Fv) -> OrientationAngles:
    """Calculate the ABangle orientation measures of an Fv.

    Args:
        fv: The Fv, numbered in the Chothia or IMGT scheme.

    Returns:
        The angles and distance.

    Raises:
        ValueError: If the Fv is numbered in another scheme, or either domain
            has CA atoms at fewer than three core positions.
    """
    if fv.scheme not in _SCHEMES:
        raise ValueError(
            f"ABangle needs Chothia or IMGT numbering, not {fv.scheme.value}"
        )
    heavy, light = _frame(fv.heavy), _frame(fv.light)

    separation = heavy[0] - light[0]
    c = _unit(separation)
    h1, h2 = _unit(heavy[1] - heavy[0]), _unit(heavy[2] - heavy[0])
    l1, l2 = _unit(light[1] - light[0]), _unit(light[2] - light[0])

    # Project H1 and L1 on to the plane perpendicular to C, and measure the
    # torsion between them, signed by its direction about C.
    n_x = np.cross(l1, c)
    n_y = np.cross(c, n_x)
    l1_projected = _unit(np.array([0.0, l1 @ n_x, l1 @ n_y]))
    h1_projected = _unit(np.array([0.0, h1 @ n_x, h1 @ n_y]))
    hl = _degrees(l1_projected @ h1_projected)
    if np.cross(l1_projected, h1_projected)[0] < 0:
        hl = -hl

    return OrientationAngles(
        hl=hl,
        hc1=_degrees(h1 @ -c),
        lc1=_degrees(l1 @ c),
        hc2=_degrees(h2 @ -c),
        lc2=_degrees(l2 @ c),
        dc=float(np.linalg.norm(separation)),
    )
