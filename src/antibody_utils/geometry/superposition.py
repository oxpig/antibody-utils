"""Superpose antibody structures and measure their RMSD.

Coordinates are compared atom by atom at equivalent numbered positions, so the
structures being compared must be numbered in the same scheme.  Only positions
present in both structures are used.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from antibody_utils.numbering.positions import Chain, Position
from antibody_utils.regions import Definition, Region, _expand
from antibody_utils.structure.models import AntibodyChain, Fv

__all__ = [
    "Superposition",
    "orientation_rmsd",
    "region_rmsd",
    "rmsd",
    "superpose",
]


@dataclass(frozen=True, slots=True)
class Superposition:
    """A rigid-body transformation that superposes one set of points on another.

    Attributes:
        rotation: The 3 x 3 rotation matrix.
        translation: The translation, applied after the rotation.
        rmsd: The RMSD between the superposed points and the target points.
    """

    rotation: np.ndarray
    translation: np.ndarray
    rmsd: float

    def apply(self, coordinates: ArrayLike) -> np.ndarray:
        """Transform coordinates: rotate them, then translate them.

        Args:
            coordinates: A point, or an n x 3 array of points.

        Returns:
            The transformed coordinates, in the same shape.
        """
        return np.asarray(coordinates, dtype=float) @ self.rotation.T + self.translation


def _points(coordinates: ArrayLike, name: str) -> np.ndarray:
    points = np.asarray(coordinates, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"{name} must be an n x 3 array of points")
    return points


def rmsd(a: ArrayLike, b: ArrayLike) -> float:
    """Calculate the RMSD between two sets of equivalent points, as they are.

    Args:
        a: An n x 3 array of points.
        b: The equivalent n x 3 array of points.

    Returns:
        The root-mean-square deviation.

    Raises:
        ValueError: If the arrays are not n x 3 arrays of the same shape, or
            are empty.

    >>> rmsd([[0, 0, 0], [1, 0, 0]], [[0, 0, 1], [1, 0, 1]])
    1.0
    """
    a, b = _points(a, "a"), _points(b, "b")
    if a.shape != b.shape or not len(a):
        raise ValueError("a and b must hold the same, non-zero, number of points")
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def superpose(mobile: ArrayLike, target: ArrayLike) -> Superposition:
    """Find the superposition of one set of points on another (Kabsch).

    Args:
        mobile: An n x 3 array of points, to be moved.
        target: The equivalent n x 3 array of points, held fixed.

    Returns:
        The rotation and translation that minimise the RMSD between the moved
        `mobile` points and the `target` points.

    Raises:
        ValueError: If the arrays are not n x 3 arrays of the same shape, or
            hold fewer than three points.

    >>> import numpy as np
    >>> target = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> mobile = target @ [[0, -1, 0], [1, 0, 0], [0, 0, 1]] + [5, 5, 5]
    >>> superposition = superpose(mobile, target)
    >>> round(superposition.rmsd, 6)
    0.0
    >>> np.allclose(superposition.apply(mobile), target)
    True
    """
    mobile, target = _points(mobile, "mobile"), _points(target, "target")
    if mobile.shape != target.shape:
        raise ValueError("mobile and target must hold the same number of points")
    if len(mobile) < 3:
        raise ValueError("Superposition needs at least three points")
    mobile_centre, target_centre = mobile.mean(axis=0), target.mean(axis=0)
    covariance = (mobile - mobile_centre).T @ (target - target_centre)
    u, _, vt = np.linalg.svd(covariance)
    # Correct for a reflection, so that the result is a proper rotation.
    sign = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
    rotation = vt.T @ np.diag([1.0, 1.0, sign]) @ u.T
    translation = target_centre - mobile_centre @ rotation.T
    moved = mobile @ rotation.T + translation
    return Superposition(rotation, translation, rmsd(moved, target))


def _atoms(
    source: AntibodyChain | Fv,
    regions: Iterable[Region | str],
    definition: Definition | str,
    atoms: Iterable[str],
) -> dict[tuple[Chain, Position, str], np.ndarray]:
    """The coordinates of the selected atoms, keyed by chain, position and name."""
    chains = [source] if isinstance(source, AntibodyChain) else list(source)
    wanted = _expand(regions)
    names = tuple(atoms)
    found = {}
    for chain in chains:
        for (position, _, region), residue in zip(
            chain.numbered.annotate_regions(definition), chain.residues, strict=True
        ):
            if region not in wanted:
                continue
            for name in names:
                atom = residue.find_atom(name, "*")
                if atom is not None:
                    pos = atom.pos
                    found[chain.chain, position, name] = np.array([pos.x, pos.y, pos.z])
    return found


def _equivalent(
    a: AntibodyChain | Fv,
    b: AntibodyChain | Fv,
    regions: Iterable[Region | str],
    definition: Definition | str,
    atoms: Iterable[str],
) -> tuple[np.ndarray, np.ndarray]:
    if a.scheme is not b.scheme:
        raise ValueError("The structures are numbered in different schemes")
    regions, atoms = tuple(regions), tuple(atoms)
    a_atoms = _atoms(a, regions, definition, atoms)
    b_atoms = _atoms(b, regions, definition, atoms)
    common = [key for key in a_atoms if key in b_atoms]
    if not common:
        raise ValueError("The structures have no selected atoms in common")
    return (
        np.array([a_atoms[key] for key in common]),
        np.array([b_atoms[key] for key in common]),
    )


def region_rmsd(
    a: AntibodyChain | Fv,
    b: AntibodyChain | Fv,
    regions: Iterable[Region | str] = ("fv",),
    *,
    definition: Definition | str,
    atoms: Iterable[str] = ("CA",),
    superposed: bool = True,
) -> float:
    """Calculate the RMSD between antibody structures over chosen regions.

    Atoms are matched by chain type, numbered position and name, and only those
    present in both structures are compared.  Regions are assigned as by
    `antibody_utils.structure.AntibodyChain.region_residues`.

    Args:
        a: An antibody chain or Fv.
        b: Another, numbered in the same scheme.
        regions: The regions, or region groups such as `"cdrs"`, to compare.
            By default, the whole variable domain (`"fv"`).
        definition: The region definition to apply.
        atoms: The names of the atoms to compare.  By default, CA atoms only.
        superposed: Superpose `a` on `b` over the selected atoms first, as is
            usual.  If false, compare the coordinates as they are.

    Returns:
        The RMSD, in angstroms.

    Raises:
        ValueError: If the structures are numbered in different schemes, have
            no selected atoms in common, or have fewer than three when
            superposing, or a region name is not recognised.
    """
    a_points, b_points = _equivalent(a, b, regions, definition, atoms)
    if superposed:
        return superpose(a_points, b_points).rmsd
    return rmsd(a_points, b_points)


def orientation_rmsd(a: Fv, b: Fv, *, definition: Definition | str) -> float:
    """Measure how differently two Fvs orient their VH and VL domains.

    The VH frameworks of `a` and `b` are superposed, and `a`'s VL framework is
    moved with them.  Separately, the VL frameworks are superposed.  The
    result is the RMSD between the two placements of `a`'s VL framework (CA
    atoms).  This follows legacy `ABDB`'s `calculate_orientation_rmsd`.

    The measure is not quite symmetric: swapping `a` and `b` can give a
    slightly different value, so consider averaging both.

    Args:
        a: An Fv.
        b: Another Fv, numbered in the same scheme.
        definition: The region definition that delimits the frameworks.

    Returns:
        The orientation RMSD, in angstroms.

    Raises:
        ValueError: If the Fvs are numbered in different schemes, or share
            fewer than three framework CA atoms in either domain.
    """
    a_heavy, b_heavy = _equivalent(a.heavy, b.heavy, ["hframework"], definition, ["CA"])
    a_light, b_light = _equivalent(a.light, b.light, ["lframework"], definition, ["CA"])
    by_heavy = superpose(a_heavy, b_heavy).apply(a_light)
    by_light = superpose(a_light, b_light).apply(a_light)
    return rmsd(by_heavy, by_light)
