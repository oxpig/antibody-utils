import csv
import json
import math
from pathlib import Path

import gemmi
import numpy as np
import pytest

from antibody_utils.geometry import (
    OrientationAngles,
    abangle,
    orientation_rmsd,
    region_rmsd,
    rmsd,
    superpose,
)
from antibody_utils.structure import Fv, from_gemmi, read_structure, select

DATA = Path(__file__).parent / "data"
LEGACY = DATA / "legacy"

with open(LEGACY / "abangle.json") as _f:
    GOLDEN_ANGLES = json.load(_f)


def _rotation(axis, degrees):
    axis = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    x, y, z = axis
    c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
    return np.array(
        [
            [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
            [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
            [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
        ]
    )


def _legacy_fv(entry, scheme):
    heavy, light = entry["heavy"], entry["light"]
    path = LEGACY / "structures" / f"{entry['files'][scheme]}.gz"
    structure = read_structure(path, scheme, chains={heavy: "H", light: "L"})
    return Fv(structure[heavy], structure[light])


@pytest.fixture
def fv():
    return _legacy_fv(GOLDEN_ANGLES[0], "imgt")


def _moved(fv, rotation, translation):
    """A copy of an Fv's structure, rigidly moved."""
    structure = fv.heavy._structure.structure.clone()
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    p = rotation @ [atom.pos.x, atom.pos.y, atom.pos.z] + translation
                    atom.pos = gemmi.Position(*p)
    moved = from_gemmi(
        structure, fv.scheme, chains={fv.heavy.name: "H", fv.light.name: "L"}
    )
    return Fv(moved[fv.heavy.name], moved[fv.light.name])


# --- Superposition and RMSD of points ----------------------------------------


def test_superpose_recovers_a_rigid_motion():
    rng = np.random.default_rng(0)
    target = rng.normal(size=(20, 3)) * 10
    rotation = _rotation([1, 2, 3], 70)
    mobile = (target - [4, 5, 6]) @ rotation
    superposition = superpose(mobile, target)
    assert superposition.rmsd == pytest.approx(0, abs=1e-9)
    assert np.allclose(superposition.apply(mobile), target)
    assert np.allclose(superposition.rotation, rotation)
    assert np.linalg.det(superposition.rotation) == pytest.approx(1)


def test_superpose_never_reflects():
    rng = np.random.default_rng(1)
    target = rng.normal(size=(10, 3))
    mirrored = target * [1, 1, -1]
    superposition = superpose(mirrored, target)
    assert np.linalg.det(superposition.rotation) == pytest.approx(1)
    assert superposition.rmsd > 0.1


def test_superposition_applies_to_one_point():
    superposition = superpose(np.eye(3), np.eye(3) + 1)
    assert np.allclose(superposition.apply([0, 0, 0]), [1, 1, 1])


@pytest.mark.parametrize(
    ("mobile", "target", "message"),
    [
        (np.zeros((3, 2)), np.zeros((3, 2)), "n x 3"),
        (np.zeros((3, 3)), np.zeros((4, 3)), "same number"),
        (np.zeros((2, 3)), np.zeros((2, 3)), "at least three"),
    ],
)
def test_superpose_rejects_bad_input(mobile, target, message):
    with pytest.raises(ValueError, match=message):
        superpose(mobile, target)


def test_rmsd():
    assert rmsd([[0, 0, 0], [0, 0, 0]], [[3, 4, 0], [0, 0, 0]]) == pytest.approx(
        math.sqrt(12.5)
    )
    with pytest.raises(ValueError, match="same, non-zero"):
        rmsd(np.zeros((0, 3)), np.zeros((0, 3)))
    with pytest.raises(ValueError, match="same, non-zero"):
        rmsd(np.zeros((1, 3)), np.zeros((2, 3)))


# --- RMSD of antibody structures ---------------------------------------------


def test_region_rmsd_of_a_structure_with_itself(fv):
    assert region_rmsd(fv, fv, definition="imgt") == pytest.approx(0, abs=1e-6)


def test_region_rmsd_superposes(fv):
    moved = _moved(fv, _rotation([0, 1, 1], 40), np.array([10.0, -3.0, 2.0]))
    assert region_rmsd(fv, moved, definition="imgt") == pytest.approx(0, abs=1e-3)
    assert region_rmsd(fv, moved, definition="imgt", superposed=False) > 1


def test_region_rmsd_between_copies():
    entry_hl, entry_pm = GOLDEN_ANGLES[0], GOLDEN_ANGLES[1]
    hl, pm = _legacy_fv(entry_hl, "imgt"), _legacy_fv(entry_pm, "imgt")
    whole = region_rmsd(hl, pm, definition="imgt")
    framework = region_rmsd(hl, pm, ["framework"], definition="imgt")
    cdrh3 = region_rmsd(hl.heavy, pm.heavy, ["cdrh3"], definition="imgt")
    backbone = region_rmsd(hl, pm, definition="imgt", atoms=["N", "CA", "C", "O"])
    assert 0 < framework < 2
    assert 0 < whole < 2
    assert cdrh3 > 0
    assert backbone > 0


def test_region_rmsd_needs_the_same_scheme():
    entry = GOLDEN_ANGLES[0]
    with pytest.raises(ValueError, match="different schemes"):
        region_rmsd(
            _legacy_fv(entry, "imgt"), _legacy_fv(entry, "chothia"), definition="imgt"
        )


def test_region_rmsd_needs_common_atoms(fv):
    with pytest.raises(ValueError, match="no selected atoms in common"):
        region_rmsd(fv.heavy, fv.light, definition="imgt")


def test_orientation_rmsd(fv):
    assert orientation_rmsd(fv, fv, definition="imgt") == pytest.approx(0, abs=1e-6)
    moved = _moved(fv, _rotation([1, 0, 0], 90), np.array([0.0, 5.0, 0.0]))
    assert orientation_rmsd(fv, moved, definition="imgt") == pytest.approx(0, abs=1e-3)
    pm = _legacy_fv(GOLDEN_ANGLES[1], "imgt")
    there, back = (
        orientation_rmsd(fv, pm, definition="imgt"),
        orientation_rmsd(pm, fv, definition="imgt"),
    )
    # The two copies in 12e8 differ by about 8.6 degrees in HL.
    assert 1 < there < 5
    assert there == pytest.approx(back, rel=0.2)


# --- ABangle -----------------------------------------------------------------


@pytest.mark.parametrize("scheme", ["chothia", "imgt"])
@pytest.mark.parametrize(
    "entry", GOLDEN_ANGLES, ids=lambda e: f"{e['pdb']}_{e['heavy']}{e['light']}"
)
def test_abangle_matches_legacy(entry, scheme):
    angles = abangle(_legacy_fv(entry, scheme)).as_dict()
    assert angles == pytest.approx(entry["angles"], abs=1e-4)


def test_abangle_is_unchanged_by_moving_the_fv(fv):
    moved = _moved(fv, _rotation([3, -1, 2], 123), np.array([7.0, 8.0, 9.0]))
    assert abangle(moved).as_dict() == pytest.approx(abangle(fv).as_dict())


def test_abangle_names():
    angles = OrientationAngles(hl=1, hc1=2, lc1=3, hc2=4, lc2=5, dc=6)
    assert angles.as_dict() == {
        "HL": 1,
        "HC1": 2,
        "LC1": 3,
        "HC2": 4,
        "LC2": 5,
        "dc": 6,
    }


def test_abangle_needs_chothia_or_imgt():
    entry = GOLDEN_ANGLES[0]
    path = LEGACY / "structures" / f"{entry['files']['chothia']}.gz"
    structure = read_structure(path, "kabat", chains={"H": "H", "L": "L"})
    with pytest.raises(ValueError, match="Chothia or IMGT"):
        abangle(Fv(structure["H"], structure["L"]))


def test_abangle_needs_core_positions(fv, tmp_path):
    path = tmp_path / "cdrs.pdb"
    select(fv, ["cdrs"], definition="imgt").write_pdb(str(path))
    structure = read_structure(path, "imgt", chains={"H": "H", "L": "L"})
    with pytest.raises(ValueError, match="ABangle core positions"):
        abangle(Fv(structure["H"], structure["L"]))


def test_abangle_on_sabdab_files_matches_legacy():
    with open(DATA / "sabdab" / "SAbDab_summary.csv", newline="") as f:
        rows = {(r["PDB"][-4:], r["Hchain"], r["Lchain"]): r for r in csv.DictReader(f)}
    for entry in GOLDEN_ANGLES:
        heavy, light = entry["heavy"], entry["light"]
        row = rows[entry["pdb"], heavy, light]
        [path] = (DATA / "sabdab").glob(f"{row['PDB']}_{heavy}_{light}_*.cif.gz")
        structure = read_structure(path, "imgt", chains={heavy: "H", light: "L"})
        angles = abangle(Fv(structure[heavy], structure[light])).as_dict()
        assert angles == pytest.approx(entry["angles"], abs=1e-4)
