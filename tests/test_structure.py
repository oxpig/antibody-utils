import json
from importlib.util import find_spec
from pathlib import Path

import gemmi
import pytest

from antibody_utils.numbering import Chain, Position, Scheme
from antibody_utils.regions import Region
from antibody_utils.structure import (
    BACKBONE_ATOMS,
    AntibodyStructure,
    Fv,
    from_gemmi,
    read_structure,
    select,
)

DATA = Path(__file__).parent / "data"
LEGACY = DATA / "legacy"

# Fv files written by the legacy fixture script: (entry, heavy, light).
LEGACY_FVS = [
    ("12e8", "H", "L"),
    ("12e8", "P", "M"),
    ("1ahw", "B", "A"),
    ("1ahw", "E", "D"),
    ("3hfm", "H", "L"),
    ("7fab", "H", "L"),
]

needs_anarcii = pytest.mark.skipif(
    find_spec("anarcii") is None, reason="needs the `numbering` extra (ANARCII)"
)


@pytest.fixture(scope="module")
def legacy_sequences():
    return json.loads((LEGACY / "numbered_sequences.json").read_text())


def _legacy_fv(entry, heavy, light, scheme):
    path = LEGACY / "structures" / f"{entry}_{heavy}{light}_{scheme}.pdb.gz"
    return read_structure(path, scheme, chains={heavy: "H", light: "L"})


@pytest.fixture
def structure():
    return _legacy_fv("12e8", "H", "L", "imgt")


@pytest.fixture
def fv(structure):
    return Fv(structure["H"], structure["L"])


def _numbering(chain):
    return [[p.number, p.insertion, r] for p, r in chain.numbered]


def _expected(legacy_sequences, name, scheme):
    return [[n, i.strip(), r] for n, i, r in legacy_sequences[name][scheme]]


def _structure(chains):
    """Build a one-model structure from `{chain name: [(position, name), ...]}`.

    Positions are numbers or strings such as `"111A"`.  Each chain needs at
    least two residues for GEMMI to treat it as a polymer.
    """
    structure = gemmi.Structure()
    structure.name = "test"
    model = gemmi.Model(1)
    for name, residues in chains.items():
        chain = gemmi.Chain(name)
        for raw_position, residue_name in residues:
            if isinstance(raw_position, int):
                position = Position(raw_position)
            else:
                position = Position.parse(raw_position)
            residue = gemmi.Residue()
            residue.name = residue_name
            residue.seqid = gemmi.SeqId(position.number, position.insertion or " ")
            atom = gemmi.Atom()
            atom.name = "O" if residue_name == "HOH" else "CA"
            residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    return structure


# --- Reading numbering from the file -----------------------------------------


@pytest.mark.parametrize("scheme", ["imgt", "chothia"])
@pytest.mark.parametrize(("entry", "heavy", "light"), LEGACY_FVS)
def test_file_numbering_matches_legacy(legacy_sequences, entry, heavy, light, scheme):
    structure = _legacy_fv(entry, heavy, light, scheme)
    assert [chain.name for chain in structure.heavy_chains] == [heavy]
    assert [chain.name for chain in structure.light_chains] == [light]
    for chain in structure:
        expected = _expected(legacy_sequences, f"{entry}_{chain.name}", scheme)
        assert _numbering(chain) == expected


def test_file_numbering_uses_only_the_variable_domain():
    residues = [(n, "ALA") for n in (-1, 0, 1, 2, 3, 128, 129)]
    structure = from_gemmi(_structure({"H": residues}), chains={"H": "H"})
    chain = structure["H"]
    assert [str(p) for p in chain.numbered.positions] == ["1", "2", "3", "128"]
    assert (chain.numbered.start, chain.numbered.end) == (2, 6)


def test_file_numbering_unknown_chain():
    with pytest.raises(KeyError, match="No chain 'X'"):
        from_gemmi(_structure({"H": [(1, "ALA")]}), chains={"X": "H"})


def test_file_numbering_without_a_domain():
    with pytest.raises(ValueError, match="no residues numbered"):
        from_gemmi(_structure({"H": [(200, "ALA")]}), chains={"H": "H"})


def test_file_numbering_out_of_order():
    residues = [(1, "ALA"), (200, "ALA"), (2, "ALA")]
    with pytest.raises(ValueError, match="outside the variable domain"):
        from_gemmi(_structure({"H": residues}), chains={"H": "H"})


def test_file_numbering_repeated_position():
    residues = [(1, "ALA"), (2, "ALA"), (1, "GLY")]
    with pytest.raises(ValueError, match="more than one residue per position"):
        from_gemmi(_structure({"H": residues}), chains={"H": "H"})


def test_sequence_skips_water_and_marks_unknown_residues():
    residues = [(1, "ALA"), (2, "ZZZ"), (3, "MSE"), (4, "HOH")]
    structure = from_gemmi(_structure({"H": residues}), chains={"H": "H"})
    assert structure["H"].numbered.sequence == "AXM"


def test_from_gemmi_copies_the_structure():
    original = _structure({"H": [(1, "ALA"), (2, "GLY")]})
    structure = from_gemmi(original, chains={"H": "H"})
    original[0]["H"][0].name = "TRP"
    assert structure["H"]["1"].name == "ALA"


# --- Antibody chains and structures ------------------------------------------


def test_chain(structure):
    heavy = structure["H"]
    assert heavy.chain is Chain.HEAVY
    assert heavy.scheme is Scheme.IMGT
    assert len(heavy) == len(heavy.numbered) == len(heavy.residues)
    assert heavy.gemmi_chain.name == "H"
    assert (
        repr(heavy)
        == "<AntibodyChain H heavy, 120 residues numbered in the imgt scheme>"
    )


def test_chain_residues_follow_the_numbering(structure):
    heavy = structure["H"]
    one_letter = "".join(
        gemmi.find_tabulated_residue(residue.name).one_letter_code.upper()
        for residue in heavy.residues
    )
    assert one_letter == heavy.numbered.sequence
    assert [position for position, _ in heavy] == list(heavy.numbered.positions)


def test_chain_lookup_by_position(structure):
    heavy = structure["H"]
    assert heavy[Position(23)].name == "CYS"
    assert heavy["23"].seqid == gemmi.SeqId(23, " ")
    assert "23" in heavy
    assert Position(23) in heavy
    assert "200" not in heavy
    assert 23 not in heavy
    with pytest.raises(KeyError):
        heavy["200"]


def test_chain_lookup_with_insertions():
    residues = [(111, "ALA"), ("111A", "GLY"), ("112A", "SER"), (112, "THR")]
    structure = from_gemmi(_structure({"H": residues}), chains={"H": "H"})
    assert structure["H"]["111A"].name == "GLY"
    assert structure["H"][Position(112, "A")].name == "SER"
    assert "111B" not in structure["H"]


def test_region_residues(structure):
    heavy = structure["H"]
    grouped = heavy.region_residues("imgt")
    assert list(grouped) == [
        Region.FWH1,
        Region.CDRH1,
        Region.FWH2,
        Region.CDRH2,
        Region.FWH3,
        Region.CDRH3,
        Region.FWH4,
    ]
    assert sum(len(residues) for residues in grouped.values()) == len(heavy)
    assert [r.seqid.num for r in grouped[Region.CDRH3]][:2] == [105, 106]


def test_structure(structure):
    assert len(structure) == 2
    assert [chain.name for chain in structure] == ["L", "H"]
    assert structure.scheme is Scheme.IMGT
    assert structure.other_chains == []
    assert len(structure.model) == 2
    with pytest.raises(KeyError):
        structure["X"]


def test_other_chains():
    structure = from_gemmi(
        _structure({"H": [(1, "ALA"), (2, "ALA")], "A": [(1, "GLY"), (2, "GLY")]}),
        chains={"H": "H"},
    )
    assert [chain.name for chain in structure.other_chains] == ["A"]


# --- Fvs ---------------------------------------------------------------------


def test_fv(structure, fv):
    assert fv.name == "HL"
    assert fv.scheme is Scheme.IMGT
    assert list(fv) == [structure["H"], structure["L"]]
    assert fv.numbered == (structure["H"].numbered, structure["L"].numbered)
    assert repr(fv) == "<Fv HL heavy=H light=L>"
    regions = fv.region_residues("north")
    assert len(regions) == 14
    assert Region.CDRH3 in regions
    assert Region.CDRL3 in regions


def test_fv_chain_types(structure):
    with pytest.raises(ValueError, match="not a heavy chain"):
        Fv(structure["L"], structure["H"])
    with pytest.raises(ValueError, match="not a light chain"):
        Fv(structure["H"], structure["H"])


def test_fv_from_different_structures(structure):
    other = _legacy_fv("12e8", "H", "L", "imgt")
    with pytest.raises(ValueError, match="different structures"):
        Fv(structure["H"], other["L"])


def test_fv_cannot_be_repaired(fv, structure):
    with pytest.raises(AttributeError):
        fv.light = structure["H"]
    with pytest.raises(AttributeError):
        structure["H"].numbered = structure["L"].numbered


def test_structure_rejects_chains_in_another_scheme(structure):
    with pytest.raises(ValueError, match="numbered in imgt, not chothia"):
        AntibodyStructure(structure.structure, [structure["H"]], Scheme.CHOTHIA, 0)


def test_fv_different_schemes(structure):
    chothia = _legacy_fv("12e8", "H", "L", "chothia")
    with pytest.raises(ValueError, match="different schemes"):
        Fv(structure["H"], chothia["L"])


# --- Selection ---------------------------------------------------------------


def test_select_fv(fv):
    copy = select(fv, definition="imgt")
    assert [chain.name for chain in copy[0]] == ["H", "L"]
    for antibody_chain, chain in zip(fv, copy[0], strict=True):
        assert [Position(r.seqid.num, r.seqid.icode) for r in chain] == list(
            antibody_chain.numbered.positions
        )
    assert copy.cell.a == pytest.approx(fv.heavy._structure.structure.cell.a)


def test_select_renumbers_in_the_scheme():
    # Renumber the file's residues from 101; the copy keeps the chain's numbering.
    structure = from_gemmi(
        _structure({"H": [(1, "ALA"), (2, "GLY")]}), chains={"H": "H"}
    )
    for residue in structure["H"].residues:
        residue.seqid = gemmi.SeqId(residue.seqid.num + 100, " ")
    copy = select(structure["H"], definition="imgt")
    assert [r.seqid.num for r in copy[0]["H"]] == [1, 2]


def test_select_regions_and_exclusions(fv):
    cdrh3 = select(fv, ["cdrh3"], definition="imgt")
    assert [chain.name for chain in cdrh3[0]] == ["H"]
    assert len(cdrh3[0]["H"]) == len(fv.heavy.region_residues("imgt")[Region.CDRH3])

    without = select(fv, exclude=["cdrh3"], definition="imgt")
    whole = select(fv, definition="imgt")
    assert without[0].count_atom_sites() + cdrh3[0].count_atom_sites() == (
        whole[0].count_atom_sites()
    )


def test_select_atoms(fv):
    copy = select(fv.heavy, ["cdrh1"], definition="imgt", atoms=BACKBONE_ATOMS)
    names = {atom.name for residue in copy[0]["H"] for atom in residue}
    assert names == BACKBONE_ATOMS


def test_select_leaves_the_original_unchanged(structure, fv):
    before = structure.structure.make_pdb_string()
    select(fv, ["cdrs"], definition="imgt", atoms=["CA"])
    assert structure.structure.make_pdb_string() == before


def test_select_round_trip(fv, tmp_path):
    path = tmp_path / "fv.pdb"
    select(fv, definition="imgt").write_pdb(str(path))
    again = read_structure(path, "imgt", chains={"H": "H", "L": "L"})
    for original, reread in zip(fv, Fv(again["H"], again["L"]), strict=True):
        assert _numbering(reread) == _numbering(original)


def test_select_unknown_region(fv):
    with pytest.raises(ValueError, match="Unknown region"):
        select(fv, ["cdrh4"], definition="imgt")


def test_select_from_different_structures(structure):
    other = _legacy_fv("12e8", "H", "L", "imgt")
    with pytest.raises(ValueError, match="different structures"):
        select([structure["H"], other["L"]], definition="imgt")


# --- Numbering with ANARCII --------------------------------------------------


@pytest.mark.numbering
@needs_anarcii
@pytest.mark.parametrize("scheme", ["imgt", "kabat", "chothia", "martin"])
def test_anarcii_numbering_matches_legacy(legacy_sequences, scheme):
    structure = read_structure(DATA / "pdb" / "12e8.pdb.gz", scheme, cpu=True)
    assert [chain.name for chain in structure] == ["L", "H", "M", "P"]
    assert [chain.name for chain in structure.heavy_chains] == ["H", "P"]
    for chain in structure:
        expected = _expected(legacy_sequences, f"12e8_{chain.name}", scheme)
        assert _numbering(chain) == expected
        one_letter = "".join(
            gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
            for r in chain.residues
        )
        assert one_letter == chain.numbered.sequence
