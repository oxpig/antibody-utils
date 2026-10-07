"""Consistency with the IMGT-numbered structure files of the current SAbDab."""

import csv
import json
from pathlib import Path

import pytest

from antibody_utils.structure import Fv, read_structure

DATA = Path(__file__).parent / "data"
SABDAB = DATA / "sabdab"

with (SABDAB / "SAbDab_summary.csv").open(newline="") as _f:
    SUMMARY = list(csv.DictReader(_f))

# SAbDab numbers the residues after the variable domain by counting on from
# its last position, so the first constant-domain residue of a light chain
# whose domain ends at L127 is numbered L128, a valid IMGT position, and is
# read as part of the domain.  To be fixed in SAbDab's structure renumbering.
FLANKING_RESIDUE = pytest.mark.xfail(
    reason="SAbDab numbers the first constant-domain residue as IMGT L128",
    strict=True,
)


def _cases(light_marks=()):
    return [
        pytest.param(row, chain, id=f"{row['INSTANCE']}-{chain}", marks=marks)
        for row in SUMMARY
        for chain, marks in (("H", ()), ("L", light_marks))
    ]


def _instance(row):
    pdb, heavy, light = row["PDB"], row["Hchain"], row["Lchain"]
    [path] = SABDAB.glob(f"{pdb}_{heavy}_{light}_*.cif.gz")
    structure = read_structure(path, "imgt", chains={heavy: "H", light: "L"})
    return Fv(structure[heavy], structure[light])


def _domain(row, chain):
    fv = _instance(row)
    return fv.heavy if chain == "H" else fv.light


@pytest.fixture(scope="module")
def legacy_sequences():
    return json.loads((DATA / "legacy" / "numbered_sequences.json").read_text())


@pytest.mark.parametrize(("row", "chain"), _cases(FLANKING_RESIDUE))
def test_numbering_matches_legacy(legacy_sequences, row, chain):
    domain = _domain(row, chain)
    name = f"{row['PDB'][-4:]}_{domain.name}"
    expected = [[n, i.strip(), r] for n, i, r in legacy_sequences[name]["imgt"]]
    assert [[p.number, p.insertion, r] for p, r in domain.numbered] == expected


@pytest.mark.parametrize(("row", "chain"), _cases(FLANKING_RESIDUE))
def test_sequence_matches_summary(row, chain):
    domain = _domain(row, chain)
    assert domain.numbered.sequence == row[f"V{chain}"]


@pytest.mark.parametrize(("row", "chain"), _cases())
def test_cdrs_match_summary(row, chain):
    sequences = {}
    for _, residue, region in _domain(row, chain).numbered.annotate_regions("imgt"):
        sequences[region] = sequences.get(region, "") + residue
    for n in (1, 2, 3):
        assert sequences[f"cdr{chain.lower()}{n}"] == row[f"CDR-{chain}{n}"]


@pytest.mark.parametrize("row", SUMMARY[:2], ids=lambda row: row["INSTANCE"])
def test_full_structure_matches_instances(row):
    heavy, light = row["Hchain"], row["Lchain"]
    full = read_structure(
        SABDAB / "pdb_000012e8_sabdab.cif.gz",
        "imgt",
        chains={"H": "H", "L": "L", "P": "H", "M": "L"},
    )
    instance = _instance(row)
    for name, domain in ((heavy, instance.heavy), (light, instance.light)):
        assert list(full[name].numbered) == list(domain.numbered)
