import json
from importlib.util import find_spec
from pathlib import Path

import pytest

from antibody_utils._data import load_toml
from antibody_utils.liabilities import Liability, find_liabilities
from antibody_utils.numbering import Chain, NumberedSequence, Position
from antibody_utils.regions import REGION_GROUPS

LEGACY = Path(__file__).parent / "data" / "legacy"


def _domain(name, chain_type, numbering, scheme="imgt"):
    return NumberedSequence(
        name=name,
        chain_type=chain_type,
        scheme=scheme,
        positions=tuple(Position(n, i) for n, i, _ in numbering),
        sequence="".join(residue for *_, residue in numbering),
        score=0.0,
        start=0,
        end=len(numbering),
    )


@pytest.fixture(scope="module")
def legacy_cases():
    with open(LEGACY / "liabilities.json") as f:
        return json.load(f)


def _key(liability: Liability):
    return (
        liability.name,
        liability.chain.value,
        liability.start,
        liability.end - 1,
        liability.positions[0],
        liability.positions[-1],
    )


def test_matches_legacy(legacy_cases):
    for case in legacy_cases:
        heavy = _domain("H", "H", case["heavy_numbering"])
        light = _domain("L", "K", case["light_numbering"])
        found = find_liabilities(heavy, light)
        expected = {
            (
                name,
                chain,
                start,
                end,
                Position.from_tuple(first),
                Position.from_tuple(last),
            )
            for name, start, end, first, last, chain, _ in case["found"]
        }
        assert {_key(liability) for liability in found} == expected, case["name"]
        assert len(found) == len(case["found"])


def test_liability_contents(legacy_cases):
    case = next(c for c in legacy_cases if c["name"] == "motifs")
    heavy = _domain("H", "H", case["heavy_numbering"])
    glycosylation = [
        liability
        for liability in find_liabilities(heavy, None)
        if liability.name.startswith("N-linked")
    ]
    assert len(glycosylation) == 1
    (site,) = glycosylation
    assert site.chain is Chain.HEAVY
    assert site.residues == heavy.sequence[site.start : site.end] == "NSS"
    assert site.positions == heavy.positions[site.start : site.end]


def test_order_follows_the_data_file(legacy_cases):
    case = next(c for c in legacy_cases if c["name"] == "motifs")
    found = find_liabilities(
        _domain("H", "H", case["heavy_numbering"]),
        _domain("L", "K", case["light_numbering"]),
    )
    names = [entry["name"] for entry in load_toml("liabilities.toml")["liability"]]
    keys = [(names.index(f.name), f.chain != Chain.HEAVY, f.start) for f in found]
    assert keys == sorted(keys)


def test_n_terminal_glutamate_needs_both_chains(legacy_cases):
    case = next(c for c in legacy_cases if c["name"] == "motifs")
    heavy = _domain("H", "H", case["heavy_numbering"])
    light = _domain("L", "K", case["light_numbering"])
    name = "N-terminal glutamate (VH and VL) (E)"
    both = [f for f in find_liabilities(heavy, light) if f.name == name]
    assert {f.chain for f in both} == {Chain.HEAVY, Chain.LIGHT}
    assert not [f for f in find_liabilities(heavy) if f.name == name]


def test_conserved_cysteines_are_ignored(legacy_cases):
    case = next(c for c in legacy_cases if c["name"] == "trastuzumab")
    heavy = _domain("H", "H", case["heavy_numbering"])
    assert "C" in heavy.sequence
    assert not [f for f in find_liabilities(heavy) if f.name.startswith("Unpaired")]


def test_rejects_wrong_chain_and_scheme(legacy_cases):
    case = legacy_cases[0]
    heavy = _domain("H", "H", case["heavy_numbering"])
    with pytest.raises(ValueError, match="Expected a light chain"):
        find_liabilities(light=heavy)
    kabat = _domain("H", "H", case["heavy_numbering"], scheme="kabat")
    with pytest.raises(ValueError, match="IMGT-numbered"):
        find_liabilities(kabat)


def test_position_sets_and_region_groups_have_distinct_names():
    positions = load_toml("liabilities.toml")["positions"]
    assert not set(REGION_GROUPS) & set(positions)


@pytest.mark.numbering
@pytest.mark.skipif(find_spec("anarcii") is None, reason="needs the `numbering` extra")
def test_numbers_raw_sequences(legacy_cases):
    case = next(c for c in legacy_cases if c["name"] == "trastuzumab")
    from_sequences = find_liabilities(case["heavy"], case["light"], cpu=True)
    from_numbering = find_liabilities(
        _domain("H", "H", case["heavy_numbering"]),
        _domain("L", "K", case["light_numbering"]),
    )
    assert [_key(f) for f in from_sequences] == [_key(f) for f in from_numbering]
