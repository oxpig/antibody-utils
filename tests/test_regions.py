import json
from pathlib import Path

import pytest

from antibody_utils._data import load_toml
from antibody_utils.numbering import Chain, Position, Scheme
from antibody_utils.regions import (
    REGION_GROUPS,
    Definition,
    Region,
    RegionSelector,
    annotate_regions,
    get_region,
)

LEGACY = Path(__file__).parent / "data" / "legacy"


def _legacy(name):
    return json.loads((LEGACY / name).read_text())


# --- Data integrity --------------------------------------------------------


@pytest.mark.parametrize("definition", list(Definition))
@pytest.mark.parametrize("chain", list(Chain))
def test_definition_tiles_imgt_positions(definition, chain):
    boundaries = load_toml("regions.toml")[definition][chain]
    c = chain.lower()
    expected_order = [
        f"fw{c}1",
        f"cdr{c}1",
        f"fw{c}2",
        f"cdr{c}2",
        f"fw{c}3",
        f"cdr{c}3",
        f"fw{c}4",
    ]
    assert list(boundaries) == expected_order
    next_start = 1
    for start, end in boundaries.values():
        assert start == next_start
        assert end >= start
        next_start = end + 1
    assert next_start == 129


def test_imgt_definition_is_the_imgt_cdrs():
    for chain in Chain:
        c = chain.lower()
        cdrs = {
            name: tuple(bounds)
            for name, bounds in load_toml("regions.toml")["imgt"][chain].items()
            if name.startswith("cdr")
        }
        assert cdrs == {
            f"cdr{c}1": (27, 38),
            f"cdr{c}2": (56, 65),
            f"cdr{c}3": (105, 117),
        }


# --- Regions and groups ----------------------------------------------------


def test_region_properties():
    assert Region.FWL4.chain is Chain.LIGHT
    assert not Region.FWL4.is_cdr
    assert {region for region in Region if region.is_cdr} == REGION_GROUPS["cdrs"]


def test_region_groups():
    assert REGION_GROUPS["fv"] == set(Region)
    assert REGION_GROUPS["vh"] == {region for region in Region if region.chain == "H"}
    assert REGION_GROUPS["framework"].isdisjoint(REGION_GROUPS["cdrs"])


# --- get_region ------------------------------------------------------------


@pytest.mark.parametrize(
    ("position", "chain", "scheme", "definition", "expected"),
    [
        ("27", "H", "imgt", "imgt", Region.CDRH1),
        ("26", "H", "imgt", "imgt", Region.FWH1),
        ("111A", "H", "imgt", "imgt", Region.CDRH3),
        ("35A", "H", "kabat", "kabat", Region.CDRH1),
        ("35B", "H", "kabat", "kabat", Region.CDRH1),
        ("35C", "H", "kabat", "kabat", Region.FWH2),
        ("33", "H", "kabat", "chothia", Region.FWH2),
        ("95", "H", "chothia", "chothia", Region.CDRH3),
        ("129", "H", "imgt", "imgt", None),
        ("114", "H", "chothia", "chothia", None),
    ],
)
def test_get_region(position, chain, scheme, definition, expected):
    assert get_region(position, chain, scheme=scheme, definition=definition) == expected


def test_get_region_accepts_enums_and_positions():
    assert (
        get_region(
            Position(100, "A"),
            Chain.HEAVY,
            scheme=Scheme.CHOTHIA,
            definition=Definition.CHOTHIA,
        )
        is Region.CDRH3
    )


def test_get_region_rejects_unknown_definition():
    with pytest.raises(ValueError, match="not a valid Definition"):
        get_region("27", "H", scheme="imgt", definition="wolfguy")


def _get_region_cases():
    fixture = _legacy("get_region.json")
    low, high = fixture["index_range"]
    for key, by_insertion in fixture["results"].items():
        scheme, chain, definition = key.split("/")
        for insertion, regions in by_insertion.items():
            for number, legacy in zip(range(low, high + 1), regions, strict=True):
                yield scheme, chain, definition, Position(number, insertion), legacy


def test_get_region_matches_legacy():
    """Every residue number and insertion, in every scheme and definition."""
    mismatches = []
    count = 0
    for scheme, chain, definition, position, legacy in _get_region_cases():
        expected = None if legacy == "?" else Region(legacy)
        got = get_region(position, chain, scheme=scheme, definition=definition)
        count += 1
        if got != expected:
            mismatches.append((scheme, chain, definition, str(position), legacy, got))
    assert count == 4 * 2 * 5 * 4 * 131
    assert not mismatches, mismatches[:10]


# --- annotate_regions ------------------------------------------------------


def test_annotate_regions_matches_legacy():
    """Real and synthetic sequences, in every scheme and definition."""
    sequences = _legacy("numbered_sequences.json")
    cases = _legacy("annotate_regions.json")
    assert len(cases) == len(sequences) * 4 * 5
    mismatches = []
    for case in cases:
        numbering = [
            (Position(number, insertion), residue)
            for number, insertion, residue in sequences[case["name"]][case["scheme"]]
        ]
        annotated = annotate_regions(
            numbering,
            case["chain"],
            scheme=case["scheme"],
            definition=case["definition"],
        )
        got = [str(region) if region else "" for *_, region in annotated]
        if got != case["regions"]:
            mismatches.append((case["name"], case["scheme"], case["definition"]))
    assert not mismatches, mismatches


def test_annotate_regions_exercises_cdrh1_heuristics():
    """The fixtures include CDR-H1s that need each heuristic to be applied."""
    sequences = _legacy("numbered_sequences.json")
    long_h1 = sequences["vh4_39_longh1"]
    chothia_h31_insertions = [i for n, i, _ in long_h1["chothia"] if n == 31 and i]
    kabat_h35_insertions = [i for n, i, _ in long_h1["kabat"] if n == 35 and i]
    assert len(chothia_h31_insertions) > 2
    assert len(kabat_h35_insertions) > 2


def test_annotate_regions_heavy_imgt_kabat():
    """Kabat's CDR-H1 is the same residues whether numbered in Kabat or IMGT."""
    sequences = _legacy("numbered_sequences.json")
    for name, entry in sequences.items():
        if entry["chain_type"] != "H":
            continue
        cdrh1 = {}
        for scheme in ("kabat", "imgt"):
            numbering = [(Position(n, i), aa) for n, i, aa in entry[scheme]]
            annotated = annotate_regions(
                numbering, "H", scheme=scheme, definition="kabat"
            )
            cdrh1[scheme] = "".join(
                aa for _, aa, region in annotated if region is Region.CDRH1
            )
        assert cdrh1["kabat"] == cdrh1["imgt"], name


def _regions(numbering, chain, scheme, definition):
    annotated = annotate_regions(numbering, chain, scheme=scheme, definition=definition)
    return [region for *_, region in annotated]


def test_annotate_regions_beyond_variable_domain():
    numbering = [("1", "E"), ("110", "D"), ("128", "S"), ("129", "A")]
    assert _regions(numbering, "H", "imgt", "imgt") == [
        Region.FWH1,
        Region.CDRH3,
        Region.FWH4,
        None,
    ]


def test_annotate_regions_accepts_gaps():
    numbering = [("1", "D"), ("27", "-"), ("28", "S"), ("39", "W")]
    assert _regions(numbering, "L", "imgt", "imgt") == [
        Region.FWL1,
        Region.CDRL1,
        Region.CDRL1,
        Region.FWL2,
    ]


def test_annotate_regions_frameworks_follow_cdrs():
    """A fragment starting after CDR 1 is labelled as if it were framework 1."""
    numbering = [("94", "R"), ("95", "D"), ("103", "W")]
    assert _regions(numbering, "H", "chothia", "chothia") == [
        Region.FWH1,
        Region.CDRH3,
        Region.FWH4,
    ]


# --- RegionSelector --------------------------------------------------------


def test_selector_groups_and_regions():
    selector = RegionSelector(["framework", "cdrh3"], scheme="imgt", definition="imgt")
    assert selector.regions == REGION_GROUPS["framework"] | {Region.CDRH3}
    assert selector.accepts("110", "H")
    assert not selector.accepts("110", "L")
    assert selector.accepts("1", "L")


def test_selector_extra_and_excluded_positions():
    selector = RegionSelector(scheme="chothia", definition="chothia")
    assert not selector.accepts("1", "H")
    selector.add_positions(["1", Position(2)], "H")
    assert selector.accepts("1", "H")
    assert not selector.accepts("1", "L")
    selector.add_regions(["vh"])
    selector.exclude_positions(["2"], "H")
    assert not selector.accepts("2", "H")
    assert selector.accepts("3", "H")


def test_selector_rejects_unknown_region():
    with pytest.raises(ValueError, match="Unknown region"):
        RegionSelector(["cdrh4"], scheme="imgt", definition="imgt")
