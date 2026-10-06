import json
import math
from pathlib import Path

import pytest

from antibody_utils.numbering import NumberedSequence, Position
from antibody_utils.regions import Region
from antibody_utils.sequence import (
    cdr_lengths,
    cdr_sequences,
    identity,
    region_sequences,
    similarity,
)

LEGACY = Path(__file__).parent / "data" / "legacy"

TRASTUZUMAB_IMGT_CDRS = {
    Region.CDRH1: "GFNIKDTY",
    Region.CDRH2: "IYPTNGYT",
    Region.CDRH3: "SRWGGDGFYAMDY",
}


def _domain(name, entry, scheme):
    numbering = entry[scheme]
    return NumberedSequence(
        name=name,
        chain_type=entry["chain_type"],
        scheme=scheme,
        positions=tuple(Position(n, i) for n, i, _ in numbering),
        sequence="".join(residue for *_, residue in numbering),
        score=0.0,
        start=0,
        end=len(numbering),
    )


@pytest.fixture(scope="module")
def legacy_sequences():
    with open(LEGACY / "numbered_sequences.json") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def trastuzumab(legacy_sequences):
    return {
        chain: _domain(name, legacy_sequences[name], "imgt")
        for chain, name in (("H", "trast_H"), ("L", "trast_L"))
    }


# --- Regions -----------------------------------------------------------------


def test_cdr_sequences(trastuzumab):
    assert cdr_sequences(trastuzumab["H"], "imgt") == TRASTUZUMAB_IMGT_CDRS


def test_cdr_lengths(trastuzumab):
    assert cdr_lengths(trastuzumab["H"], "imgt") == {
        cdr: len(seq) for cdr, seq in TRASTUZUMAB_IMGT_CDRS.items()
    }


def test_region_sequences_reassemble_the_domain(trastuzumab):
    for domain in trastuzumab.values():
        regions = region_sequences(domain, "north")
        assert len(regions) == 7
        assert "".join(regions.values()) == domain.sequence


def test_missing_cdr_is_empty():
    truncated = NumberedSequence(
        name="x",
        chain_type="K",
        scheme="imgt",
        positions=(Position(1), Position(27)),
        sequence="DQ",
        score=0.0,
        start=0,
        end=2,
    )
    assert cdr_sequences(truncated, "imgt") == {
        Region.CDRL1: "Q",
        Region.CDRL2: "",
        Region.CDRL3: "",
    }


# --- Identity and similarity: behaviour --------------------------------------


def test_identical_sequences(trastuzumab):
    pair = (trastuzumab["H"], trastuzumab["L"])
    assert identity(pair, pair, definition="imgt") == 1.0
    assert identity(pair, pair, definition="imgt", regions=["cdrh3"]) == 1.0


def test_only_like_chains_are_compared(trastuzumab):
    heavy, light = trastuzumab["H"], trastuzumab["L"]
    # No common chain: nothing to compare.
    assert identity(heavy, light, definition="imgt") == 0.0
    # A VH/VL pair against the VH alone compares just the heavy chains.
    assert identity((heavy, light), heavy, definition="imgt") == 1.0


def test_pooled_over_chains(legacy_sequences):
    d = {
        name: _domain(name, legacy_sequences[name], "imgt") for name in legacy_sequences
    }
    pooled = identity(
        (d["trast_H"], d["trast_L"]), (d["vh4_39"], d["lambda_L"]), definition="imgt"
    )
    heavy = identity(d["trast_H"], d["vh4_39"], definition="imgt")
    light = identity(d["trast_L"], d["lambda_L"], definition="imgt")
    assert min(heavy, light) < pooled < max(heavy, light)


def test_similarity_normalised_is_mean_score(trastuzumab):
    heavy = trastuzumab["H"]
    total = similarity(heavy, heavy, definition="imgt")
    mean = similarity(heavy, heavy, definition="imgt", normalise=True)
    assert math.isclose(total / len(heavy), mean)


def test_two_domains_of_one_chain_rejected(trastuzumab):
    heavy = trastuzumab["H"]
    with pytest.raises(ValueError, match="More than one heavy chain"):
        identity((heavy, heavy), heavy, definition="imgt")


def test_different_schemes_rejected(legacy_sequences):
    imgt = _domain("trast_H", legacy_sequences["trast_H"], "imgt")
    kabat = _domain("trast_H", legacy_sequences["trast_H"], "kabat")
    with pytest.raises(ValueError, match="different schemes"):
        identity(imgt, kabat, definition="kabat")


def test_unscored_residues_are_skipped():
    def domain(sequence):
        return NumberedSequence(
            name="x",
            chain_type="H",
            scheme="imgt",
            positions=(Position(1), Position(2)),
            sequence=sequence,
            score=0.0,
            start=0,
            end=2,
        )

    # "U" (selenocysteine) has no BLOSUM62 score; W/W scores 11.
    assert similarity(domain("WU"), domain("WU"), definition="imgt") == 11.0


# --- Identity and similarity: legacy equivalence -----------------------------


def test_matches_legacy(legacy_sequences):
    """Identity and similarity for every legacy fixture pair and region."""
    with open(LEGACY / "identity_similarity.json") as f:
        cases = json.load(f)
    mismatches = []
    for case in cases:
        a = _domain(case["a"], legacy_sequences[case["a"]], case["scheme"])
        b = _domain(case["b"], legacy_sequences[case["b"]], case["scheme"])
        kwargs = {"definition": case["definition"], "regions": case["region"]}
        got = {
            "identity": identity(a, b, **kwargs),
            "similarity": similarity(a, b, **kwargs),
            "similarity_normalised": similarity(a, b, normalise=True, **kwargs),
        }
        for key, value in got.items():
            # Legacy similarity(region=None) raised NameError, so has no value;
            # here it covers the whole variable domain, as identity does.
            expected = case[key]
            if expected is None:
                expected = got[key] if case["region"] is None else None
            if expected is None or not math.isclose(value, expected):
                mismatches.append((case["a"], case["b"], case["region"], key, value))
    assert not mismatches, mismatches[:10]


def test_similarity_without_regions_is_the_whole_fv(trastuzumab, legacy_sequences):
    other = _domain("vh4_39", legacy_sequences["vh4_39"], "imgt")
    heavy = trastuzumab["H"]
    assert similarity(heavy, other, definition="imgt") == similarity(
        heavy, other, definition="imgt", regions=["fv"]
    )
