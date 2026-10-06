import dataclasses
import json
import logging
import subprocess
import sys
from importlib.util import find_spec
from itertools import islice
from pathlib import Path

import pytest

from antibody_utils.numbering import (
    Chain,
    NumberedSequence,
    NumberingError,
    Position,
    Scheme,
    iter_number,
    number,
    number_sequence,
)
from antibody_utils.numbering import anarcii as anarcii_module
from antibody_utils.regions import Region

LEGACY = Path(__file__).parent / "data" / "legacy"

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISAD"
    "TSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTL"
    "TISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
NOT_AN_ANTIBODY = "MKTAYIAKQRQISFVKSHFSRQ"

needs_anarcii = pytest.mark.skipif(
    find_spec("anarcii") is None, reason="needs the `numbering` extra (ANARCII)"
)


# --- Without ANARCII ---------------------------------------------------------


def test_importing_does_not_import_anarcii_or_torch():
    code = (
        "import sys, antibody_utils, antibody_utils.numbering, antibody_utils.regions,"
        " antibody_utils.structure, antibody_utils.geometry;"
        "print('anarcii' in sys.modules, 'torch' in sys.modules)"
    )
    result = subprocess.run(  # noqa: S603 (fixed command, this interpreter)
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.split() == ["False", "False"]


@pytest.mark.parametrize(
    "module",
    [
        "antibody_utils",
        "antibody_utils.numbering",
        "antibody_utils.numbering.anarcii",
        "antibody_utils.numbering.positions",
        "antibody_utils.regions",
        "antibody_utils.structure",
        "antibody_utils.structure.models",
        "antibody_utils.structure.parser",
        "antibody_utils.structure.selection",
        "antibody_utils.geometry",
        "antibody_utils.geometry.orientation",
        "antibody_utils.geometry.superposition",
    ],
)
def test_each_module_imports_first(module):
    """No circular imports, whichever module a fresh interpreter imports first."""
    subprocess.run([sys.executable, "-c", f"import {module}"], check=True)  # noqa: S603


def test_missing_anarcii_gives_install_hint(monkeypatch):
    monkeypatch.setitem(sys.modules, "anarcii", None)
    anarcii_module._model.cache_clear()
    try:
        with pytest.raises(ImportError, match=r"antibody-utils\[numbering\]"):
            number_sequence(TRASTUZUMAB_VH)
    finally:
        anarcii_module._model.cache_clear()


@pytest.mark.parametrize("sequence", ["EVQL-VES", "EVQL/DIQM", "EVQL VES", ""])
def test_rejects_non_letters_before_loading_anarcii(sequence, monkeypatch):
    monkeypatch.setattr(anarcii_module, "_model", None)  # must not be reached
    with pytest.raises(ValueError, match="one-letter amino-acid codes"):
        number({"x": sequence})


def _numbered(chain_type="H"):
    return NumberedSequence(
        name="x",
        chain_type=chain_type,
        scheme=Scheme.IMGT,
        positions=(Position(1), Position(27), Position(111, "A")),
        sequence="EGD",
        score=30.0,
        start=0,
        end=3,
    )


def test_numbered_sequence_needs_one_residue_per_position():
    with pytest.raises(ValueError, match="3 positions but 2 residues"):
        NumberedSequence(
            name="x",
            chain_type="H",
            scheme=Scheme.IMGT,
            positions=(Position(1), Position(2), Position(3)),
            sequence="EV",
            score=30.0,
            start=0,
            end=2,
        )


def test_numbered_sequence():
    numbered = _numbered()
    assert numbered.chain is Chain.HEAVY
    assert numbered.sequence == "EGD"
    assert len(numbered) == 3
    assert list(numbered)[2] == (Position(111, "A"), "D")
    assert numbered.numbering == list(numbered)
    assert [region for *_, region in numbered.annotate_regions("imgt")] == [
        Region.FWH1,
        Region.CDRH1,
        Region.CDRH3,
    ]


@pytest.mark.parametrize("chain_type", ["K", "L"])
def test_kappa_and_lambda_are_light(chain_type):
    assert _numbered(chain_type).chain is Chain.LIGHT


class FakeAnarcii:
    """Stands in for ANARCII: numbers "EVQL…" sequences, fails on anything else."""

    def __init__(self):
        self.batches = []

    def run(self, model, sequences, scheme):
        self.batches.append(list(sequences))
        return {
            name: (
                {
                    "numbering": [((i + 1, " "), aa) for i, aa in enumerate(seq)]
                    + [((len(seq) + 1, " "), "-")],
                    "chain_type": "H",
                    "score": 25.0,
                    "query_start": 0,
                    "query_end": len(seq) - 1,
                    "error": None,
                }
                if seq.startswith("EVQL")
                else {"numbering": None, "error": "not an antibody"}
            )
            for name, seq in sequences.items()
        }


@pytest.fixture
def fake_anarcii(monkeypatch):
    fake = FakeAnarcii()
    monkeypatch.setattr(anarcii_module, "_model", lambda cpu: object())
    monkeypatch.setattr(anarcii_module, "_run", fake.run)
    return fake


def test_numbers_in_batches(fake_anarcii):
    numbered = number([f"EVQL{'A' * i}" for i in range(5)], batch_size=2)
    assert fake_anarcii.batches == [["0", "1"], ["2", "3"], ["4"]]
    assert [len(n) for n in numbered.values()] == [4, 5, 6, 7, 8]


def test_iter_number_is_lazy(fake_anarcii):
    consumed = []

    def sequences():
        for i in range(10):
            consumed.append(i)
            yield "EVQLVES"

    results = iter_number(sequences(), batch_size=3)
    assert consumed == []
    next(results)
    assert consumed == [0, 1, 2]
    assert fake_anarcii.batches == [["0", "1", "2"]]


def test_invalid_sequence_raises_when_its_batch_is_reached(fake_anarcii):
    results = iter_number(["EVQL", "EVQL", "EVQL-VES"], batch_size=2)
    assert [name for name, _ in islice(results, 2)] == ["0", "1"]
    with pytest.raises(ValueError, match="'2' must contain only"):
        next(results)


def test_failures_in_a_batch_map_to_none(fake_anarcii, caplog):
    with caplog.at_level(logging.WARNING):
        numbered = number({"good": "EVQLVES", "bad": "MKTAYIA"})
    assert numbered["bad"] is None
    assert numbered["good"].sequence == "EVQLVES"
    assert "Could not number 'bad': not an antibody" in caplog.text


def test_positions_are_shared_between_sequences(fake_anarcii):
    numbered = number(["EVQLVES", "EVQLVKS"])
    first, second = numbered["0"].positions, numbered["1"].positions
    assert all(a is b for a, b in zip(first, second, strict=True))


def test_equivalent_positions_are_one_object():
    position = anarcii_module._position(111, " ")
    assert position == Position(111)
    for residue_number, insertion in [(111, ""), ("111", " "), ("111", "")]:
        assert anarcii_module._position(residue_number, insertion) is position
    assert anarcii_module._position(111, "A") is not position


def test_shared_positions_cannot_be_modified(fake_anarcii):
    """Sharing `Position`s between sequences relies on all of this."""
    numbered = number(["EVQLVES", "EVQLVKS"])
    first, second = numbered["0"], numbered["1"]
    position = first.positions[0]
    assert position is second.positions[0]
    with pytest.raises(dataclasses.FrozenInstanceError):
        position.insertion = "A"  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        first.positions = ()  # type: ignore[misc]
    assert isinstance(first.positions, tuple)
    # Changing a position means making a new one, leaving the shared one alone.
    assert dataclasses.replace(position, insertion="A") is not position
    assert second.positions[0] == Position(1)


@pytest.mark.parametrize(
    "sequences",
    [
        (seq for seq in ["EVQLVES"]),
        iter(["EVQLVES"]),
        {"EVQLVES"},
    ],
    ids=["generator", "iterator", "set"],
)
def test_number_rejects_streams(sequences, fake_anarcii):
    with pytest.raises(TypeError, match="use iter_number"):
        number(sequences)
    assert fake_anarcii.batches == []


@pytest.mark.parametrize("sequences", [["EVQLVES"], ("EVQLVES",), {"x": "EVQLVES"}])
def test_number_accepts_collections(sequences, fake_anarcii):
    assert len(number(sequences)) == 1


def test_iter_number_accepts_streams(fake_anarcii):
    results = dict(iter_number(seq for seq in ["EVQLVES", "EVQLVKS"]))
    assert list(results) == ["0", "1"]


@pytest.mark.parametrize("function", [number, iter_number])
def test_single_string_is_rejected(function, fake_anarcii):
    with pytest.raises(TypeError, match="use number_sequence"):
        dict(function("EVQLVES"))


@pytest.mark.parametrize("batch_size", [0, 100_001])
def test_batch_size_limits(batch_size, fake_anarcii):
    with pytest.raises(ValueError, match="batch_size"):
        number(["EVQL"], batch_size=batch_size)


# --- With ANARCII ------------------------------------------------------------


@pytest.fixture(scope="module")
def legacy_sequences():
    with open(LEGACY / "numbered_sequences.json") as f:
        sequences = json.load(f)
    # Hand-made numbering with no natural sequence behind it.
    del sequences["trast_L_l54_insertions"]
    return sequences


# Legacy ANARCI numbers 7FAB's lambda chain, which has an unusually short
# CDR-L2, differently from ANARCII's conversion of its (identical) IMGT
# numbering to the Kabat-like schemes: L51/L59 rather than L52/L59A.
KNOWN_DIFFERENCES = {("7fab_L", scheme) for scheme in ("kabat", "chothia", "martin")}


@pytest.mark.numbering
@needs_anarcii
@pytest.mark.parametrize("scheme", list(Scheme))
def test_matches_legacy_anarci(scheme, legacy_sequences):
    sequences = {
        name: "".join(residue for *_, residue in entry["imgt"])
        for name, entry in legacy_sequences.items()
    }
    numbered = number(sequences, scheme, cpu=True)
    assert list(numbered) == list(sequences)
    for name, entry in legacy_sequences.items():
        expected = [(Position(n, i), aa) for n, i, aa in entry[scheme]]
        got = list(numbered[name])
        assert numbered[name].chain == entry["chain_type"]
        if (name, scheme) in KNOWN_DIFFERENCES:
            assert got != expected, f"{name} now matches legacy; update the test"
        else:
            assert got == expected, name


@pytest.mark.numbering
@needs_anarcii
def test_trastuzumab():
    numbered = number({"VH": TRASTUZUMAB_VH, "VL": TRASTUZUMAB_VL}, cpu=True)
    heavy, light = numbered["VH"], numbered["VL"]
    assert (heavy.chain_type, light.chain_type) == ("H", "K")
    assert heavy.sequence == TRASTUZUMAB_VH
    assert light.chain is Chain.LIGHT
    cdrl3 = "".join(
        aa for _, aa, r in light.annotate_regions("imgt") if r is Region.CDRL3
    )
    assert cdrl3 == "QQHYTTPPT"


@pytest.mark.numbering
@needs_anarcii
def test_domain_bounds_exclude_flanking_residues():
    sequence = "AAA" + TRASTUZUMAB_VH + "ASTKGPSVF"
    numbered = number_sequence(sequence, cpu=True)
    assert sequence[numbered.start : numbered.end] == TRASTUZUMAB_VH
    assert numbered.sequence == TRASTUZUMAB_VH


@pytest.mark.numbering
@needs_anarcii
def test_iterable_input_is_named_by_index():
    numbered = number([TRASTUZUMAB_VH, TRASTUZUMAB_VL], "chothia", cpu=True)
    assert list(numbered) == ["0", "1"]
    assert numbered["0"].scheme is Scheme.CHOTHIA


@pytest.mark.numbering
@needs_anarcii
def test_failure_maps_to_none_with_warning(caplog):
    with caplog.at_level(logging.WARNING, logger="antibody_utils.numbering.anarcii"):
        numbered = number({"VH": TRASTUZUMAB_VH, "junk": NOT_AN_ANTIBODY}, cpu=True)
    assert numbered["junk"] is None
    assert numbered["VH"] is not None
    assert "Could not number 'junk'" in caplog.text


@pytest.mark.numbering
@needs_anarcii
def test_number_sequence_raises_on_failure():
    with pytest.raises(NumberingError, match="Could not number"):
        number_sequence(NOT_AN_ANTIBODY, cpu=True)
