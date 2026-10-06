"""Tests of PyIgClassify2 support.

No PyIgClassify2 data is distributed with antibody-utils, so these tests use
synthetic data in the same format.  The tests marked `pyigclassify2` run on a
licensed copy of the real data, if `ANTIBODY_UTILS_PYIGCLASSIFY2_DATA` points
to one, and are skipped otherwise.
"""

import gzip
import json
import os
import sys
import tomllib
from pathlib import Path

import numpy as np
import pytest

from antibody_utils import canonicals
from antibody_utils.canonicals import (
    ClusterAssignment,
    LicenceNotDeclaredError,
    assign,
    assign_cdr,
    declare_licence,
    load_pyigclassify2,
)
from antibody_utils.numbering import NumberedSequence, Position

LEGACY = Path(__file__).parent / "data" / "legacy"
REAL_DATA = os.environ.get(canonicals.DATA_ENV)

COLUMNS = ("model", "pdb", "chain", "cdr", "cdr_ordered", "cluster", "cdr_len")
COLUMNS += ("cdr_seq", "rama")

# Synthetic CDRs: (pdb, chain, cdr, ordered, cluster, sequence), optionally
# followed by a length that disagrees with the sequence.
ROWS = [
    ("1AAA", "H", "H1", "TRUE", "H1-5-1", "AAAAA"),
    ("1AAB", "H", "H1", "TRUE", "H1-5-1", "AAAAC"),
    ("1AAC", "H", "H1", "TRUE", "H1-5-1", "AAAAC"),  # a repeat
    ("1AAD", "H", "H1", "TRUE", "H1-5-2", "WWWWW"),
    ("1AAE", "H", "H1", "TRUE", "H1-5-2", "WWWWY"),
    ("1AAF", "H", "H1", "TRUE", "H1-5-*", "WWWYY"),  # unclustered
    ("1AAG", "H", "H1", "FALSE", "H1-5-1", "YYYYY"),  # disordered
    ("1AAH", "H", "H1", "TRUE", "H1-5-1", "AAAA", 5),  # wrong length
    ("1AAI", "H", "H1", "TRUE", "H1-5-1", "AAUAA"),  # unknown residue
    ("1AAJ", "H", "H1", "TRUE", "H1-6-1", "AAAAAA"),
    # One sequence seen in two clusters, and once unclustered.
    ("1AAK", "L", "L1", "TRUE", "L1-3-1", "DEF"),
    ("1AAL", "L", "L1", "TRUE", "L1-3-2", "DEF"),
    ("1AAM", "L", "L1", "TRUE", "L1-3-2", "DEF"),
    ("1AAN", "L", "L1", "TRUE", "L1-3-*", "DEF"),
]


def write_data(path, rows=ROWS):
    lines = ["  ".join(COLUMNS)]
    for pdb, chain, cdr, ordered, cluster, sequence, *length in rows:
        model = f"1{chain}_{pdb}{chain}_1"
        length = length[0] if length else len(sequence)
        fields = (model, pdb, chain, cdr, ordered, cluster, str(length))
        lines.append("  ".join((*fields, sequence, "B" * len(sequence))))
    text = "\n".join(lines) + "\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    return path


@pytest.fixture(autouse=True)
def isolated(monkeypatch, tmp_path):
    """Use a temporary configuration directory, and no declared licence."""
    import platformdirs

    monkeypatch.setattr(
        platformdirs, "user_config_path", lambda name: tmp_path / "config" / name
    )
    monkeypatch.delenv(canonicals.LICENCE_ENV, raising=False)
    monkeypatch.delenv(canonicals.DATA_ENV, raising=False)
    monkeypatch.setattr(canonicals, "_interactive", lambda: False)
    canonicals._load.cache_clear()
    yield
    canonicals._load.cache_clear()


@pytest.fixture
def declared(monkeypatch):
    monkeypatch.setenv(canonicals.LICENCE_ENV, "declared")


@pytest.fixture
def data_path(tmp_path):
    return write_data(tmp_path / "pyig_cdr_data.txt.gz")


@pytest.fixture
def data(declared, data_path):
    return load_pyigclassify2(data_path)


# --- The licence declaration -------------------------------------------------


def test_undeclared_in_a_script(data_path):
    with pytest.raises(LicenceNotDeclaredError, match="declare_licence") as error:
        load_pyigclassify2(data_path)
    assert canonicals.LICENCE_ENV in str(error.value)
    assert "dunbrack.fccc.edu/lab/PyIgClassify2_lic" in str(error.value)


def test_declared_by_environment(declared, data_path):
    assert len(load_pyigclassify2(data_path)) > 0


@pytest.mark.parametrize("value", ["", "yes", "1"])
def test_environment_needs_declared(monkeypatch, data_path, value):
    monkeypatch.setenv(canonicals.LICENCE_ENV, value)
    with pytest.raises(LicenceNotDeclaredError):
        load_pyigclassify2(data_path)


def test_declare_licence(tmp_path, data_path):
    path = declare_licence()
    assert path == tmp_path / "config" / "antibody-utils" / "licences.toml"
    recorded = tomllib.loads(path.read_text())["pyigclassify2"]
    assert set(recorded) == {"declared", "antibody_utils_version"}
    assert len(load_pyigclassify2(data_path)) > 0


def test_unreadable_declaration(tmp_path, data_path):
    path = tmp_path / "config" / "antibody-utils" / "licences.toml"
    path.parent.mkdir(parents=True)
    path.write_text("not = [toml")
    with pytest.raises(LicenceNotDeclaredError):
        load_pyigclassify2(data_path)


def test_prompt_yes(monkeypatch, capsys, data_path):
    monkeypatch.setattr(canonicals, "_interactive", lambda: True)
    monkeypatch.setattr("builtins.input", lambda prompt: "y")
    assert len(load_pyigclassify2(data_path)) > 0
    output = capsys.readouterr().out
    assert "research use only" in output
    assert "Recorded in" in output
    # Declared once, so no second prompt.
    monkeypatch.setattr("builtins.input", pytest.fail)
    canonicals._load.cache_clear()
    load_pyigclassify2(data_path)


@pytest.mark.parametrize("answer", ["", "n", "no", "maybe"])
def test_prompt_declined(monkeypatch, tmp_path, data_path, answer):
    monkeypatch.setattr(canonicals, "_interactive", lambda: True)
    monkeypatch.setattr("builtins.input", lambda prompt: answer)
    with pytest.raises(LicenceNotDeclaredError, match="No PyIgClassify2 licence"):
        load_pyigclassify2(data_path)
    assert not (tmp_path / "config").exists()


def test_missing_platformdirs(monkeypatch, data_path):
    monkeypatch.setitem(sys.modules, "platformdirs", None)
    with pytest.raises(ImportError, match=r"antibody-utils\[canonicals\]"):
        load_pyigclassify2(data_path)


# --- Loading the data --------------------------------------------------------


def test_no_data_path(declared):
    with pytest.raises(FileNotFoundError, match="PyIgClassify2_download"):
        load_pyigclassify2()


def test_data_path_from_environment(declared, monkeypatch, data_path):
    monkeypatch.setenv(canonicals.DATA_ENV, str(data_path))
    assert load_pyigclassify2().path == data_path


def test_missing_file(declared, tmp_path):
    with pytest.raises(FileNotFoundError, match="No PyIgClassify2 data"):
        load_pyigclassify2(tmp_path / "missing.txt.gz")


def test_wrong_format(declared, tmp_path):
    path = tmp_path / "other.txt"
    path.write_text("a b c\n1 2 3\n")
    with pytest.raises(ValueError, match="not PyIgClassify2"):
        load_pyigclassify2(path)


def test_no_cdrs(declared, tmp_path):
    path = write_data(tmp_path / "empty.txt", rows=[])
    with pytest.raises(ValueError, match="No CDRs"):
        load_pyigclassify2(path)


def test_uncompressed_file(declared, tmp_path):
    path = write_data(tmp_path / "pyig_cdr_data.txt")
    assert len(load_pyigclassify2(path)) == len(load_pyigclassify2(path))


def test_loaded_once(declared, data_path):
    assert load_pyigclassify2(data_path) is load_pyigclassify2(data_path)


def test_contents(data):
    # Disordered CDRs, wrong lengths and unknown residues are left out, and
    # repeats are counted once.
    groups = data._groups
    assert groups["H1", 5].sequences == ("AAAAA", "AAAAC", "WWWWW", "WWWWY", "WWWYY")
    assert groups["H1", 5].clusters == ("H1-5-1", "H1-5-1", "H1-5-2", "H1-5-2", None)
    assert groups["H1", 5].examples[1] == ("1AABH", "1AACH")
    assert len(data) == 7
    assert data.clusters("H1") == ["H1-5-1", "H1-5-2", "H1-6-1"]
    assert data.clusters("cdrl1") == ["L1-3-2"]
    assert repr(data).endswith(": 7 distinct CDR sequences>")


def test_majority_cluster():
    assert canonicals._most_common({"a": 2, "b": 1}) == "a"
    # Ties favour a cluster over none, then the first by name.
    assert canonicals._most_common({None: 1, "b": 1}) == "b"
    assert canonicals._most_common({"b": 1, "a": 1}) == "a"


def test_unknown_cdr(data):
    with pytest.raises(ValueError, match="Unknown CDR 'H5'"):
        data.clusters("H5")


# --- Assignment --------------------------------------------------------------


def test_assign_cdr(data):
    assignment = assign_cdr("WWWWF", "H1", data, neighbours=3)
    assert assignment.cluster == "H1-5-2"
    assert assignment.confidence == pytest.approx(2 / 3)
    assert [n.sequence for n in assignment.neighbours] == ["WWWWY", "WWWWW", "WWWYY"]
    assert [n.cluster for n in assignment.neighbours] == ["H1-5-2", "H1-5-2", None]
    scores = [n.score for n in assignment.neighbours]
    assert scores == sorted(scores, reverse=True)
    assert assignment.neighbours[0].examples == ("1AAEH",)


def test_assign_cdr_unclustered(data):
    assignment = assign_cdr("WWWYY", "H1", data, neighbours=1)
    assert assignment.cluster is None
    assert assignment.confidence == 1.0


@pytest.mark.parametrize(("query", "expected"), [("AAA", "H2-3-1"), ("AAC", "H2-3-2")])
def test_assign_cdr_tie_goes_to_the_nearest(declared, tmp_path, query, expected):
    rows = [
        ("1AAA", "H", "H2", "TRUE", "H2-3-1", "AAA"),
        ("1AAB", "H", "H2", "TRUE", "H2-3-2", "CCC"),
    ]
    data = load_pyigclassify2(write_data(tmp_path / "tie.txt", rows))
    assignment = assign_cdr(query, "H2", data, neighbours=2)
    assert assignment.confidence == 0.5
    assert assignment.cluster == assignment.neighbours[0].cluster == expected


def test_assign_cdr_default_neighbours(data):
    assert len(assign_cdr("AAAAA", "H1", data).neighbours) == 5


def test_assign_cdr_unknown_length(data):
    assert assign_cdr("AAA", "H1", data) == ClusterAssignment(
        "H1", "AAA", None, 0.0, ()
    )


def test_assign_cdr_accepts_region_names_and_lower_case(data):
    assert assign_cdr("wwwww", "cdrh1", data).cluster == "H1-5-2"


@pytest.mark.parametrize("sequence", ["", "AA-AA", "AAUAA"])
def test_assign_cdr_bad_sequence(data, sequence):
    with pytest.raises(ValueError, match="Not a sequence"):
        assign_cdr(sequence, "H1", data)


def test_assign_cdr_bad_neighbours(data):
    with pytest.raises(ValueError, match="at least 1"):
        assign_cdr("AAAAA", "H1", data, neighbours=-1)


def test_assign_cdr_loads_the_data(declared, monkeypatch, data_path):
    monkeypatch.setenv(canonicals.DATA_ENV, str(data_path))
    assert assign_cdr("AAAAA", "H1").cluster == "H1-5-1"


def test_assign_cdr_needs_a_licence(monkeypatch, data_path):
    monkeypatch.setenv(canonicals.DATA_ENV, str(data_path))
    with pytest.raises(LicenceNotDeclaredError):
        assign_cdr("AAAAA", "H1")


@pytest.fixture(scope="module")
def legacy_sequences():
    with open(LEGACY / "numbered_sequences.json") as f:
        return json.load(f)


def _domain(legacy_sequences, name, scheme):
    entry = legacy_sequences[name]
    return NumberedSequence(
        name=name,
        chain_type=entry["chain_type"],
        scheme=scheme,
        positions=tuple(Position(n, i) for n, i, _ in entry[scheme]),
        sequence="".join(residue for *_, residue in entry[scheme]),
        score=0.0,
        start=0,
        end=len(entry[scheme]),
    )


# 12e8's CDRs: North for 1-3, and IMGT 80-87 for 4.
CDRS_12E8 = {
    "H1": "TASGFNIKDYYIH",
    "H2": "WIDPEIGDTE",
    "H3": "NAGHDYDRGRFPY",
    "H4": "ADTSSNTA",
    "L1": "KASQNVGTAVA",
    "L2": "YSASNRYT",
    "L3": "QQYSSYPLT",
    "L4": "GSGTDF",
}


@pytest.mark.parametrize("scheme", ["imgt", "kabat", "chothia", "martin"])
def test_assign_extracts_the_cdrs(legacy_sequences, declared, tmp_path, scheme):
    rows = [
        ("9ZZZ", cdr[0], cdr, "TRUE", f"{cdr}-{len(seq)}-1", seq)
        for cdr, seq in CDRS_12E8.items()
    ]
    data = load_pyigclassify2(write_data(tmp_path / "cdrs.txt", rows))
    light = _domain(legacy_sequences, "12e8_L", scheme)
    heavy = _domain(legacy_sequences, "12e8_H", scheme)
    assignments = assign([light, heavy], data)
    assert list(assignments) == list(CDRS_12E8)
    for cdr, assignment in assignments.items():
        assert assignment.sequence == CDRS_12E8[cdr]
        assert assignment.cluster == f"{cdr}-{len(CDRS_12E8[cdr])}-1"
    assert list(assign(heavy, data)) == ["H1", "H2", "H3", "H4"]


def test_assign_rejects_two_heavy_chains(legacy_sequences, data):
    heavy = _domain(legacy_sequences, "12e8_H", "imgt")
    with pytest.raises(ValueError, match="same chain"):
        assign([heavy, heavy], data)


def test_assign_leaves_out_missing_cdrs(data):
    domain = NumberedSequence(
        name="x",
        chain_type="H",
        scheme="imgt",
        positions=tuple(Position(n) for n in range(27, 32)),
        sequence="AAAAA",
        score=0.0,
        start=0,
        end=5,
    )
    assert list(assign(domain, data)) == ["H1"]


# --- The real data -----------------------------------------------------------

needs_real_data = pytest.mark.skipif(
    not REAL_DATA, reason=f"set {canonicals.DATA_ENV} to a licensed copy"
)


@pytest.fixture(scope="module")
def real_data():
    previous = os.environ.get(canonicals.LICENCE_ENV)
    os.environ[canonicals.LICENCE_ENV] = "declared"
    try:
        return canonicals.PyIgClassify2(
            Path(REAL_DATA), canonicals._read(Path(REAL_DATA))
        )
    finally:
        if previous is None:
            del os.environ[canonicals.LICENCE_ENV]
        else:
            os.environ[canonicals.LICENCE_ENV] = previous


@pytest.mark.pyigclassify2
@needs_real_data
def test_real_data_has_every_cdr(real_data):
    for cdr in canonicals.CDRS:
        assert real_data.clusters(cdr)


@pytest.mark.pyigclassify2
@needs_real_data
@pytest.mark.parametrize(
    ("cdr", "minimum"),
    [("H1", 0.8), ("H2", 0.8), ("H4", 0.9), ("L1", 0.8), ("L2", 0.9), ("L3", 0.75)],
)
def test_real_data_leave_one_out(real_data, cdr, minimum):
    """Assign each clustered sequence from the others of its CDR and length."""
    correct = total = 0
    for (name, length), group in real_data._groups.items():
        if name != cdr:
            continue
        for i, cluster in enumerate(group.clusters):
            keep = np.arange(len(group.sequences)) != i
            if cluster is None or not keep.any():
                continue
            others = canonicals._Group(
                sequences=tuple(np.array(group.sequences)[keep]),
                encoded=group.encoded[keep],
                clusters=tuple(c for j, c in enumerate(group.clusters) if keep[j]),
                examples=tuple(e for j, e in enumerate(group.examples) if keep[j]),
            )
            held_out = canonicals.PyIgClassify2(
                real_data.path, {(name, length): others}
            )
            assigned = assign_cdr(group.sequences[i], cdr, held_out)
            correct += assigned.cluster == cluster
            total += 1
    assert correct / total >= minimum
