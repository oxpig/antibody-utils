"""Assign CDRs to the canonical clusters of PyIgClassify2.

PyIgClassify2 (Kelow et al. 2022) clusters the conformations of antibody CDRs
in the PDB.  Its data is licensed by Fox Chase Cancer Center (FCCC) for
research use (academic and government licences) or under a separate
commercial licence, and may not be redistributed, so none of it is included
in antibody-utils.  To use this module:

1. Obtain a licence and the data from FCCC
   (https://dunbrack.fccc.edu/lab/PyIgClassify2_lic), and download
   `pyig_cdr_data.txt.gz`.
2. Declare that you hold a licence: answer the prompt that appears the first
   time the data is loaded in a terminal, or call `declare_licence`, or set
   `ANTIBODY_UTILS_PYIGCLASSIFY2_LICENCE=declared`.
3. Point `load_pyigclassify2` at the data file, or set
   `ANTIBODY_UTILS_PYIGCLASSIFY2_DATA` to its path.

The data is read into memory and never written elsewhere.  Please cite
PyIgClassify2 in any publication of results that use it.

PyIgClassify2 assigns clusters from the backbone dihedral angles of CDRs in
structures.  This module assigns a sequence instead: it finds the known CDRs
of the same type and length whose sequences are most similar (by BLOSUM62),
and takes the most common cluster among them.  That is an approximation.  In
leave-one-out tests on PyIgClassify2's 2022 data, it gave the same cluster (or
the same lack of one) for 72-88% of distinct CDR sequences, depending on the
CDR.  For sequences that PyIgClassify2 places in a cluster, it found that
cluster for 82-99% of them, except for CDR-H3, where it found only 23%: most
H3 loops are unclustered, so the neighbours of a clustered one usually are
too.

CDRs 1-3 are those of the North definition (`North, Lehmann and Dunbrack 2011
<https://doi.org/10.1016/j.jmb.2010.10.030>`__).  CDR4 is the loop between
strands D and E of framework 3, at IMGT positions 80-87.

This module needs the `canonicals` extra:
`pip install "antibody-utils[canonicals]"`.
"""

from __future__ import annotations

import gzip
import os
import sys
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import UTC, datetime
from functools import cache
from pathlib import Path

import numpy as np

from antibody_utils._data import load_blosum62, load_toml
from antibody_utils.numbering.anarcii import NumberedSequence
from antibody_utils.numbering.positions import Chain, Scheme
from antibody_utils.regions import _imgt_positions

__all__ = [
    "CDRS",
    "ClusterAssignment",
    "LicenceNotDeclaredError",
    "Neighbour",
    "PyIgClassify2",
    "assign",
    "assign_cdr",
    "declare_licence",
    "load_pyigclassify2",
]

#: The CDRs that PyIgClassify2 clusters.
CDRS = ("H1", "H2", "H3", "H4", "L1", "L2", "L3", "L4")

LICENCE_ENV = "ANTIBODY_UTILS_PYIGCLASSIFY2_LICENCE"
DATA_ENV = "ANTIBODY_UTILS_PYIGCLASSIFY2_DATA"

_SETTINGS = "canonicals.toml"
_LICENCES_FILE = "licences.toml"


def _settings() -> dict:
    return load_toml(_SETTINGS)


# --- The licence declaration -------------------------------------------------


class LicenceNotDeclaredError(RuntimeError):
    """No PyIgClassify2 licence has been declared."""


def _licences_path() -> Path:
    try:
        from platformdirs import user_config_path
    except ImportError as error:
        raise ImportError(
            "PyIgClassify2 support needs platformdirs, which is not installed.  "
            'Install the canonicals extra: pip install "antibody-utils[canonicals]"'
        ) from error
    return user_config_path("antibody-utils") / _LICENCES_FILE


def _terms() -> str:
    settings = _settings()["pyigclassify2"]
    return (
        "PyIgClassify2 is licensed by Fox Chase Cancer Center (FCCC).  In "
        "summary, its academic/government licence permits research use only, "
        "is non-transferable, and forbids redistributing the data or copying it "
        "except as authorised use requires; commercial use needs a separate "
        "licence from FCCC.  Publications of results must acknowledge "
        "PyIgClassify2's authors.  The full terms are at\n"
        f"  {settings['licence_url']}\n"
        f"Please cite: {settings['citation']}"
    )


def declare_licence() -> Path:
    """Record that you, or your institution, hold a PyIgClassify2 licence.

    Only declare this if a licence from FCCC covers your use: an
    academic/government licence for research, or a commercial licence.  The
    declaration is recorded in your user configuration directory, so it is
    made once per user and machine.

    Returns:
        The file that the declaration was recorded in.
    """
    path = _licences_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    from antibody_utils import __version__

    path.write_text(
        "# Licences declared for antibody-utils.\n"
        "[pyigclassify2]\n"
        f'declared = "{datetime.now(UTC).isoformat(timespec="seconds")}"\n'
        f'antibody_utils_version = "{__version__}"\n'
    )
    return path


def _declared() -> bool:
    if os.environ.get(LICENCE_ENV, "").strip().lower() == "declared":
        return True
    path = _licences_path()
    if not path.is_file():
        return False
    import tomllib

    try:
        return "pyigclassify2" in tomllib.loads(path.read_text())
    except tomllib.TOMLDecodeError:
        return False


def _interactive() -> bool:
    return sys.stdin.isatty() and sys.stdout.isatty()


def _require_licence() -> None:
    """Check for a licence declaration, prompting for one in a terminal.

    Raises:
        LicenceNotDeclaredError: If no licence is declared, and none is
            declared at the prompt (or there is no terminal to prompt in).
    """
    if _declared():
        return
    if _interactive():
        print(_terms())
        answer = input(
            "\nDo you, or does your institution, hold a PyIgClassify2 licence "
            "from FCCC that covers this use? [y/N] "
        )
        if answer.strip().lower() in ("y", "yes"):
            path = declare_licence()
            print(f"Recorded in {path}.")
            return
        raise LicenceNotDeclaredError("No PyIgClassify2 licence was declared.")
    raise LicenceNotDeclaredError(
        "Using PyIgClassify2 needs a licence from FCCC.  If you, or your "
        "institution, hold one that covers this use, declare it by running "
        "antibody_utils.canonicals.declare_licence() once, or by setting "
        f"{LICENCE_ENV}=declared.\n\n{_terms()}"
    )


# --- The data ----------------------------------------------------------------

# Columns of `pyig_cdr_data.txt` that are read.  The file is space-delimited;
# its last two columns can run together, so only earlier columns are used.
_COLUMNS = ("pdb", "chain", "cdr", "cdr_ordered", "cluster", "cdr_len", "cdr_seq")


@dataclass(frozen=True, slots=True)
class _Group:
    """The distinct known sequences of one CDR and length."""

    sequences: tuple[str, ...]
    encoded: np.ndarray  # one row of residue indices per sequence
    clusters: tuple[str | None, ...]  # `None` where PyIgClassify2 has none
    examples: tuple[tuple[str, ...], ...]  # e.g. ("12E8H",) for each sequence


class PyIgClassify2:
    """PyIgClassify2's CDR clusters, read from a licensed copy of its data.

    Create one with `load_pyigclassify2`.
    """

    __slots__ = ("_groups", "path")

    def __init__(self, path: Path, groups: dict[tuple[str, int], _Group]) -> None:
        """Hold the data read from `path`."""
        #: The data file.
        self.path = path
        self._groups = groups

    def __repr__(self) -> str:
        """Describe the data."""
        return f"<PyIgClassify2 {self.path}: {len(self)} distinct CDR sequences>"

    def __len__(self) -> int:
        """The number of distinct CDR sequences."""
        return sum(len(group.sequences) for group in self._groups.values())

    def clusters(self, cdr: str) -> list[str]:
        """List the clusters of a CDR.

        Args:
            cdr: The CDR, such as `"H1"`.

        Returns:
            The names of its clusters, such as `"H1-13-1"`, sorted.
        """
        cdr = _cdr_name(cdr)
        return sorted(
            {
                cluster
                for (name, _), group in self._groups.items()
                if name == cdr
                for cluster in group.clusters
                if cluster is not None
            }
        )


def _cdr_name(cdr: str) -> str:
    name = str(cdr).upper().removeprefix("CDR")
    if name not in CDRS:
        raise ValueError(f"Unknown CDR {cdr!r}; expected one of {', '.join(CDRS)}")
    return name


def _residue_index() -> dict[str, int]:
    return {residue: i for i, (residue, _) in enumerate(_blosum_order())}


@cache
def _blosum_order() -> tuple[tuple[str, int], ...]:
    residues = sorted({a for a, _ in load_blosum62()})
    return tuple((residue, i) for i, residue in enumerate(residues))


@cache
def _blosum_array() -> np.ndarray:
    blosum62 = load_blosum62()
    residues = [residue for residue, _ in _blosum_order()]
    return np.array(
        [[blosum62[a, b] for b in residues] for a in residues], dtype=np.int16
    )


def _read(path: Path) -> dict[tuple[str, int], _Group]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        names = f.readline().split()
        try:
            columns = [names.index(column) for column in _COLUMNS]
        except ValueError:
            raise ValueError(
                f"{path} is not PyIgClassify2's pyig_cdr_data.txt: it lacks "
                f"some of the columns {', '.join(_COLUMNS)}"
            ) from None
        index = _residue_index()
        # For each CDR and length: each distinct sequence's clusters and
        # examples, counted over every structure it occurs in.
        found: dict[tuple[str, int], dict[str, tuple[Counter, set[str]]]] = {}
        for line in f:
            fields = line.split()
            if len(fields) <= max(columns):
                continue
            pdb, chain, cdr, ordered, cluster, length, sequence = (
                fields[i] for i in columns
            )
            # Only fully ordered CDRs have reliable clusters.
            if ordered != "TRUE" or cdr not in CDRS:
                continue
            if not length.isdigit() or len(sequence) != int(length):
                continue
            if any(residue not in index for residue in sequence):
                continue
            clusters, examples = found.setdefault((cdr, len(sequence)), {}).setdefault(
                sequence, (Counter(), set())
            )
            # Starred clusters (such as "H3-13-*") mean no cluster.
            clusters[None if cluster.endswith("*") else cluster] += 1
            examples.add(pdb + chain)

    groups = {}
    for key, by_sequence in found.items():
        sequences = tuple(sorted(by_sequence))
        groups[key] = _Group(
            sequences=sequences,
            encoded=np.array(
                [[index[r] for r in sequence] for sequence in sequences],
                dtype=np.int8,
            ),
            clusters=tuple(
                _most_common(by_sequence[sequence][0]) for sequence in sequences
            ),
            examples=tuple(
                tuple(sorted(by_sequence[sequence][1])) for sequence in sequences
            ),
        )
    if not groups:
        raise ValueError(f"No CDRs found in {path}")
    return groups


def _most_common(counts: Counter) -> str | None:
    # Break ties in favour of a cluster over none, then by name.
    return min(counts, key=lambda c: (-counts[c], c is None, c or ""))


def load_pyigclassify2(path: str | os.PathLike[str] | None = None) -> PyIgClassify2:
    """Load a licensed copy of PyIgClassify2's CDR data.

    The data is read once per file and process, and is kept in memory only.
    A PyIgClassify2 licence must have been declared (see the module
    documentation); in a terminal, you are asked the first time.

    Args:
        path: The path to `pyig_cdr_data.txt.gz` (or the uncompressed
            `pyig_cdr_data.txt`).  By default, the path in the environment
            variable `ANTIBODY_UTILS_PYIGCLASSIFY2_DATA`.

    Returns:
        The data.

    Raises:
        LicenceNotDeclaredError: If no licence is declared.
        FileNotFoundError: If no path is given and none is set in the
            environment, or the file does not exist.
        ValueError: If the file is not in the expected format.
    """
    _require_licence()
    if path is None:
        path = os.environ.get(DATA_ENV)
        if not path:
            settings = _settings()["pyigclassify2"]
            raise FileNotFoundError(
                f"Give the path to PyIgClassify2's {settings['data_file']}, or set "
                f"{DATA_ENV} to it.  Licence holders can download it from "
                f"{settings['download_url']}"
            )
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"No PyIgClassify2 data at {path}")
    return _load(path)


@cache
def _load(path: Path) -> PyIgClassify2:
    return PyIgClassify2(path, _read(path))


# --- Assignment --------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Neighbour:
    """A known CDR similar to the one being assigned.

    Attributes:
        sequence: Its sequence.
        cluster: Its PyIgClassify2 cluster, or `None` if it has none.
        score: Its BLOSUM62 similarity to the CDR being assigned.
        examples: Structures with this CDR sequence, as PDB code and chain,
            such as `"12E8H"`.
    """

    sequence: str
    cluster: str | None
    score: int
    examples: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ClusterAssignment:
    """The cluster assigned to a CDR sequence.

    Attributes:
        cdr: The CDR, such as `"H1"`.
        sequence: The CDR sequence.
        cluster: The assigned cluster, such as `"H1-13-1"`, or `None` if the
            most similar known CDRs are mostly unclustered, or no known CDR
            has this length.
        confidence: The fraction of the neighbours in that cluster (or with no
            cluster), from 0 to 1.  0 if there are no neighbours.
        neighbours: The most similar known CDRs of the same length, most
            similar first.
    """

    cdr: str
    sequence: str
    cluster: str | None
    confidence: float
    neighbours: tuple[Neighbour, ...]


def assign_cdr(
    sequence: str,
    cdr: str,
    data: PyIgClassify2 | None = None,
    *,
    neighbours: int | None = None,
) -> ClusterAssignment:
    """Assign one CDR sequence to a PyIgClassify2 cluster.

    The `neighbours` known CDRs of the same type and length with the highest
    BLOSUM62 similarity vote, and the most common cluster wins.  Ties go to
    the cluster of the most similar neighbour.

    Args:
        sequence: The CDR sequence, by the North definition for CDRs 1-3, or
            IMGT positions 80-87 for CDR4.
        cdr: The CDR: `"H1"`-`"H4"` or `"L1"`-`"L4"` (or the region names
            `"cdrh1"` and so on).
        data: The PyIgClassify2 data.  By default, loaded with
            `load_pyigclassify2()`.
        neighbours: How many similar CDRs vote (by default, 5).

    Returns:
        The assignment.

    Raises:
        ValueError: If `cdr` or `sequence` is not recognised.
        LicenceNotDeclaredError: If `data` is not given and no licence is
            declared.
        FileNotFoundError: If `data` is not given and the data cannot be
            found.
    """
    cdr = _cdr_name(cdr)
    sequence = sequence.upper()
    index = _residue_index()
    if not sequence or any(residue not in index for residue in sequence):
        raise ValueError(f"Not a sequence of one-letter residue codes: {sequence!r}")
    k = neighbours or _settings()["assignment"]["neighbours"]
    if k < 1:
        raise ValueError("neighbours must be at least 1")
    if data is None:
        data = load_pyigclassify2()

    group = data._groups.get((cdr, len(sequence)))
    if group is None:
        return ClusterAssignment(cdr, sequence, None, 0.0, ())
    query = np.array([index[r] for r in sequence], dtype=np.intp)
    scores = _blosum_array()[query, group.encoded].sum(axis=1)
    # Most similar first; ties in sequence order, for reproducibility.
    nearest = np.argsort(-scores, kind="stable")[:k]
    found = tuple(
        Neighbour(
            sequence=group.sequences[i],
            cluster=group.clusters[i],
            score=int(scores[i]),
            examples=group.examples[i],
        )
        for i in nearest
    )
    votes = Counter(neighbour.cluster for neighbour in found)
    top = max(votes.values())
    cluster = next(n.cluster for n in found if votes[n.cluster] == top)
    return ClusterAssignment(cdr, sequence, cluster, top / len(found), found)


def _cdr4(numbered: NumberedSequence) -> str:
    chain = numbered.chain
    first, last = _settings()["cdr4"][chain.value]
    to_imgt = _imgt_positions(Scheme(numbered.scheme), chain)
    return "".join(
        residue
        for position, residue in numbered
        if first <= to_imgt.get(position.number, 0) <= last
    )


def assign(
    numbered: NumberedSequence | Iterable[NumberedSequence],
    data: PyIgClassify2 | None = None,
    *,
    neighbours: int | None = None,
) -> dict[str, ClusterAssignment]:
    """Assign the CDRs of numbered domains to PyIgClassify2 clusters.

    Args:
        numbered: A numbered domain, or several (such as a VH and its VL), in
            any scheme.
        data: The PyIgClassify2 data.  By default, loaded with
            `load_pyigclassify2()`.
        neighbours: How many similar CDRs vote (by default, 5).

    Returns:
        The assignment of each CDR present, such as `"H1"`, in order.  CDRs
        with no residues are left out.

    Raises:
        ValueError: If two domains are of the same chain.
        LicenceNotDeclaredError: If `data` is not given and no licence is
            declared.
        FileNotFoundError: If `data` is not given and the data cannot be
            found.
    """
    domains = [numbered] if isinstance(numbered, NumberedSequence) else list(numbered)
    chains = [domain.chain for domain in domains]
    if len(set(chains)) != len(chains):
        raise ValueError("More than one domain of the same chain given")
    if data is None:
        data = load_pyigclassify2()

    assignments = {}
    for domain in sorted(domains, key=lambda d: d.chain is not Chain.HEAVY):
        c = domain.chain.value
        sequences = dict.fromkeys((f"{c}{n}" for n in (1, 2, 3, 4)), "")
        for _, residue, region in domain.annotate_regions("north"):
            if region is not None and region.is_cdr:
                sequences[f"{c}{region.value[-1]}"] += residue
        sequences[f"{c}4"] = _cdr4(domain)
        for cdr, sequence in sequences.items():
            if sequence:
                assignments[cdr] = assign_cdr(
                    sequence, cdr, data, neighbours=neighbours
                )
    return assignments
