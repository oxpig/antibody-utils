# Contributing to antibody-utils

Thank you for helping.  This guide covers setting up a development environment,
the conventions the code follows, and how changes are reviewed and released.

## Setting up

You need [uv](https://docs.astral.sh/uv/) and Git.

```console
git clone https://github.com/oxpig/antibody-utils.git
cd antibody-utils
uv run --no-project python scripts/fetch_blosum62.py
uv sync
uv run pre-commit install
```

The BLOSUM62 substitution matrix is not tracked in Git.
`scripts/fetch_blosum62.py` downloads it from the NCBI and checks it against a
pinned SHA-256 checksum.  Run it once per clone; the package raises a clear
error if the matrix is missing.

`uv sync` installs the package in editable mode, with the `dev` dependency
group: pytest, Sybil, pre-commit and Bump My Version.

### ANARCII

Most tests run without ANARCII.  The tests marked `numbering`, and the
documentation examples that number sequences, are skipped unless ANARCII is
installed.  ANARCII depends on PyTorch, whose default Linux wheels bundle
several gigabytes of CUDA libraries.  If your development platform lacks a GPU,
first install the CPU-only build into the environment:

```console
uv pip install --torch-backend=cpu -e ".[numbering]" --group dev
```

A later `uv sync` removes it again; use `uv sync --inexact` to keep it.

### pre-commit hooks

`uv run pre-commit install` runs Ruff (linting and formatting) and some basic
checks before each commit.  Pull requests are also checked by
[pre-commit.ci](https://pre-commit.ci/), which pushes any automatic fixes to
the branch.

Git cannot install hooks for you when you clone a repository.  To have
pre-commit hooks installed automatically in every repository you clone or
initialise from now on, install pre-commit as a tool and configure a Git
template directory, once per machine:

```console
uv tool install pre-commit
git config --global init.templateDir ~/.git-template
pre-commit init-templatedir ~/.git-template
```

Repositories without a `.pre-commit-config.yaml` are unaffected.  Existing
clones still need `pre-commit install`.

## Running the tests

```console
uv run pytest
```

The tests include:

- the `>>>` examples in docstrings, collected by pytest's `--doctest-modules`;
- the ```` ```pycon ```` examples in the README and the documentation,
  collected by [Sybil](https://sybil.readthedocs.io/) (see `conftest.py`).
  ```` ```python ```` blocks are illustrative and are not run.

So keep examples runnable, and update their output when behaviour changes.
In the documentation, examples that need ANARCII go between
`% skip: start if(not HAS_ANARCII, reason="needs ANARCII")` and
`% skip: end`.  Documentation examples run in a temporary directory that holds
the example structure files listed in `EXAMPLE_FILES` in `conftest.py`.

Some tests compare results with reference outputs from SAbDab's legacy `ABDB`
module, in `tests/data/legacy/`; its README describes them.  Others read
structure files from the PDB and from SAbDab, in `tests/data/pdb/` and
`tests/data/sabdab/`.

To build the documentation as CI does, treating warnings as errors:

```console
uv run --group docs sphinx-build -W docs docs/_build/html
```

## Conventions

- **Python**: 3.11 or later.  The core dependencies are only NumPy and GEMMI;
  anything heavier belongs in an optional extra, imported when first used.
- **Style**: Ruff checks and formats the code; see `[tool.ruff]` in
  `pyproject.toml`.
- **Docstrings**: [Google style](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings),
  for every public module, class and function, with `Args:`, `Returns:` and
  `Raises:` sections as needed.  Docstrings are reStructuredText, so write a
  hyperlink as `` `text <https://…>`__ ``.  Cite published methods with their
  DOI.
- **Data**: configuration and reference data live in `src/antibody_utils/data/`
  as TOML (hand-maintained) or JSON (generated), not as Python literals.

## Proposing a change

Make your change on a branch from `main`, with tests, and open a pull request.

`main` is protected: changes reach it through pull requests, which are squash
merged.  A pull request can be merged once the **All checks** job passes; it
depends on every other CI job (the test matrix, the tests with ANARCII,
pre-commit and the documentation build).

Describe the change itself in the pull request, relative to `main`.

## Releasing

Maintainers, who can push to `main` directly, release from an up-to-date
`main`:

```console
uv run bump-my-version bump minor   # or major, or patch
git push --follow-tags
```

bump-my-version updates the version in `pyproject.toml` and `uv.lock`, commits,
and creates a signed `v*` tag.  Pushing the tag runs the publish workflow,
which builds the distributions, publishes them, and creates a GitHub release
with Sigstore signatures.  Until 0.1.0, releases go to TestPyPI and are marked
as pre-releases.
