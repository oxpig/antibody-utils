# antibody-utils

Utilities for antibody sequence and structure analysis: numbering, region
definitions, sequence comparison, VH/VL orientation and more.

antibody-utils is developed by the
[Oxford Protein Informatics Group](https://opig.stats.ox.ac.uk/).

## Installation

With [uv](https://docs.astral.sh/uv/), add antibody-utils to your project:

```console
uv add antibody-utils
```

Or with pip:

```console
pip install antibody-utils
```

To number raw sequences and structures with
[ANARCII](https://github.com/oxpig/ANARCII), install the `numbering` extra:

```console
uv add "antibody-utils[numbering]"
```

```console
pip install "antibody-utils[numbering]"
```

> [!TIP]
> ANARCII depends on PyTorch.  On Linux, the default PyTorch wheels from PyPI
> bundle CUDA libraries, a download of several gigabytes.  To get the build that
> suits your machine (CPU-only, or your CUDA or ROCm version):
>
> - With uv's pip interface, let uv detect your hardware:
>
>   ```console
>   uv pip install --torch-backend=auto "antibody-utils[numbering]"
>   ```
>
>   `uv add` has no such flag; in a uv project, point `torch` at the right
>   PyTorch index in your `pyproject.toml` instead, as described in
>   [uv's PyTorch guide](https://docs.astral.sh/uv/guides/integration/pytorch/).
>
> - pip cannot detect your hardware, so install PyTorch first from the
>   [PyTorch index](https://pytorch.org/get-started/locally/) for your platform,
>   e.g. CPU-only:
>
>   ```console
>   pip install torch --index-url https://download.pytorch.org/whl/cpu
>   pip install "antibody-utils[numbering]"
>   ```

## Documentation

See the [documentation](https://antibody-utils.readthedocs.io) for a user guide
and API reference.

## Development

The BLOSUM62 substitution matrix is fetched from NCBI and verified against a
pinned checksum, rather than tracked in this repository.  Fetch it once before
running the tests:

```console
uv run --no-project python scripts/fetch_blosum62.py
uv run pytest
```

The tests include the `>>>` examples in docstrings, and the ```` ```pycon ````
examples in the README and documentation, so keep them runnable.

Install the [pre-commit](https://pre-commit.com/) hooks with
`uv run pre-commit install`.  Pull requests are also checked by
[pre-commit.ci](https://pre-commit.ci/), which pushes any automatic fixes.

> [!TIP]
> Git cannot install hooks for you when you clone a repository.  To have
> pre-commit hooks installed automatically in every repository you clone or
> initialise from now on, install pre-commit as a tool and configure a Git
> template directory, once per machine:
>
> ```console
> uv tool install pre-commit
> git config --global init.templateDir ~/.git-template
> pre-commit init-templatedir ~/.git-template
> ```
>
> Repositories without a `.pre-commit-config.yaml` are unaffected.  Existing
> clones still need `pre-commit install`.

## Licence

BSD 3-Clause; see [LICENCE](LICENCE).
