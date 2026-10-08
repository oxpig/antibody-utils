# antibody-utils

[![Tests](https://github.com/oxpig/antibody-utils/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/oxpig/antibody-utils/actions/workflows/tests.yml)
[![Documentation](https://readthedocs.org/projects/antibody-utils/badge/?version=latest)](https://antibody-utils.readthedocs.io/en/latest/)
[![pre-commit.ci](https://results.pre-commit.ci/badge/github/oxpig/antibody-utils/main.svg)](https://results.pre-commit.ci/latest/github/oxpig/antibody-utils/main)
[![Coverage](https://codecov.io/gh/oxpig/antibody-utils/graph/badge.svg)](https://codecov.io/gh/oxpig/antibody-utils)
[![PyPI](https://img.shields.io/pypi/v/antibody-utils)](https://pypi.org/project/antibody-utils/)
[![Python versions](https://img.shields.io/pypi/pyversions/antibody-utils)](https://pypi.org/project/antibody-utils/)
[![Downloads](https://static.pepy.tech/badge/antibody-utils/month)](https://pepy.tech/projects/antibody-utils)

[![Licence: BSD-3-Clause](https://img.shields.io/badge/licence-BSD--3--Clause-blue)](LICENCE)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)

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

## Citing

If you use antibody-utils in published work, please cite it and the methods
it implements; see [Citing](https://antibody-utils.readthedocs.io/en/latest/citing.html)
in the documentation.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for how to set up a development
environment, run the tests and propose changes.

## Licence

BSD 3-Clause; see [LICENCE](LICENCE).
