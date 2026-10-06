# Installation

antibody-utils requires Python 3.11 or later.

## With uv

Add antibody-utils to a [uv](https://docs.astral.sh/uv/) project:

```console
uv add antibody-utils
```

## With pip

```console
pip install antibody-utils
```

## Optional extras

`numbering`
: Number raw sequences and structures with [ANARCII](https://github.com/oxpig/ANARCII).
  This pulls in PyTorch, which is a large download.

  ```console
  uv add "antibody-utils[numbering]"
  ```

  ```console
  pip install "antibody-utils[numbering]"
  ```

  :::{tip}
  On Linux, the default PyTorch wheels from PyPI bundle CUDA libraries, a
  download of several gigabytes.  To get the build that suits your machine
  (CPU-only, or your CUDA or ROCm version):

  - With uv's pip interface, let uv detect your hardware:

    ```console
    uv pip install --torch-backend=auto "antibody-utils[numbering]"
    ```

    `uv add` has no such flag; in a uv project, point `torch` at the right
    PyTorch index in your `pyproject.toml` instead, as described in
    [uv's PyTorch guide](https://docs.astral.sh/uv/guides/integration/pytorch/).

  - pip cannot detect your hardware, so install PyTorch first from the
    [PyTorch index](https://pytorch.org/get-started/locally/) for your
    platform, e.g. CPU-only:

    ```console
    pip install torch --index-url https://download.pytorch.org/whl/cpu
    pip install "antibody-utils[numbering]"
    ```

  :::

`canonicals`
: Assign CDRs to the canonical clusters of
  [PyIgClassify2](https://dunbrack.fccc.edu/lab/PyIgClassify2_lic)
  (Kelow et al. 2022).  PyIgClassify2's data is licensed by Fox Chase Cancer
  Center and is not included: you need a licence from them, and your own copy
  of the data.  See `antibody_utils.canonicals` for how to declare your
  licence and load the data.

  ```console
  uv add "antibody-utils[canonicals]"
  ```

  ```console
  pip install "antibody-utils[canonicals]"
  ```

## Installing for development

Clone the repository, fetch the BLOSUM62 matrix (not tracked in git) and
install the development environment with [uv](https://docs.astral.sh/uv/):

```console
git clone https://github.com/oxpig/antibody-utils.git
cd antibody-utils
uv run --no-project python scripts/fetch_blosum62.py
uv sync
uv run pytest
```
