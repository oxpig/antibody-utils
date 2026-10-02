"""Sphinx configuration for antibody-utils."""

from importlib.metadata import version as _version

project = "antibody-utils"
author = "Oxford Protein Informatics Group"
copyright = "2026, University of Oxford"
release = _version("antibody-utils")
version = ".".join(release.split(".")[:2])

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]

# Google-style docstrings only.
napoleon_google_docstring = True
napoleon_numpy_docstring = False

autosummary_generate = True
autodoc_default_options = {"members": True, "show-inheritance": True}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "gemmi": ("https://gemmi.readthedocs.io/en/latest", None),
}

myst_enable_extensions = ["colon_fence", "deflist"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build"]

html_theme = "furo"
html_title = f"antibody-utils {release}"
