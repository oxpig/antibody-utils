# Changelog

All notable changes to antibody-utils are recorded here.  The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).  Until 1.0.0, minor
releases may change the API.

## Unreleased

## 0.1.1 (2026-10-09)

## 0.1.0 (2026-10-08)

The first release.  antibody-utils replaces the reusable parts of `ABDB`, the
Python package behind the original SAbDab.  See
[Migrating from ABDB](https://antibody-utils.readthedocs.io/en/latest/migrating.html).

### Added

- `antibody_utils.numbering`: residue positions (`Position`), numbering schemes
  (`Scheme`: IMGT, Kabat, Chothia and Martin) and chain types (`Chain`).
  `number`, `number_sequence` and `iter_number` number sequences with ANARCII,
  from the optional `numbering` extra.
- `antibody_utils.regions`: the CDR and framework regions of the IMGT, Kabat,
  Chothia, Contact and North definitions, in any of the numbering schemes
  (`get_region`, `annotate_regions` and `RegionSelector`).
- `antibody_utils.sequence`: region and CDR sequences and lengths, and the
  identity and BLOSUM62 similarity of numbered sequences.
- `antibody_utils.liabilities`: sequence-liability motifs
  (`find_liabilities`).
- `antibody_utils.structure`: read structures with GEMMI, numbering their
  antibody chains with ANARCII or from the file's own numbering
  (`read_structure` and `from_gemmi`); pair chains as an `Fv`; and copy
  regions and atoms into a new structure (`select`).
- `antibody_utils.geometry`: Kabsch superposition and RMSD, by region and
  for VH–VL orientation (`superpose`, `rmsd`, `region_rmsd` and
  `orientation_rmsd`), and the ABangle VH–VL orientation measures (`abangle`).
