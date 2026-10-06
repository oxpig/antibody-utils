# SAbDab structure files

IMGT-numbered structure files from the current version of
[SAbDab](https://sabdab.opig.stats.ox.ac.uk), the Structural Antibody
Database, downloaded on 6 October 2026 for the same PDB entries as the legacy
fixtures in `../legacy/`.  The tests check that `antibody-utils` reads the
numbering of these files consistently with legacy `ABDB` and with SAbDab's own
summary.

| File | Contents |
| --- | --- |
| `pdb_{entry}_{heavy}_{light}_ab.cif.gz`, `…_abag.cif.gz` | One antibody instance, with its antigen (`abag`) if it has one.  Gzipped, as downloaded. |
| `pdb_000012e8_sabdab.cif.gz` | The whole of entry 12e8, with both instances. |
| `SAbDab_summary.csv` | SAbDab's summary rows for the instances, merged from the downloads.  `VH`, `VL` and the `CDR-…` columns are IMGT. |

SAbDab numbers the residues after the variable domain by counting on from its
last position.  Where a light-chain domain ends at L127, the first residue of
the constant domain is therefore numbered L128, which is a valid IMGT position,
so these files' light chains read as one residue longer than SAbDab's `VL`.
The tests that this affects are marked as expected failures until SAbDab's
renumbering is fixed.

If you use SAbDab, please cite Dunbar J, Krawczyk K, Leem J, Baker T, Fuchs A,
Georges G, Shi J and Deane CM, "SAbDab: the structural antibody database",
*Nucleic Acids Research* 42:D1140–D1146 (2014),
https://doi.org/10.1093/nar/gkt1043.
