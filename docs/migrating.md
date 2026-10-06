# Migrating from ABDB

antibody-utils replaces the reusable parts of `ABDB`, the closed-source
Python package behind the original SAbDab.  `ABDB` needed a per-host
configuration file and a local copy of the SAbDab database to import at all;
antibody-utils needs neither.

## What changes

- **No database.** The parts of `ABDB` that served SAbDab's database
  (`Database_interface`, `ABDB_updater`, `config` and the command-line tool)
  are not included.  To find structures, use
  [SAbDab](https://sabdab.opig.stats.ox.ac.uk) itself.
- **Numbering with ANARCII.** Sequences and structures are numbered with
  [ANARCII](https://github.com/oxpig/ANARCII), from the optional `numbering`
  extra.  ABnum, IMGT/DomainGapAlign, MUSCLE-based numbering and the Wolfguy
  scheme are gone.
- **IMGT by default.** `ABDB` defaulted to the Chothia scheme and definition.
  antibody-utils numbers in IMGT by default, and functions that assign regions
  take the region definition as a required keyword argument, `definition=`.
- **GEMMI, not Biopython.** Structures are read with GEMMI and wrapped, not
  subclassed, so residues and atoms are GEMMI objects.
- **Positions.** Residue positions are `Position(number, insertion)` objects,
  where no insertion is `""` rather than `" "`.  `Position.from_tuple` and
  `Position.to_tuple` convert from and to `ABDB`'s `(number, " ")` tuples.
- **Explicit pairing.** Heavy and light chains are not paired automatically,
  and antigens are not identified.  Pair chains with `Fv(heavy, light)`.

## Equivalent names

Names on the right are relative to `antibody_utils`.

| `ABDB` | antibody-utils |
| --- | --- |
| **Regions** | |
| `AB_Utils.region_definitions.get_region` | `regions.get_region(position, chain, *, scheme, definition)` |
| `AB_Utils.region_definitions.annotate_regions` | `regions.annotate_regions`, or `NumberedSequence.annotate_regions(definition)` |
| `AB_Utils.region_definitions.Accept` | `regions.RegionSelector` |
| `AB_Utils.regions_tuples` and the other region tables | `regions.REGION_GROUPS`, and `get_region` |
| `AB_Utils.tuple_interpret` | `numbering.Position.parse`, `Position.from_tuple` |
| **Numbering** | |
| `Annotate.annotate`, `Annotate.anarci` | `numbering.number_sequence`, `number`, `iter_number` |
| `AB_Utils.sequence_liabilities.annotate_sequences` | `numbering.number` |
| **Sequences** | |
| `AB_Utils.identity`, `fab_identity` | `sequence.identity(a, b, *, definition, regions)` |
| `AB_Utils.similarity` | `sequence.similarity(a, b, *, definition, regions, normalise)` |
| `ABchain.get_CDRs`, `get_sequence` (as sequences) | `sequence.cdr_sequences`, `region_sequences` |
| `AB_Utils.sequence_liabilities.get_liabilities` | `liabilities.find_liabilities(heavy, light)` |
| **Structures** | |
| `AbPDB.AntibodyParser().get_antibody_structure(id, file)` | `structure.read_structure(path, scheme)` |
| `AntibodyStructure.get_fabs()` | `structure.Fv(heavy, light)`, paired by you |
| `AntibodyStructure.get_abchains()`, `get_chains()` | `AntibodyStructure.chains`, `heavy_chains`, `light_chains` |
| `Fab.get_VH()`, `Fab.get_VL()` | `Fv.heavy`, `Fv.light` |
| `ABchain.get_fragments()`, `get_CDRs()`, `residue.region` | `AntibodyChain.region_residues(definition)` |
| `ABchain.get_sequence()` | `AntibodyChain.numbered` |
| `AbPDB.Select.fv_only`, `CDRH3`, `fv_only_no_CDRH3` | `structure.select(source, regions, *, definition, exclude)` |
| `AbPDB.Select.backbone`, `fv_only_backbone` | `structure.select(..., atoms=structure.BACKBONE_ATOMS)` |
| **Geometry** | |
| `AB_Utils.superimpose`, `get_equivalent_arrays` | `geometry.superpose`, `geometry.rmsd` |
| `AB_Utils.calculate_region_rmsd` | `geometry.region_rmsd` |
| `AB_Utils.calculate_orientation_rmsd` | `geometry.orientation_rmsd` |
| `ABangle.abangle().calculate_angles(fab)` | `geometry.abangle(fv)` |

## Not included

- Antigen identification and classification (`get_antigens`, `get_antigen`,
  `is_bound`), and automatic pairing of chains into Fabs and scFvs.
- Header analysis (`is_engineered`, `is_scfv`, hapten detection).
- North canonical classes (`AB_Utils.canonicals.assign_canonical`); for
  sequence-based canonical classes, see
  [SCALOP](https://github.com/oxpig/SCALOP).
- The `interface` selection, TMalign superposition, ABangle's angle
  predictor and plots, and `Visualise`.
