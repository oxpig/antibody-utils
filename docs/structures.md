# Working with structures

`antibody_utils.structure` reads structure files with
[GEMMI](https://gemmi.readthedocs.io/), finds the antibody chains, and maps
the numbered positions of each variable domain to residues of the structure.

## Reading a structure

There are two ways to number a structure's antibody chains.

**With ANARCII** (the `numbering` extra), every protein chain is numbered,
and the chains that ANARCII recognises are the antibody chains:

% skip: start if(not HAS_ANARCII, reason="needs ANARCII")

```pycon
>>> from antibody_utils.structure import read_structure
>>> structure = read_structure("12e8.pdb.gz", "imgt")
>>> [(chain.name, chain.chain.name) for chain in structure]
[('L', 'LIGHT'), ('H', 'HEAVY'), ('M', 'LIGHT'), ('P', 'HEAVY')]

```

The structure itself is not renumbered: its residues keep the numbers in the
file, and each antibody chain maps the numbered positions to them.  In the
PDB's file for 12E8, the residue at IMGT position H111 is numbered 99:

```pycon
>>> structure["H"]["111"].seqid.num
99

```

% skip: end

**From the file's own numbering**, if it is already numbered, as
[SAbDab](https://sabdab.opig.stats.ox.ac.uk)'s files are in the IMGT scheme.
Name the antibody chains and their types, and give the scheme:

```pycon
>>> from antibody_utils.structure import read_structure
>>> structure = read_structure(
...     "pdb_000012e8_H_L_ab.cif.gz", "imgt", chains={"H": "H", "L": "L"}
... )
>>> structure
<AntibodyStructure 12E8 antibody chains: L, H>

```

The residues whose numbers fall in the scheme's variable domain (1–128 for
IMGT) are taken as the domain.  The numbering is trusted as it is: a residue
outside the domain that the file numbers within that range, such as the first
residue of a constant domain numbered on from the end of a short domain, is
read as part of the domain.

Either way, only the first variable domain on each chain is numbered, and
only one model is read (`model=0` by default).  An already-loaded
`gemmi.Structure` can be used with `from_gemmi` instead.

## Antibody chains

An `AntibodyChain` holds the numbered domain as a `NumberedSequence`, so the
sequence tools work on it directly, and gives access to its residues as
`gemmi.Residue` objects:

```pycon
>>> heavy = structure["H"]
>>> heavy
<AntibodyChain H heavy, 120 residues numbered in the imgt scheme>
>>> heavy.numbered.sequence[:10]
'EVQLQQSGAE'
>>> heavy["23"].name
'CYS'
>>> for position, residue in list(heavy)[:3]:
...     print(position, residue.name)
1 GLU
2 VAL
3 GLN

```

`region_residues` groups the residues by region:

```pycon
>>> regions = heavy.region_residues("imgt")
>>> {str(region): len(residues) for region, residues in regions.items()}
{'fwh1': 25, 'cdrh1': 8, 'fwh2': 17, 'cdrh2': 8, 'fwh3': 38, 'cdrh3': 13, 'fwh4': 11}

```

The whole GEMMI chain, including any residues outside the variable domain,
ligands and water, is `heavy.gemmi_chain`, and the whole structure is
`structure.structure`.  Treat them as read-only: adding or removing chains or
residues invalidates the antibody chains' mapping.

## Fvs

An `Fv` pairs a heavy and a light chain from the same structure.  Chains are
never paired automatically, so choose the pair yourself, for example from the
structure's metadata or SAbDab's summary file:

```pycon
>>> from antibody_utils.structure import Fv
>>> fv = Fv(structure["H"], structure["L"])
>>> fv.name
'HL'

```

`fv.numbered` gives the two numbered domains, to compare with
`antibody_utils.sequence.identity` and `similarity`.  For the geometry of an
Fv, see {doc}`orientation`.

## Selecting and saving

`select` copies regions and atoms of antibody chains into a new
`gemmi.Structure`, numbered in their scheme, and leaves the original
unchanged.  For example, the backbone of the Fv without CDR-H3:

```pycon
>>> from antibody_utils.structure import BACKBONE_ATOMS, select
>>> copy = select(fv, exclude=["cdrh3"], definition="imgt", atoms=BACKBONE_ATOMS)
>>> len(copy[0]["H"])
107
>>> copy.write_pdb("12e8_fv_backbone.pdb")

```

Write mmCIF with GEMMI's `copy.make_mmcif_document().write_file(...)`.
