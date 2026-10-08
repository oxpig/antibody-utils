# Quick start

This page takes an antibody from sequence to structure: numbering its
variable domains, extracting CDRs, comparing sequences, finding sequence
liabilities, and reading a structure to measure its VH/VL orientation.

Numbering raw sequences needs ANARCII, from the `numbering` extra (see
{doc}`installation`).  Everything else works with the core install, given
sequences or structures that are already numbered.

## Number a sequence

`number_sequence` numbers one variable-domain sequence, in the IMGT scheme by
default.  Here are trastuzumab's VH and VL:

% skip: start if(not HAS_ANARCII, reason="needs ANARCII")

```pycon
>>> from antibody_utils.numbering import number_sequence
>>> vh = (
...     "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
...     "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
... )
>>> vl = (
...     "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
...     "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
... )
>>> heavy = number_sequence(vh)
>>> light = number_sequence(vl)
>>> heavy.chain_type, light.chain_type, len(heavy)
('H', 'K', 120)
>>> heavy.numbering[:3]
[(Position(number=1, insertion=''), 'E'), (Position(number=2, insertion=''), 'V'), (Position(number=3, insertion=''), 'Q')]

```

The result is a `NumberedSequence`: the position of each residue, the chain
type (`H`, or `K` or `L` for kappa and lambda light chains) and where the
domain lies in the input.  To number many sequences at once, use `number`, or
`iter_number` for a stream too large to hold in memory.

## Extract the CDRs

A *numbering scheme* gives each residue a position; a *region definition*
says which positions form each CDR.  Any definition can be applied to any
scheme (see {doc}`numbering-and-regions`):

```pycon
>>> from antibody_utils.sequence import cdr_sequences
>>> for cdr, sequence in cdr_sequences(heavy, "imgt").items():
...     print(cdr, sequence)
cdrh1 GFNIKDTY
cdrh2 IYPTNGYT
cdrh3 SRWGGDGFYAMDY
>>> for cdr, sequence in cdr_sequences(heavy, "north").items():
...     print(cdr, sequence)
cdrh1 AASGFNIKDTYIH
cdrh2 RIYPTNGYTR
cdrh3 SRWGGDGFYAMDY

```

## Find sequence liabilities

`find_liabilities` looks for motifs that may affect developability, such as
deamidation and isomerisation sites, in the regions where they matter:

```pycon
>>> from antibody_utils.liabilities import find_liabilities
>>> for liability in find_liabilities(heavy, light):
...     if "deamidation" in liability.name:
...         print(liability.chain, liability.residues, liability.positions)
H NG (Position(number=62, insertion=''), Position(number=63, insertion=''))
L NT (Position(number=36, insertion=''), Position(number=37, insertion=''))

```

## Read a structure

`read_structure` reads a PDB or mmCIF file with
[GEMMI](https://gemmi.readthedocs.io/), and numbers each protein chain that is
an antibody chain.  The PDB entry 12E8 holds two copies of the 2E8 Fab:

```pycon
>>> from antibody_utils.structure import Fv, read_structure
>>> structure = read_structure("12e8.pdb.gz")
>>> structure
<AntibodyStructure 12e8 antibody chains: L, H, M, P>
>>> [chain.name for chain in structure.heavy_chains]
['H', 'P']

```

Heavy and light chains are not paired automatically.  Pair them yourself, from
what you know of the structure:

```pycon
>>> fv = Fv(structure["H"], structure["L"])
>>> fv
<Fv HL heavy=H light=L>

```

## Compare sequences

`identity` and `similarity` compare numbered domains position by position, so
they need no alignment.  A VH/VL pair is compared with another chain by chain:

```pycon
>>> from antibody_utils.sequence import identity
>>> round(identity((heavy, light), fv.numbered, definition="imgt"), 3)
0.634
>>> round(identity((heavy, light), fv.numbered, definition="imgt", regions=["cdrs"]), 3)
0.489

```

% skip: end

## Use a numbered structure file

Structure files from [SAbDab](https://sabdab.opig.stats.ox.ac.uk) are already
numbered in the IMGT scheme.  Name the antibody chains, and their numbering is
used as it is, without renumbering with ANARCII:

```pycon
>>> from antibody_utils.structure import Fv, read_structure
>>> sabdab = read_structure("pdb_000012e8_H_L_ab.cif.gz", chains={"H": "H", "L": "L"})
>>> fv = Fv(sabdab["H"], sabdab["L"])
>>> fv.heavy
<AntibodyChain H heavy, 120 residues numbered in the imgt scheme>

```

## Measure the VH/VL orientation

`abangle` measures how the VH and VL domains are oriented relative to each
other (see {doc}`orientation`):

```pycon
>>> from antibody_utils.geometry import abangle
>>> {name: round(value, 1) for name, value in abangle(fv).as_dict().items()}
{'HL': -71.0, 'HC1': 69.2, 'LC1': 118.6, 'HC2': 106.5, 'LC2': 74.6, 'dc': 16.0}

```

## Save part of a structure

`select` copies chosen regions of antibody chains into a new GEMMI structure,
numbered in their scheme, ready to write out:

```pycon
>>> from antibody_utils.structure import select
>>> cdrs = select(fv, ["cdrs"], definition="imgt")
>>> [(chain.name, len(chain)) for chain in cdrs[0]]
[('H', 29), ('L', 18)]
>>> cdrs.write_pdb("12e8_cdrs.pdb")

```
