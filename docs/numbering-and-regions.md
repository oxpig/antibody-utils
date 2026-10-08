# Numbering schemes and region definitions

Two separate choices decide which residues of an antibody variable domain
form each CDR:

- A **numbering scheme** gives each residue a *position*, so that equivalent
  residues in different antibodies share a number.  Positions with an
  insertion code, such as 111A, hold residues that some antibodies have and
  others lack.
- A **region definition** says which positions belong to each
  complementarity-determining region (CDR) and framework region.

Several definitions were first published alongside a scheme of the same name,
which makes the two easy to confuse, but they are independent: antibody-utils
applies any definition to sequences numbered in any scheme.

## Numbering schemes

| Scheme | Reference |
| --- | --- |
| `imgt` | Lefranc M-P *et al.*, [IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains](https://doi.org/10.1016/S0145-305X(02)00039-3), *Dev Comp Immunol* 27:55–77 (2003) |
| `kabat` | Kabat EA *et al.*, *Sequences of Proteins of Immunological Interest*, 5th edition, NIH publication 91-3242 (1991) |
| `chothia` | Chothia C and Lesk AM, [Canonical structures for the hypervariable regions of immunoglobulins](https://doi.org/10.1016/0022-2836(87)90412-8), *J Mol Biol* 196:901–917 (1987) |
| `martin` | Abhinandan KR and Martin ACR, [Analysis and improvements to Kabat and structurally correct numbering of antibody variable domains](https://doi.org/10.1016/j.molimm.2008.05.022), *Mol Immunol* 45:3832–3839 (2008) |

ANARCII numbers sequences in all four (see {doc}`quickstart`).  IMGT is the
default throughout antibody-utils, and it is the scheme of SAbDab's structure
files.

A position is a number and an optional insertion code:

```pycon
>>> from antibody_utils.numbering import Position, parse_chain_position
>>> Position.parse("111A")
Position(number=111, insertion='A')
>>> str(Position(111, "A"))
'111A'
>>> parse_chain_position("H100A")
(<Chain.HEAVY: 'H'>, Position(number=100, insertion='A'))

```

## Region definitions

| Definition | Reference |
| --- | --- |
| `imgt` | Lefranc *et al.* (2003), above |
| `kabat` | Kabat *et al.* (1991), above |
| `chothia` | Chothia and Lesk (1987), above |
| `contact` | MacCallum RM, Martin ACR and Thornton JM, [Antibody-antigen interactions: contact analysis and binding site topography](https://doi.org/10.1006/jmbi.1996.0548), *J Mol Biol* 262:732–745 (1996) |
| `north` | North B, Lehmann A and Dunbrack RL Jr, [A new clustering of antibody CDR loop conformations](https://doi.org/10.1016/j.jmb.2010.10.030), *J Mol Biol* 406:228–256 (2011) |

antibody-utils stores each definition as ranges of IMGT positions, and maps
positions in the other schemes on to IMGT positions.  In IMGT positions, the
CDRs are:

| Definition | CDR-H1 | CDR-H2 | CDR-H3 | CDR-L1 | CDR-L2 | CDR-L3 |
| --- | --- | --- | --- | --- | --- | --- |
| `imgt` | 27–38 | 56–65 | 105–117 | 27–38 | 56–65 | 105–117 |
| `kabat` | 36–40 | 55–74 | 107–117 | 24–40 | 56–69 | 105–117 |
| `chothia` | 27–37 | 57–64 | 107–117 | 24–40 | 56–69 | 105–117 |
| `contact` | 31–40 | 52–66 | 105–116 | 36–42 | 52–68 | 105–116 |
| `north` | 24–40 | 55–66 | 105–117 | 24–40 | 55–69 | 105–117 |

The definitions published in another scheme are exact in IMGT positions
except where the schemes place insertions differently, around CDR-H1 and
CDR-L2.  There, `annotate_regions` uses the rest of the sequence to give the
CDR that the definition gives in its own scheme.

## Assigning regions

`annotate_regions` assigns every residue of a numbered domain to its region,
and is what the rest of antibody-utils uses.  `NumberedSequence` has it as a
method, and `antibody_utils.sequence.region_sequences` and `cdr_sequences`
build on it.

`get_region` looks up one position on its own.  It can't see insertions
elsewhere in the sequence, so near CDR-H1 and CDR-L2 it can differ from
`annotate_regions` when the definition and the scheme differ.  The same
Chothia position falls in different regions under different definitions:

```pycon
>>> from antibody_utils.regions import get_region
>>> for definition in ["imgt", "kabat", "chothia", "contact", "north"]:
...     print(definition, get_region("35", "H", scheme="chothia", definition=definition))
imgt fwh2
kabat cdrh1
chothia fwh2
contact cdrh1
north cdrh1

```

Regions are named `fwh1`–`fwh4`, `cdrh1`–`cdrh3`, `fwl1`–`fwl4` and
`cdrl1`–`cdrl3`.  Wherever antibody-utils takes region names, it also takes
these groups:

| Group | Regions |
| --- | --- |
| `hcdrs`, `lcdrs`, `cdrs` | the CDRs of the heavy chain, light chain or both |
| `hframework`, `lframework`, `framework` | the framework regions likewise |
| `vh`, `vl`, `fv` | the whole heavy domain, light domain, or both |

`RegionSelector` selects positions by region, and can add or exclude
individual positions:

```pycon
>>> from antibody_utils.regions import RegionSelector
>>> selector = RegionSelector(["cdrs"], scheme="imgt", definition="north")
>>> selector.accepts("111A", "H"), selector.accepts("50", "H")
(True, False)

```
