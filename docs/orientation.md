# VH/VL orientation and RMSD

`antibody_utils.geometry` compares antibody structures: it superposes them,
measures their RMSD over chosen regions, and describes how an Fv's VH and VL
domains are oriented relative to each other.

The examples on this page use SAbDab's IMGT-numbered file for the PDB entry
12E8, which holds two copies of the 2E8 Fab: chains H and L, and chains P and
M.

```pycon
>>> from antibody_utils.structure import Fv, read_structure
>>> structure = read_structure(
...     "pdb_000012e8_sabdab.cif.gz",
...     chains={"H": "H", "L": "L", "P": "H", "M": "L"},
... )
>>> first = Fv(structure["H"], structure["L"])
>>> second = Fv(structure["P"], structure["M"])

```

## ABangle

ABangle describes the orientation of the VH and VL domains with five angles
and a distance (Dunbar J, Fuchs A, Shi J and Deane CM,
[ABangle: characterising the VH–VL orientation in antibodies](https://doi.org/10.1093/protein/gzt020),
*Protein Eng Des Sel* 26:611–620, 2013).

Each domain is given a coordinate frame by superposing a consensus domain on
to the CA atoms of a *core set* of framework positions, which vary least in
structure.  That places, on each domain, a fixed point and two vectors in a
plane fitted to the domain: **H1** and **H2** on VH, and **L1** and **L2** on
VL.  The vector **C** joins the two points.  Then:

| Measure | Meaning |
| --- | --- |
| HL | The torsion angle between **H1** and **L1**, looking along **C**. |
| HC1, LC1 | The angles between **C** and **H1**, and **C** and **L1**: how each domain bends towards or away from the other. |
| HC2, LC2 | The angles between **C** and **H2**, and **C** and **L2**: how each domain twists. |
| dc | The length of **C**, in ångströms: how far apart the domains sit. |

```pycon
>>> from antibody_utils.geometry import abangle
>>> angles = abangle(first)
>>> round(angles.hl, 1), round(angles.dc, 1)
(-71.0, 16.0)
>>> {name: round(value, 1) for name, value in abangle(second).as_dict().items()}
{'HL': -62.4, 'HC1': 74.1, 'LC1': 118.8, 'HC2': 113.2, 'LC2': 81.6, 'dc': 15.9}

```

The two copies of the Fab in this crystal differ in HL by almost 9°.

ABangle's core sets are defined in Chothia numbering, and have exact
equivalents in IMGT numbering, so `abangle` needs an Fv numbered in either of
those schemes.

## RMSD

`region_rmsd` superposes two antibody chains or Fvs over chosen regions and
atoms, and returns the RMSD in ångströms.  Atoms are matched by chain type,
numbered position and atom name, so both structures must be numbered in the
same scheme, and only positions present in both are compared:

```pycon
>>> from antibody_utils.geometry import region_rmsd
>>> round(region_rmsd(first, second, definition="imgt"), 2)
0.93
>>> round(region_rmsd(first.heavy, second.heavy, ["cdrh3"], definition="imgt"), 2)
0.28
>>> round(
...     region_rmsd(
...         first, second, ["cdrs"], definition="north", atoms=["N", "CA", "C", "O"]
...     ),
...     2,
... )
0.66

```

By default, CA atoms are compared.  Pass `superposed=False` to compare the
coordinates as they are, for structures that are already superposed.

`orientation_rmsd` measures how differently two Fvs place the VL domain
relative to the VH domain.  It superposes the VH frameworks, and separately
the VL frameworks, and gives the RMSD between the two resulting placements of
the first Fv's VL framework.  It isn't quite symmetric, so consider averaging
both directions:

```pycon
>>> from antibody_utils.geometry import orientation_rmsd
>>> round(orientation_rmsd(first, second, definition="imgt"), 2)
3.14
>>> round(orientation_rmsd(second, first, definition="imgt"), 2)
3.13

```

## Superposing coordinates

`superpose` and `rmsd` work on any arrays of equivalent points.  `superpose`
finds the rotation and translation that best fit one set of points on to
another (the Kabsch algorithm), and never returns a reflection:

```pycon
>>> import numpy as np
>>> from antibody_utils.geometry import superpose
>>> target = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
>>> mobile = target @ [[0, -1, 0], [1, 0, 0], [0, 0, 1]] + [5, 5, 5]
>>> superposition = superpose(mobile, target)
>>> np.allclose(superposition.apply(mobile), target)
True

```
