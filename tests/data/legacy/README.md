# Legacy reference outputs

These fixtures record the output of SAbDab's legacy `ABDB` module, which
`antibody-utils` replaces, in its Python 3 version as last changed on
14 February 2025 (SAbDab commit fd76334).  The tests check that the new
implementation reproduces them exactly.  They were generated once, by running
the legacy code; that code is not part of this repository.

| File | Contents |
| --- | --- |
| `get_region.json` | `get_region` for every residue number 0–130 and insertion code (none, A, B, C), for each scheme (IMGT, Kabat, Chothia, Martin), chain (H, L) and definition (IMGT, Kabat, Chothia, Contact, North).  `"?"` means the legacy code found no region. |
| `numbered_sequences.json` | Variable-domain sequences numbered by legacy ANARCI in each scheme, as `[number, insertion, residue]` triples.  They come from PDB entries 12e8, 1ahw, 3hfm and 7fab, plus trastuzumab, an IGHV4-39-like heavy chain (with a synthetic variant whose CDR-H1 has more than two insertions), a lambda light chain and a VHH.  `trast_L_l54_insertions` is trastuzumab's VL with two residues inserted on L54 by hand in the Kabat-like schemes (its IMGT numbering is unchanged), to exercise the CDR-L2 heuristic. |
| `annotate_regions.json` | `annotate_regions` for each of those sequences in each scheme and definition. |
| `identity_similarity.json` | `identity` and `similarity` (total and normalised) for every pair of like-chain sequences, in Chothia numbering with the Chothia and North definitions and in IMGT numbering with the IMGT definition, over the whole Fv and over the CDRs, the framework, CDR-H3 and CDR-L1.  Legacy `similarity` without a region raised `NameError`; those values are `null`. |
| `liabilities.json` | Sequence liabilities found in three VH/VL pairs (trastuzumab, an IGHV4-39-like VH with trastuzumab's VL, and a synthetic pair rich in liability motifs), with the IMGT numbering that legacy ANARCI gave them. |
| `structures/` | The VH/VL pairs of 12e8, 1ahw, 3hfm and 7fab, cut from the PDB files and renumbered with legacy ANARCI in Chothia and IMGT numbering, keeping only the numbered residues.  Named `{entry}_{heavy}{light}_{scheme}.pdb.gz`.  Their numbering is that of `numbered_sequences.json`. |
