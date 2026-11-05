# λN-boxB structure panels

Figures 3E and 6G: two views of the λN peptide bound to the boxB hairpin.

| Figure | Panel | View |
| --- | --- | --- |
| [boxb_lambdan_structure](./figures/boxb_lambdan_structure.pdf) | 3E | The complex, with the boxB bases the mutagenesis windows randomise |
| [boxb_lambdan_structure_dms_positions](./figures/boxb_lambdan_structure_dms_positions.pdf) | 6G | The modelled λN helix through residue 22, with the arginines the deep mutational scan resolves |

## Coordinates

The structure is the λN-boxB complex determined by
[Legault et al. 1998](https://doi.org/10.1016/S0092-8674(00)81579-2). It has no Protein Data Bank
deposit, so the coordinates were obtained from its authors and are not redistributed here. Readers who
want them should contact that paper's authors.

[structure_source.csv](./annotations/structure_source.csv) records the chain assignment and a residue
offset against PDB entry [1QFQ](https://www.rcsb.org/structure/1QFQ). It also records the
SHA-256 of the coordinate file these panels were rendered from, so anyone who obtains the
coordinates from the authors can confirm they hold the same file. 1QFQ is a *different* structure,
of a different peptide, by a different group (Schärpf et al. 2000, N36-boxB); it is recorded only so
residue numbers can be mapped between the two, and it is not a substitute for the coordinates these
panels use. Displayed RNA numbering follows the manuscript.

Because the coordinates are not in this repository, the two panels are provided as rendered PDFs and
96-DPI PNGs rather than as a workflow.

## Annotations

| File | Contents |
| --- | --- |
| [dms_highlights.csv](./annotations/dms_highlights.csv) | The two boxB bases labelled in Figure 6G, with their chain, residue and anchor atom |
| [dms_position_summary.csv](./annotations/dms_position_summary.csv) | Mean missense editing per λN position at each measured recorder site, which is the selection behind the residues Figure 6G highlights |
| [structure_source.csv](./annotations/structure_source.csv) | Provenance of the coordinates |

The position summary is derived from the [λN deep mutational scan](../ln_dms_sequencing/README.md) at
the five-read linkage threshold its figures use.
