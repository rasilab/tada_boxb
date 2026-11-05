# BoxB library, purified enzymes

Figures 2-4 and Supplementary 1-7: editing of a synthesised reporter library by purified TadA–λN,
unfused TadA, pAG–TadA and TadA-GFPNb, across enzyme concentration, temperature and time.

## Design

Each reporter carries a 19-nucleotide boxB hairpin, a 17-nucleotide recorder cassette holding eight
measured adenosines, and a 10-nucleotide barcode that declares which mutagenesis window the reporter
randomises and where the recorder sits relative to boxB. Six windows tile boxB; four recorder
placements vary the boxB-to-recorder distance.

| File | Contents |
| --- | --- |
| [sample_info.csv](./annotations/sample_info.csv) | The 20 reactions: enzyme, concentration, temperature and incubation |
| [sample_annotations.csv](./annotations/sample_annotations.csv) | Read layout of each library: barcode and UMI positions |
| [barcode_annotations.csv](./annotations/barcode_annotations.csv) | Every construct barcode, its randomised window, recorder placement and boxB coordinates |
| [recorder_positions.csv](./annotations/recorder_positions.csv) | The eight measured recorder adenosines, their sequence contexts and count columns |
| [boxb_structure_regions.csv](./annotations/boxb_structure_regions.csv) | The stem and loop spans of boxB, and its RNA sequence |
| [boxb_wildtype_mutant_list.csv](./annotations/boxb_wildtype_mutant_list.csv) | The stem-disrupting mutants and their matched wild-type inserts |
| [figure_outputs.csv](./annotations/figure_outputs.csv) | Every figure this analysis writes, under its PDF filename |
| [source_data_panels.csv](./annotations/source_data_panels.csv) | The manuscript panel each figure is exported as |

The remaining annotation CSVs select the reactions, sites, sequences and layout of one figure family
each. [annotations/README.md](./annotations/README.md) describes the columns.

## Measurement

- Reads are assigned to a construct by exact barcode match, their invariant flanking sequences are
  checked, and a reverse-primer handle appearing before its designed position marks a truncated
  amplicon and is rejected. Reads with an ambiguous base in the UMI or the variable insert are
  discarded. Molecules are deduplicated on the seven-nucleotide UMI within each sample, barcode and
  insert.
- Recorder bases below Phred 30 are written as N rather than discarding the molecule. A site's
  editing percentage is its edited molecule count over the molecules callable at that site;
  whole-recorder edit-count classes use molecules callable at all eight adenosines.
- Recorder positions are reported in the RNA orientation, so an A-to-I edit appears as T-to-C on the
  sequenced complementary strand.
- A boxB sequence's effect pools edited and callable counts across every retained placement and
  window representing the same full sequence, then divides, so a sequence shared between windows has
  one estimate. Mutant fractions are divided by wild-type fractions pooled over the same contributing
  barcodes and shown as a percentage of wild type. No coverage cutoff, context correction or
  pseudocount is applied, and measured zeros are left off logarithmic axes and counted in the log.
- Wild-type and mutant comparisons of the recorder sites and placements average the four mutants
  sharing a wild-type barcode first, then compare with that wild type using paired two-sided t-tests
  across eight matched window and placement groups. The per-adenosine boxB figure instead uses an
  unpaired two-sided Welch t-test over exactly the constructs it draws, 15 wild type against 32
  mutant per position, so its stars test the plotted means and construct counts; the plotting script
  asserts that they do. No test is corrected for multiple comparisons, and means and standard errors
  describe variation among constructs.
- The two spacer5_0 libraries at windows 1-3 and 7-10 are excluded for read-quality defects, leaving
  22 distinct wild-type constructs for the in vitro summaries.

## Run

From the repository root:

```bash
./run_everything.sh --only boxb_in_vitro
```

```bash
./run_everything.sh --plots --only boxb_in_vitro
```

Counting reads splits each library by barcode and writes per-construct summaries into the ignored
data directory; it needs a compute allocation. Adding a trailing rule name, for example
`plot_recorder_contexts`, runs one figure family. Figures are in [figures](./figures) and the
per-panel data in [source_data](../../source_data).
