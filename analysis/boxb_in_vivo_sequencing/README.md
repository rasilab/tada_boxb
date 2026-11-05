# BoxB library, in cells

Figure 5 and Supplementary 8: editing of the same boxB reporter library in human cells, by
TadA–λN recruited to boxB, by TadA-GFPNb, and in a no-editor control.

## Libraries

| Sample | Editor |
| --- | --- |
| `i79_p21` | TadA–λN, recruited to boxB |
| `i79_p22` | TadA-GFPNb |
| `i79_p23` | No editor |

TadA-GFPNb and the no-editor library are carried through counting; the no-editor library is drawn
beside TadA–λN in the Supplementary 8A edit-count panel.

| File | Contents |
| --- | --- |
| [sample_info.csv](./annotations/sample_info.csv) | Experimental metadata of each library |
| [sample_annotations.csv](./annotations/sample_annotations.csv) | Read layout: barcode and UMI positions, and the FASTQ of each stream |
| [oligo_design.csv](./annotations/oligo_design.csv) | The 32 design oligos every construct coordinate is derived from |
| [barcode_annotations.csv](./annotations/barcode_annotations.csv) | Derived construct coordinates; regenerate with [derive_barcode_annotations.py](./scripts/derive_barcode_annotations.py) rather than editing |
| [recorder_backgrounds.csv](./annotations/recorder_backgrounds.csv) | The two fixed recorder cassettes and the invariant base beside boxB |
| [recorder_positions.csv](./annotations/recorder_positions.csv) | The eight measured recorder adenosines |
| [sample_comparisons.csv](./annotations/sample_comparisons.csv) | Which libraries a figure family may put side by side |
| [excluded_constructs.csv](./annotations/excluded_constructs.csv) | Constructs excluded from the analysis; empty means none |
| [figure_outputs.csv](./annotations/figure_outputs.csv) | Every figure this analysis writes |

## Measurement

- One single-pass counter reads each library's reporter and UMI streams together, requiring matching
  read identifiers and valid FASTQ records, and writes recorder counts, boxB counts and attrition in
  one go. Only R1 and R2 are read, which are the streams deposited in SRA.
- Reads are assigned by exact barcode match; the reverse-primer handle is checked at its designed
  position, since the amplicon ends there and the 79-base read stops before it, so an earlier handle
  marks a truncated molecule. Reads with an ambiguous base in the UMI or genotype are rejected and
  counted.
- Phred 30 is required at every measured recorder base and every boxB position; a sub-threshold base
  is written as N so it leaves only its own site's numerator and denominator.
- Recorder columns keep the read orientation, `pos_<i>_C` and `num_<n>_C`, matching the in vitro
  definitions. BoxB columns use sense coordinates and sense bases, `boxb_<position>_<base>`, so an
  A→G edit inside boxB appears as `boxb_<position>_G`.
- A boxB sequence's effect is calculated once, pooling counts across the placements and windows that
  encode the same full sequence, and matched to the wild-type construct sharing its barcode. Measured
  zeros are left off logarithmic axes and counted in the log.
- The two recorder backgrounds are kept apart: the G-depleted cassette changes A2's own sequence
  context, so its editing is not comparable site by site with the standard background.
- No condition was sequenced twice, so means and standard errors describe variation among
  constructs. Tests are paired across matched window and placement groups, without multiple-testing
  correction.

Shared calculations and the visual vocabulary are loaded from the in vitro analysis through
[in_vitro_reuse.py](./scripts/in_vitro_reuse.py), so both experiments take each measurement's
denominator from one definition.

## Run

From the repository root:

```bash
./run_everything.sh --only boxb_in_vivo
```

```bash
./run_everything.sh --plots --only boxb_in_vivo
```

Supplementary 8B compares these libraries with the purified-enzyme reaction, so it reads count tables
from [boxb_in_vitro_sequencing](../boxb_in_vitro_sequencing). Figures are in [figures](./figures) and
the per-panel data in [source_data](../../source_data).
