# λN deep mutational scan

Figure 6: editing of a reporter library in which every codon of the λN peptide is mutated, measured
at two recorder adenosines and linked to each variant through a barcode.

## Design

Two sequencing libraries are counted. A genomic-DNA library links each barcode to the λN variant it
carries; an RNA library counts recorder editing per barcode. Variant identity therefore comes from
the linkage library and editing from the RNA library, joined on the barcode.

| File | Contents |
| --- | --- |
| [dataset_manifest.csv](./annotations/dataset_manifest.csv) | The two libraries, their role, the run each came from, and which RNA library the figures use |
| [read_files.csv](./annotations/read_files.csv) | The FASTQ each library is counted from |
| [read_layouts.csv](./annotations/read_layouts.csv) | Where the barcode, insert and recorder sit in each read, and the constant features checked |
| [insert_annotations.csv](./annotations/insert_annotations.csv) | Every designed insert, its codon position, codon and amino acid; regenerate with [generate_insert_annotations.py](./scripts/generate_insert_annotations.py) |
| [measured_adenosines.csv](./annotations/measured_adenosines.csv) | The recorder adenosines measured, their read positions and sequence contexts |
| [amino_acid_groups.csv](./annotations/amino_acid_groups.csv) | The amino-acid grouping and order the figures draw |

## Measurement

- Phred 30 is required for barcode and insert identity. A measured recorder base below threshold is
  masked as N independently, so a site's callable count is its total reads minus its own N count.
  Edited and callable counts are pooled before dividing.
- Constant regions are validated with editing-aware matching, since a reference T may read as C.
  Features marked `required` filter reads; features marked `diagnostic` only report observations.
  Rejection accounting and count conservation are checked.
- Barcode-to-variant assignment requires three designed reads and 90% designed-insert purity. Every
  final figure additionally requires five designed reads per barcode. Off-design support is reported
  separately, and barcodes one substitution from another barcode are excluded.
- Codon coverage and class distributions retain every distinct DNA insert, with the canonical wild
  type once. Amino-acid figures pool codons and barcodes per substitution, pool synonymous variants
  into one shared wild-type reference, and pool stop codons per position.
- A variant or a position is tested against the pooled synonymous wild type on the barcode rather
  than on the read: wild-type barcodes disagree beyond what their read counts imply, so a read-level
  test rejects wild type against itself. Variances are inflated by barcode-level dispersion estimated
  from wild type, with degrees of freedom equal to the number of contributing barcodes minus one.
- Uncovered rates are left undefined. Measured zeros are included on linear axes and omitted, without
  pseudocounts, on logarithmic axes. The two measured sites share reads and are not independent
  replicates.

Statistics quoted outside a figure are computed by
[compute_wildtype_statistics.py](./scripts/compute_wildtype_statistics.py) and written to
[wildtype_statistics.csv](./tables/wildtype_statistics.csv), one row per statistic.

## Run

From the repository root:

```bash
./run_everything.sh --only ln_dms
```

```bash
./run_everything.sh --plots --only ln_dms
```

Figures are in [figures](./figures) and the per-panel data in [source_data](../../source_data).
