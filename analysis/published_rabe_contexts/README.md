# Sequence contexts of published rABE editing

Supplementary 3B: the flanking-base preference of rABE editing in published RNA sequencing, for
comparison with the recorder contexts measured here.

Fifteen samples from [Lin et al. 2023](https://doi.org/10.1038/s41592-023-02046-z), covering
PUM1–rABE, PUM2–rABE, unfused rABE, rbFOX2–rABE and an uninduced control, three replicates each.

## Inputs

| File | Contents |
| --- | --- |
| [sample_annotations.csv](./annotations/sample_annotations.csv) | The 15 samples: experiment accession, condition, replicate and display order |
| [run_annotations.csv](./annotations/run_annotations.csv) | ENA run accessions, FASTQ URLs and MD5 checksums |
| [reference_annotations.csv](./annotations/reference_annotations.csv) | The genome assembly and splice annotation URLs |
| [containers.csv](./annotations/containers.csv) | The pinned alignment and analysis images |

Everything needed to run the workflow is in these files; the reads are fetched from ENA by
accession. [fetch_run_metadata.py](./scripts/fetch_run_metadata.py) refreshes the run annotation and
is optional.

## Measurement

- Reads are aligned with HISAT-3N in A-to-G conversion mode against the unmasked GRCh38 primary
  assembly with Ensembl 108 splice annotations, using forward RNA-strand specification and a repeat
  limit of 1,000. Read 1 only.
- A site is retained if it has at least one converted read and more than ten reads in total.
- Contexts are oriented to the edited strand, and contexts containing an N are excluded.
- Each context's value is the mean site-level editing fraction, pooling retained sites across the
  three replicates of a condition.

Choices are declared at the top of [run_analysis.smk](./scripts/run_analysis.smk).

## Run

From the repository root:

```bash
./run_everything.sh --only published_rabe
```

Alignment needs a compute allocation, a genome index and substantial scratch space. The figure is in
[figures](./figures), the per-context table in
[rabe_rna_sequence_context.csv](./tables/rabe_rna_sequence_context.csv), and the panel data in
[figure_s3b.csv](../../source_data/figure_s3b.csv).
