# Low-throughput reporter editing

Figure 1C: the fraction of reporter RNA edited by TadA–λN and by unfused TadA, with and without a
boxB hairpin in the reporter.

## Inputs

| File | Contents |
| --- | --- |
| [tada_lambdaN_replicates_tada.csv](./source_data/tada_lambdaN_replicates_tada.csv) | Read counts per sample: total reads, reads with and without indels, and reads carrying A→G substitutions |
| [figure_1c_samples.csv](./annotations/figure_1c_samples.csv) | The four samples drawn, their enzyme and boxB state, and the label each carries in the exported data |
| [edit_classes.csv](./annotations/edit_classes.csv) | The two edit classes drawn, and the count column each reads |

## Measurement

- Each bar is the fraction of indel-free reads carrying at least one or at least two A→G
  substitutions. The two classes are cumulative and overlapping.
- Error bars are 95% percentile intervals from 100,000 binomial draws, with sample size equal to the
  indel-free read count and probability equal to the observed fraction, under the NumPy seed fixed at
  the top of the plotting script.
- Annotations above each pair give the ratio of edited-read fractions with boxB to without, and a
  two-sided [Fisher exact test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html)
  on reads meeting versus not meeting each threshold. Stars are \*q < 0.05, \*\*q < 0.01,
  \*\*\*q < 0.001, with no star above 0.05. Every numerator, denominator, percentage and interval is
  in [figure_1c.csv](../../source_data/figure_1c.csv), and the tests are in
  [figure_1c_statistics.csv](./tables/figure_1c_statistics.csv).

## Run

From the repository root:

```bash
./run_everything.sh --plots --only nanopore_reporter_editing
```

Outputs are an editable SVG, a vector PDF and a 96-DPI PNG in [figures](./figures). Reporter colours
are taken from the Figure 2 palette in
[recorder_edit_summary_panels.csv](../boxb_in_vitro_sequencing/annotations/recorder_edit_summary_panels.csv),
so the two figures agree.
