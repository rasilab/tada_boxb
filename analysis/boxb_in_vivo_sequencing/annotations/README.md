What the annotation files in this directory mean. These files carry no in-file comments, so that
GitHub renders them as tables.

### sample_annotations.csv

Sequencing metadata for the three in vivo boxB RNA libraries of the original batch. One row per
sequenced library. Experimental metadata lives in sample_info.csv; the two are joined on sample_id.

`read1_fastq` carries the reporter; `barcode_offset` is the zero-based offset of the construct
barcode inside it. `umi_fastq`/`umi_start`/`umi_length` name the stream and zero-based slice that
holds the molecular identifier. Only these two streams are read, because they are the two deposited
in SRA (data/sra_accessions.csv); the i7 index and reverse reads were sequenced but are not used, so
the counts reproduce from the public runs.

`read1_length` is the cycle count, verified against the raw files in raw_fastq_manifest.csv.

### sample_info.csv

Experimental metadata for the three in vivo boxB RNA libraries. Sequencing metadata is in
sample_annotations.csv; the two are joined on sample_id.

Each row is one sequenced library prepared from one dish of cells, so summaries describe construct
variation rather than biological replication.

`enzyme` is the key every figure selects on; `enzyme_label` is its presentation name.
`recruitment` says whether the editor was targeted to boxB. Panel order is a figure choice and
lives in sample_comparisons.csv, not here.

### sample_comparisons.csv

The libraries each figure family draws, grouped into the comparisons that may be put side by side.

`role` says what a library contributes: `editor` is a boxB-targeted editor and `background` is the
no-editor library. Figure selections name their libraries explicitly. `panel_order` fixes the
left-to-right or top-to-bottom order within the comparison.

### barcode_annotations.csv

Derived from oligo_design.csv by scripts/derive_barcode_annotations.py. Do not edit by hand: rerun
the deriving script instead.

Coordinates are zero-based and relative to the first base of the construct barcode in read 1.
`variable_start` is the randomized window, `target_start` the measured 17-base recorder,
`boxb_start` the 19-base boxB region, and `handle_start` the first base of the reverse primer, which
sits one base past the end of the amplicon.

`target_dist` counts the designed spacer bases between boxB and the recorder, after the invariant
base that always sits beside boxB (`boxb_adjacent_base`), and matches the in vitro annotation of the
same constructs. `designed_spacer_label` is the number in the oligo name, which states that distance
only for the 5′ placements.

`upstream_seq` and `downstream_seq` are ten-base flank patterns, with a wildcard wherever a
randomized base falls inside the flank.

### recorder_backgrounds.csv

The two fixed recorder cassettes of the in vivo boxB library, in read (antisense) orientation.

`recorder_read_sequence` is the measured 17-mer. `boxb_adjacent_base` is the invariant base that
always separates boxB from the designed spacer, so the declared placement distance counts the
designed bases between boxB and the recorder after it.

### recorder_positions.csv

The eight measured recorder adenosines.

`read_index` is the zero-based offset inside the 17-base recorder read sequence; `editing_column` is
the counter column that holds its edited molecules. `context` is the sense-strand trinucleotide in
the standard background only: the G-depleted background changes sense position 3 and therefore A2's
context, which the helper derives per background rather than reading from this file.

### recorder_joint_events.csv

Joint recorder events the counter tabulates directly, rather than as a product of marginals.

`sites` names the recorder positions, most 5′ first; the counter writes the full four-state table
over the molecules callable at every named site, so a joint estimate and each marginal can be
compared on exactly one set of molecules.

### boxb_substitution_measurements.csv

The recorder measurements the substitution figures report.

`sites` names one position for a marginal measurement and several for a joint one, which the counter
tabulates directly. `comparison_sets` scopes a measurement to the comparisons that draw it.

### boxb_substitution_correlations.csv

Compare direct site editing percentages between two recorder sites in the original TadA–λN library.

Site identifiers resolve through boxb_substitution_measurements.csv. `variant_set` names the
sequences a panel draws: `single` the 57 single substitutions, `designed` every mutant the
combinatorial library designs, singles included. Wild type is drawn in neither.

### boxb_wildtype_mutant_list.csv

The predefined stem-disrupting boxB mutants and their matched wild-type inserts, taken unchanged
from the in vitro analysis so the two experiments test the same sequences. A row's `insert` is in
read (antisense) orientation. These are stem-disrupting mutants, never a boxB deletion or no-boxB
control.

Source: analysis/boxb_in_vitro_sequencing/annotations/boxb_wildtype_mutant_list.csv

### boxb_in_vitro_comparison.csv

Per-variant agreement between the purified-enzyme reaction and the cellular library.

Both axes carry one quantity built one way: each designed boxB sequence's editing as a percentage of
its matched wild type, taken as the geometric mean of the named marginal sites. Normalising to
matched wild type removes the activity difference between a two-hour reaction and a cellular library
and leaves the sequence determinants, which is what the two experiments are compared on.

`variant_sets` names the panels a figure draws, in order: `designed` every designed sequence pooled,
`single` the 57 single substitutions, and `stem` and `loop` the mutants whose randomized window does
not and does overlap the loop. Panels share one pair of axes and one marker, so their correlations
are read against each other rather than each against its own scale.

### boxb_effect_heatmaps.csv

Pool edited and callable counts at each site, then average the site editing percentages equally.
Values are direct editing percentages, without wild-type normalization.

### boxb_mutation_distance.csv

Direct editing percentages: pool constructs within each site, then average sites equally. Select
designed sequences only; position zero shows distinct wild-type constructs as in the in vitro
violin.

### boxb_read_coverage.csv

Coverage of the same library and standard recorder background as the editing figures. Pool reads
over all four placements and every window encoding a given full boxB sequence.

### recorder_site_panels.csv

TadA–λN only, averaging distinct constructs over all four placements.

### recorder_edit_count_panels.csv

One figure per row. The main-figure row draws TadA–λN alone; the supplementary row draws it beside
the no-editor library on one shared scale, which is the only way the hundredfold difference between
a recruited editor and the background reads off the bars themselves.

`fold_labels` names the libraries whose bars carry fold changes and significance stars. The
background panel carries none: a ratio between two background rates is not a measurement.

`star_offset` holds the fold label, and the stars below it, in points from the top of the taller bar
and its error bar. It is negative where tall bars leave the label room against them and positive
where a background panel's error bars reach the label instead.

### excluded_constructs.csv

Constructs excluded from the in vivo analysis, established from these libraries rather than copied
from the in vitro exclusion list. Every figure reports its results with and without them. An empty
file means nothing is excluded.

### figure_outputs.csv

Every figure this analysis writes, registered under its PDF filename. Each figure is also written as
a 96-DPI PNG of the same stem.

`source` is the name the drawing script passes to save_figure, and no source or filename may repeat.
`comparison_set` names the libraries a figure draws, and `position_id` the recorder site a family
measures at.

### source_data_panels.csv

Each manuscript panel whose plotted values are exported to the repository-root source_data/
directory, and the registered figure it is drawn from. `source` and `filename` must name a row of
figure_outputs.csv, so a panel cannot be exported from a figure this analysis does not register.
