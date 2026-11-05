What the columns of the annotation files in this directory mean.

### sample_annotations.csv

sample_name: descriptive name of library
sample_id: short identifier for the library, used as the key in every other annotation
barcode1_read: illumina read where barcode is located
barcode1_start: which nucleotide in the read to start extracting barcode from
barcode1_length: how many nucleotides is the barcode
umi_read: illumina read where umi is located
umi_start: which nucleotide in read to start extracting umi from
umi_length: how many nucleotides is the umi
linkage_ref_barcode_column: which column in barcode_annotations.csv to match barcodes to

### barcode_annotations.csv

barcode: sequence of barcode
reverse_complement: reverse complement of barcode sequence for matching reads
oligo_name: descriptive name of oligo
n_variants: how many variants in the oligo library with this barcode
variable_loc: which nucleotide in the read does the variable (NNN...) region start
variable_length: how many nucleotides is the variable region
target_loc: which nucleotide in the read does the a-rich target region start
target_length: how many nucleotides is the target region
target_pos_to_boxb: is the target region located 5' or 3' of the boxb hairpin
target_dist: how many nucleotides away from boxb hairpin is the a-rich target region
variable_type: what type of variable region is this- an a-rich variable region or boxb variable region
variable_subpos: for boxb variable libraries, what positions in the hairpin are varied
g_depleted: whether the spacer regions use the G-depleted sequence, which avoids G-quadruplex secondary structure
### boxb_recruitment_comparison.csv

The direct-recruitment reaction and the split-recruitment reactions it is read against, one row per panel.

source: registered source name of the figure, matching `source` in figure_outputs.csv
panel_order: left-to-right panel position, pAG–TadA first as in every figure comparing the two
y_sample_id: direct-recruitment reaction on the y axis of every panel, matching boxb_reaction_sets.csv
x_sample_id: split-recruitment reaction on this panel's x axis, matching boxb_reaction_sets.csv
measurement: recorder measurement both axes carry, as `selected_measurement` resolves it

### boxb_hairpin_site_figures.csv

The reactions the per-adenosine boxB figure draws, one row per panel. These files carry no in-file
comments, since the annotation loader reads them as plain CSV.

figure: registered source name of the figure this panel belongs to, matching `source` in figure_outputs.csv
enzyme: the editor whose palette name and colour this panel carries, `none` for the shared control panel
sample_id: reaction drawn in this panel, matching sample_info.csv
tada_type: editor in this reaction, `none` for the control panel every figure opens on
tada_conc: TadA concentration in this reaction
condition: temperature and incubation time of this reaction
panel_label: this panel's own heading where a figure is one row of panels, and its row's name where a figure is a grid
panel_row: top-to-bottom position of this panel, numbered from 0; the no-TadA control sits at row 0, column 0
panel_column: left-to-right position within that row, numbered from 0; a row holding fewer panels than the widest is centred
