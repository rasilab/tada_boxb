"""Shared CSV loading, construct matching, and recorder calculations; no file exports."""

import itertools
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, ttest_rel

ANALYSIS_DIRECTORY = "analysis/boxb_in_vitro_sequencing"
ANNOTATION_DIRECTORY = f"{ANALYSIS_DIRECTORY}/annotations"
COUNT_DIRECTORY = f"{ANALYSIS_DIRECTORY}/data/summary_stats_combined"
LOOP_COUNT_DIRECTORY = f"{ANALYSIS_DIRECTORY}/data/boxb_stats_combined"
LOOP_EDIT_BINS = {"frac_1edit": (1,), "frac_2edit": (2, 3)}
LIBRARY_TYPE = "boxb"
BACKGROUND = "no"
LAYOUT_PATTERN = r"^(spacer[35]_\d+)"
EXPECTED_PAIRS = 32
EXPECTED_WT = 22
EXPECTED_WT_WINDOWS = 6
EXPECTED_PAIRED_WINDOWS = 8
EXPECTED_DISTANCES = [0, 10, 20, 30]
EXPECTED_SITE_MUT = 8
EXPECTED_SET_ENZYMES = 2
EXCLUDED_LAYOUT = "spacer5_0"
EXCLUDED_WINDOWS = ("1_3", "7_10")
BOXB_LENGTH = 19
# The counter measures these read offsets from the hairpin origin, as its BOXB_MEASURED_POSITIONS,
# and writes them as pos_0 onwards in the order listed here: the three loop adenosines, then the one
# that closes the stem. Its edit-count bins stay over the loop offsets alone.
BOXB_READ_OFFSETS = [7, 9, 10, 6]
DNA_BASES = "ACGT"
RNA_BASES = "ACGU"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
READ_TO_RNA = dict(zip("ACGT", "UGCA"))
SAMPLE_KEYS = ["sample_id", "tada_type", "tada_conc", "condition", "enzyme_label", "panel_order"]
GEOMETRY_KEYS = ["target_pos_to_boxb", "target_dist"]
EDIT_TYPES = ["frac_1edit", "frac_2edit"]
# The `tada_type` of a reaction carrying no TadA, which every per-adenosine figure opens on.
CONTROL_ENZYME = "none"
ORIENTATION_GROUPS = {"5′": [0, 10], "3′": [20, 30]}
SITE_STARS = [(.001, "***"), (.01, "**"), (.05, "*")]
MEASUREMENT_SEPARATOR = ";"
MEAN_PREFIX = "mean"
# How a measurement naming several recorder sites turns them into one number. A mean is derived from
# the sites' own log₂ ratios rather than counted from reads, so it adds no count of its own.
COMBINE_KINDS = ["site", MEAN_PREFIX]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def annotation(name):
    return pd.read_csv(f"{ANNOTATION_DIRECTORY}/{name}.csv", dtype={"variable_subpos": str})


def enzyme_metadata():
    return annotation("recorder_context_selection")[SAMPLE_KEYS].drop_duplicates().sort_values("panel_order")


def boxb_reaction_sets():
    """The boxB reactions, grouped into the pairs each figure draws, with the scale each pair shares.

    A set is two reactions drawn as two panels, named for whatever distinguishes them: the enzyme
    concentration for the TadA-λN against TadA series, the editor for the pAG/GFP-nanobody pair, whose
    two panels sit at different concentrations. Sets that belong in one figure share a `scale_group`
    and therefore one axis; sets in different groups are scaled independently, since a figure assembled
    from one set should not be stretched by a set it will never sit beside.
    """
    sets = annotation("boxb_reaction_sets")
    # A panel is titled by whatever distinguishes it within its set. That is the enzyme in every
    # set but the time course, whose two panels are one enzyme at two incubation times, so a set
    # may name its panels itself and falls back to the enzyme name when it does not.
    sets["panel_title"] = sets.panel_title if "panel_title" in sets else np.nan
    reactions = annotation("sample_info")
    require(sets[["reaction_set", "scale_group", *SAMPLE_KEYS]].notna().all().all() and not sets.sample_id.duplicated().any(), "Invalid boxB reaction sets")
    require(len(sets.merge(reactions, on=SAMPLE_KEYS[:4])) == len(sets), "boxB reaction sets disagree with the reaction annotation")
    require(sets.groupby("reaction_set").size().eq(EXPECTED_SET_ENZYMES).all(), "Every reaction set must carry two reactions")
    require(sets.groupby("reaction_set").scale_group.nunique().eq(1).all(), "A reaction set must sit in one scale group")
    require(sets.groupby("reaction_set").panel_order.apply(lambda order: sorted(order) == list(range(EXPECTED_SET_ENZYMES))).all(), "Panel order must number the panels of each set")
    require(sets.assign(named=sets.panel_title.fillna(sets.tada_type)).groupby("reaction_set").named.nunique().eq(EXPECTED_SET_ENZYMES).all(), "The panels of a set must be told apart by their titles")
    require(sets.groupby("tada_type")[["enzyme_label"]].nunique().eq(1).all().all(), "Conflicting enzyme presentation across reaction sets")
    return sets


def reaction_set(set_id):
    """One set's reactions, in panel order, in the sample columns the pooled effect calculation reads."""
    sets = boxb_reaction_sets()
    selected = sets.loc[sets.reaction_set.eq(set_id), SAMPLE_KEYS + ["panel_title"]].sort_values("panel_order")
    require(len(selected) == EXPECTED_SET_ENZYMES, f"Unannotated boxB reaction set {set_id}")
    return selected


def reaction_set_ids(scale_group=None):
    """Every declared set, or those of one scale group, in the order the annotation lists them."""
    sets = boxb_reaction_sets()
    if scale_group is not None:
        sets = sets.loc[sets.scale_group.eq(scale_group)]
    return sets.reaction_set.drop_duplicates().tolist()


def scale_groups():
    """Each scale group and the sets that share its axis."""
    sets = boxb_reaction_sets()
    return {group: frame.reaction_set.drop_duplicates().tolist() for group, frame in sets.groupby("scale_group")}


def coverage_of(column):
    """Each measurement's denominator. Recorder bases below the counter's quality floor are N, so a
    single site drops only its own uncallable calls, and the edit-count bins take the molecules whose
    whole recorder is callable."""
    if column.startswith("pos_"):
        return f"coverage_{column}"
    return "num_0_N"


def counts_for(sample_id, columns, directory=COUNT_DIRECTORY):
    sites = [column for column in columns if column.startswith("pos_")]
    uncallable = [f"pos_{column.split('_')[1]}_N" for column in sites]
    requested = list(dict.fromkeys(["sample_id", "barcode", "insert", "umi_counts", "num_0_N", *columns, *uncallable]))
    data = pd.read_csv(f"{directory}/{sample_id}.csv.gz", usecols=requested)
    require(data.sample_id.eq(sample_id).all(), f"Unexpected sample identifiers in {sample_id}")
    require(not data.duplicated(["sample_id", "barcode", "insert"]).any(), f"Duplicate constructs in {sample_id}")
    for column, missing in zip(sites, uncallable):
        data[f"coverage_{column}"] = data.umi_counts - data[missing]
    return data


def valid_counts(data, columns, positive_coverage=True):
    """Counts must be finite integers lying inside the molecules each measurement is defined on.
    Pass positive_coverage=False for constructs the analysis has not yet narrowed to the ones it
    divides: a construct too thin to carry one fully callable molecule is not yet a zero divisor."""
    denominators = list(dict.fromkeys(coverage_of(column) for column in columns))
    matrix = data[["umi_counts", *columns, *denominators]].to_numpy()
    require(np.isfinite(matrix).all() and (matrix == np.floor(matrix)).all(), "Counts must be finite integers")
    require(data.umi_counts.gt(0).all(), "Construct UMI coverage must be positive")
    require(data[denominators].ge(0).all().all() and data[denominators].le(data.umi_counts, axis=0).all().all(), "Callable coverage must lie within the molecule count")
    if positive_coverage:
        require(data[denominators].gt(0).all().all(), "Callable coverage must be positive where it is a denominator")
    for column in columns:
        require((data[column].ge(0) & data[column].le(data[coverage_of(column)])).all(), f"Edited counts exceed callable coverage for {column}")


def barcode_metadata():
    data = annotation("barcode_annotations")
    data = data.loc[data.variable_type.eq(LIBRARY_TYPE) & data.g_depleted.eq(BACKGROUND)].copy()
    data = data.drop(columns="barcode").rename(columns={"reverse_complement": "barcode"})
    data["layout"] = data.oligo_name.str.extract(LAYOUT_PATTERN, expand=False)
    require(not data.barcode.duplicated().any() and data.layout.notna().all(), "Invalid barcode annotations")
    return data


def matched_stems(samples, columns, directory=COUNT_DIRECTORY):
    """WT is pooled over all six windows; only the stem windows carry paired mutants."""
    barcodes = barcode_metadata()[["barcode", "variable_subpos", "layout", *GEOMETRY_KEYS]]
    listed = annotation("boxb_wildtype_mutant_list")[["variable_subpos", "mutant_bases", "insert_type", "insert"]]
    require(not listed.duplicated().any(), "Duplicate construct annotations")
    require(listed.loc[listed.mutant_bases.isna(), "insert_type"].eq("wt").all(), "Unpaired constructs must be wild type")
    frames = []
    for sample_id in samples.sample_id:
        data = counts_for(sample_id, columns, directory=directory).merge(barcodes, on="barcode", validate="many_to_one")
        data = data.merge(listed, on=["variable_subpos", "insert"], validate="many_to_many")
        data = data.loc[~(data.layout.eq(EXCLUDED_LAYOUT) & data.variable_subpos.isin(EXCLUDED_WINDOWS))].copy()
        valid_counts(data, columns)
        data["window_id"] = data.layout + "__" + data.variable_subpos
        data["pair_id"] = data.window_id + "__" + data.mutant_bases
        pairs = data.loc[data.pair_id.notna()].pivot(index="pair_id", columns="insert_type", values="umi_counts")
        require(pairs.shape == (EXPECTED_PAIRS, 2) and pairs.notna().all().all(), f"Incomplete stem pairs in {sample_id}")
        unique = data.drop_duplicates(["barcode", "insert"])
        require(unique.insert_type.value_counts().to_dict() == {"mut": EXPECTED_PAIRS, "wt": EXPECTED_WT}, "Unexpected unique WT/MUT counts")
        require(unique.loc[unique.insert_type.eq("wt"), "variable_subpos"].nunique() == EXPECTED_WT_WINDOWS, "Unexpected wild-type windows")
        frames.append(data)
    return pd.concat(frames, ignore_index=True).merge(samples, on="sample_id", validate="many_to_one")


def recorder_site_tests(data, plotted=False):
    """Paired construct tests per recorder site, reported without multiple-testing correction.

    With `plotted`, the WT/MUT test instead compares exactly the constructs a figure draws, every
    distinct WT against every distinct MUT, by an unpaired two-sided Welch t-test, so the stars test
    the same means and the same construct counts the points and error bars show.

    WT/MUT tests average mutants sharing a barcode before pairing each mean with its
    wild type. Each of the eight window-placement groups contributes once. The construct-pair
    P values are reported beside the barcode-mean ones.
    Orientation (5′ = mean of 0 and 10 nt versus 3′ = mean of 20 and 30 nt): pair identical window/insert
    constructs that carry all four placements, separately for WT and MUT.
    Neither constructs nor placements are biological replicates.
    """
    keys = ["sample_id", "position_id"]
    unique = data.drop_duplicates(keys + ["barcode", "insert"])
    require(unique.percent_edited.between(0, 100).all(), "Invalid site editing percentages")
    matched = unique.loc[unique.mutant_bases.notna()]
    paired = matched.groupby(keys + ["barcode", "insert_type"]).percent_edited.mean().unstack("insert_type").reset_index()
    # Each mutant takes the wild type sharing its barcode, which fixes its window and placement, so the
    # construct pairing holds the four mutants of a barcode against one reference rather than averaging.
    references = matched.loc[matched.insert_type.eq("wt")].groupby(keys + ["barcode"]).percent_edited.mean()
    by_construct = matched.loc[matched.insert_type.eq("mut")].rename(columns={"percent_edited": "mut"})
    by_construct["wt"] = by_construct.set_index(keys + ["barcode"]).index.map(references)
    tests = []
    for values, group in paired.groupby(keys):
        require(group[["wt", "mut"]].notna().all().all(), "Incomplete WT/MUT site pairing")
        constructs = by_construct
        for key, value in zip(keys, values):
            constructs = constructs.loc[constructs[key].eq(value)]
        require(len(constructs) == EXPECTED_PAIRS and constructs[["wt", "mut"]].notna().all().all(),
                "Incomplete WT/MUT construct pairing")
        both = {"construct_pairs": (constructs.wt, constructs.mut), "barcode_means": (group.wt, group.mut)}
        p_values = {f"p_value_{name}": ttest_rel(*values_pair).pvalue for name, values_pair in both.items()}
        first, second = both["barcode_means"]
        if plotted:
            drawn = unique.loc[unique.sample_id.eq(values[0]) & unique.position_id.eq(values[1])]
            first, second = (drawn.loc[drawn.insert_type.eq(variant), "percent_edited"] for variant in ["wt", "mut"])
        tests.append((dict(zip(keys, values), **p_values), "WT/MUT", first, second))
    (near_name, near), (far_name, far) = ORIENTATION_GROUPS.items()
    for values, group in unique.groupby(keys + ["insert_type"]):
        wide = group.pivot(index=["variable_subpos", "insert"], columns="target_dist", values="percent_edited")
        require(set(near + far) <= set(wide.columns), "Missing requested distance")
        common = wide.dropna(subset=near + far)
        tests.append((dict(zip(keys + ["insert_type"], values), n_unmatched_constructs=len(wide) - len(common)),
                      f"{values[-1].upper()} {near_name}/{far_name}", common[near].mean(axis=1), common[far].mean(axis=1)))
    records = []
    for metadata, comparison, first, second in tests:
        unpaired = plotted and comparison == "WT/MUT"
        require(min(len(first), len(second)) >= 2, "Insufficient matched constructs for a site comparison")
        result = ttest_ind(first, second, equal_var=False) if unpaired else ttest_rel(first, second)
        require(np.isfinite(result.pvalue), "Non-finite recorder-site test")
        records.append(dict(metadata, comparison=comparison, n_pairs=np.nan if unpaired else len(first),
                            **({"n_first": len(first), "n_second": len(second)} if plotted else {}),
                            first_mean_percent=first.mean(), second_mean_percent=second.mean(),
                            difference_percentage_points=first.mean() - second.mean(),
                            statistic=result.statistic, p_value=result.pvalue))
    statistics = pd.DataFrame(records)
    statistics["stars"] = np.select([statistics.p_value.lt(cutoff) for cutoff, _ in SITE_STARS],
                                    [stars for _, stars in SITE_STARS], default="ns")
    return statistics


def boxb_hairpin():
    """The boxB sense sequence and its annotated regions, checked to tile the hairpin end to end."""
    regions = annotation("boxb_structure_regions").sort_values("start", ignore_index=True)
    sequence = "".join(regions.rna_sequence)
    require(len(sequence) == BOXB_LENGTH and regions.start.iloc[0] == 1 and regions.end.iloc[-1] == BOXB_LENGTH
            and regions.start.iloc[1:].eq(regions.end.iloc[:-1].to_numpy() + 1).all(), "boxB regions must tile the hairpin")
    require((regions.end - regions.start + 1).eq(regions.rna_sequence.str.len()).all(), "A boxB region's span and sequence disagree")
    return sequence, regions


def loop_region():
    """The annotated loop, whose bases every loop measurement reads."""
    _, regions = boxb_hairpin()
    loop = regions.loc[regions.region.eq("loop")]
    require(len(loop) == 1, "Expected one annotated loop")
    return loop.iloc[0]


def hairpin_positions():
    """Every adenosine of the hairpin, its RNA context and region, and the counter's column for it.

    A read offset r carries boxB sense position BOXB_LENGTH - r, the origin rule the window
    coordinates follow, so the counter's offsets fix which column measures which adenosine. Checked
    against reads: in window 10-13 the dominant base at pos_0, pos_1 and pos_3 tracks the designed
    base at sense 12, 10 and 13 in every construct, while pos_2, at sense 9, stays A throughout.
    """
    sequence, regions = boxb_hairpin()
    records = []
    for column_index, offset in enumerate(BOXB_READ_OFFSETS):
        position = BOXB_LENGTH - offset
        require(1 < position < BOXB_LENGTH and sequence[position - 1] == "A", "A measured column holds no interior adenosine")
        region = regions.loc[regions.start.le(position) & regions.end.ge(position)]
        require(len(region) == 1, "A measured adenosine falls outside the annotated regions")
        records.append(dict(position_id=f"A{position}", boxb_position=position, region=region.iloc[0].region,
                            region_label=region.iloc[0].label, context=sequence[position - 2:position + 1],
                            editing_column=f"pos_{column_index}_C"))
    positions = pd.DataFrame(records).sort_values("boxb_position", ignore_index=True)
    positions["position_index"] = positions.index + 1
    require(set(positions.boxb_position) == {index + 1 for index, base in enumerate(sequence) if base == "A"},
            "The measured columns miss an adenosine of the hairpin")
    return positions


def measured_windows(windows, positions):
    """The randomized windows a per-adenosine measurement can still be read in.

    A window covering the loop or a measured adenosine puts its own designed base where a measurement
    reads, so there an edit cannot be told from a substitution. Position 13 falls only in window
    10-13, which the loop already excludes, so measuring it costs no construct.
    """
    bounds = windows.str.split("_", expand=True).astype(int)
    loop = loop_region()
    covered = bounds[0].le(loop.end) & bounds[1].ge(loop.start)
    for position in positions.boxb_position:
        covered |= bounds[0].le(position) & bounds[1].ge(position)
    return ~covered


def hairpin_site_figures():
    """The reactions each per-adenosine figure draws, one row per panel, from the versioned selection.

    Each figure opens on a reaction carrying no TadA, which is the measurement's floor rather than a
    formality: the wild-type/mutant contrast at the AAG-context adenosines is reproduced in full with
    no enzyme present, so a series starting at the lowest enzyme concentration reads a property of the
    mutant constructs as editing. `panel_row` and `panel_column` place each panel, so one figure is a
    row of concentrations and another a grid of editors over concentrations; a row holding fewer panels
    than the widest is centred, which is where a no-TadA reaction two columns share sits.
    """
    figures = annotation("boxb_hairpin_site_figures")
    require(figures.notna().all().all() and not figures.duplicated(["figure", "sample_id"]).any(),
            "Invalid hairpin site figure selection")
    reactions = annotation("sample_info")
    panels = figures.merge(reactions, on=SAMPLE_KEYS[:4], validate="many_to_one")
    require(len(panels) == len(figures), "A hairpin site panel names an unannotated reaction")
    for figure, group in panels.groupby("figure"):
        controls = group.loc[group.tada_type.eq(CONTROL_ENZYME)]
        require(len(controls) == 1 and controls.panel_row.eq(0).all() and controls.panel_column.eq(0).all(),
                f"{figure} needs one no-TadA panel, first in its grid")
        editors = group.loc[~group.tada_type.eq(CONTROL_ENZYME)]
        require(editors.tada_type.eq(editors.enzyme).all(), f"{figure} gives a panel another editor's colour")
        require(controls.enzyme.isin(set(editors.enzyme) | {CONTROL_ENZYME}).all(),
                f"{figure} control panel names an editor the figure does not draw")
        require(sorted(group.panel_row.unique()) == list(range(group.panel_row.max() + 1)),
                f"{figure} panel rows are not contiguous")
        require(editors.groupby("panel_column").enzyme.nunique().eq(1).all(),
                f"{figure} draws two editors in one column")
        for row, across in group.groupby("panel_row"):
            require(sorted(across.panel_column) == list(range(len(across))),
                    f"{figure} row {row} columns are not contiguous")
            require(row == 0 or across.panel_label.nunique() == 1, f"{figure} row {row} carries two labels")
    return panels.sort_values(["figure", "panel_row", "panel_column"], ignore_index=True)


def hairpin_site_data(reactions):
    """Per-adenosine boxB editing, with the four recorder placements pooled as replicates.

    Every adenosine of the hairpin is measured, the three in the loop and the one closing the stem.
    Windows reaching a measured adenosine are dropped, as they are from the loop edit counts: there
    the designed base is the base a measurement reads, so a substitution cannot be told from an edit.
    Pooling makes each construct-placement one observation of its adenosine. Placements are not
    biological replicates, so `recorder_site_tests` keeps the 5'/3' contrast that says how far
    pooling them is fair. One call draws one figure, so a reaction shared between figures, such as the
    no-enzyme control, is summarised within each rather than once across all of them.
    """
    positions = hairpin_positions()
    require(not reactions.sample_id.duplicated().any(), "A hairpin site figure repeats a reaction")
    constructs = matched_stems(reactions, positions.editing_column.tolist(), directory=LOOP_COUNT_DIRECTORY)
    constructs = constructs.loc[measured_windows(constructs.variable_subpos, positions)].copy()
    require(sorted(constructs.target_dist.unique()) == EXPECTED_DISTANCES, "Invalid recorder placements")
    unique = constructs.drop_duplicates(["sample_id", "barcode", "insert"])
    # Each adenosine carries its own callable denominator, so build the long frame one site at a time.
    pieces = []
    for column in positions.editing_column:
        piece = unique.drop(columns=positions.editing_column.tolist()).copy()
        piece["editing_column"] = column
        piece["edited_count"] = unique[column].to_numpy()
        piece["callable_umis"] = unique[coverage_of(column)].to_numpy()
        pieces.append(piece)
    data = pd.concat(pieces, ignore_index=True).merge(positions, on="editing_column", validate="many_to_one")
    data["percent_edited"] = 100 * data.edited_count / data.callable_umis
    keys = ["sample_id", "panel_label", "panel_row", "panel_column", "tada_type", "tada_conc",
            "position_id", "context", "position_index", "insert_type"]
    summary = data.groupby(keys).agg(mean_percent=("percent_edited", "mean"), se_percent=("percent_edited", "sem"),
                                     n=("percent_edited", "size"), edited_umis=("edited_count", "sum"),
                                     callable_umis=("callable_umis", "sum")).reset_index()
    barcodes = barcode_metadata()
    retained = measured_windows(barcodes.variable_subpos, positions)
    retained &= ~(barcodes.layout.eq(EXCLUDED_LAYOUT) & barcodes.variable_subpos.isin(EXCLUDED_WINDOWS))
    expected = np.where(summary.insert_type.eq("wt"), int(retained.sum()), EXPECTED_PAIRS)
    require(len(summary) == 2 * len(positions) * len(reactions) and summary.n.eq(expected).all(), "Unexpected hairpin construct counts")
    require(summary.mean_percent.between(0, 100).all() and np.isfinite(summary.se_percent).all(), "Invalid hairpin site summaries")
    return summary, positions, data


def comparison_definitions():
    # A panel names the figure it belongs to, so the conditions of one enzyme are drawn as the rows of
    # one figure rather than as a figure each; a panel drawn on its own names itself.
    panels = annotation("recorder_edit_summary_panels").fillna({"row_label": ""})
    sites = annotation("recorder_edit_summary_sites")
    definitions = []
    for row in panels.to_dict("records"):
        if "edit_count" in row["measurement_modes"].split(";"):
            definitions.append(dict(row, measurement_mode="edit_count"))
    for row in sites.to_dict("records"):
        panel = panels.loc[panels.name.eq(row["comparison"])]
        require(len(panel) == 1, "Site comparison has no unique panel definition")
        definitions.append(dict(panel.iloc[0], **row, name=row["output_name"], measurement_mode="site", width=row["site_width"], height=row["site_height"]))
    return definitions


def comparison_data():
    definitions = comparison_definitions()
    samples = annotation("sample_info")
    require(not samples.sample_id.duplicated().any(), "Duplicate reaction annotations")
    panel_samples = {}
    for definition in definitions:
        selected = samples.loc[samples.short_key.isin(definition["short_keys"].split(";")) & samples.condition.isin(definition["conditions"].split(";")) & samples.tada_conc.isin(definition["concentrations"].split(";"))]
        require(len(selected) == definition["expected_samples"], f"Unexpected reactions for {definition['name']}")
        require(set(selected.tada_type) == set(definition["enzymes"].split(";")) and set(selected[definition["column"]]) == set(definition["levels"].split(";")), "Incomplete comparison conditions")
        panel_samples[definition["name"]] = selected.sample_id.tolist()
    positions = annotation("recorder_positions")
    needed_sites = positions.loc[positions.position_id.isin(annotation("recorder_edit_summary_sites").position_id)]
    columns = [*[f"num_{i}_C" for i in range(1, 8)], *needed_sites.editing_column]
    ids = set(itertools.chain.from_iterable(panel_samples.values()))
    constructs = matched_stems(samples.loc[samples.sample_id.isin(ids)], columns)
    require(constructs[[f"num_{i}_C" for i in range(1, 8)]].sum(axis=1).le(constructs.num_0_N).all(), "Edit-count bins exceed callable coverage")
    # Each measurement carries its own numerator and denominator beside its fraction, so the drawn
    # percentage can state the molecules behind it without a second pass over the counts.
    constructs["edited_frac_1edit"] = constructs.num_1_C
    constructs["edited_frac_2edit"] = constructs[[f"num_{i}_C" for i in range(2, 8)]].sum(axis=1)
    constructs["callable_frac_1edit"] = constructs.num_0_N
    constructs["callable_frac_2edit"] = constructs.num_0_N
    constructs["frac_1edit"] = constructs.edited_frac_1edit / constructs.callable_frac_1edit
    constructs["frac_2edit"] = constructs.edited_frac_2edit / constructs.callable_frac_2edit
    counted = [f"{prefix}{edit}" for edit in EDIT_TYPES + ["frac_site"] for prefix in ["edited_", "callable_"]]
    results = {}
    for definition in definitions:
        data = constructs.loc[constructs.sample_id.isin(panel_samples[definition["name"]])].copy()
        edits = EDIT_TYPES
        if definition["measurement_mode"] == "site":
            site = positions.loc[positions.position_id.eq(definition["position_id"])]
            require(len(site) == 1 and site.iloc[0].context == definition["dna_context"].replace("T", "U"), "Invalid site/context selection")
            data["edited_frac_site"] = data[site.iloc[0].editing_column]
            data["callable_frac_site"] = data[coverage_of(site.iloc[0].editing_column)]
            data["frac_site"] = data.edited_frac_site / data.callable_frac_site
            edits = ["frac_site"]
        identifiers = [c for c in data.columns if c not in edits and c not in counted]
        long = data.melt(id_vars=identifiers, value_vars=edits, var_name="edit_type", value_name="fraction_edited")
        # Melting the numerators and denominators on their own keeps each one beside its own fraction:
        # melt stacks the chosen columns in order, each block holding the frame's own row order.
        long["edited_umis"] = data.melt(value_vars=[f"edited_{edit}" for edit in edits]).value.to_numpy()
        long["callable_umis"] = data.melt(value_vars=[f"callable_{edit}" for edit in edits]).value.to_numpy()
        data = long
        keys = list(dict.fromkeys(["tada_type", "tada_conc", definition["column"], "edit_type"]))
        unique = data.drop_duplicates(["sample_id", "barcode", "insert", "edit_type"])
        summary = unique.groupby(keys + ["insert_type"], sort=True).agg(
            n=("fraction_edited", "size"), mean=("fraction_edited", "mean"), se=("fraction_edited", "sem"),
            edited_umis=("edited_umis", "sum"), callable_umis=("callable_umis", "sum")).reset_index()
        summary[["mean", "se"]] *= 100
        matched = data.loc[data.pair_id.notna()]
        paired = matched.pivot(index=keys + ["pair_id"], columns="insert_type", values="fraction_edited").reset_index()
        windows = matched.groupby(keys + ["window_id", "insert_type"]).fraction_edited.mean().unstack("insert_type").reset_index()
        tests = []
        for values, group in paired.groupby(keys):
            window = windows
            for key, value in zip(keys, values):
                window = window.loc[window[key].eq(value)]
            require(len(group) == EXPECTED_PAIRS and len(window) == EXPECTED_PAIRED_WINDOWS and group[["wt", "mut"]].notna().all().all(), "Incomplete statistical pairing")
            tests.append(dict(zip(keys, values), p_value=ttest_rel(window.wt, window.mut).pvalue, p_value_window_mean=ttest_rel(window.wt, window.mut).pvalue, p_value_construct_pairs=ttest_rel(group.wt, group.mut).pvalue, n_pairs=len(window)))
        tests = pd.DataFrame(tests)
        for source in ["p_value", "p_value_window_mean"]:
            require(np.isfinite(tests[source]).all(), "Non-finite paired test")
        tests["significance"] = np.select([tests.p_value.lt(.001), tests.p_value.lt(.01), tests.p_value.lt(.05)], ["***", "**", "*"], default="ns")
        statistics = summary.pivot(index=keys, columns="insert_type", values=["mean", "se", "n"])
        statistics.columns = [f"{variant}_{value}" for value, variant in statistics.columns]
        statistics = statistics.reset_index().merge(tests, on=keys, validate="one_to_one")
        statistics["n_construct_pairs"] = EXPECTED_PAIRS
        statistics["n_wt_constructs"] = statistics.wt_n.astype(int)
        statistics["difference_percentage_points"] = statistics.wt_mean - statistics.mut_mean
        statistics["fold_change"] = statistics.wt_mean / statistics.mut_mean
        statistics = statistics.rename(columns={f"{variant}_{value}": f"{variant}_{value}_percent" for variant in ["wt", "mut"] for value in ["mean", "se"]})
        require(np.isfinite(statistics.select_dtypes("number")).all().all(), "Non-finite comparison summary")
        results[definition["name"]] = (summary, statistics)
    return definitions, results


def distance_data():
    metadata = enzyme_metadata()
    positions = annotation("recorder_positions").copy()
    positions["position_index"] = positions.position_id.str[1:].astype(int).rank(method="dense").astype(int)
    require(len(positions) == 8 and not positions.isna().any().any(), "Invalid recorder positions")
    concentrations = annotation("recorder_edit_summary_panels").iloc[0].concentrations.split(";")
    samples = annotation("sample_info")
    selected = samples.loc[samples.tada_type.isin(metadata.tada_type) & samples.tada_conc.isin(concentrations) & samples.condition.eq(metadata.iloc[0].condition)]
    require(len(selected) == len(concentrations) * len(metadata), "Incomplete distance reactions")
    constructs = matched_stems(selected, positions.editing_column.tolist())
    require(sorted(constructs.target_dist.unique()) == EXPECTED_DISTANCES, "Invalid recorder placements")
    unique = constructs.drop_duplicates(["sample_id", "barcode", "insert"])
    # Each site carries its own callable denominator, so build the long frame one site at a time.
    pieces = []
    for column in positions.editing_column:
        piece = unique.drop(columns=positions.editing_column.tolist()).copy()
        piece["editing_column"] = column
        piece["edited_count"] = unique[column].to_numpy()
        piece["callable_umis"] = unique[coverage_of(column)].to_numpy()
        pieces.append(piece)
    data = pd.concat(pieces, ignore_index=True).merge(positions, on="editing_column", validate="many_to_one")
    data["percent_edited"] = 100 * data.edited_count / data.callable_umis
    keys = ["tada_type", "tada_conc", "target_dist", "position_id", "context", "position_index", "insert_type"]
    summary = data.groupby(keys).agg(mean_percent=("percent_edited", "mean"), se_percent=("percent_edited", "sem"),
                                     n=("percent_edited", "size"), edited_umis=("edited_count", "sum"),
                                     callable_umis=("callable_umis", "sum")).reset_index()
    excluded = barcode_metadata().loc[barcode_metadata().layout.eq(EXCLUDED_LAYOUT), "target_dist"].unique()
    require(len(excluded) == 1, "Excluded layout spans multiple placements")
    expected_wt = np.where(summary.target_dist.eq(excluded[0]), EXPECTED_WT_WINDOWS - len(EXCLUDED_WINDOWS), EXPECTED_WT_WINDOWS)
    expected = np.where(summary.insert_type.eq("wt"), expected_wt, EXPECTED_SITE_MUT)
    require(len(summary) == 128 * len(concentrations) and summary.n.eq(expected).all(), "Unexpected distance construct counts")
    return summary, metadata, positions.sort_values("position_index"), data


def context_data():
    """Pool geometries and backgrounds independently for each selected recorder site."""
    selection = annotation("recorder_context_heatmap_selection")
    sites = annotation("recorder_context_selection").drop(columns=GEOMETRY_KEYS)
    positions = annotation("recorder_positions").rename(columns={"context": "reference_context"})
    choices = selection.merge(sites, on=["sample_id", "position_id"], how="left", validate="many_to_one").merge(positions, on="position_id", how="left", validate="many_to_one")
    required = [c for c in choices if c not in ["point_color", "marker"]]
    site_count = selection.position_id.nunique()
    sample_count = selection.sample_id.nunique()
    expected_choices = sample_count * site_count * len(EXPECTED_DISTANCES)
    require(site_count > 0 and len(choices) == expected_choices and not choices[required].isna().any().any() and not selection.duplicated().any(), "Incomplete context selections")
    require(set(selection.sample_id) == set(enzyme_metadata().sample_id) and selection.groupby(["sample_id", "position_id"]).size().eq(len(EXPECTED_DISTANCES)).all(), "Incomplete context reactions or placements")
    samples = annotation("sample_info").set_index("sample_id")
    barcodes = annotation("barcode_annotations")
    geometries = barcodes.loc[barcodes.variable_type.eq(LIBRARY_TYPE), GEOMETRY_KEYS + ["target_length"]].drop_duplicates()
    require(len(geometries) == 4 and geometries.target_length.nunique() == 1, "Invalid context geometries")
    backgrounds = {"".join(bases) for bases in itertools.product(DNA_BASES, repeat=3)}
    frames = []
    for sample_id, sample_choices in choices.groupby("sample_id"):
        sample = samples.loc[sample_id]
        require(all(sample_choices[key].eq(sample[key]).all() for key in ["tada_type", "tada_conc", "condition"]), "Context reaction annotation mismatch")
        counts = counts_for(sample_id, sample_choices.editing_column.unique())
        for choice in sample_choices.to_dict("records"):
            barcode = barcodes
            for key in ["variable_type", "variable_subpos", *GEOMETRY_KEYS]:
                barcode = barcode.loc[barcode[key].eq(choice[key])]
            require(len(barcode) == 1 and barcode.iloc[0].target_length == geometries.iloc[0].target_length, "Context choice needs one matching barcode")
            data = counts.loc[counts.barcode.eq(barcode.iloc[0].reverse_complement)].copy()
            column = choice["editing_column"]
            valid_counts(data, [column])
            require(len(data) > 0 and not data["insert"].duplicated().any() and data["insert"].str.fullmatch(r"[ACGT]{5}").all(), "Invalid context inserts")
            flank = [int(choice["insert_fiveprime_index"]) - 1, int(choice["insert_threeprime_index"]) - 1]
            require(len(set(flank)) == 2 and min(flank) >= 0 and max(flank) < 5, "Invalid context flank indices")
            other = [i for i in range(5) if i not in flank]
            data["rna_context"] = data["insert"].str[flank[0]].map(READ_TO_RNA) + choice["reference_context"][1] + data["insert"].str[flank[1]].map(READ_TO_RNA)
            data["background"] = data["insert"].map(lambda value: "".join(value[i] for i in other))
            require(data.rna_context.nunique() == 16 and all(len(group) == 64 and set(group.background) == backgrounds for _, group in data.groupby("rna_context")), "Incomplete context backgrounds")
            data["callable_umis"] = data[coverage_of(column)]
            pooled = data.groupby("rna_context").agg(n_variants=("insert", "size"), total_umis=("callable_umis", "sum"), total_edited_umis=(column, "sum"), minimum_umis=("callable_umis", "min"), maximum_umis=("callable_umis", "max")).reset_index()
            for key, value in choice.items():
                pooled[key] = value
            pooled["geometry"] = f"{choice['target_pos_to_boxb']}prime_{choice['target_dist']}nt"
            frames.append(pooled)
    site_summary = pd.concat(frames, ignore_index=True)
    keys = ["sample_id", "tada_type", "tada_conc", "condition", "enzyme_label", "panel_order", "position_id", "rna_context"]
    summary = site_summary.groupby(keys).agg(n_sites=("position_id", "nunique"), n_geometries=("geometry", "nunique"), n_site_variants=("n_variants", "sum"), total_umis=("total_umis", "sum"), total_edited_umis=("total_edited_umis", "sum")).reset_index().sort_values(["panel_order", "rna_context"])
    require(len(site_summary) == expected_choices * 16 and len(summary) == sample_count * site_count * 16 and summary.n_sites.eq(1).all() and summary.n_geometries.eq(len(EXPECTED_DISTANCES)).all() and summary.n_site_variants.eq(len(EXPECTED_DISTANCES) * 64).all(), "Incomplete site-specific context tiles")
    summary["fiveprime"] = summary.rna_context.str[0]
    summary["threeprime"] = summary.rna_context.str[2]
    summary["pooled_fraction_central_a_edited"] = summary.total_edited_umis / summary.total_umis
    return summary, site_summary


def drawn_at(selections, set_id):
    """Keep the selections a reaction set draws, from the sets each one declares.

    A measurement is scoped by set because editing falls with enzyme, and a measurement whose rarest
    sequences reach zero edited molecules there carries no ratio to plot. Dropping it from that set is
    a declared selection, not a silent gap: the sets it keeps are in its own row.
    """
    declared = selections.reaction_sets.str.split(";")
    require(declared.explode().isin(boxb_reaction_sets().reaction_set.tolist()).all(), "Unannotated reaction set in a boxB selection")
    return selections.loc[declared.map(lambda sets: set_id in sets)].reset_index(drop=True)


def substitution_measurements(set_id=None):
    measurements = annotation("boxb_substitution_measurements")
    positions = annotation("recorder_positions").set_index("position_id", verify_integrity=True)
    require(not measurements.position_id.duplicated().any() and measurements[["position_id", "sites", "label", "reaction_sets", "combine_sites"]].notna().all().all(), "Invalid substitution measurements")
    require(measurements.combine_sites.isin(COMBINE_KINDS).all(), "Every measurement must say how its sites combine")
    for index, measurement in measurements.iterrows():
        members = measurement.sites.split(MEASUREMENT_SEPARATOR)
        require(len(set(members)) == len(members) and set(members).issubset(positions.index), "Invalid recorder sites")
        require((len(members) == 1) == (measurement.combine_sites == "site"), "Only a single-site measurement combines as one site")
        require(len(members) > 1 or measurement.position_id == members[0], "Single-site identifier mismatch")
        measurements.loc[index, "context"] = " & ".join(positions.loc[members, "context"])
        measurements.loc[index, "axis_label"] = (f"Mean {'/'.join(members)} editing" if measurement.combine_sites == MEAN_PREFIX
                                                 else f"Editing at {measurement.label}")
        measurements.loc[index, "editing_column"] = (positions.loc[members[0], "editing_column"] if measurement.combine_sites == "site"
                                                     else measurement.position_id)
    if set_id is None:
        return measurements
    selected = drawn_at(measurements, set_id)
    require(len(selected) > 0, f"No substitution measurement is drawn at {set_id}")
    return selected


def site_measurement(position_id):
    """One single-site measurement, in the columns the pooled effect calculation reads."""
    positions = annotation("recorder_positions").set_index("position_id", verify_integrity=True)
    require(position_id in positions.index, f"Unannotated recorder site {position_id}")
    site = positions.loc[position_id]
    return pd.DataFrame([{"position_id": position_id, "sites": position_id, "label": position_id,
                          "context": site.context, "editing_column": site.editing_column, "combine_sites": "site",
                          "axis_label": f"Editing at {position_id}"}])


def mean_measurement(sites):
    """A derived measurement: the geometric mean of several single-site effects.

    The sites are two adenosines of one recorder, read from the same molecules, so one sequence's
    effect is measured twice rather than once. Averaging them halves the counting noise without
    averaging away anything the sites disagree on structurally: across the 19 boxB positions the two
    report the same picture, at r = 0.94 to 0.99 on the position means the heat map and logo draw.
    """
    members = sites.split(MEASUREMENT_SEPARATOR)
    require(len(members) > 1 and len(set(members)) == len(members), "A site average needs distinct sites")
    positions = annotation("recorder_positions").set_index("position_id", verify_integrity=True)
    require(set(members).issubset(positions.index), "Unannotated recorder site in a site average")
    position_id = f"{MEAN_PREFIX}_{'_'.join(members)}"
    return pd.DataFrame([{"position_id": position_id, "sites": sites, "label": f"{'/'.join(members)} mean",
                          "context": " & ".join(positions.loc[members, "context"]),
                          "editing_column": position_id, "combine_sites": "mean",
                          "axis_label": f"Mean {'/'.join(members)} editing"}])


def selected_measurement(position_id):
    """Resolve the measurement a figure's selection names: one recorder site, or a mean of several.

    A mean is named for the sites it averages, `mean_A8_A2`, so one annotation column still carries the
    whole choice and a figure moves between a single site and a mean by that cell alone.
    """
    if position_id.startswith(f"{MEAN_PREFIX}_"):
        return mean_measurement(MEASUREMENT_SEPARATOR.join(position_id[len(MEAN_PREFIX) + 1:].split("_")))
    return site_measurement(position_id)


def site_average(summary, measurement):
    """Average several single-site effects of one sequence into the measurement that derives them.

    Every effect in this analysis is a ratio to matched wild type carried on a log₂ axis, so the mean
    that belongs on that axis is the geometric one: the plain mean of the log₂ ratios. The sites are
    read from the same molecules, so their molecule counts are not independent and are never added.
    Each support figure is carried as the weaker of the sites instead, which keeps a cell honest about
    the thinnest evidence under it, and a sequence with no edited molecule at either site keeps the
    −∞ that its own site already carries rather than being averaged into a finite value.
    """
    members = measurement.sites.split(MEASUREMENT_SEPARATOR)
    chosen = summary.loc[summary.position_id.isin(members)]
    keys = SAMPLE_KEYS + ["hairpin", "n_mutations"]
    require(chosen.groupby(keys).position_id.nunique().eq(len(members)).all(), "A site average is missing one of its sites")
    # The sites are measured on one construct set, so everything describing that set must already agree.
    shared = ["n_ratio_geometries", "n_ratio_windows", "n_ratio_constructs", "ratio_windows"]
    require(chosen.groupby(keys)[shared].nunique().eq(1).all().all(), "The sites of an average rest on different constructs")
    averaged = chosen.groupby(keys, as_index=False).agg(
        log2_pooled_ratio_to_matched_wt=("log2_pooled_ratio_to_matched_wt", "mean"),
        pooled_umis=("pooled_umis", "min"), pooled_edited_umis=("pooled_edited_umis", "min"),
        pooled_matched_wt_umis=("pooled_matched_wt_umis", "min"),
        pooled_matched_wt_edited_umis=("pooled_matched_wt_edited_umis", "min"),
        minimum_umis=("minimum_umis", "min"),
        pooled_fraction_edited=("pooled_fraction_edited", "mean"),
        pooled_matched_wt_fraction_edited=("pooled_matched_wt_fraction_edited", "mean"),
        **{column: (column, "first") for column in shared})
    averaged["position_id"] = measurement.position_id
    averaged["context"] = measurement.context
    averaged["editing_column"] = measurement.editing_column
    averaged["pooled_ratio_to_matched_wt"] = 2. ** averaged.log2_pooled_ratio_to_matched_wt
    averaged["percent_of_matched_wt"] = 100 * averaged.pooled_ratio_to_matched_wt
    require(np.isfinite(averaged.log2_pooled_ratio_to_matched_wt).eq(averaged.pooled_edited_umis.gt(0)).all(), "Only a sequence unedited at a site may have no mean")
    require(averaged.loc[averaged.n_mutations.eq(0), "pooled_ratio_to_matched_wt"].eq(1).all(), "Wild type is not the unit reference of a site average")
    return averaged[summary.columns]


def boxb_wildtype():
    """The wild-type boxB sequence, as RNA and as read DNA, with its annotated stem-loop regions."""
    regions = annotation("boxb_structure_regions").sort_values("start")
    rna = "".join(regions.rna_sequence)
    require(regions.iloc[0].start == 1 and regions.iloc[-1].end == len(rna) and regions.rna_sequence.str.fullmatch(r"[ACGU]+").all() and regions.rna_sequence.str.len().eq(regions.end - regions.start + 1).all() and np.array_equal(regions.start.to_numpy()[1:], regions.end.to_numpy()[:-1] + 1), "Invalid WT structure annotation")
    return regions, rna, rna.replace("U", "T")


def sequence_effects(samples, measurements):
    """One matched-control effect per reaction, measurement, and complete boxB sequence.

    Every construct is reconstructed into its full boxB sequence and matched to the wild-type
    construct carrying the same barcode, which fixes its randomized window and recorder placement.
    Edited and callable counts are summed over every retained representation, and the matched
    wild-type counts over exactly those barcode records, before one division. A sequence that two
    windows carry therefore receives one pooled estimate rather than one estimate per window, so the
    single-substitution and combinatorial figures cannot disagree about the same sequence.
    """
    positions = annotation("recorder_positions").set_index("position_id", verify_integrity=True)
    require("combine_sites" in measurements and measurements.combine_sites.isin(COMBINE_KINDS).all(), "Every measurement must say how its sites combine")
    # A mean is derived from sites measured independently, so bring in any of its own the caller did not
    # ask for, and drop them again once the mean is taken. Only the counted measurements read counts.
    derived = measurements.loc[measurements.combine_sites.eq(MEAN_PREFIX)]
    counted = measurements.loc[~measurements.combine_sites.eq(MEAN_PREFIX)]
    extra = [site for sites in derived.sites for site in sites.split(MEASUREMENT_SEPARATOR) if site not in set(counted.position_id)]
    if extra:
        counted = pd.concat([counted, *[site_measurement(site) for site in dict.fromkeys(extra)]], ignore_index=True)
    require(len(counted) > 0, "A measurement set must count at least one site")
    measurements, requested = counted, measurements.position_id.tolist()
    members = list(dict.fromkeys(MEASUREMENT_SEPARATOR.join(measurements.sites).split(MEASUREMENT_SEPARATOR)))
    site_columns = positions.loc[members, "editing_column"].tolist()
    _, _, wildtype = boxb_wildtype()
    barcodes = barcode_metadata()[["barcode", "oligo_name", "variable_subpos", "variable_length", "layout", *GEOMETRY_KEYS]]
    require(len(barcodes[GEOMETRY_KEYS].drop_duplicates()) == len(EXPECTED_DISTANCES), "Invalid boxB recorder placements")
    reactions = annotation("sample_info")
    require(len(samples) > 0 and len(samples.merge(reactions, on=SAMPLE_KEYS[:4])) == len(samples), "Reaction annotation mismatch")
    frames = []
    for sample in samples.to_dict("records"):
        data = counts_for(sample["sample_id"], site_columns).merge(barcodes, on="barcode", validate="many_to_one")
        bounds = data.variable_subpos.str.split("_", expand=True).astype(int)
        require(bounds[0].ge(1).all() and bounds[1].le(len(wildtype)).all() and data.variable_length.eq(bounds[1] - bounds[0] + 1).all() and data["insert"].str.len().eq(bounds[1] - bounds[0] + 1).all(), "Invalid boxB window coordinates")
        # The insert is sequenced on the antisense strand; the sense window is its reverse complement.
        data["hairpin"] = [wildtype[:start - 1] + insert.translate(COMPLEMENT)[::-1] + wildtype[end:] for start, end, insert in zip(bounds[0], bounds[1], data["insert"])]
        data["n_mutations"] = data.hairpin.map(lambda sequence: sum(a != b for a, b in zip(sequence, wildtype)))
        data["included_in_ratio"] = ~(data.layout.eq(EXCLUDED_LAYOUT) & data.variable_subpos.isin(EXCLUDED_WINDOWS))
        require((~data.included_in_ratio & data.n_mutations.eq(0)).sum() == len(EXCLUDED_WINDOWS), "Failed WT exclusion changed")
        require(not data.duplicated(["layout", "variable_subpos", "hairpin"]).any(), "Duplicate sequences within one window and placement")
        # Pooled sums divide, so a construct too thin to carry one callable molecule is held to
        # count consistency here rather than to being a positive denominator on its own.
        valid_counts(data, measurements.editing_column.tolist(), positive_coverage=False)
        # Each measurement carries its own callable denominator, so build the long frame per site.
        pieces = []
        for column in measurements.editing_column:
            piece = data.drop(columns=measurements.editing_column.tolist()).copy()
            piece["editing_column"] = column
            piece["edited_umis"] = data[column].to_numpy()
            piece["callable_umis"] = data[coverage_of(column)].to_numpy()
            pieces.append(piece)
        data = pd.concat(pieces, ignore_index=True).merge(measurements[["position_id", "context", "editing_column"]], on="editing_column", validate="many_to_one")
        for key, value in sample.items():
            data[key] = value
        frames.append(data)
    constructs = pd.concat(frames, ignore_index=True)
    wt_keys = ["sample_id", "position_id", "barcode"]
    controls = constructs.loc[constructs.n_mutations.eq(0), wt_keys + ["callable_umis", "edited_umis"]].rename(columns={"callable_umis": "matched_wt_umis", "edited_umis": "matched_wt_edited_umis"})
    constructs = constructs.merge(controls, on=wt_keys, how="left", validate="many_to_one")
    require(constructs[["matched_wt_umis", "matched_wt_edited_umis"]].notna().all().all() and constructs.loc[constructs.included_in_ratio, "matched_wt_edited_umis"].gt(0).all(), "Missing matched WT denominator")
    keys = SAMPLE_KEYS + ["position_id", "context", "editing_column", "hairpin", "n_mutations"]
    summary = constructs.loc[constructs.included_in_ratio].groupby(keys).agg(
        pooled_umis=("callable_umis", "sum"), pooled_edited_umis=("edited_umis", "sum"),
        pooled_matched_wt_umis=("matched_wt_umis", "sum"), pooled_matched_wt_edited_umis=("matched_wt_edited_umis", "sum"),
        n_ratio_geometries=("target_dist", "nunique"), n_ratio_windows=("variable_subpos", "nunique"),
        n_ratio_constructs=("insert", "size"), minimum_umis=("callable_umis", "min"),
        ratio_windows=("variable_subpos", lambda windows: ";".join(sorted(set(windows))))).reset_index()
    summary["pooled_fraction_edited"] = summary.pooled_edited_umis / summary.pooled_umis
    summary["pooled_matched_wt_fraction_edited"] = summary.pooled_matched_wt_edited_umis / summary.pooled_matched_wt_umis
    summary["pooled_ratio_to_matched_wt"] = summary.pooled_fraction_edited / summary.pooled_matched_wt_fraction_edited
    summary["percent_of_matched_wt"] = 100 * summary.pooled_ratio_to_matched_wt
    with np.errstate(divide="ignore"):
        # A covered sequence with no edited molecule keeps its −∞ effect rather than a pseudocount.
        summary["log2_pooled_ratio_to_matched_wt"] = np.log2(summary.pooled_ratio_to_matched_wt)
    require(summary.pooled_umis.gt(0).all() and summary.pooled_matched_wt_edited_umis.gt(0).all(), "Pooled denominators must be positive")
    require(np.isfinite(summary.log2_pooled_ratio_to_matched_wt).eq(summary.pooled_edited_umis.gt(0)).all(), "Only zero-edit sequences may be non-finite")
    require(summary.loc[summary.n_mutations.eq(0), "pooled_ratio_to_matched_wt"].eq(1).all(), "Wild type is not the unit reference")
    if len(derived):
        summary = pd.concat([summary, *[site_average(summary, row) for row in derived.itertuples()]], ignore_index=True)
        summary = summary.loc[summary.position_id.isin(requested)].reset_index(drop=True)
        require(set(summary.position_id) == set(requested), "A requested measurement was lost in deriving a site average")
    return summary, constructs


def multimutant_window_data(set_id):
    """Lay the shared pooled effect of each complete boxB sequence onto the six designed window grids."""
    selection = annotation("boxb_multimutant_heatmap_selection").sort_values("window_order")
    selection["excluded_layouts"] = selection.excluded_layouts.fillna("")
    require(len(selection) == EXPECTED_WT_WINDOWS and not selection.variable_subpos.duplicated().any() and not selection.window_order.duplicated().any(), "Invalid combinatorial window selections")
    require(selection.tada_type.nunique() == 1 and selection.position_id.nunique() == 1, "Combinatorial heat maps use one enzyme and one recorder site")
    reactions = reaction_set(set_id)
    sample = reactions.loc[reactions.tada_type.eq(selection.iloc[0].tada_type)]
    require(len(sample) == 1, "The combinatorial enzyme is not one reaction of its set")
    measurement = selected_measurement(selection.iloc[0].position_id)
    _, _, wildtype = boxb_wildtype()
    require(set(selection.variable_subpos) == set(barcode_metadata().variable_subpos), "Selections must cover every annotated boxB window")
    effects, constructs = sequence_effects(sample, measurement)
    effects = effects.loc[effects.position_id.eq(measurement.iloc[0].position_id)]
    frames = []
    for choice in selection.to_dict("records"):
        rows = [int(position) for position in str(choice["row_positions"]).split(";")]
        columns = [int(position) for position in str(choice["column_positions"]).split(";")]
        length = choice["window_end"] - choice["window_start"] + 1
        require(rows + columns == list(range(choice["window_start"], choice["window_end"] + 1)), "Axis positions must reconstruct the window in order")
        require(choice["variable_subpos"] == f"{choice['window_start']}_{choice['window_end']}", "Window coordinates disagree with the barcode annotation")
        # The window's own exclusion is the shared failed-library rule, not a separate selection.
        excluded = {layout for layout in choice["excluded_layouts"].split(";") if layout}
        require(excluded == ({EXCLUDED_LAYOUT} if choice["variable_subpos"] in EXCLUDED_WINDOWS else set()), "Window exclusions disagree with the failed-library exclusion")
        reference = wildtype[choice["window_start"] - 1:choice["window_end"]]
        designed = pd.DataFrame({"window_sequence": ["".join(bases) for bases in itertools.product(DNA_BASES, repeat=length)]})
        designed["hairpin"] = wildtype[:choice["window_start"] - 1] + designed.window_sequence + wildtype[choice["window_end"]:]
        designed["wildtype_sequence"] = reference
        designed["n_substitutions"] = [sum(a != b for a, b in zip(sequence, reference)) for sequence in designed.window_sequence]
        designed["is_wildtype"] = designed.n_substitutions.eq(0)
        designed["row_label"] = designed.window_sequence.str[:len(rows)].str.replace("T", "U")
        designed["column_label"] = designed.window_sequence.str[len(rows):].str.replace("T", "U")
        # Window depth is coverage metadata; the plotted estimate rests on every window carrying the sequence.
        window = constructs.loc[constructs.variable_subpos.eq(choice["variable_subpos"]) & constructs.included_in_ratio]
        coverage = window.groupby("hairpin").agg(window_umis=("callable_umis", "sum"), window_edited_umis=("edited_umis", "sum"), window_placements=("target_dist", "nunique"), window_constructs=("insert", "size"), window_minimum_umis=("callable_umis", "min")).reset_index()
        designed = designed.merge(coverage, on="hairpin", how="left", validate="one_to_one").merge(effects, on="hairpin", how="left", validate="one_to_one")
        require(designed.notna().all().all(), "Incomplete window design or missing pooled effect")
        for key, value in choice.items():
            designed[key] = value
        frames.append(designed)
    cells = pd.concat(frames, ignore_index=True).sort_values(["window_order", "row_label", "column_label"])
    require(len(cells) == sum(len(DNA_BASES) ** (row.window_end - row.window_start + 1) for row in selection.itertuples()), "Unexpected combinatorial design size")
    require(not cells.duplicated(["variable_subpos", "window_sequence"]).any() and not cells.duplicated(["variable_subpos", "row_label", "column_label"]).any(), "Window cells are not unique")
    require(cells.n_substitutions.eq(cells.n_mutations).all(), "Window substitution counts disagree with the full sequence")
    require(cells.is_wildtype.sum() == len(selection) and cells.loc[cells.is_wildtype, "log2_pooled_ratio_to_matched_wt"].eq(0).all(), "Wild-type cells are not the zero reference")
    require(np.isfinite(cells.log2_pooled_ratio_to_matched_wt).eq(cells.pooled_edited_umis.gt(0)).all(), "Only zero-edit observations may be non-finite")
    # Both representations of a sequence two windows carry must read the same pooled estimate.
    require(cells.groupby("hairpin").pooled_ratio_to_matched_wt.nunique().eq(1).all(), "A shared sequence carries more than one estimate")
    return cells, selection, sample.iloc[0], measurement.iloc[0]


def designed_boxb_sequences():
    """Every complete boxB sequence the combinatorial library designs, with its distance from wild type.

    `sequence_effects` returns every sequence a reaction observes, which includes a handful outside the
    design, so an analysis over the library must select the designed set rather than take what was seen.
    Those strays also differ between reactions, and selecting here keeps every reaction on one sequence set.
    """
    selection = annotation("boxb_multimutant_heatmap_selection").sort_values("window_order")
    require(len(selection) == EXPECTED_WT_WINDOWS, "Invalid combinatorial window selections")
    _, _, wildtype = boxb_wildtype()
    frames = []
    for choice in selection.to_dict("records"):
        length = choice["window_end"] - choice["window_start"] + 1
        designed = pd.DataFrame({"window_sequence": ["".join(bases) for bases in itertools.product(DNA_BASES, repeat=length)]})
        designed["hairpin"] = wildtype[:choice["window_start"] - 1] + designed.window_sequence + wildtype[choice["window_end"]:]
        frames.append(designed)
    cells = pd.concat(frames, ignore_index=True)
    cells["n_substitutions"] = [sum(a != b for a, b in zip(sequence, wildtype)) for sequence in cells.hairpin]
    sequences = cells[["hairpin", "n_substitutions"]].drop_duplicates()
    require(not sequences.hairpin.duplicated().any(), "A designed sequence carries two distances")
    require(sequences.n_substitutions.eq(0).sum() == 1, "Wild type is not one designed sequence")
    return sequences


def designed_variant_effects(set_id, position_ids):
    """Every designed boxB sequence's pooled effect at the given recorder sites, for one set's reactions.

    These are the same pooled estimates the combinatorial window maps draw, selected over the whole
    designed library rather than over its single substitutions, so a point placed here is the number a
    map cell shows. A sequence with no edited molecule keeps its −∞ effect rather than a pseudocount,
    so the figure that draws it can count and report it instead of losing it here.
    """
    metadata = reaction_set(set_id)
    require(len(position_ids) > 0 and len(set(position_ids)) == len(position_ids), "Designed-variant effects need distinct recorder sites")
    measurements = pd.concat([site_measurement(position_id) for position_id in position_ids], ignore_index=True)
    sequences = designed_boxb_sequences()
    effects, _ = sequence_effects(metadata, measurements)
    data = effects.merge(sequences, on="hairpin", validate="many_to_one")
    require(len(data) == len(metadata) * len(measurements) * len(sequences), "Every reaction and measurement must carry every designed sequence")
    require(data.n_substitutions.eq(data.n_mutations).all(), "Designed distance disagrees with the full sequence")
    require(np.isfinite(data.log2_pooled_ratio_to_matched_wt).eq(data.pooled_edited_umis.gt(0)).all(), "Only zero-edit sequences may be non-finite")
    require(data.loc[data.n_substitutions.eq(0), "pooled_ratio_to_matched_wt"].eq(1).all(), "Wild type is not the unit reference")
    return data


def mutation_distance_data(set_id):
    """Each designed boxB sequence's pooled effect for both enzymes, keyed by its distance from wild type.

    One window is randomized per construct, so a sequence's distance is the number of substitutions in
    that window. The effects are the same pooled estimates the heat maps draw, so a point in a violin is
    the number a map cell shows.
    """
    selection = annotation("boxb_multimutant_heatmap_selection")
    require(selection.position_id.nunique() == 1, "Mutation-distance violins use one recorder site")
    measurement = selected_measurement(selection.iloc[0].position_id)
    reactions = reaction_set(set_id)
    sequences = designed_boxb_sequences()
    effects, constructs = sequence_effects(reactions, measurement)
    effects = effects.loc[effects.position_id.eq(measurement.iloc[0].position_id)]
    distances = effects.merge(sequences, on="hairpin", validate="many_to_one")
    require(len(distances) == len(sequences) * len(reactions), "Every reaction must carry every designed sequence")
    require(distances.n_substitutions.eq(distances.n_mutations).all(), "Designed distance disagrees with the full sequence")
    # A sequence with no edited molecule has no logarithm. It stays in the frame, so the figure can
    # count and report it, rather than being dropped here where the count would be lost.
    require(np.isfinite(distances.log2_pooled_ratio_to_matched_wt).eq(distances.pooled_edited_umis.gt(0)).all(), "Only zero-edit sequences may be non-finite")
    require(distances.loc[distances.n_mutations.eq(0), "pooled_ratio_to_matched_wt"].eq(1).all(), "Wild type is not the unit reference")
    # The wild-type controls, one construct per window and placement, each against the pooled wild type
    # the mutant ratios are taken against. A mutant is matched to the wild type sharing its barcode, so
    # its own ratio is one by construction and carries no spread; these controls carry the spread instead,
    # and show how far apart the wild type reads across the library the mutants are measured in.
    controls = constructs.loc[constructs.n_mutations.eq(0) & constructs.included_in_ratio].copy()
    sites = measurement.iloc[0].sites.split(MEASUREMENT_SEPARATOR)
    # Each site carries its own copy of every control construct, so count them per site.
    require(controls.groupby(["sample_id", "editing_column"]).size().eq(EXPECTED_WT).all()
            and len(controls) == EXPECTED_WT * len(reactions) * len(sites), "Unexpected wild-type control set")
    require(controls.callable_umis.gt(0).all() and controls.edited_umis.gt(0).all(), "A wild-type control carries no callable or no edited molecule")
    controls["fraction_edited"] = controls.edited_umis / controls.callable_umis
    # Refer each site to its own pooled wild type before averaging the sites, exactly as the sequence
    # effects are built: one ratio per site, then the mean of the log₂ ratios.
    keys = ["sample_id", "editing_column"]
    controls["pooled_wt_fraction_edited"] = (controls.groupby(keys).edited_umis.transform("sum")
                                             / controls.groupby(keys).callable_umis.transform("sum"))
    controls["ratio_to_pooled_wt"] = controls.fraction_edited / controls.pooled_wt_fraction_edited
    controls["log2_ratio_to_pooled_wt"] = np.log2(controls.ratio_to_pooled_wt)
    require(np.isfinite(controls.log2_ratio_to_pooled_wt).all(), "A wild-type control has no ratio to its pooled mean")
    if len(sites) > 1:
        construct_keys = SAMPLE_KEYS + ["barcode", "insert"]
        controls = controls.groupby(construct_keys, as_index=False).agg(
            n_sites=("editing_column", "nunique"), log2_ratio_to_pooled_wt=("log2_ratio_to_pooled_wt", "mean"))
        require(controls.n_sites.eq(len(sites)).all() and len(controls) == EXPECTED_WT * len(reactions), "A wild-type control is missing one of its sites")
        controls["ratio_to_pooled_wt"] = 2. ** controls.log2_ratio_to_pooled_wt
    return distances, controls, measurement.iloc[0]


def single_mutant_data(set_id):
    """Select the wild type and the 57 single substitutions from the shared pooled sequence effects."""
    metadata = reaction_set(set_id)
    measurements = substitution_measurements(set_id)
    regions, wildtype_rna, wildtype = boxb_wildtype()
    grid = pd.DataFrame([dict(boxb_position=i + 1, nucleotide=base.replace("T", "U"), wildtype_base=reference.replace("T", "U"), is_wildtype=base == reference, hairpin=wildtype[:i] + base + wildtype[i + 1:]) for i, reference in enumerate(wildtype) for base in DNA_BASES])
    expected_sequences = set(grid.hairpin)
    effects, constructs = sequence_effects(metadata, measurements)
    singles = constructs.loc[constructs.n_mutations.le(1)]
    require(all(set(group.hairpin) == expected_sequences for _, group in singles.groupby(GEOMETRY_KEYS)), "Incomplete single-substitution library")
    summary = effects.loc[effects.n_mutations.le(1)].copy()
    # A substitution keeps four placements unless every window carrying it loses spacer5_0. Position
    # 10 lies in both windows 7-10 and 10-13, so it keeps four through its 10-13 representation.
    windows_per_hairpin = singles.loc[singles.n_mutations.eq(1)].groupby("hairpin").variable_subpos.agg(set)
    affected = {hairpin for hairpin, windows in windows_per_hairpin.items() if windows <= set(EXCLUDED_WINDOWS)}
    require(len(summary) == len(metadata) * len(measurements) * len(expected_sequences) and summary.n_ratio_geometries.eq(np.where(summary.hairpin.isin(affected), 3, 4)).all(), "Single-substitution support changed")
    tiles = summary.merge(grid, on="hairpin", validate="many_to_many").sort_values(["position_id", "panel_order", "boxb_position", "nucleotide"])
    tiles["wildtype_rna"] = wildtype_rna
    require(len(tiles) == len(metadata) * len(measurements) * len(wildtype) * len(DNA_BASES) and not tiles.duplicated(["sample_id", "position_id", "boxb_position", "nucleotide"]).any() and np.isfinite(tiles.log2_pooled_ratio_to_matched_wt).all() and tiles.loc[tiles.is_wildtype, "log2_pooled_ratio_to_matched_wt"].eq(0).all(), "Invalid single-substitution tiles")
    return tiles, constructs, regions
