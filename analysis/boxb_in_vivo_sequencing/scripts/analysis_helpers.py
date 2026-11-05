"""Shared CSV loading, construct matching, and recorder calculations for the in vivo boxB libraries.

The pooled effect contract is the in vitro one: every construct is reconstructed into its full boxB
sequence, matched to the wild-type construct carrying the same barcode, which fixes its randomized
window and recorder placement, and the edited and callable counts, and the matched wild-type counts, are
summed over exactly the retained representations before one division. A sequence two windows carry
therefore receives one pooled estimate, so no two figures can disagree about it. `coverage_of` is
imported from the in vitro helper rather than restated, so both analyses take each measurement's
denominator from one definition; nothing else is imported from it, and that module is left unchanged so
its figures cannot move.

No condition was sequenced twice, so summaries describe construct variation rather than biological
replication. The two recorder backgrounds are kept apart, because the G-depleted recorder changes
sense position 3 and so changes A2's own sequence context. No file is exported here.
"""

import itertools

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel

from in_vitro_reuse import helpers as in_vitro

ANALYSIS_DIRECTORY = "analysis/boxb_in_vivo_sequencing"
IN_VITRO_DIRECTORY = "analysis/boxb_in_vitro_sequencing"
ANNOTATION_DIRECTORY = f"{ANALYSIS_DIRECTORY}/annotations"
PRODUCTION_STAGE = "quality"
COUNT_DIRECTORY = f"{ANALYSIS_DIRECTORY}/data/recorder_counts/{PRODUCTION_STAGE}"
QC_DIRECTORY = f"{ANALYSIS_DIRECTORY}/data/qc"
# One definition of each measurement's denominator serves both analyses.
coverage_of = in_vitro.coverage_of

STANDARD_BACKGROUND = "no"
DNA_BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
READ_TO_RNA = dict(zip("ACGT", "UGCA"))
SAMPLE_KEYS = ["sample_id", "batch", "enzyme", "enzyme_label", "recruitment"]
GEOMETRY_KEYS = ["target_pos_to_boxb", "target_dist"]
EXPECTED_DISTANCES = [0, 10, 20, 30]
# The same thresholds the in vitro site figures star by, restated rather than imported, so the two
# analyses read the same way while this module keeps taking only `coverage_of` from that one.
SITE_STARS = [(.001, "***"), (.01, "**"), (.05, "*")]
EXPECTED_SINGLE_SUBSTITUTIONS = 57
EXPECTED_DESIGNED_SEQUENCES = 760


def require(condition, message):
    if not condition:
        raise ValueError(message)


def annotation(name):
    return pd.read_csv(f"{ANNOTATION_DIRECTORY}/{name}.csv", comment="#", dtype={"variable_subpos": str})


def samples():
    """Every library, with its experimental and sequencing metadata joined on sample_id."""
    experimental = annotation("sample_info")
    sequencing = annotation("sample_annotations")
    require(not experimental.sample_id.duplicated().any() and not sequencing.sample_id.duplicated().any(), "Duplicate library annotations")
    merged = experimental.merge(sequencing, on="sample_id", validate="one_to_one", suffixes=("", "_sequencing"))
    require(merged.batch.eq(merged.batch_sequencing).all(), "The two library annotations disagree on batch")
    require(merged[SAMPLE_KEYS].notna().all().all(), "Incomplete library metadata")
    return merged


def comparison_sets():
    """The libraries each figure family draws, grouped into the comparisons that may sit side by side."""
    sets = annotation("sample_comparisons")
    known = samples()
    require(not sets.sample_id.duplicated().any(), "A library belongs to two comparisons")
    require(set(sets.sample_id) <= set(known.sample_id), "A comparison names an unannotated library")
    require(sets.groupby("comparison_set").scale_group.nunique().eq(1).all(), "A comparison spans two scale groups")
    require(sets.groupby("comparison_set").panel_order.apply(lambda order: sorted(order) == list(range(len(order)))).all(),
            "Panel order must number a comparison's panels")
    merged = sets.merge(known[["sample_id", "batch"]], on="sample_id", validate="one_to_one")
    require(merged.groupby("comparison_set").batch.nunique().eq(1).all(), "A comparison must not mix batches")
    return sets


def comparison_samples(set_id):
    """One comparison's libraries, in panel order, with the metadata the effect calculation reads."""
    sets = comparison_sets()
    selected = sets.loc[sets.comparison_set.eq(set_id)]
    require(len(selected) > 0, f"Unannotated comparison {set_id}")
    merged = selected.merge(samples(), on="sample_id", validate="one_to_one")
    return merged.sort_values("panel_order")[["comparison_set", "scale_group", "role", "panel_order", *SAMPLE_KEYS]].reset_index(drop=True)


def comparison_set_ids(scale_group=None):
    sets = comparison_sets()
    if scale_group is not None:
        sets = sets.loc[sets.scale_group.eq(scale_group)]
    return sets.comparison_set.drop_duplicates().tolist()


def excluded_constructs():
    """The in vivo exclusion list, established from these libraries rather than copied from in vitro."""
    excluded = annotation("excluded_constructs")
    require(list(excluded.columns) == ["layout", "variable_subpos", "g_depleted", "reason", "source"], "Unexpected exclusion columns")
    if len(excluded):
        require(excluded.notna().all().all(), "An exclusion must name its reason and source")
        known = barcode_metadata(background=None)[["layout", "variable_subpos", "g_depleted"]]
        require(len(excluded.merge(known, on=["layout", "variable_subpos", "g_depleted"])) == len(excluded),
                "An exclusion names a construct the design does not have")
    return excluded


def barcode_metadata(background=STANDARD_BACKGROUND):
    """The boxB constructs, by default the standard recorder background alone."""
    data = annotation("barcode_annotations")
    require(data.variable_type.eq("boxb").all(), "The in vivo pool is entirely boxB constructs")
    if background is not None:
        data = data.loc[data.g_depleted.eq(background)].copy()
    data = data.drop(columns="barcode").rename(columns={"reverse_complement": "barcode"})
    require(not data.barcode.duplicated().any() and data.layout.notna().all(), "Invalid barcode annotations")
    return data


def recorder_positions(g_depleted=STANDARD_BACKGROUND):
    """The eight measured sites, with each site's sense context derived from its own background.

    The G-depleted recorder differs from the standard one at sense position 3, which is A2's 3′
    neighbour, so a context cannot be copied between backgrounds.
    """
    positions = annotation("recorder_positions").sort_values("read_index").reset_index(drop=True)
    backgrounds = annotation("recorder_backgrounds").set_index("g_depleted", verify_integrity=True)
    recorder = backgrounds.loc[g_depleted, "recorder_read_sequence"]
    sense = recorder.translate(COMPLEMENT)[::-1].replace("T", "U")
    require(positions.editing_column.tolist() == [f"pos_{index}_C" for index in positions.index],
            "Recorder positions must number the pos_<i> columns in read order")
    sense_positions = len(recorder) - positions.read_index
    require(positions.position_id.eq("A" + sense_positions.astype(str)).all(), "A site identifier disagrees with its read offset")
    require(sense_positions.map(lambda position: sense[position - 1]).eq("A").all(), "A measured recorder site is not an adenosine")
    positions["sense_position"] = sense_positions
    positions["derived_context"] = [sense[position - 2:position + 1] for position in sense_positions]
    if g_depleted == STANDARD_BACKGROUND:
        require(positions.derived_context.eq(positions.context).all(), "A derived context disagrees with the annotation")
    return positions


def site_measurement(position_id):
    positions = recorder_positions().set_index("position_id", verify_integrity=True)
    require(position_id in positions.index, f"Unannotated recorder site {position_id}")
    site = positions.loc[position_id]
    return pd.DataFrame([{"position_id": position_id, "sites": position_id, "label": position_id,
                          "context": site.derived_context, "editing_column": site.editing_column}])


def drawn_at(selections, set_id):
    """Keep the selections one comparison draws, from the comparisons each one declares."""
    declared = selections.comparison_sets.str.split(";")
    require(declared.explode().isin(comparison_sets().comparison_set.tolist()).all(), "Unannotated comparison in a selection")
    return selections.loc[declared.map(lambda sets: set_id in sets)].reset_index(drop=True)


def substitution_measurements(set_id=None):
    """The marginal and joint recorder measurements the substitution figures report."""
    measurements = annotation("boxb_substitution_measurements")
    positions = recorder_positions().set_index("position_id", verify_integrity=True)
    require(not measurements.position_id.duplicated().any() and measurements.notna().all().all(), "Invalid substitution measurements")
    for index, measurement in measurements.iterrows():
        members = measurement.sites.split(";")
        require(len(set(members)) == len(members) and set(members) <= set(positions.index), "Invalid joint recorder sites")
        require(len(members) > 1 or measurement.position_id == members[0], "Single-site identifier mismatch")
        measurements.loc[index, "context"] = " & ".join(positions.loc[members, "derived_context"])
        # A joint event is counted directly; its numerator is the state where every named site is edited.
        measurements.loc[index, "editing_column"] = (positions.loc[members[0], "editing_column"] if len(members) == 1
                                                     else f"joint_{measurement.position_id}_state_{'G' * len(members)}")
    if set_id is None:
        return measurements
    selected = drawn_at(measurements, set_id)
    require(len(selected) > 0, f"No substitution measurement is drawn at {set_id}")
    return selected


def coverage_column(column):
    """Each measurement's denominator. The joint states carry their own callable count."""
    if column.startswith("joint_"):
        return column if column.endswith("_callable") else f"{column.rsplit('_state_', 1)[0]}_callable"
    return coverage_of(column)


def counts_for(sample_id, columns, directory=COUNT_DIRECTORY):
    """One library's counts, with the callable denominator each requested measurement needs."""
    sites = [column for column in columns if column.startswith("pos_")]
    uncallable = [f"pos_{column.split('_')[1]}_N" for column in sites]
    joint = [coverage_column(column) for column in columns if column.startswith("joint_")]
    requested = list(dict.fromkeys(["sample_id", "barcode", "oligo_name", "insert", "umi_counts", "read_counts",
                                    "num_0_N", *columns, *uncallable, *joint]))
    data = pd.read_csv(f"{directory}/{sample_id}.csv.gz", usecols=requested)
    require(data.sample_id.eq(sample_id).all(), f"Unexpected library identifiers in {sample_id}")
    require(not data.duplicated(["sample_id", "barcode", "insert"]).any(), f"Duplicate constructs in {sample_id}")
    for column, missing in zip(sites, uncallable):
        data[f"coverage_{column}"] = data.umi_counts - data[missing]
    return data


def valid_counts(data, columns, positive_coverage=True):
    """Counts must be finite integers lying inside the molecules each measurement is defined on."""
    denominators = list(dict.fromkeys(coverage_column(column) for column in columns))
    matrix = data[["umi_counts", *columns, *denominators]].to_numpy()
    require(np.isfinite(matrix).all() and (matrix == np.floor(matrix)).all(), "Counts must be finite integers")
    require(data.umi_counts.gt(0).all(), "Construct UMI coverage must be positive")
    require(data[denominators].ge(0).all().all() and data[denominators].le(data.umi_counts, axis=0).all().all(),
            "Callable coverage must lie within the molecule count")
    if positive_coverage:
        require(data[denominators].gt(0).all().all(), "Callable coverage must be positive where it is a denominator")
    for column in columns:
        require((data[column].ge(0) & data[column].le(data[coverage_column(column)])).all(),
                f"Edited counts exceed callable coverage for {column}")


def boxb_wildtype():
    """The wild-type boxB sequence, as RNA and as read DNA, with its annotated stem-loop regions."""
    regions = annotation("boxb_structure_regions").sort_values("start")
    rna = "".join(regions.rna_sequence)
    require(regions.iloc[0].start == 1 and regions.iloc[-1].end == len(rna)
            and regions.rna_sequence.str.fullmatch(r"[ACGU]+").all()
            and regions.rna_sequence.str.len().eq(regions.end - regions.start + 1).all()
            and np.array_equal(regions.start.to_numpy()[1:], regions.end.to_numpy()[:-1] + 1), "Invalid WT structure annotation")
    return regions, rna, rna.replace("U", "T")


def designed_grid(background=STANDARD_BACKGROUND):
    """Every complete boxB sequence the library designs, with its distance from wild type.

    The grid is enumerated from the design before any count is read, so a designed sequence that no
    molecule reached is a measured absence rather than a row that silently never appears.
    """
    windows = barcode_metadata(background)[["variable_subpos", "window_start", "window_end"]].drop_duplicates()
    _, _, wildtype = boxb_wildtype()
    frames = []
    for window in windows.itertuples():
        length = window.window_end - window.window_start + 1
        designed = pd.DataFrame({"window_sequence": ["".join(bases) for bases in itertools.product(DNA_BASES, repeat=length)]})
        designed["variable_subpos"] = window.variable_subpos
        designed["window_start"], designed["window_end"] = window.window_start, window.window_end
        designed["hairpin"] = wildtype[:window.window_start - 1] + designed.window_sequence + wildtype[window.window_end:]
        frames.append(designed)
    cells = pd.concat(frames, ignore_index=True)
    cells["n_mutations"] = [sum(a != b for a, b in zip(sequence, wildtype)) for sequence in cells.hairpin]
    require(not cells.duplicated(["variable_subpos", "window_sequence"]).any(), "Duplicate designed window cells")
    return cells


def designed_boxb_sequences(background=STANDARD_BACKGROUND):
    """The distinct designed sequences, which is what an analysis over the library must select."""
    cells = designed_grid(background)
    sequences = cells[["hairpin", "n_mutations"]].drop_duplicates()
    require(not sequences.hairpin.duplicated().any(), "A designed sequence carries two distances")
    require(sequences.n_mutations.eq(0).sum() == 1, "Wild type is not one designed sequence")
    if background == STANDARD_BACKGROUND:
        require(len(sequences) == EXPECTED_DESIGNED_SEQUENCES, "Unexpected designed sequence count")
        require(sequences.n_mutations.eq(1).sum() == EXPECTED_SINGLE_SUBSTITUTIONS, "Unexpected single-substitution count")
    return sequences


def designed_construct_grid(background=None):
    """Every construct and genotype the library designs, before any count is read.

    One row per barcode and designed insert, which is the unit a read is binned into, so a genotype no
    molecule reached appears here as a measured absence rather than as a row that never appears. Both
    recorder backgrounds are included by default, since library-level coverage is a property of what was
    sequenced rather than of the subset a given figure analyses.
    """
    barcodes = barcode_metadata(background)
    frames = []
    for construct in barcodes.itertuples():
        inserts = ["".join(bases) for bases in itertools.product(DNA_BASES, repeat=int(construct.variable_length))]
        require(len(inserts) == int(construct.n_variants), f"Designed genotype count disagrees for {construct.oligo_name}")
        frames.append(pd.DataFrame({"barcode": construct.barcode, "oligo_name": construct.oligo_name,
                                    "variable_subpos": construct.variable_subpos, "layout": construct.layout,
                                    "g_depleted": construct.g_depleted, "insert": inserts}))
    grid = pd.concat(frames, ignore_index=True)
    require(not grid.duplicated(["barcode", "insert"]).any(), "Duplicate designed construct genotypes")
    return grid


def construct_measurements(sample_metadata, measurements, background=STANDARD_BACKGROUND, apply_exclusions=True):
    """One long row per library, construct, observed insert, and measurement, with its own denominator.

    `included_in_ratio` marks the constructs the pooled estimate sums over, so a figure can report its
    results with and without the flagged placements rather than losing them here.
    """
    members = list(dict.fromkeys(";".join(measurements.sites).split(";")))
    positions = recorder_positions(background).set_index("position_id", verify_integrity=True)
    site_columns = positions.loc[members, "editing_column"].tolist()
    requested = list(dict.fromkeys(site_columns + measurements.editing_column.tolist()))
    _, _, wildtype = boxb_wildtype()
    barcodes = barcode_metadata(background)[["barcode", "oligo_name", "variable_subpos", "variable_length", "window_start",
                                            "window_end", "layout", "g_depleted", *GEOMETRY_KEYS]]
    require(len(barcodes[GEOMETRY_KEYS].drop_duplicates()) == len(EXPECTED_DISTANCES), "Invalid boxB recorder placements")
    excluded = excluded_constructs() if apply_exclusions else excluded_constructs().iloc[:0]
    known = samples()
    require(len(sample_metadata) > 0 and len(sample_metadata.merge(known, on=SAMPLE_KEYS)) == len(sample_metadata),
            "Library annotation mismatch")

    frames = []
    for sample in sample_metadata.to_dict("records"):
        data = counts_for(sample["sample_id"], requested).merge(barcodes, on=["barcode", "oligo_name"], validate="many_to_one")
        require(data.variable_length.eq(data.window_end - data.window_start + 1).all()
                and data["insert"].str.len().eq(data.variable_length).all(), "Invalid boxB window coordinates")
        # The insert is sequenced on the antisense strand; the sense window is its reverse complement.
        data["hairpin"] = [wildtype[:start - 1] + insert.translate(COMPLEMENT)[::-1] + wildtype[end:]
                           for start, end, insert in zip(data.window_start, data.window_end, data["insert"])]
        data["n_mutations"] = data.hairpin.map(lambda sequence: sum(a != b for a, b in zip(sequence, wildtype)))
        flagged = data.merge(excluded[["layout", "variable_subpos", "g_depleted"]].assign(flagged=True),
                             on=["layout", "variable_subpos", "g_depleted"], how="left", validate="many_to_one")
        data["included_in_ratio"] = flagged.flagged.isna().to_numpy()
        require(not data.duplicated(["layout", "variable_subpos", "hairpin"]).any(),
                "Duplicate sequences within one window and placement")
        for event in measurements.loc[measurements.sites.str.contains(";")].itertuples():
            marginals = data[positions.loc[event.sites.split(";"), "editing_column"]]
            # A joint event needs every one of its sites callable, while a marginal keeps molecules that
            # are uncallable elsewhere, so only the upper bound follows from the marginals.
            require(data[event.editing_column].le(marginals.min(axis=1)).all(), "Joint edits exceed their marginals")
        # Pooled sums divide, so a construct too thin to carry one callable molecule is held to count
        # consistency here rather than to being a positive denominator on its own.
        valid_counts(data, measurements.editing_column.tolist(), positive_coverage=False)
        pieces = []
        for column in measurements.editing_column:
            piece = data.drop(columns=requested).copy()
            piece["editing_column"] = column
            piece["edited_umis"] = data[column].to_numpy()
            piece["callable_umis"] = data[coverage_column(column)].to_numpy()
            pieces.append(piece)
        long = pd.concat(pieces, ignore_index=True).merge(
            measurements[["position_id", "context", "label", "editing_column"]], on="editing_column", validate="many_to_one")
        for key, value in sample.items():
            long[key] = value
        frames.append(long)
    return pd.concat(frames, ignore_index=True)


def sequence_effects(sample_metadata, measurements, background=STANDARD_BACKGROUND, apply_exclusions=True):
    """One matched-control effect per library, measurement, and complete boxB sequence."""
    constructs = construct_measurements(sample_metadata, measurements, background, apply_exclusions)
    wt_keys = ["sample_id", "position_id", "barcode"]
    controls = constructs.loc[constructs.n_mutations.eq(0), wt_keys + ["callable_umis", "edited_umis"]].rename(
        columns={"callable_umis": "matched_wt_umis", "edited_umis": "matched_wt_edited_umis"})
    constructs = constructs.merge(controls, on=wt_keys, how="left", validate="many_to_one")
    require(constructs[["matched_wt_umis", "matched_wt_edited_umis"]].notna().all().all(), "Missing matched WT denominator")
    keys = list(sample_metadata.columns) + ["position_id", "context", "label", "editing_column", "hairpin", "n_mutations"]
    retained = constructs.loc[constructs.included_in_ratio]
    summary = retained.groupby(keys, sort=False).agg(
        pooled_umis=("callable_umis", "sum"), pooled_edited_umis=("edited_umis", "sum"),
        pooled_matched_wt_umis=("matched_wt_umis", "sum"), pooled_matched_wt_edited_umis=("matched_wt_edited_umis", "sum"),
        n_ratio_geometries=("target_dist", "nunique"), n_ratio_windows=("variable_subpos", "nunique"),
        n_ratio_constructs=("insert", "size"), minimum_umis=("callable_umis", "min"),
        ratio_windows=("variable_subpos", lambda windows: ";".join(sorted(set(windows))))).reset_index()
    summary["pooled_fraction_edited"] = summary.pooled_edited_umis / summary.pooled_umis
    summary["pooled_matched_wt_fraction_edited"] = summary.pooled_matched_wt_edited_umis / summary.pooled_matched_wt_umis
    summary["pooled_ratio_to_matched_wt"] = summary.pooled_fraction_edited / summary.pooled_matched_wt_fraction_edited
    summary["percent_of_matched_wt"] = 100 * summary.pooled_ratio_to_matched_wt
    with np.errstate(divide="ignore", invalid="ignore"):
        # A covered sequence with no edited molecule keeps its −∞ effect rather than a pseudocount, and a
        # sequence whose matched wild type carries no edited molecule has no ratio at all.
        summary["log2_pooled_ratio_to_matched_wt"] = np.log2(summary.pooled_ratio_to_matched_wt)
    require(summary.pooled_umis.gt(0).all(), "Pooled callable coverage must be positive")
    require(summary.loc[summary.n_mutations.eq(0) & summary.pooled_matched_wt_edited_umis.gt(0),
                        "pooled_ratio_to_matched_wt"].eq(1).all(), "Wild type is not the unit reference")
    # Both representations of a sequence two windows carry must read one pooled estimate.
    shared = summary.groupby(["sample_id", "position_id", "hairpin"]).pooled_ratio_to_matched_wt.nunique()
    require(shared.le(1).all(), "A shared sequence carries more than one estimate")
    return summary, constructs


def matched_stem_pairs(sample_metadata, measurements, background=STANDARD_BACKGROUND):
    """The predefined stem-disrupting mutants and their matched wild types, as the in vitro set defines them.

    The expected construct counts are derived from the design and the exclusion list rather than fixed, so
    excluding a placement changes what is required instead of breaking the check. Every retained barcode
    carries one wild type, and each listed mutant appears once per retained placement of its own window.
    """
    listed = annotation("boxb_wildtype_mutant_list")[["variable_subpos", "mutant_bases", "insert_type", "insert"]]
    require(not listed.duplicated().any(), "Duplicate construct annotations")
    require(listed.loc[listed.mutant_bases.isna(), "insert_type"].eq("wt").all(), "Unpaired constructs must be wild type")
    constructs = construct_measurements(sample_metadata, measurements, background)
    retained = constructs.loc[constructs.included_in_ratio].drop_duplicates(["layout", "variable_subpos"])
    expected_wt = len(retained)
    mutants = listed.loc[listed.insert_type.eq("mut")].merge(
        retained[["layout", "variable_subpos"]], on="variable_subpos", validate="many_to_many")
    expected_pairs = len(mutants)
    expected_windows = retained.variable_subpos.nunique()
    data = constructs.loc[constructs.included_in_ratio].merge(listed, on=["variable_subpos", "insert"], validate="many_to_many")
    data["window_id"] = data.layout + "__" + data.variable_subpos
    data["pair_id"] = data.window_id + "__" + data.mutant_bases
    data["fraction_edited"] = data.edited_umis / data.callable_umis
    for (sample_id, position_id), group in data.groupby(["sample_id", "position_id"]):
        unique = group.drop_duplicates(["barcode", "insert"])
        counts = unique.insert_type.value_counts().to_dict()
        require(counts.get("mut") == expected_pairs and counts.get("wt") == expected_wt,
                f"Unexpected WT/MUT construct counts in {sample_id} at {position_id}: {counts}, "
                f"expected {expected_pairs} mutants and {expected_wt} wild types")
        require(unique.loc[unique.insert_type.eq("wt"), "variable_subpos"].nunique() == expected_windows,
                f"Unexpected wild-type windows in {sample_id} at {position_id}")
    return data


def edit_count_fractions(sample_metadata, background=STANDARD_BACKGROUND):
    """One- and two-or-more-edit fractions per construct and genotype, over the fully callable molecules.

    An edit-count bin counts edits across all eight recorder sites, so its denominator is the molecules
    whose whole recorder is callable, `num_0_N`. The `2+ edits` label includes the eight-edit class; the
    bins are summed rather than taken as one minus the zero-edit fraction, so nothing silently omits it.
    """
    bins = [f"num_{number}_C" for number in range(1, 9)]
    barcodes = barcode_metadata(background)[["barcode", "oligo_name", "variable_subpos", "layout", *GEOMETRY_KEYS]]
    frames = []
    for sample in sample_metadata.to_dict("records"):
        data = counts_for(sample["sample_id"], bins).merge(barcodes, on=["barcode", "oligo_name"], validate="many_to_one")
        require(data[bins].sum(axis=1).le(data.num_0_N).all(), "Edit-count bins exceed callable coverage")
        data["frac_1edit"] = np.where(data.num_0_N.gt(0), data.num_1_C / data.num_0_N, np.nan)
        data["frac_2edit"] = np.where(data.num_0_N.gt(0), data[bins[1:]].sum(axis=1) / data.num_0_N, np.nan)
        data["edited_1"] = data.num_1_C
        data["edited_2plus"] = data[bins[1:]].sum(axis=1)
        for key, value in sample.items():
            data[key] = value
        frames.append(data)
    return pd.concat(frames, ignore_index=True)


def boxb_self_editing(sample_metadata, background=STANDARD_BACKGROUND):
    """Editing of boxB's own constant adenosines, per library, construct, and sense position.

    A true A→G edit inside boxB is indistinguishable from a designed A→G substitution when the genotype
    is read from the same molecule, so this is the measurement that says how far genotype-level claims
    can be trusted. Only the positions a construct holds constant are measurable in that construct.
    """
    _, wildtype_rna, _ = boxb_wildtype()
    barcodes = barcode_metadata(background)[["barcode", "oligo_name", "variable_subpos", "window_start", "window_end",
                                             "layout", "g_depleted", *GEOMETRY_KEYS]]
    adenosines = [position for position, base in enumerate(wildtype_rna, start=1) if base == "A"]
    require(len(adenosines) > 0, "The wild-type boxB carries no adenosine")
    frames = []
    for sample in sample_metadata.to_dict("records"):
        columns = [f"boxb_{position}_{base}" for position in adenosines for base in "ACGTN"]
        data = counts_for(sample["sample_id"], []).merge(
            pd.read_csv(f"{COUNT_DIRECTORY}/{sample['sample_id']}.csv.gz", usecols=["barcode", "insert", *columns]),
            on=["barcode", "insert"], validate="one_to_one").merge(barcodes, on=["barcode", "oligo_name"], validate="many_to_one")
        for position in adenosines:
            constant = ~data.window_start.le(position) | ~data.window_end.ge(position)
            observed = data[[f"boxb_{position}_{base}" for base in "ACGTN"]].sum(axis=1)
            require(observed.loc[constant].eq(data.umi_counts.loc[constant]).all(),
                    f"boxB position {position} base counts do not sum to the molecule count where it is constant")
            require(observed.loc[~constant].eq(0).all(), f"boxB position {position} was measured inside a randomized window")
            piece = data.loc[constant, ["barcode", "oligo_name", "insert", "umi_counts", "variable_subpos", "layout",
                                        "g_depleted", *GEOMETRY_KEYS]].copy()
            piece["sense_position"] = position
            piece["callable_umis"] = data.loc[constant, [f"boxb_{position}_{base}" for base in "ACGT"]].sum(axis=1)
            piece["edited_umis"] = data.loc[constant, f"boxb_{position}_G"]
            piece["other_base_umis"] = data.loc[constant, [f"boxb_{position}_{base}" for base in "CT"]].sum(axis=1)
            for key, value in sample.items():
                piece[key] = value
            frames.append(piece)
    editing = pd.concat(frames, ignore_index=True)
    editing["fraction_edited"] = np.where(editing.callable_umis.gt(0), editing.edited_umis / editing.callable_umis, np.nan)
    regions = annotation("boxb_structure_regions")
    editing = editing.merge(regions[["region", "label", "start", "end"]].rename(columns={"label": "region_label"}).assign(key=1),
                            how="cross").query("start <= sense_position <= end").drop(columns=["key", "start", "end"])
    return editing


