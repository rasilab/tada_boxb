"""All configured recorder comparisons, with distinct-construct means and paired tests."""

from analysis_helpers import comparison_data, annotation, require
from plot_style import plt, palettes, MUTANT_COLORS, FONT_SIZE, save_pdf, place_shared_ylabel
from matplotlib.text import Annotation
from math import isfinite

# Side margins are held in inches rather than as a fraction of the page, so a figure narrowed to its
# bars keeps the room its tick numbers, y label, and right-hand strip label actually need.
SIDE_MARGINS_IN = {"edit_count": (.54, .27), "site": (.592, .222)}
# A figure that stacks rows of panels sets its bands in inches as well, so the panels take the height
# the annotation says they need rather than a share of the page. The band above each row holds that
# row's condition names and the air that separates one comparison from the next.

STACKED_TOP_IN = .18
STACKED_TITLE_IN = .32
STACKED_TITLE_PAD_PT = 1.5
STACKED_BOTTOM_IN = .40
# A stacked site figure names only its WT/MUT pairs under the bottom row, not edit groups as well.
STACKED_SITE_BOTTOM_IN = .22
STACKED_WSPACE = .32
# The edit names sit a fixed distance below a stacked panel, not a share of its height, so shortening
# the panels does not walk them into the pair names above.
STACKED_EDIT_LABEL_PT = 18
# Panels that share a y axis need no room between them for a second set of tick numbers.
SHARED_WSPACE = .12
LABEL_CLEARANCE_PT = .5
SHARED_GROUP_MARGIN = .5
# A shared axis holds its outermost bar inside the panel by this much, measured in inches like the
# bar widths and spacings themselves.
BAR_CLEARANCE_IN = .02
FOLD_GAP_PT = 1.5
FOLD_DECIMAL_DIGITS = 1
ANNOTATION_TOP_PADDING_PT = 1
EDIT_LABEL_Y = -.32
COMPONENT_KEY_HEIGHT = .65
COMPONENT_SYMBOLS = {False: "−", True: "+"}
VARIANT_ORDER = ["wt", "mut"]
EDIT_LABELS = {"frac_1edit": "1 edit", "frac_2edit": "2+ edits"}

definitions, prepared = comparison_data()
labels, light_colors, dark_colors = palettes()
samples = annotation("sample_info")
components = annotation("recorder_reaction_components")
figure_names = list(dict.fromkeys(row["figure_group"] for row in definitions))


def fitted_top(figure, ax, renderer):
    """The upper limit that clears every annotation of one panel, with no fixed headroom."""
    height_pt = ax.get_window_extent(renderer).height * 72 / figure.dpi
    upper = ax.get_ylim()[1] / (1 - ANNOTATION_TOP_PADDING_PT / height_pt)
    for text in ax.texts:
        if isinstance(text, Annotation):
            text_height_pt = text.get_window_extent(renderer).height * 72 / figure.dpi
            reserve_pt = text.xyann[1] + text_height_pt + ANNOTATION_TOP_PADDING_PT
            upper = max(upper, text.xy[1] / (1 - reserve_pt / height_pt))
    return upper


for figure_name in figure_names:
    rows = [row for row in definitions if row["figure_group"] == figure_name]
    site_mode = rows[0]["measurement_mode"] == "site"
    # A shared-axis row draws every condition of one comparison on one axis, whether it measures one
    # adenosine or the edit-count groups; the panel layout gives each condition an axis of its own.
    single_axis = rows[0]["layout"] in ["shared_axis", "components"]
    component_key = rows[0]["layout"] == "components"
    enzyme_order = rows[0]["enzymes"].split(";")
    # One y axis per figure draws the panels against each other rather than each against itself; only
    # the leftmost column then needs tick numbers, and a stacked figure needs the label only once.
    shared_y = bool(rows[0]["shared_y"])
    nrows = len(rows) if len(rows) > 1 else (1 if single_axis or rows[0]["column"] == "tada_type" else len(enzyme_order))
    ncols = 1 if single_axis else len(rows[0]["levels"].split(";"))
    # A stacked figure names the measurement once beside the panels, whether or not they share a scale:
    # the label names the quantity, which every row carries, rather than the scale each row is drawn on.
    # Rows of panels, whether they carry edit counts or one adenosine on a shared axis, set their
    # bands in inches and name the measurement once beside them.
    stacked = len(rows) > 1 and (not site_mode or single_axis)
    shared_label = stacked
    width = rows[0]["width"]
    bottom_band = STACKED_SITE_BOTTOM_IN if site_mode else STACKED_BOTTOM_IN
    height = STACKED_TOP_IN + sum(STACKED_TITLE_IN + row["height"] for row in rows) + bottom_band if stacked \
        else sum(row["height"] for row in rows) + (.30 if len(rows) > 1 else 0) + (COMPONENT_KEY_HEIGHT if component_key else 0)
    left_margin, right_margin = SIDE_MARGINS_IN["site" if site_mode else "edit_count"]
    vertical = dict(top=1 - (STACKED_TOP_IN + STACKED_TITLE_IN) / height, bottom=bottom_band / height, hspace=STACKED_TITLE_IN / rows[0]["height"], wspace=STACKED_WSPACE) if stacked \
        else dict(top=.87 if len(rows) > 1 else .84, bottom=.16 if nrows > 1 else (.36 if component_key else .29), hspace=.72, wspace=SHARED_WSPACE if shared_y else .38)
    gap_scale = rows[0]["panel_gap_scale"]
    if gap_scale != 1:
        # Remove the released gaps from the canvas, preserving the physical panel and margin sizes.
        panel_width = (width - left_margin - right_margin) / (ncols + (ncols - 1) * vertical["wspace"])
        panel_height = height * (vertical["top"] - vertical["bottom"]) / (nrows + (nrows - 1) * vertical["hspace"])
        top_band, bottom_band = height * (1 - vertical["top"]), height * vertical["bottom"]
        vertical["wspace"] *= gap_scale
        vertical["hspace"] *= gap_scale
        width = left_margin + right_margin + panel_width * (ncols + (ncols - 1) * vertical["wspace"])
        height = top_band + bottom_band + panel_height * (nrows + (nrows - 1) * vertical["hspace"])
        vertical.update(top=1 - top_band / height, bottom=bottom_band / height)
    figure, axes = plt.subplots(nrows, ncols, figsize=(width, height), squeeze=False)
    figure.subplots_adjust(left=left_margin / width, right=1 - right_margin / width, **vertical)
    if len(rows) > 1:
        figure.suptitle(labels[enzyme_order[0]], y=1 - .04 / height if stacked else .97)
    for row_index in range(nrows):
        definition = rows[row_index] if len(rows) > 1 else rows[0]
        summary, statistics = prepared[definition["name"]]
        levels = definition["levels"].split(";")
        condition_labels = definition["labels"].split(";")
        edits = ["frac_site"] if site_mode else ["frac_1edit", "frac_2edit"]
        # Every row of a stacked edit-count figure repeats the same pairs and the same edit groups, so
        # only the bottom row names them.
        bottom_row = row_index == nrows - 1
        name_pairs = bottom_row if stacked or nrows > 1 else True
        enzyme = None if single_axis or definition["column"] == "tada_type" else enzyme_order[row_index] if len(rows) == 1 else enzyme_order[0]
        row_data = summary if enzyme is None else summary.loc[summary.tada_type.eq(enzyme)]
        row_stats = statistics if enzyme is None else statistics.loc[statistics.tada_type.eq(enzyme)]
        for column_index, ax in enumerate(axes[row_index]):
            level = levels[column_index]
            data = row_data if single_axis else row_data.loc[row_data[definition["column"]].eq(level)]
            tests = row_stats if single_axis else row_stats.loc[row_stats[definition["column"]].eq(level)]
            panel_top = (data["mean"] + data.se).max()
            ticks, tick_labels = [], []
            group_levels = levels if single_axis else [level]
            group_spacing = definition["group_spacing"]
            panel_width = figure.get_figwidth() * ax.get_position().width
            if single_axis:
                # The margin has to clear the outermost bar rather than hold a fixed number of axis
                # units: a condition carrying two edit groups reaches further from its centre than one
                # carrying a single pair, and a fixed margin cut the first and last bars in half. The
                # margin and the unit scale each set the other, so they are solved together and the
                # clearance stays physical. A comparison whose bars already fit keeps the former margin.
                span = (len(levels) - 1) * group_spacing
                reach = (definition["edit_spacing"] if len(edits) > 1 else 0) / 2 \
                    + definition["pair_spacing"] / 2 + definition["bar_width"] / 2 + BAR_CLEARANCE_IN
                require(2 * reach < panel_width, f"{figure_name} is narrower than the bars of one condition")
                margin = max(SHARED_GROUP_MARGIN, (span / (1 - 2 * reach / panel_width) - span) / 2)
                x_limits = (-margin, span + margin)
            else:
                x_limits = (-.65, (len(edits) - 1) * definition["edit_gap"] + .65)
            units_per_inch = (x_limits[1] - x_limits[0]) / panel_width
            bar_width = definition["bar_width"] * units_per_inch
            pair_offset = definition["pair_spacing"] * units_per_inch / 2 if single_axis else definition["pair_offset"]
            for group_index, group_level in enumerate(group_levels):
                group = data.loc[data[definition["column"]].eq(group_level)]
                group_tests = tests.loc[tests[definition["column"]].eq(group_level)]
                # On a shared axis the edit groups of one condition sit a declared distance apart in
                # inches, measured from the centre of their condition; a site row draws one group and
                # keeps that centre.
                edit_spacing = definition["edit_spacing"] * units_per_inch if len(edits) > 1 else 0
                require(not single_axis or isfinite(edit_spacing), f"{figure_name} needs a declared edit spacing")
                for edit_index, edit in enumerate(edits):
                    offset = (edit_index - (len(edits) - 1) / 2) * edit_spacing
                    center = group_index * group_spacing + offset if single_axis else edit_index * definition["edit_gap"]
                    points = group.loc[group.edit_type.eq(edit)].set_index("insert_type")
                    test = group_tests.loc[group_tests.edit_type.eq(edit)].iloc[0]
                    for variant_index, variant in enumerate(VARIANT_ORDER):
                        point = points.loc[variant]
                        x = center + (variant_index * 2 - 1) * pair_offset
                        palette = light_colors if edit == "frac_1edit" else dark_colors
                        color = palette[point.tada_type] if variant == "wt" else MUTANT_COLORS["light" if edit == "frac_1edit" else "dark"]
                        ax.bar(x, point["mean"], width=bar_width, color=color, edgecolor="black", linewidth=.2)
                        ax.errorbar(x, point["mean"], yerr=point.se, color="black", fmt="none", linewidth=.35, capsize=1.5)
                        ticks.append(x)
                        tick_labels.append(variant.upper())
                    top = (points["mean"] + points.se).max()
                    hidden = group_level in definition["hide_annotations"].split(";")
                    if not hidden:
                        # The offset is measured to the text box, whose bottom sits a line's descent
                        # below the asterisks themselves, so a negative offset is what puts the stars
                        # against the bar; the fold label keeps its place above them either way.
                        star_offset_pt = definition["star_offset"]
                        fold_offset_pt = star_offset_pt
                        if test.significance != "ns":
                            ax.annotate(test.significance, (center, top), xytext=(0, star_offset_pt), textcoords="offset points", ha="center", va="bottom")
                            fold_offset_pt += FONT_SIZE + FOLD_GAP_PT
                        ax.annotate(f"{test.fold_change:.{FOLD_DECIMAL_DIGITS}f}x", (center, top), xytext=(0, fold_offset_pt), textcoords="offset points", ha="center", va="bottom")
                    if not site_mode and name_pairs:
                        edit_label_y = -STACKED_EDIT_LABEL_PT / 72 / definition["height"] if stacked else EDIT_LABEL_Y
                        ax.text(center, edit_label_y, EDIT_LABELS[edit], ha="center", transform=ax.get_xaxis_transform())
                if single_axis and not component_key:
                    ax.text(group_index * group_spacing, 1.05, condition_labels[group_index], ha="center", va="bottom", transform=ax.get_xaxis_transform())
            ax.set_ylim(0, panel_top)
            ax.set_xlim(x_limits)
            ax.set_xticks(ticks, tick_labels if name_pairs else [])
            ax.tick_params(axis="x", length=0, labelsize=definition["pair_label_size"])
            # Panels that share one scale carry one drawn axis, at the left; the rest keep neither its
            # numbers nor its spine, which is what lets them sit close together.
            if shared_y and column_index > 0:
                ax.tick_params(axis="y", left=False, labelleft=False)
                ax.spines["left"].set_visible(False)
            if not single_axis and (row_index == 0 or len(rows) > 1):
                ax.set_title(condition_labels[column_index], pad=STACKED_TITLE_PAD_PT if stacked else 8)
            if column_index == 0 and not shared_label:
                ax.set_ylabel(f"% Edited Adenosine ({definition['position_id']})" if site_mode else "% Edited RNA")
            if column_index == ncols - 1 and (stacked or not single_axis):
                strip = definition["row_label"] if len(rows) > 1 else labels[enzyme] if enzyme is not None else ""
                ax.text(1.04, .5, strip, rotation=270, va="center", ha="left", transform=ax.transAxes)
            if component_key:
                ax.set_xticklabels([])
                selected = samples.loc[samples.short_key.isin(definition["short_keys"].split(";")) & samples.condition.isin(definition["conditions"].split(";")) & samples.tada_conc.isin(definition["concentrations"].split(";"))].set_index(definition["column"])
                key_rows = [(definition["variant_row_label"], None)] + [(component.component, component) for component in components.itertuples()]
                for key_index, (label, component) in enumerate(key_rows):
                    y = -.12 - key_index * .17
                    ax.text(x_limits[0] - .07, y, label, ha="right", va="center", transform=ax.get_xaxis_transform())
                    for group_index, group_level in enumerate(levels):
                        for variant_index, variant in enumerate(VARIANT_ORDER):
                            text = variant.upper() if component is None else COMPONENT_SYMBOLS[str(selected.loc[group_level, component.annotation_column]) != str(component.absent_value)]
                            ax.text(group_index * group_spacing + (variant_index * 2 - 1) * pair_offset, y, text, ha="center", va="center", transform=ax.get_xaxis_transform())
    # Fit every panel to its bars and actual annotation heights, without fixed headroom. Panels that
    # share a y axis all take the deepest fit, so one height means one measurement across the figure.
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    tops = [fitted_top(figure, ax, renderer) for ax in axes.flat]
    for ax, top in zip(axes.flat, tops):
        ax.set_ylim(0, max(tops) if shared_y else top)
    # A pair tightened past the width of its own names runs WT into MUT, so measure the drawn names
    # and fail here rather than ship a figure whose labels touch.
    for ax in axes.flat:
        drawn = [label.get_window_extent(renderer) for label in ax.get_xticklabels() if label.get_text()]
        for left, right in zip(drawn, drawn[1:]):
            require(right.x0 - left.x1 >= LABEL_CLEARANCE_PT * figure.dpi / 72, f"Bar names overlap in {figure_name}")
    if shared_label:
        # Panels on their own scales carry tick numbers of their own widths, so the label is placed
        # against whichever left-hand panel reaches furthest out.
        widest = min(axes[:, 0], key=lambda ax: min(number.get_window_extent(renderer).x0 for number in ax.get_yticklabels() if number.get_text()))
        place_shared_ylabel(figure, widest, f"% Edited Adenosine ({rows[0]['position_id']})" if site_mode else "% Edited RNA")
    save_pdf(figure, figure_name, png=bool(rows[0]["png_companion"]))
