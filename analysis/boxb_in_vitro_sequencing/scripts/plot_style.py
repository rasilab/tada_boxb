"""Shared manuscript styling and vector PDF output."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from analysis_helpers import annotation, enzyme_metadata, require

FONT = "Helvetica"
FONT_SIZE = 7
LABEL_GAP = .09
AXIS_COLOR = "#999999"
AXIS_WIDTH = 0.25
MUTANT_COLORS = {"light": "#cccccc", "dark": "#888888"}
DISTANCE_SHADES = ["#555555", "#888888", "#aaaaaa", "#cccccc"]
ENZYME_DISTANCE_SHADES = [["#76266B", "#AA4499", "#CD83BB", "#E2B8D6"], ["#064D23", "#117733", "#68AF84", "#ADD5BC"]]
NUCLEOTIDE_COLORS = {"A": "#117733", "C": "#0072B2", "G": "#E69F00", "U": "#C03030"}
FIGURE_DIRECTORY = "analysis/boxb_in_vitro_sequencing/figures"
OUTPUTS = annotation("figure_outputs")
PDF_METADATA = {"Creator": "BoxB Python analysis", "CreationDate": None, "ModDate": None}
PNG_DPI = 96
COMPACT_PADDING = .04
CONCENTRATION_CORRELATION_STYLE = {
    "size": (3.1, 2.85), "point_size": 5, "marker": "o", "alpha": .65, "wt_size": 24,
    "density_points": 256, "density_width": .16, "density_gap": .04,
    "density_alpha": .5, "reference_color": "#999999", "reference_width": .4,
    "lowest_tick": -4, "tick_step": 2, "axis_padding": .12, "ylabel_gap": .06,
    "stat_xy": (.04, .97), "row_label_x": 1.4,
}

font_manager.findfont(FONT, fallback_to_default=False)
plt.rcParams.update({"font.family": FONT, "font.size": FONT_SIZE, "figure.titlesize": FONT_SIZE,
                    "figure.labelsize": FONT_SIZE, "axes.titlesize": FONT_SIZE,
                    "axes.labelsize": FONT_SIZE, "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
                    "legend.fontsize": FONT_SIZE, "legend.title_fontsize": FONT_SIZE,
                    "axes.edgecolor": AXIS_COLOR, "axes.linewidth": AXIS_WIDTH,
                    "xtick.major.width": AXIS_WIDTH, "ytick.major.width": AXIS_WIDTH,
                    "xtick.major.size": 2, "ytick.major.size": 2, "axes.spines.top": False,
                    "axes.spines.right": False, "pdf.fonttype": 42})


def palettes():
    panels = annotation("recorder_edit_summary_panels")
    labels, light, dark = {}, {}, {}
    for row in panels.to_dict("records"):
        enzymes = row["enzymes"].split(";")
        for field, target in [("enzyme_labels", labels), ("light_colors", light), ("dark_colors", dark)]:
            values = row[field].split(";")
            require(len(values) == len(enzymes), "Incomplete enzyme presentation annotation")
            for enzyme, value in zip(enzymes, values):
                if field != "enzyme_labels":
                    require(enzyme not in target or target[enzyme] == value, "Conflicting enzyme colors")
                target.setdefault(enzyme, value)
    return labels, light, dark


def distance_palettes():
    return dict(zip(enzyme_metadata().tada_type, ENZYME_DISTANCE_SHADES))


def place_shared_ylabel(figure, ax, label, gap=LABEL_GAP):
    """Set a figure-wide y label the standard gap from the tick numbers it names, not at the page edge.

    `supylabel` sits at a fixed fraction of the page, so the same label drifts away from its ticks as
    soon as a figure is resized, and two figures of the set then carry it at different distances. This
    measures the drawn tick numbers of `ax` and places the label `gap` inches to their left, defaulting
    to the gap the substitution heat map and logo panels carry.
    Measuring needs the resolved layout, so the layout engine is frozen first and the panels stay put.
    A label broken over two lines centres them on each other and sits on the middle of the panels it
    names, rather than on the middle of a page whose top and bottom bands are not the same height.
    """
    text = figure.supylabel(label, multialignment="center")
    figure.canvas.draw()
    figure.set_layout_engine("none")
    renderer = figure.canvas.get_renderer()
    ticks = [number.get_window_extent(renderer) for number in ax.get_yticklabels() if number.get_text()]
    require(len(ticks) > 0, "The labelled axis carries no tick numbers to place its label against")
    shift = min(box.x0 for box in ticks) - gap * figure.dpi - text.get_window_extent(renderer).x1
    text.set_x(text.get_position()[0] + shift / figure.bbox.x1)
    spans = [panel.get_position() for panel in figure.axes]
    text.set_y((min(span.y0 for span in spans) + max(span.y1 for span in spans)) / 2)
    require(text.get_window_extent(renderer).x0 >= 0, "The y label does not fit the page")
    return text


def centre_shared_xlabel(figure, text):
    """Put a figure-wide x label on the middle of the panels it names, not the middle of the page.

    `supxlabel` sits at a fixed fraction of the page, and the page is not symmetric about its panels:
    the y label and its tick numbers take the left margin, so a page-centred label reads to the left of
    the axes it names. This measures the drawn panels and centres the label on them, as
    `place_shared_ylabel` does vertically, and is called once the layout engine is frozen.
    """
    spans = [panel.get_position() for panel in figure.axes]
    text.set_x((min(span.x0 for span in spans) + max(span.x1 for span in spans)) / 2)
    return text


def logo_letter(ax, base, x, bottom, width, height):
    """Draw one vector sequence-logo glyph scaled to its stacked height."""
    path = TextPath((0, 0), base, size=1, prop=font_manager.FontProperties(family=FONT, weight="bold"))
    box = path.get_extents()
    require(box.width > 0 and box.height > 0, f"Empty logo glyph for {base}")
    transform = Affine2D().translate(-box.x0, -box.y0).scale(width / box.width, height / box.height).translate(x - width / 2, bottom)
    ax.add_patch(PathPatch(transform.transform_path(path), facecolor=NUCLEOTIDE_COLORS[base], edgecolor="none"))


def save_pdf(figure, name, transparent=False, png=False, compact=False):
    os.makedirs(FIGURE_DIRECTORY, exist_ok=True)
    output = OUTPUTS.loc[OUTPUTS.source.eq(name), "filename"]
    require(len(output) == 1, f"Expected one registered PDF output for {name}")
    filename = f"{FIGURE_DIRECTORY}/{output.item()}"
    bounds = {"bbox_inches": "tight", "pad_inches": COMPACT_PADDING} if compact else {}
    figure.savefig(filename, transparent=transparent, metadata=PDF_METADATA, **bounds)
    if png:
        png_filename = f"{filename[:-4]}.png"
        figure.savefig(png_filename, transparent=transparent, dpi=PNG_DPI, **bounds)
        print(f"Saved {png_filename}; PDF/PNG size ratio={os.path.getsize(filename) / os.path.getsize(png_filename):.2f}")
    plt.close(figure)
    print(f"Saved {filename}")

