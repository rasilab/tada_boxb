"""Manuscript styling and figure output for the in vivo boxB figures.

The visual vocabulary is the in vitro one, imported rather than restated, so the two analyses carry the
same font, sizes, label gap, axis colours, nucleotide colours, and logo glyphs, and so a change to that
vocabulary reaches both. Only what must differ is defined here: the figure directory, this analysis's
output manifest, and the enzyme palette, which has to name the no-editor background the in vitro
figures do not carry. The in vitro module is imported, never edited, so its figures cannot move.

Every figure is written as a vector PDF and as a 96-DPI PNG, so a panel can be read on its own as
well as assembled into the manuscript figure.
"""

import os

from analysis_helpers import annotation, require, samples
from in_vitro_reuse import style as in_vitro_style

ANALYSIS_DIRECTORY = "analysis/boxb_in_vivo_sequencing"
# The whole visual vocabulary comes from the in vitro style module, so one change reaches both analyses.
AXIS_COLOR = in_vitro_style.AXIS_COLOR
AXIS_WIDTH = in_vitro_style.AXIS_WIDTH
DISTANCE_SHADES = in_vitro_style.DISTANCE_SHADES
FONT = in_vitro_style.FONT
FONT_SIZE = in_vitro_style.FONT_SIZE
LABEL_GAP = in_vitro_style.LABEL_GAP
MUTANT_COLORS = in_vitro_style.MUTANT_COLORS
NUCLEOTIDE_COLORS = in_vitro_style.NUCLEOTIDE_COLORS
PDF_METADATA = in_vitro_style.PDF_METADATA
centre_shared_xlabel = in_vitro_style.centre_shared_xlabel
logo_letter = in_vitro_style.logo_letter
place_shared_ylabel = in_vitro_style.place_shared_ylabel
plt = in_vitro_style.plt

FIGURE_DIRECTORY = f"{ANALYSIS_DIRECTORY}/figures"
PNG_DPI = 96
# Editors keep the in vitro colours so the same enzyme reads the same way across the two experiments.
# The no-editor background is new here and takes its own entry. Every enzyme the library annotation
# carries must appear exactly once, which is what palettes() checks.
ENZYME_COLORS = {"tada_ln": "#AA4499", "tada_gfpnb": "#E69F00", "none": "#777777"}
ENZYME_LIGHT = {"tada_ln": "#E2B8D6", "tada_gfpnb": "#F3CF80", "none": "#CCCCCC"}


def palettes():
    """Each enzyme's presentation name and colours, taken from the library annotation and this palette.

    A figure takes the name from here rather than from any other column, so every panel of the set names
    an enzyme the same way.
    """
    labels = samples().drop_duplicates("enzyme").set_index("enzyme").enzyme_label.to_dict()
    require(set(labels) == set(ENZYME_COLORS) == set(ENZYME_LIGHT), "The enzyme palette and the library annotation disagree")
    return labels, ENZYME_LIGHT, ENZYME_COLORS


def figure_outputs(analysis):
    """This analysis's registered outputs, which is the only place a filename is decided."""
    outputs = annotation("figure_outputs")
    require(not outputs.duplicated(["source"]).any() and not outputs.duplicated(["filename"]).any(),
            "The output manifest repeats a source or a filename")
    selected = outputs.loc[outputs.analysis.eq(analysis)]
    require(len(selected) > 0, f"No registered output for {analysis}")
    return selected


def save_figure(figure, name, transparent=False):
    """Write one registered figure as a vector PDF and a 96-DPI PNG."""
    os.makedirs(FIGURE_DIRECTORY, exist_ok=True)
    outputs = annotation("figure_outputs")
    output = outputs.loc[outputs.source.eq(name), "filename"]
    require(len(output) == 1, f"Expected one registered output for {name}")
    stem = f"{FIGURE_DIRECTORY}/{output.item()}"
    require(stem.endswith(".pdf"), "Register every figure under its PDF filename")
    figure.savefig(stem, transparent=transparent, metadata=PDF_METADATA)
    figure.savefig(stem.replace(".pdf", ".png"), transparent=transparent, dpi=PNG_DPI)
    plt.close(figure)
    print(f"Saved {stem} and its PNG")
