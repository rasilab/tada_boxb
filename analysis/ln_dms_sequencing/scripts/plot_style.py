"""Shared Helvetica typography, display settings, and PNG/PDF export for the seven retained figures."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from analysis_helpers import ANALYSIS_DIRECTORY

FONT = "Helvetica"
FONT_SIZE = 7
LABEL_GAP = .09
AXIS_COLOR = "#999999"
AXIS_WIDTH = 0.25
PNG_DPI = 96
FIGURE_DIRECTORY = f"{ANALYSIS_DIRECTORY}/figures"
PDF_METADATA = {"Creator": "λN DMS Python analysis", "CreationDate": None, "ModDate": None}

font_manager.findfont(FONT, fallback_to_default=False)
plt.rcParams.update({"font.family": FONT, "font.size": FONT_SIZE, "figure.titlesize": FONT_SIZE,
                     "figure.labelsize": FONT_SIZE, "axes.titlesize": FONT_SIZE,
                     "axes.labelsize": FONT_SIZE, "xtick.labelsize": FONT_SIZE,
                     "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
                     "legend.title_fontsize": FONT_SIZE, "axes.edgecolor": AXIS_COLOR,
                     "axes.linewidth": AXIS_WIDTH, "xtick.major.width": AXIS_WIDTH,
                     "ytick.major.width": AXIS_WIDTH, "xtick.major.size": 2,
                     "ytick.major.size": 2, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})


def save_figure(figure, name):
    """Write the named figure as vector PDF and 96-DPI PNG."""
    os.makedirs(FIGURE_DIRECTORY, exist_ok=True)
    stem = f"{FIGURE_DIRECTORY}/{name}"
    figure.savefig(f"{stem}.pdf", transparent=False, metadata=PDF_METADATA)
    figure.savefig(f"{stem}.png", transparent=False, dpi=PNG_DPI)
    plt.close(figure)
    print(f"Saved {stem}.pdf and {stem}.png")
