"""Load the in vitro helper and style modules under their own names, so both analyses can be reused.

The two analyses each carry an `analysis_helpers` and a `plot_style`, so a plain import from the in vivo
scripts would either find the wrong one or re-enter the module that is still being imported. Loading the
in vitro files explicitly, under distinct module names, lets the in vivo scripts take the shared
denominator definition and the whole visual vocabulary from the in vitro modules without editing them,
which is what keeps the in vitro figures from moving when this analysis changes.

The in vitro style module resolves `analysis_helpers` by name, so that name points at the in vitro helper
for the duration of its import and is then restored to whatever the importing script had.
"""

import importlib.util
import sys

IN_VITRO_SCRIPTS = "analysis/boxb_in_vitro_sequencing/scripts"


def _load(alias, filename):
    if alias in sys.modules:
        return sys.modules[alias]
    spec = importlib.util.spec_from_file_location(alias, f"{IN_VITRO_SCRIPTS}/{filename}")
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {filename} from {IN_VITRO_SCRIPTS}; run from the repository root")
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias] = module
    spec.loader.exec_module(module)
    return module


helpers = _load("in_vitro_analysis_helpers", "analysis_helpers.py")

_previous = sys.modules.get("analysis_helpers")
sys.modules["analysis_helpers"] = helpers
try:
    style = _load("in_vitro_plot_style", "plot_style.py")
finally:
    if _previous is None:
        sys.modules.pop("analysis_helpers", None)
    else:
        sys.modules["analysis_helpers"] = _previous
