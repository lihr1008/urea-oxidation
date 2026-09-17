from pathlib import Path
import sys

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
ML_ROOT = next(parent for parent in SCRIPT_DIR.parents if parent.name == "MachineLearning")
if str(ML_ROOT) not in sys.path:
    sys.path.insert(0, str(ML_ROOT))

from composition_ilr_pca import save_ilr_pca_outputs

data = pd.read_excel(SCRIPT_DIR / 'analysis.xlsx')
save_ilr_pca_outputs(
    data,
    SCRIPT_DIR,
    ML_ROOT / "data" / "experiment.xlsx",
    target_column=data.columns[-1],
    target_label="potential",
)
