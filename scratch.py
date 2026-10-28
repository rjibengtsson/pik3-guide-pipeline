from __future__ import annotations

import pandas as pd
import numpy as np
import os, sys


DEFAULT_MIN_EFFICACY = 0.9
DEFAULT_MERGE_GAP = 0
DEFAULT_MIN_WINDOW = 1500
DEFAULT_MAX_WINDOW = 2000
DEFAULT_SIZE_STEP = 250
DEFAULT_STEP = 50


def find_hotspot_regions(predictions: pd.DataFrame,
                         min_efficacy: float = DEFAULT_MIN_EFFICACY,
                         merge_gap: int = DEFAULT_MERGE_GAP):

    qualifying = predictions[predictions["predicted_efficacy"] >= min_efficacy]

    for transcript, group in qualifying.groupby("transcript_id", sort=True):
        rows = group.sort_values(["start", "end"]).itertuples(index=False)
        for row in rows:
            print(row.start)
            exit(0)



def main():
    prediction_file = os.path.join(os.path.dirname(__file__), "results/pik3/pik3-predictions.csv")
    data = pd.read_csv(prediction_file)
    find_hotspot_regions(data)

if __name__ == "__main__":
    main()