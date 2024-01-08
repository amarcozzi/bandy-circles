"""
This file generates the grid search for the 1D wind speed experiment. It outputs a parameters.csv file with the following columns:
    - experiment_id
    - wind_speed
"""

import numpy as np
import pandas as pd
from pathlib import Path

NUMBER_SAMPLES = 100
WIND_SPEED_MIN = 0.5  # m/s
WIND_SPEED_MAX = 5.0  # m/s

current_file_path = Path(__file__).parent


def main():
    # Generate the wind speed values
    wind_speeds = np.linspace(WIND_SPEED_MIN, WIND_SPEED_MAX, NUMBER_SAMPLES)
    wind_speeds = np.round(wind_speeds, 4)

    # Create the dataframe
    df = pd.DataFrame(
        {"experiment_id": np.arange(NUMBER_SAMPLES), "wind_speed": wind_speeds}
    )

    # Write the dataframe to a CSV file
    df.to_csv(current_file_path / "parameters.csv", index=False)


if __name__ == "__main__":
    main()
