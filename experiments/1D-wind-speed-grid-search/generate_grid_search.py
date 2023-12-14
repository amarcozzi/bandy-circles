"""
This file generates the grid search for the 1D wind speed experiment. It outputs a CSV file with the following columns:
    - experiment_id
    - wind_speed
"""

import numpy as np
import pandas as pd

NUMBER_SAMPLES = 100
WIND_SPEED_MIN = 0.5  # m/s
WIND_SPEED_MAX = 5.0  # m/s


def main():
    # Generate the wind speed values
    wind_speeds = np.linspace(WIND_SPEED_MIN, WIND_SPEED_MAX, NUMBER_SAMPLES)

    # Create the dataframe
    df = pd.DataFrame(
        {"experiment_id": np.arange(NUMBER_SAMPLES), "wind_speed": wind_speeds}
    )

    # Write the dataframe to a CSV file
    df.to_csv("grid_search.csv", index=False)


if __name__ == "__main__":
    main()
