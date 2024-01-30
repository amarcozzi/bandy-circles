"""
This file generates the grid search for the 1D wind speed experiment. It outputs a parameters.csv file with the following columns:
    - experiment_id
    - wind_speed
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import product

NUMBER_SAMPLES = 100
WIND_SPEED_MIN = 0.5  # m/s
WIND_SPEED_MAX = 5.0  # m/s

current_file_path = Path(__file__).parent


def main():
    # Load the congiguration file
    with open(current_file_path / "config.json") as f:
        config = json.load(f)

    # Generate the wind speed values
    wind_speeds = np.linspace(
        config["wind_speed"]["min"],
        config["wind_speed"]["max"],
        config["wind_speed"]["number_of_samples"],
    )
    wind_speeds = np.round(wind_speeds, 4)

    # Define a dictionary to hold the parameters for each simulation
    data_list = []
    simulation_id = 0
    for resolution in config["resolution"]:
        for fuel_model in config["fuel_model"]:
            for fuel_load_category in config["fuel_load_category"]:
                # Generate fuel load values for the control and treatment for each fuel load category
                control_fuel_load = config["fuel_load"][fuel_load_category]
                treatment_fuel_loads = np.linspace(
                    control_fuel_load * config["fuel_load_ratios"]["min"],
                    control_fuel_load * config["fuel_load_ratios"]["max"],
                    config["fuel_load_ratios"]["number_of_samples"],
                )
                for treatment_fuel_load in treatment_fuel_loads:
                    for wind_speed in wind_speeds:
                        data_list.append(
                            {
                                "simulation_id": simulation_id,
                                "wind_speed": wind_speed,
                                "fuel_load_category": fuel_load_category,
                                "control_fuel_load": control_fuel_load,
                                "treatment_fuel_load": treatment_fuel_load,
                                "control_fuel_height": config["control_fuel_height"],
                                "treatment_fuel_height": config[
                                    "treatment_fuel_height"
                                ],
                                "control_fuel_moisture_content": config[
                                    "control_fuel_moisture_content"
                                ],
                                "treatment_fuel_moisture_content": config[
                                    "treatment_fuel_moisture_content"
                                ],
                                "control_sav": config["control_sav"],
                                "treatment_sav": config["treatment_sav"],
                                "circle_radius": config["circle_radius"],
                                "fuel_model": fuel_model,
                                "resolution": resolution,
                            }
                        )
                        simulation_id += 1

    # Create the dataframe
    df = pd.DataFrame(data_list)

    # Write the dataframe to a CSV file
    df.to_csv(current_file_path / "parameters.csv", index=False)


if __name__ == "__main__":
    main()
