import json
import pandas as pd
from pathlib import Path
from itertools import product

current_file_path = Path(__file__).parent


def main():
    # Load the configuration file
    with open(current_file_path / "config.json") as f:
        config = json.load(f)

    # Define a list to hold the parameters for each simulation
    data_list = []
    simulation_id = 0

    # Generate all combinations of parameters
    for (
        fuel_model,
        circle_radius,
        control_fuel_height,
        treatment_fuel_height,
        control_fuel_moisture_content,
        treatment_fuel_moisture_content,
        control_sav,
        treatment_sav,
        fuel_load,
        wind_speed,
    ) in product(
        config["fuel_model"],
        config["circle_radius"],
        config["control_fuel_height"],
        config["treatment_fuel_height"],
        config["control_fuel_moisture_content"],
        config["treatment_fuel_moisture_content"],
        config["control_sav"],
        config["treatment_sav"],
        config["fuel_load"],
        config["wind_speed"],
    ):
        data_list.append(
            {
                "simulation_id": simulation_id,
                "fuel_model": fuel_model,
                "circle_radius": circle_radius,
                "control_fuel_height": control_fuel_height,
                "treatment_fuel_height": treatment_fuel_height,
                "control_fuel_moisture_content": control_fuel_moisture_content,
                "treatment_fuel_moisture_content": treatment_fuel_moisture_content,
                "control_sav": control_sav,
                "treatment_sav": treatment_sav,
                "fuel_load": fuel_load,
                "wind_speed": wind_speed,
            }
        )
        simulation_id += 1

    # Create the dataframe
    df = pd.DataFrame(data_list)

    # Write the dataframe to a CSV file
    df.to_csv(current_file_path / "parameters.csv", index=False)

    print(f"Generated {len(data_list)} parameter combinations.")


if __name__ == "__main__":
    main()
