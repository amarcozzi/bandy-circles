from generate_bdf_files import generate_bdf_files

import sys
import json
import shutil
import pandas as pd
from pathlib import Path
from string import Template

CURRENT_DIR = Path(__file__).parent
TEMPLATE_DIR = CURRENT_DIR / "templates"


def main(experiment_path, simulation_id):
    if isinstance(experiment_path, str):
        experiment_path = Path(experiment_path)

    # Read the configuration JSON file
    with open(experiment_path / "config.json") as f:
        config = json.load(f)
    resolution = config["resolution"]
    fuel_model = config["fuel_model"]

    # Try and read the circle, control, and treatment data from the config file.
    # If it doesn't exist, data will come from the parameters CSV file.
    control_fuel_height = config.get("control_fuel_height")
    control_fuel_bulk_density = config.get("control_fuel_bulk_density")
    treatment_fuel_height = config.get("treatment_fuel_height")
    treatment_fuel_bulk_density = config.get("treatment_fuel_bulk_density")
    circle_radius = config.get("circle_radius")

    # Read data from the parameters CSV file
    parameters_df = pd.read_csv(experiment_path / "parameters.csv")
    wind_speed = float(parameters_df.iloc[int(simulation_id)]["wind_speed"])

    # Create the template for the simulation file
    template_dict = {}
    template_dict["wind_speed"] = wind_speed
    template_dict["title"] = f"Experiment iteration: {simulation_id}"
    template_dict["chid"] = f"out_{simulation_id}"

    # Write the FDS input file using the template
    template_path = TEMPLATE_DIR / f"{fuel_model}_template_{resolution}.fds"
    output_path = experiment_path / "output"
    simulation_path = output_path / f"{simulation_id}"
    fds_input_file_path = simulation_path / f"input_{simulation_id}.fds"
    with open(template_path, "r") as f:
        template = Template(f.read())
        fds_input = template.substitute(template_dict)
        with open(fds_input_file_path, "w") as f:
            f.write(fds_input)

    # Copy vegetation_model.txt from templates to the simulation folder
    vegetation_model_template_path = TEMPLATE_DIR / "vegetation_model.txt"
    vegetation_model_out_path = simulation_path / "vegetation_model.txt"
    shutil.copy(vegetation_model_template_path, vegetation_model_out_path)

    # If the fuel model is pfm, write the bulk density files (.bdf)
    dx = 0.1 if resolution == "fine" else 0.2
    dy = 0.1 if resolution == "fine" else 0.2
    if fuel_model == "pfm":
        generate_bdf_files(
            out_path=simulation_path,
            circle_radius=circle_radius,
            control_fuel_height=control_fuel_height,
            control_fuel_bulk_density=control_fuel_bulk_density,
            treatment_fuel_height=treatment_fuel_height,
            treatment_fuel_bulk_density=treatment_fuel_bulk_density,
            dx=dx,
            dy=dy,
        )


if __name__ == "__main__":
    """
    Arguments:
    0: Name of the script
    1: Path to the experiment folder
    2: Simulation ID
    """
    try:
        main(Path(sys.argv[1]), sys.argv[2])
    except IndexError:
        main(
            "/home/anthony/Work/UM/bandy-circles/experiments/1D-wind-speed-grid-search",
            0,
        )
