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

    # Read data from the parameters CSV file
    parameters_df = pd.read_csv(experiment_path / "parameters.csv")
    simulation_id = int(simulation_id)
    sim_params = parameters_df.iloc[simulation_id]
    wind_speed = sim_params["wind_speed"]
    control_fuel_height = sim_params["control_fuel_height"]
    control_fuel_load = sim_params["control_fuel_load"]
    control_fuel_moisture_content = sim_params["control_fuel_moisture_content"]
    control_fuel_sav = sim_params["control_sav"]
    treatment_fuel_height = sim_params["treatment_fuel_height"]
    treatment_fuel_load = sim_params["treatment_fuel_load"]
    treatment_fuel_moisture_content = sim_params["treatment_fuel_moisture_content"]
    treatment_fuel_sav = sim_params["treatment_sav"]
    circle_radius = sim_params["circle_radius"]
    fuel_model = sim_params["fuel_model"]
    resolution = sim_params["resolution"]

    # Get the correct number of cells for the resolution
    i = 180 if resolution == "fine" else 90
    j = 450 if resolution == "fine" else 225
    k = 200 if resolution == "fine" else 100
    ijk_lines = f"{i}, {j}, {k}"

    # Create the template for the simulation file
    template_dict = {}
    template_dict["wind_speed"] = wind_speed
    template_dict["control_fuel_height"] = control_fuel_height
    template_dict["control_fuel_load"] = control_fuel_load
    template_dict["control_fuel_moisture_content"] = control_fuel_moisture_content
    template_dict["control_sav"] = control_fuel_sav
    template_dict["treatment_fuel_height"] = treatment_fuel_height
    template_dict["treatment_fuel_load"] = treatment_fuel_load
    template_dict["treatment_fuel_moisture_content"] = treatment_fuel_moisture_content
    template_dict["treatment_sav"] = treatment_fuel_sav
    template_dict["circle_radius"] = circle_radius
    template_dict["fuel_model"] = fuel_model
    template_dict["ijk_lines"] = ijk_lines
    template_dict["title"] = f"Experiment iteration: {simulation_id}"
    template_dict["chid"] = f"out_{simulation_id}"

    # Calculate the mass per volume for the fuel models if using BFM
    if fuel_model == "bfm":
        treatment_mass_per_volume = treatment_fuel_load / treatment_fuel_height
        control_mass_per_volume = control_fuel_load / control_fuel_height
        template_dict["treatment_mass_per_volume"] = treatment_mass_per_volume
        template_dict["control_mass_per_volume"] = control_mass_per_volume

    # Write the FDS input file using the template
    template_path = TEMPLATE_DIR / f"{fuel_model}_template.fds"
    output_path = experiment_path / "output"
    simulation_path = output_path / f"simulation_{simulation_id}"
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
            control_fuel_load=control_fuel_load,
            treatment_fuel_height=treatment_fuel_height,
            treatment_fuel_load=treatment_fuel_load,
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
            "/home/anthony/Work/UM/bandy-circles/experiments/1D-grid-search-coarse",
            0,
        )
