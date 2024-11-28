from generate_bdf_files import generate_bdf_files

import sys
import shutil
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from string import Template

CURRENT_DIR = Path(__file__).parent
TEMPLATE_DIR = CURRENT_DIR / "templates"


def preprocess_simulation(experiment_path, parameters_df, simulation_id):
    simulation_id = int(simulation_id)
    sim_params = parameters_df.loc[
        parameters_df["simulation_id"] == simulation_id
    ].iloc[0]

    # Extract parameters
    wind_speed = sim_params["wind_speed"]
    fuel_load = sim_params["fuel_load"]
    control_fuel_height = sim_params["control_fuel_height"]
    control_fuel_moisture_content = sim_params["control_fuel_moisture_content"]
    control_sav = sim_params["control_sav"]
    treatment_fuel_height = sim_params["treatment_fuel_height"]
    treatment_fuel_moisture_content = sim_params["treatment_fuel_moisture_content"]
    treatment_sav = sim_params["treatment_sav"]
    circle_radius = sim_params["circle_radius"]
    fuel_model = sim_params["fuel_model"]

    # Create the template for the simulation file
    template_dict = {
        "wind_speed": wind_speed,
        "control_fuel_height": control_fuel_height,
        "control_fuel_load": fuel_load,
        "control_fuel_moisture_content": control_fuel_moisture_content,
        "control_sav": control_sav,
        "treatment_fuel_height": treatment_fuel_height,
        "treatment_fuel_load": fuel_load,
        "treatment_fuel_moisture_content": treatment_fuel_moisture_content,
        "treatment_sav": treatment_sav,
        "circle_radius": circle_radius,
        "fuel_model": fuel_model,
        "title": f"Experiment iteration: {simulation_id}",
        "chid": f"out_{simulation_id}",
    }

    # Calculate the mass per volume for the fuel models if using BFM
    if fuel_model == "bfm":
        treatment_mass_per_volume = fuel_load / treatment_fuel_height
        control_mass_per_volume = fuel_load / control_fuel_height
        template_dict["treatment_mass_per_volume"] = treatment_mass_per_volume
        template_dict["control_mass_per_volume"] = control_mass_per_volume

    # Write the FDS input file using the template
    template_path = TEMPLATE_DIR / f"{fuel_model}_template.fds"
    all_simulations_path = experiment_path / "simulations"
    if not all_simulations_path.exists():
        all_simulations_path.mkdir()
    simulation_path = all_simulations_path / f"simulation_{simulation_id}"
    if not simulation_path.exists():
        simulation_path.mkdir()
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
    dx = 0.1
    dy = 0.1
    xmin = -10
    xmax = 10
    ymin = -12
    ymax = 12
    if fuel_model == "pfm":
        generate_bdf_files(
            out_path=simulation_path,
            circle_radius=circle_radius,
            control_fuel_height=control_fuel_height,
            control_fuel_load=fuel_load,
            treatment_fuel_height=treatment_fuel_height,
            treatment_fuel_load=fuel_load,
            dx=dx,
            dy=dy,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )

    # Add a slurm template to the simulation folder
    slurm_template_path = TEMPLATE_DIR / "submit_template.slurm"
    slurm_path = simulation_path / "submit.slurm"
    slurm_data = {
        "job_id": simulation_id,
        "job_name": f"simulation_{simulation_id}",
        "job_name_output": f"simulation_{simulation_id}.out",
        "job_name_error": f"simulation_{simulation_id}.err",
        "fds_input_file_name": f"input_{simulation_id}.fds",
    }
    with open(slurm_template_path, "r") as f:
        template = Template(f.read())
        template.delimiter = "%"
        slurm = template.substitute(slurm_data)
        with open(slurm_path, "w") as f:
            f.write(slurm)


def main(experiment_path, simulation_id=None):
    if isinstance(experiment_path, str):
        experiment_path = Path(experiment_path)

    # Read data from the parameters CSV file
    parameters_df = pd.read_csv(experiment_path / "parameters.csv")

    if simulation_id:
        preprocess_simulation(experiment_path, parameters_df, simulation_id)
    else:
        for row in tqdm(parameters_df.iterrows()):
            simulation_id = row[0]
            preprocess_simulation(experiment_path, parameters_df, simulation_id)


if __name__ == "__main__":
    """
    Arguments:
    0: Name of the script
    1: Path to the experiment folder
    2: Simulation ID
    """
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main("/Users/anthony/Work/UM/bandy-circles/experiments/grid-search-bfm")
