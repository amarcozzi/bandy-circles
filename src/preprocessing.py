import sys
import pandas as pd
from pathlib import Path
from string import Template

CURRENT_DIR = Path(__file__).parent
TEMPLATE_DIR = CURRENT_DIR / "templates"


def main(experiment_path, experiment_id):
    print("Preprocessing")

    # Read the grid search CSV file
    grid_search = pd.read_csv(experiment_path / "grid_search.csv")

    # Get the wind speed from the grid search using the experiment ID
    wind_speed = float(grid_search.iloc[int(experiment_id)]["wind_speed"])

    # Create the template for the experiment file
    template_dict = {}
    template_dict["wind_speed"] = wind_speed
    template_dict["experiment_id"] = experiment_id
    template_dict["title"] = f"Experiment iteration: {experiment_id}"
    template_dict["chid"] = f"out_{experiment_id}"

    # Write the FDS input file using the template
    template_path = TEMPLATE_DIR / "pfm_template.fds"
    output_path = experiment_path / "output"
    fds_input_path = output_path / f"input_{experiment_id}.fds"
    with open(template_path, "r") as f:
        template = Template(f.read())
        fds_input = template.substitute(template_dict)
        with open(fds_input_path, "w") as f:
            f.write(fds_input)


if __name__ == "__main__":
    """
    Arguments:
    0: Name of the script
    1: Path to the experiment folder
    2: Experiment ID
    """
    main(Path(sys.argv[1]), sys.argv[2])
