import sys
import numpy as np
import xarray as xr
import pandas as pd
import warnings
import logging
from tqdm import tqdm
from pathlib import Path
from fdsreader import Simulation
from typing import Optional
from shutil import make_archive, rmtree

# Suppress warnings
warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)


def create_particle_data_array(sim: Simulation, sim_dir: Path):
    """
    Process FDS particle data into a zarr array. If the zarr file already exists,
    load it instead of reprocessing.
    """
    # Define grid
    resolution = 0.1
    x_coords = np.arange(-10 + resolution / 2, 10, resolution)
    y_coords = np.arange(-12 + resolution / 2, 12, resolution)
    z_coords = np.arange(0 + resolution / 2, 0.6, resolution)

    x_edges = np.arange(-10, 10 + resolution, resolution)
    y_edges = np.arange(-12, 12 + resolution, resolution)
    z_edges = np.arange(0, 0.6 + resolution, resolution)

    # Initialize dataset
    dataset = xr.Dataset(
        coords={
            "time": sim.particles.times,
            "z": z_coords,
            "y": y_coords,
            "x": x_coords,
        }
    )

    # Get particle point data positions
    # Get particle point positions
    treatment_particles = sim.particles[0]
    treatment_particles_positions = treatment_particles.positions
    control_particles = sim.particles[1]
    control_particles_positions = control_particles.positions

    # Create a single histogram across all time steps for treatment particles
    all_treatment_positions = np.vstack(treatment_particles_positions)
    treatment_hist, _ = np.histogramdd(
        all_treatment_positions,
        bins=(x_edges, y_edges, z_edges),
    )

    # Mark cells as TRUE where treatment particles exist at any time
    # Add treatment flag to dataset (3D: x, y, z)
    treatment_flag = treatment_hist > 0
    dataset["TREATMENT"] = (("x", "y", "z"), treatment_flag)

    for quantity_object in treatment_particles.quantities:
        quantity = quantity_object.name
        treatment_particle_quantity = treatment_particles.data[quantity]
        control_particle_quantity = control_particles.data[quantity]

        # Initialize data array for this quantity
        gridded_data = np.zeros(
            (len(dataset.time), len(dataset.z), len(dataset.y), len(dataset.x))
        )

        # For each timestep
        for t in range(len(sim.particles.times)):
            # Combine all particle positions
            treatment_positions = treatment_particles_positions[t]
            control_positions = control_particles_positions[t]
            positions = np.concatenate(
                [treatment_positions, control_positions]
            )  # Shape: (n_particles, 3)

            # Combine all particle quantities
            treatment_quantities = treatment_particle_quantity[t]
            control_quantities = control_particle_quantity[t]
            values = np.concatenate(
                [treatment_quantities, control_quantities]
            )  # Shape: (n_particles,)

            # Bin the particle data into 3D histogram
            hist, _ = np.histogramdd(
                (positions[:, 2], positions[:, 1], positions[:, 0]),
                bins=(z_edges, y_edges, x_edges),
                weights=values,
            )

            # Count number of particles in each cell
            counts, _ = np.histogramdd(
                (positions[:, 2], positions[:, 1], positions[:, 0]),
                bins=(z_edges, y_edges, x_edges),
            )

            # Avoid division by zero by masking where counts is 0
            mask = counts > 0
            hist[mask] /= counts[mask]  # Average value per cell

            gridded_data[t] = hist

        # Add to dataset
        dataset[quantity] = xr.DataArray(
            data=gridded_data,
            dims=("time", "z", "y", "x"),
            coords={
                "time": dataset.time,
                "z": dataset.z,
                "y": dataset.y,
                "x": dataset.x,
            },
            name=quantity,
        )

    # Remove data before t=0 seconds
    dataset = dataset.sel(time=dataset.time > 0)

    # Reverse y-axis to match array coordinates
    dataset = dataset.isel(y=slice(None, None, -1))

    dataset.to_zarr(sim_dir / "particle_data.zarr", mode="w")

    # Create a zip archive of the zarr file
    make_archive(sim_dir / "particle_data", "zip", sim_dir / "particle_data.zarr")

    # Remove the zarr directory
    rmtree(sim_dir / "particle_data.zarr")

    return dataset


def process_simulation(sim_dir: Path):
    """Process a single simulation's binary data to zarr format"""
    try:
        sim = Simulation(str(sim_dir))
        create_particle_data_array(sim, sim_dir)
        return True
    except Exception as e:
        print(f"Error processing {sim_dir}: {str(e)}")
        return False


def process_experiment(experiment_dir: Path, simulation_id: Optional[int] = None):
    """Process binary FDS data to zarr format for an entire experiment or single simulation"""
    simulations_dir = experiment_dir / "simulations"

    # Get all simulation directories
    if simulation_id is not None:
        sim_dirs = [simulations_dir / f"simulation_{simulation_id}"]
        if not sim_dirs[0].exists():
            raise ValueError(
                f"Simulation {simulation_id} not found in {simulations_dir}"
            )
    else:
        sim_dirs = sorted(simulations_dir.glob("simulation_*"))
        if not sim_dirs:
            raise ValueError(f"No simulation directories found in {simulations_dir}")

    # Process each simulation
    for sim_dir in tqdm(sim_dirs, desc="Processing simulations"):
        process_simulation(sim_dir)
        break


def main():
    # Set default test directory if no arguments provided
    default_dir = Path("/Users/anthony/Work/UM/bandy-circles/experiments/grid-search")

    if len(sys.argv) == 1:
        # No arguments - use default directory and process all simulations
        process_experiment(default_dir)
    elif len(sys.argv) == 2:
        # One argument - use provided directory and process all simulations
        process_experiment(Path(sys.argv[1]))
    elif len(sys.argv) == 3:
        # Two arguments - use provided directory and specific simulation
        process_experiment(Path(sys.argv[1]), int(sys.argv[2]))
    else:
        print(
            "Usage: python process_binary_data_to_zarr.py [experiment_dir] [simulation_id]"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
