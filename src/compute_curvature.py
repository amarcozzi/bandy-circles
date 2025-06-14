"""
compute_curvature.py
"""

import sys
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
from typing import Optional, List, Tuple
from xarray import DataArray
from multiprocessing import Pool, cpu_count
from functools import partial
from postprocess import (
    get_particle_data_array,
    get_polynomial,
    get_active_fire_array,
    get_fire_line,
    get_fire_line_ideal,
    find_fire_intersection_time,
)
from plotting import plot_fire_front_evolution, plot_fire_front_with_polynomial_fit, plot_combined_fire_front_analysis, \
    plot_polynomials

NUM_WORKERS = 7

# TARGET_SIM_ID = 267
# TARGET_SIM_ID = 125
TARGET_SIM_ID = None
QUANTITY = "heat_flux"

R_C_MAX = 2.7


def get_curvature(polynomial_coeffs: list) -> list:
    # Get curvatures (2 * quadratic coefficient) for all polynomials
    return [2 * coeffs[0] for coeffs in polynomial_coeffs]


def process_simulation(args) -> tuple:
    """
    Process a single simulation and return its ID and results.
    Modified to accept a single argument tuple for Pool.map compatibility.
    """
    sim_directory, sim_params = args
    sim_id = int(sim_params["simulation_id"])

    try:
        # Get particle data over time as a dataset
        data = get_particle_data_array(sim_directory)

        data["ACTIVE FIRE"] = (data["PARTICLE TOTAL HEAT FLUX"].mean("z") < -20).compute()
        # rolling_heat_flux = data["PARTICLE TOTAL HEAT FLUX"].rolling(time=5, center=True).mean().compute()
        # data["ACTIVE FIRE"] = rolling_heat_flux.mean("z") < -20

        # particle_temp_rolling_avg = data["PARTICLE TEMPERATURE"].mean("z").rolling(time=5, center=True).mean().compute()
        # data["ACTIVE FIRE"] = particle_temp_rolling_avg > 200
        # data["ACTIVE FIRE"] = (data["PARTICLE TEMPERATURE"].mean("z") > 350).compute()

        # mass = data["PARTICLE MASS"].compute()
        # weighting_factor = data["PARTICLE WEIGHTING FACTOR"].compute()
        # weighting_factor_adj = DataArray(
        #     data=np.where(data["TREATMENT"], weighting_factor.max(), weighting_factor.median()),
        #     dims=["x", "y", "z", ])
        # mass_adj = mass * weighting_factor_adj
        # mass_adj_sum = mass_adj.sum(dim="z")
        # start_mass = mass_adj_sum.isel(time=0).copy() * 0 + mass_adj_sum.max()
        # mass_loss_percent = (start_mass - mass_adj_sum) / start_mass * 100
        # data["ACTIVE FIRE"] = mass_loss_percent > 1

        # starting_mass = mass_xy.isel(time=0).copy() * 0 + 0.00356416
        # mass_loss_percent = (starting_mass - mass_xy) / starting_mass * 100

        # On average, at what time does the fireline reach y=0 to the left of the circle out to x=-10?
        circle_radius = sim_params["circle_radius"]
        left_points = np.arange(-10, -circle_radius - 0.25, 0.1)
        left_times = []
        left_indices = []
        for x in left_points:
            time, index = find_fire_intersection_time(data["ACTIVE FIRE"], x, 0)
            left_times.append(time)
            left_indices.append(index)
        left_time = float(np.mean(left_times))
        left_index = int(np.mean(left_indices) + 0.5)

        # On average, at what time does the fireline reach y=0 to the right of the circle out to x=10?
        right_points = np.arange(circle_radius + 0.25, 10, 0.1)
        right_times = []
        right_indices = []
        for x in right_points:
            time, index = find_fire_intersection_time(data["ACTIVE FIRE"], x, 0)
            right_times.append(time)
            right_indices.append(index)
        right_time = float(np.mean(right_times))
        right_index = int(np.mean(right_indices) + 0.5)

        if not left_time or not right_time:
            return sim_id, {"curvature": 0}

        # Get the firelines at the time range of interest
        # target_index = (left_index + right_index) // 2
        first_index = max(left_index, right_index)
        last_index = max(left_index, right_index)
        # fire_line = get_fire_line(
        #     data["ACTIVE FIRE"],
        #     first_index,
        #     last_index,
        #     -circle_radius,
        #     circle_radius,
        # )
        fire_line = get_fire_line_ideal(data["ACTIVE FIRE"], first_index, last_index, circle_radius,
                                        -circle_radius - circle_radius, circle_radius + circle_radius)

        # Get actual time values in seconds for plotting
        # Extract from data's time coordinate if available
        times_in_seconds = []
        for i in range(first_index, last_index + 1):
            time = data["time"].isel(time=i).values
            times_in_seconds.append(time)

        # plot_fire_front_evolution(fire_line, circle_radius, first_index)

        # Get the polynomial fit to the fireline
        polynomial = get_polynomial(fire_line)

        # plot_fire_front_with_polynomial_fit(fire_line, polynomial, circle_radius, first_index)

        # Create the combined visualization
        fig_combined = plot_combined_fire_front_analysis(
            fire_line=fire_line, polynomial_coeffs=polynomial, circle_radius=circle_radius, first_index=first_index,
            simulation_id=sim_id, times_in_seconds=times_in_seconds, wind_speed=sim_params["wind_speed"],
            treatment_height=sim_params["treatment_fuel_height"]
        )

        # Save the figure for later reference
        fig_combined.savefig(f"{sim_directory}/ideal_fire_front_analysis_{QUANTITY}.png", dpi=300, bbox_inches='tight')
        plt.close(fig_combined)  # Close the figure to free memory

        # Get the curvature of the fireline
        curvature = get_curvature(polynomial)
        avg_curvature = np.mean(curvature)

        scaled_curvature = avg_curvature * sim_params["circle_radius"] / R_C_MAX

        return sim_id, {"curvature": avg_curvature, "scaled_curvature": scaled_curvature}

    except Exception as e:
        print(f"Error processing simulation {sim_id}: {str(e)}")
        return sim_id, {"curvature": None}


def postprocess(experiment_directory: str | Path):
    """
    Parallel processing version of the postprocessing function.
    """
    if isinstance(experiment_directory, str):
        experiment_directory = Path(experiment_directory)
    if not experiment_directory.is_dir():
        raise ValueError(f"Experiment directory {experiment_directory} does not exist.")

    simulations_directory = experiment_directory / "simulations"
    if isinstance(simulations_directory, str):
        simulations_directory = Path(simulations_directory)
    if not simulations_directory.is_dir():
        raise ValueError(f"Outputs directory {simulations_directory} does not exist.")

    # Load the parameters.csv file into a dataframe
    inputs_df = pd.read_csv(experiment_directory / "parameters.csv")

    # Prepare arguments for parallel processing
    process_args = []
    for _, row in inputs_df.iterrows():
        sim_params = row.to_dict()
        sim_id = int(sim_params["simulation_id"])

        if TARGET_SIM_ID and sim_id != TARGET_SIM_ID:
            continue

        sim_directory = simulations_directory / f"simulation_{sim_id}"
        process_args.append((sim_directory, sim_params))

    # Use all available CPU cores except one
    print(f"Processing {len(process_args)} simulations using {NUM_WORKERS} processes...")

    # Process simulations in parallel
    outputs_dict = {}
    with Pool(processes=NUM_WORKERS) as pool:
        # Use tqdm to show progress
        results = list(tqdm(
            pool.imap(process_simulation, process_args),
            total=len(process_args),
            desc="Processing simulations"
        ))

    # Collect results
    for sim_id, sim_outputs in results:
        outputs_dict[sim_id] = sim_outputs

    # Create output dataframe
    outputs_df = pd.DataFrame.from_dict(outputs_dict, orient="index")

    # Merge with inputs_df using simulation_id as the index
    inputs_df.set_index("simulation_id", inplace=True)
    outputs_df.index.name = "simulation_id"
    final_df = pd.merge(inputs_df, outputs_df, left_index=True, right_index=True)

    # Save to CSV
    if not TARGET_SIM_ID:
        final_df.to_csv(experiment_directory / f"ideal_curvatures_{QUANTITY}.csv")
        print(f"Processing complete. Results saved to ideal_curvatures_{QUANTITY}.csv")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        postprocess(*sys.argv[1:])

    else:
        postprocess(
            "/Volumes/T7 Shield/bandy-circles/grid-search-combined",
        )
