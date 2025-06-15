import numpy as np
import pandas as pd
from pathlib import Path
from postprocess import (
    get_particle_data_array,
    get_polynomial,
    get_fire_line,
    find_fire_intersection_time,
)
import matplotlib.pyplot as plt
from typing import Optional, List
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap

plt.style.use(['science', 'bright'])
plt.rcParams['font.size'] = '24'

# TARGET_SIM_ID = 480  # Left figure
TARGET_SIM_ID = 484  # Right figure


def plot_polynomials(
        polynomial_coeffs: List[List[float]],
        times: List[float],
        circle_radius: float,
        evaluation_index=None,
        output_path: Optional[str] = None,
):
    """
    Plot polynomial curves colored by their curvature (2nd derivative).
    Colors transition from red (negative) through light gray (zero) to blue (positive).

    Parameters
    ----------
    polynomial_coeffs : List[List[float]]
        List of polynomial coefficients for each time step [a, b, c] where y = ax^2 + bx + c
    times : List[float]
        List of times corresponding to each polynomial. Only times[start_index:end_index+1]
        will be considered for the plot title and time-based analysis
    circle_radius : float
        Radius of the circle to plot for reference
    output_path : Optional[str]
        Path to save the plot. If None, displays the plot instead.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create x values for polynomial evaluation
    x = np.linspace(-circle_radius - 1, circle_radius + 1, 200)

    # Plot reference circle
    circle = plt.Circle(
        (0, 0), circle_radius, fill=False, color="black", linestyle="--", alpha=0.5, linewidth=1.5
    )
    ax.add_artist(circle)

    # Get curvatures (2 * quadratic coefficient) for all polynomials
    curvatures = [2 * circle_radius * coeffs[0] for coeffs in polynomial_coeffs]

    # Create a symmetric normalization centered at 0
    # max_abs_curvature = max(abs(min(curvatures)), abs(max(curvatures)))
    max_abs_curvature = 0.7  # Set a fixed range for curvature normalization
    norm = Normalize(-max_abs_curvature, max_abs_curvature)

    # min_curvature = min(curvatures)
    # norm = Normalize(min_curvature, 0)

    # Create custom colormap with light gray center
    colors = [
        (0.2, 0.2, 0.8),  # blue
        (0.85, 0.85, 0.85),  # light gray
        (0.8, 0.2, 0.2),  # red
    ]
    n_bins = 256
    custom_cmap = LinearSegmentedColormap.from_list("custom", colors, N=n_bins)

    # Plot each polynomial
    for i, coeffs in enumerate(polynomial_coeffs):
        # Evaluate polynomial
        y = coeffs[0] * x ** 2 + coeffs[1] * x + coeffs[2]

        # Plot with curvature-based color
        if evaluation_index and i == evaluation_index:
            line = ax.plot(x, y, color="green", alpha=0.7, linewidth=10)
        else:
            line = ax.plot(x, y, alpha=0.7, linewidth=3)
            plt.setp(line, color=custom_cmap(norm(curvatures[i])))

    # Add colorbar with centered ticks
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    cbar = plt.colorbar(
        sm, label="$\kappa^*$", ax=ax
    )

    # Set colorbar ticks to show the symmetry
    tick_locations = np.linspace(-max_abs_curvature, max_abs_curvature, 5)
    cbar.set_ticks(tick_locations)
    cbar.set_ticklabels([f"{val:.2f}" for val in tick_locations])

    # Set axis labels and title
    ax.set_xlim(-3, 3)
    ax.set_xlabel("x (m)")
    ax.set_ylim(-3, 3)
    ax.set_ylabel("y (m)")

    # Set equal aspect ratio
    ax.set_aspect("equal")

    # Add grid
    ax.grid(True, alpha=0.3)

    # Set limits with some padding
    padding = 1
    ax.set_xlim(-circle_radius - padding, circle_radius + padding)
    ax.set_ylim(-circle_radius - padding, circle_radius + padding)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=1600)
        plt.close()
    else:
        plt.show()


def main():
    experiment_dir = Path("/Volumes/T7 Shield/bandy-circles/grid-search-combined")
    df = pd.read_csv(experiment_dir / "parameters.csv")
    params = df[df["simulation_id"] == TARGET_SIM_ID]

    sim_path = experiment_dir / "simulations" / f"simulation_{TARGET_SIM_ID}"
    data = get_particle_data_array(sim_path)

    data["ACTIVE FIRE"] = (data["PARTICLE TOTAL HEAT FLUX"].mean("z") < -20).compute()

    # On average, at what time does the fireline reach y=0 to the left of the circle out to x=-5?
    circle_radius = params["circle_radius"].values[0]
    left_points = np.arange(-10, -circle_radius - 0.25, 0.1)
    left_times = []
    left_indices = []
    for x in left_points:
        time, index = find_fire_intersection_time(data["ACTIVE FIRE"], x, 0)
        left_times.append(time)
        left_indices.append(index)
    left_time = float(np.median(left_times))
    left_index = int(np.median(left_indices) + 0.5)

    # On average, at what time does the fireline reach y=0 to the right of the circle out to x=5?
    right_points = np.arange(circle_radius + 0.25, 10, 0.1)
    right_times = []
    right_indices = []
    for x in right_points:
        time, index = find_fire_intersection_time(data["ACTIVE FIRE"], x, 0)
        right_times.append(time)
        right_indices.append(index)
    right_time = float(np.median(right_times))
    right_index = int(np.median(right_indices) + 0.5)

    if not left_time or not right_time:
        return TARGET_SIM_ID, {"curvature": 0}

    # Get the firelines at the time range of interest
    # target_index = (left_index + right_index) // 2
    first_index = max(left_index, right_index) - 50
    last_index = max(left_index, right_index) + 50
    fire_line = get_fire_line(
        data["ACTIVE FIRE"],
        first_index,
        last_index,
        -circle_radius - 2,
        circle_radius + 2,
    )

    # Get actual time values in seconds for plotting
    # Extract from data's time coordinate if available
    times_in_seconds = []
    evaluation_index = 0
    j = 0
    for i in range(first_index, last_index):
        time = data["time"].isel(time=i).values
        times_in_seconds.append(time)
        if i == max(left_index, right_index):
            evaluation_index = j
        j += 1

    polynomial = get_polynomial(fire_line)

    n = 4  # Sample every n time steps to reduce the number of polynomials plotted
    polynomial = polynomial[::n]
    times_in_seconds = times_in_seconds[::n]
    plot_polynomials(polynomial,
                     times=times_in_seconds,
                     circle_radius=circle_radius,
                     evaluation_index=50 // n,
                     # output_path=None
                     output_path=sim_path / f"sim_{TARGET_SIM_ID}_polynomials.pdf"
                     )


if __name__ == "__main__":
    main()
