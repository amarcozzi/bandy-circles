import numpy as np
import pandas as pd
from pathlib import Path
from postprocess import (
    get_particle_data_array,
    get_polynomial,
    get_active_fire_array,
    get_fire_line,
    find_fire_intersection_time,
)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as pe

TARGET_SIM_ID = 122


def plot_fire_front_evolution(fire_line, circle_radius=None, simulation_id=None,
                              times_in_seconds=None, x_extent=4, colormap='plasma',
                              show_plot=True):
    """
    Plot the evolution of a fire front over time with cleaner isochrones.

    Parameters:
    -----------
    fire_line : list of arrays
        List where each element is an array of (x, y) coordinates for a timestep
    circle_radius : float, optional
        Radius of the circular obstacle, if applicable
    simulation_id : int, optional
        Simulation ID to include in title if provided
    times_in_seconds : list, optional
        Actual time values in seconds corresponding to each timestep
    x_extent : float, optional
        The x-axis extent for the plot
    colormap : str, optional
        The colormap to use for the visualization
    show_plot : bool, optional
        Whether to call plt.show() at the end
    """
    # Create figure and axis with a specific figure size for better quality
    fig, ax = plt.subplots(figsize=(10, 9))

    # Create a custom colormap that transitions more visibly
    cmap = cm.get_cmap(colormap, len(fire_line))

    # For legend
    legend_handles = []
    legend_labels = []

    # Plot each timestep with a different color
    for i, timestep in enumerate(fire_line):
        # Skip if there are no points
        if not timestep:
            continue

        # Get color for this timestep
        color = cmap(i / max(1, len(fire_line) - 1))

        # Extract x and y coordinates
        x_points = [point[0] for point in timestep]
        y_points = [point[1] for point in timestep]

        # Determine time label based on whether seconds are provided
        if times_in_seconds and i < len(times_in_seconds):
            time_label = f'{times_in_seconds[i]:.1f}s'
        else:
            time_label = f't={i}'

        # Sort points by x-coordinate for smoother lines
        sorted_indices = np.argsort(x_points)
        x_sorted = [x_points[j] for j in sorted_indices]
        y_sorted = [y_points[j] for j in sorted_indices]

        # Plot the isochrone as a line with enhanced styling
        line = ax.plot(x_sorted, y_sorted, '-', color=color, linewidth=2.5,
                       alpha=0.9, solid_capstyle='round',
                       path_effects=[pe.Stroke(linewidth=3.5, foreground='black', alpha=0.2),
                                     pe.Normal()])[0]

        # Plot points with small markers
        scatter = ax.scatter(x_points, y_points, color=color, s=20, alpha=0.7,
                             edgecolor='white', linewidth=0.5)

        # Add to legend
        legend_handles.append(line)
        legend_labels.append(time_label)

    # Add the circular obstacle if radius is provided
    if circle_radius:
        circle = plt.Circle((0, 0), circle_radius, fill=True,
                            color='gray', alpha=0.3, edgecolor='black', linewidth=1)
        ax.add_patch(circle)

        # Create a circle for the legend
        from matplotlib.patches import Circle
        circle_legend = Circle((0, 0), 5, color='gray', alpha=0.3, edgecolor='black')
        legend_handles.append(circle_legend)
        legend_labels.append('Treatment')

    # Set axis limits to the specified extent
    ax.set_xlim(-x_extent, x_extent)
    ax.set_ylim(-4.5, 4.5)  # Keeping y-axis similar to original

    # Add labels and title with enhanced styling
    ax.set_xlabel('X (m)', fontsize=12)
    ax.set_ylabel('Y (m)', fontsize=12)

    # title = 'Evolution of Fire Front Over Time'
    # if simulation_id is not None:
    #     title += f' - Simulation {simulation_id}'
    # ax.set_title(title, fontsize=14, fontweight='bold', pad=15)

    # Add a more visually appealing legend to the right of the plot
    ax.legend(legend_handles, legend_labels, loc='upper left',
              bbox_to_anchor=(1., 1), fontsize=11, title_fontsize=11,
              framealpha=0.7, edgecolor='gray')

    # Add grid for better reference - using a more subtle styling
    ax.grid(True, linestyle='--', alpha=0.4, color='gray')

    # Set equal aspect ratio
    ax.set_aspect('equal')

    # # Set background color to enhance visibility
    # ax.set_facecolor('#f8f8f8')

    # Add a box around the plot
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.2)
        spine.set_color('gray')

    fig.tight_layout()

    if show_plot:
        plt.show()

    return fig, ax


def main():
    experiment_dir = Path("/Volumes/T7 Shield/bandy-circles/grid-search")
    df = pd.read_csv(experiment_dir / "parameters.csv")
    params = df[df["simulation_id"] == TARGET_SIM_ID]

    sim_path = experiment_dir / "simulations" / f"simulation_{TARGET_SIM_ID}"
    data = get_particle_data_array(sim_path)

    # rolling_average = data["PARTICLE TOTAL HEAT FLUX"].rolling(time=10, center=True).mean("z").compute()
    # data["ACTIVE FIRE"] = rolling_average < -50
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
    first_index = max(left_index, right_index) - 80
    last_index = max(left_index, right_index) + 90
    fire_line = get_fire_line(
        data["ACTIVE FIRE"],
        first_index,
        last_index,
        -circle_radius,
        circle_radius,
    )
    fire_line_extended = get_fire_line(
        data["ACTIVE FIRE"],
        first_index,
        last_index,
        -circle_radius - 8,
        circle_radius + 8,
    )

    # Get actual time values in seconds for plotting
    # Extract from data's time coordinate if available
    times_in_seconds = []
    for i in range(first_index, last_index + 1):
        time = data["time"].isel(time=i).values
        times_in_seconds.append(time)

    # Skip every n time steps for plotting
    n = 18
    fire_line = fire_line[::n]
    fire_line_extended = fire_line_extended[::n]
    times_in_seconds = times_in_seconds[::n]

    plot_fire_front_evolution(fire_line_extended, circle_radius, first_index, x_extent=4,
                              times_in_seconds=times_in_seconds)

    # polynomial = get_polynomial(fire_line)
    # plot_fire_front_with_polynomial_fit(fire_line, polynomial, circle_radius, first_index,
    #                                     times_in_seconds=times_in_seconds)

    print()


# plot_fire_front_evolution()


if __name__ == "__main__":
    main()
