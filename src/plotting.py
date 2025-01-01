import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import animation
from typing import Optional, List, Tuple
from xarray.core.dataarray import DataArray
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.gridspec import GridSpec


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import animation
from typing import Optional, List, Tuple
from xarray.core.dataarray import DataArray
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.gridspec import GridSpec


def create_dataset_movie(
    dataset: xr.Dataset,
    quantity: str,
    output_file: str,
    fps: int = 10,
    dpi: int = 100,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    color_map: str = "inferno",
    start_time: Optional[float] = None,
    end_time: Optional[float] = None,
    z_level: Optional[float] = None,
    verbose: bool = False,
    circle_radius: Optional[float] = None,
):
    # Time selection
    if start_time is not None:
        dataset = dataset.sel(time=slice(start_time, None))
    if end_time is not None:
        dataset = dataset.sel(time=slice(None, end_time))

    # Data preparation
    if z_level is None:
        data = dataset.mean(dim="z").values
        z_label = "z-averaged"
    else:
        data = dataset.sel(z=z_level, method="nearest").values
        z_label = f"z = {z_level:.2f} m"

    # Set up the figure with better proportions and spacing
    fig = plt.figure(figsize=(16, 9))
    gs = GridSpec(1, 2, width_ratios=[6, 1], wspace=0.05)
    ax_main = fig.add_subplot(gs[0])
    ax_colorbar = fig.add_subplot(gs[1])

    # Determine value range with robust outlier handling
    if vmin is None:
        vmin = np.nanpercentile(data, 1)
    if vmax is None:
        vmax = np.nanpercentile(data, 99)

    # Create main plot with improved visualization
    im = ax_main.imshow(
        data[0],
        cmap=plt.get_cmap(color_map),
        vmin=vmin,
        vmax=vmax,
        extent=[dataset.x.min(), dataset.x.max(), dataset.y.min(), dataset.y.max()],
        animated=True,
        interpolation="bilinear",
    )

    # Add circle if radius is provided
    if circle_radius is not None:
        circle = plt.Circle(
            (0, 0), circle_radius, fill=False, color="red", linestyle="--", alpha=0.8
        )
        ax_main.add_artist(circle)

    # Enhanced colorbar
    cbar = fig.colorbar(im, cax=ax_colorbar)
    cbar.set_label(f"{quantity}\n({z_label})", fontsize=12, labelpad=10)
    cbar.ax.tick_params(labelsize=10)

    # Improved axes styling
    ax_main.set_xlabel("x [m]", fontsize=12, labelpad=8)
    ax_main.set_ylabel("y [m]", fontsize=12, labelpad=8)
    ax_main.tick_params(labelsize=10)
    ax_main.set_aspect("equal")
    ax_main.grid(True, linestyle="--", alpha=0.3)

    # Enhanced title with background
    title_text = f"t = {dataset.time[0].values:.2f} s"
    title = ax_main.set_title(
        title_text,
        pad=20,
        fontsize=14,
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
    )

    def update(frame):
        im.set_array(data[frame])
        title.set_text(f"t = {dataset.time[frame].values:.2f} s")
        return im, title

    # Create animation with progress bar
    if verbose:
        from tqdm import tqdm

        frames = tqdm(range(len(dataset.time)))
    else:
        frames = range(len(dataset.time))

    anim = animation.FuncAnimation(
        fig, update, frames=frames, interval=1000 / fps, blit=True
    )

    if verbose:
        print(f"Saving movie to {output_file}...")

    # Create writer with all parameters specified at creation time
    writer = animation.FFMpegWriter(
        fps=fps,
        metadata=dict(title=f"{quantity} Evolution", artist="Scientific Visualization"),
        extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"],
    )

    # Save with only the writer and dpi specified
    anim.save(output_file, writer=writer, dpi=dpi)

    plt.close(fig)

    if verbose:
        print(f"Movie saved successfully as {output_file}")


def plot_fire_lines_with_curvature(
    fire_line: List[List[Tuple[float, float]]],
    curvatures: Tuple[float, ...],
    circle_radius: float,
):
    """
    Plot all fire lines with their color representing the curvature at that index.

    Parameters:
    -----------
    fire_line : List[List[Tuple[float, float]]]
        A list of fire line coordinates for each time step.
    curvatures : Tuple[float, ...]
        A tuple of curvature values corresponding to each fire line.
    circle_radius : float
        The radius of the circle centered at (0, 0).
    """
    plt.figure(figsize=(12, 8))

    # Normalize curvatures for colormap
    norm = Normalize(min(curvatures), max(curvatures))

    for coords, curvature in zip(fire_line, curvatures):
        x, y = zip(*coords)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, cmap="viridis", norm=norm)
        lc.set_array(np.full(len(segments), curvature))
        line = plt.gca().add_collection(lc)

    plt.colorbar(line, label="Curvature")

    # Plot the circle
    circle = plt.Circle((0, 0), circle_radius, fill=False, color="r")
    plt.gca().add_artist(circle)

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Fire Lines Colored by Curvature")
    plt.axis("equal")
    plt.grid(True)

    # Set the limits of the plot
    all_x = [x for line in fire_line for x, _ in line]
    all_y = [y for line in fire_line for _, y in line]
    plt.xlim(min(all_x), max(all_x))
    plt.ylim(min(all_y), max(all_y))

    plt.show()


def plot_min_curvature_fire_line(
    fire_line: List[List[Tuple[float, float]]],
    curvatures: Tuple[float, ...],
    active_fire_data: DataArray,
    circle_radius: float,
):
    """
    Plot the fire line with minimum curvature over the ACTIVE FIRE grid at that time.

    Parameters:
    -----------
    fire_line : List[List[Tuple[float, float]]]
        A list of fire line coordinates for each time step.
    curvatures : Tuple[float, ...]
        A tuple of curvature values corresponding to each fire line.
    active_fire_data : DataArray
        The ACTIVE FIRE data array.
    circle_radius : float
        The radius of the circle centered at (0, 0).
    """
    min_curvature_index = np.argmax(curvatures)
    min_curvature_fire_line = fire_line[min_curvature_index]

    plt.figure(figsize=(12, 8))

    # Plot ACTIVE FIRE grid
    active_fire_at_time = active_fire_data.isel(time=min_curvature_index)
    plt.imshow(
        active_fire_at_time,
        extent=[
            active_fire_data.x.min(),
            active_fire_data.x.max(),
            active_fire_data.y.min(),
            active_fire_data.y.max(),
        ],
        cmap="Reds",
        alpha=0.7,
    )

    # Plot fire line
    x, y = zip(*min_curvature_fire_line)
    plt.plot(x, y, color="blue", linewidth=2, label="Fire Line")

    # Plot the circle
    circle = plt.Circle((0, 0), circle_radius, fill=False, color="r")
    plt.gca().add_artist(circle)

    plt.colorbar(label="ACTIVE FIRE")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"Fire Line at Minimum Curvature (Time Step: {min_curvature_index})")
    plt.legend()
    plt.axis("equal")
    plt.grid(True)
    plt.show()


def plot_polynomials(
    polynomial_coeffs: List[List[float]],
    times: List[float],
    start_index: int,
    end_index: int,
    circle_radius: float,
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
    start_index : int
        Starting index for the time range to consider
    end_index : int
        Ending index for the time range to consider
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
        (0, 0), circle_radius, fill=False, color="black", linestyle="--", alpha=0.5
    )
    ax.add_artist(circle)

    # Get curvatures (2 * quadratic coefficient) for all polynomials
    curvatures = [2 * coeffs[0] for coeffs in polynomial_coeffs]

    # Create a symmetric normalization centered at 0
    max_abs_curvature = max(abs(min(curvatures)), abs(max(curvatures)))
    norm = Normalize(-max_abs_curvature, max_abs_curvature)

    # Create custom colormap with light gray center
    colors = [
        (0.8, 0.2, 0.2),  # red
        (0.85, 0.85, 0.85),  # light gray
        (0.2, 0.2, 0.8),
    ]  # blue
    n_bins = 256
    custom_cmap = LinearSegmentedColormap.from_list("custom", colors, N=n_bins)

    # Plot each polynomial
    for i, coeffs in enumerate(polynomial_coeffs):
        # Evaluate polynomial
        y = coeffs[0] * x**2 + coeffs[1] * x + coeffs[2]

        # Plot with curvature-based color
        line = ax.plot(x, y, alpha=0.7, linewidth=2)
        plt.setp(line, color=custom_cmap(norm(curvatures[i])))

    # Add colorbar with centered ticks
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    cbar = plt.colorbar(
        sm, label="Curvature (2nd derivative)", cax=fig.add_axes([0.92, 0.1, 0.02, 0.8])
    )

    # Set colorbar ticks to show the symmetry
    tick_locations = np.linspace(-max_abs_curvature, max_abs_curvature, 5)
    cbar.set_ticks(tick_locations)
    cbar.set_ticklabels([f"{val:.2f}" for val in tick_locations])

    # Set axis labels and title
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    plot_time_range = f"{times[start_index]:.1f}s to {times[end_index]:.1f}s"
    ax.set_title(f"Fire Front Evolution with Curvature Coloring\n{plot_time_range}")

    # Set equal aspect ratio
    ax.set_aspect("equal")

    # Add grid
    ax.grid(True, alpha=0.3)

    # Set limits with some padding
    padding = 1
    ax.set_xlim(-circle_radius - padding, circle_radius + padding)
    ax.set_ylim(-circle_radius - padding, circle_radius + padding)

    if output_path:
        plt.savefig(output_path)
        plt.close()
    else:
        plt.show()
