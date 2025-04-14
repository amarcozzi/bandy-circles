import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path


def plot_curvature_heatmap(df, radius, ax=None, x_max=None, y_max=None):
    """
    Plot a curvature heatmap for a specific radius.

    Parameters:
    -----------
    df : DataFrame
        The dataset containing curvature data
    radius : float
        The circle radius to filter data for
    ax : matplotlib Axes, optional
        The axis to plot on. If None, uses the current axis
    x_max : float, optional
        Maximum wind speed to include in the plot
    y_max : float, optional
        Maximum treatment fuel height to include in the plot

    Returns:
    --------
    fig : Figure
        The figure object containing the heatmap
    """
    # Read and filter data
    df_filtered = df[df["circle_radius"] == radius].copy()

    # Apply x_max and y_max filters if provided
    if x_max is not None:
        df_filtered = df_filtered[df_filtered["wind_speed"] <= x_max]
    if y_max is not None:
        df_filtered = df_filtered[df_filtered["treatment_fuel_height"] <= y_max]

    # Create pivot table for heatmap
    pivot = df_filtered.pivot(
        index="treatment_fuel_height", columns="wind_speed", values="curvature"
    )

    # Create figure if axis not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    else:
        fig = ax.figure

    # Create heatmap
    sns.heatmap(
        pivot.sort_index(ascending=False),
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Curvature"},
        xticklabels=2,
        yticklabels=2,
        ax=ax
    )

    # Customize plot
    ax.set_title(f"Curvature Heatmap (Circle Radius = {radius}m)")
    ax.set_xlabel("Wind Speed (m/s)")
    ax.set_ylabel("Treatment Fuel Height (m)")

    return fig


def combined_heatmap_figure(df, radii, x_max=None, y_max=None):
    """
    Create a combined figure with multiple heatmaps as subplots.

    Parameters:
    -----------
    df : DataFrame
        The dataset containing curvature data
    radii : list
        List of circle radii to create heatmaps for
    x_max : float, optional
        Maximum wind speed to include in the plots
    y_max : float, optional
        Maximum treatment fuel height to include in the plots

    Returns:
    --------
    fig : Figure
        The combined figure with all heatmaps
    individual_figs : list
        List of individual figure objects for each heatmap
    """
    # Calculate grid dimensions
    n = len(radii)
    ncols = min(3, n)  # Max 3 columns
    nrows = (n + ncols - 1) // ncols  # Ceiling division

    # Create combined figure
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))

    # Flatten axes array for easy indexing
    if n > 1:
        axes = axes.flatten()
    else:
        axes = [axes]

    # Hide unused axes
    for i in range(n, len(axes)):
        axes[i].set_visible(False)

    # Create individual figures for each heatmap
    individual_figs = []

    # Plot heatmaps
    for i, radius in enumerate(radii):
        # Create individual figure for this radius
        ind_fig = plot_curvature_heatmap(df, radius, x_max=x_max, y_max=y_max)
        individual_figs.append(ind_fig)

        # Plot on combined figure
        plot_curvature_heatmap(df, radius, ax=axes[i], x_max=x_max, y_max=y_max)

    plt.tight_layout()
    return fig, individual_figs


def main():
    data_path = Path("/Volumes/T7 Shield/bandy-circles/grid-search-high")
    quantity = "heat_flux"
    df = pd.read_csv(data_path / f"curvatures_{quantity}.csv")

    # Define radii
    radii = [0.9, 1.35, 1.8, 2.25, 2.7]

    # Optional: Set max values for x and y axes
    x_max = None  # Maximum wind speed to include (e.g., 5.0)
    y_max = None  # Maximum treatment fuel height to include (e.g., 1.5)

    # Create combined figure and individual figures
    combined_fig, individual_figs = combined_heatmap_figure(df, radii, x_max=x_max, y_max=y_max)

    # Save combined figure
    combined_fig.savefig(data_path / f"combined_curvature_heatmaps_{quantity}.png", dpi=300, bbox_inches='tight')

    # Save individual figures if needed
    for i, (radius, fig) in enumerate(zip(radii, individual_figs)):
        fig.savefig(data_path / f"curvature_heatmap_{radius}_{quantity}.png", dpi=300, bbox_inches='tight')
        plt.close(fig)  # Close individual figure

    # Show combined figure
    plt.show()


if __name__ == "__main__":
    main()
