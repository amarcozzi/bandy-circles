import pandas as pd
import numpy as np
from pathlib import Path
import plotly.graph_objects as go
from scipy.interpolate import griddata


def plot_curvature_isosurface(df, threshold=0, interpolate=True):
    """
    Plot a 3D isosurface showing regions where curvature meets a threshold.

    Parameters:
    -----------
    df : DataFrame
        The dataset containing curvature data
    threshold : float or list
        The curvature threshold(s) to visualize
    interpolate : bool
        Whether to interpolate the data to a finer grid
    """
    # Extract unique values for each dimension
    wind_speeds = sorted(df["wind_speed"].unique())
    treatment_heights = sorted(df["treatment_fuel_height"].unique())[:-1]  # Exclude last value
    circle_radii = sorted(df["circle_radius"].unique())

    # Create 3D grid for visualization
    if interpolate:
        # Create fine regular grid
        x_fine = np.linspace(min(wind_speeds), max(wind_speeds), 50)
        y_fine = np.linspace(min(treatment_heights), max(treatment_heights), 50)
        z_fine = np.linspace(min(circle_radii), max(circle_radii), 50)
        X, Y, Z = np.meshgrid(x_fine, y_fine, z_fine)

        # Interpolate curvature values onto fine grid
        points = df[['wind_speed', 'treatment_fuel_height', 'circle_radius']].values
        values = df['curvature'].values
        curvature_interp = griddata(points, values, (X, Y, Z), method='linear')
    else:
        # Use original data points
        X, Y, Z = np.meshgrid(wind_speeds, treatment_heights, circle_radii)

        # Create a 3D grid of curvature values
        curvature_grid = np.zeros(X.shape)
        for i, ws in enumerate(wind_speeds):
            for j, th in enumerate(treatment_heights):
                for k, cr in enumerate(circle_radii):
                    mask = (df['wind_speed'] == ws) & \
                           (df['treatment_fuel_height'] == th) & \
                           (df['circle_radius'] == cr)
                    if mask.any():
                        curvature_grid[j, i, k] = df.loc[mask, 'curvature'].values[0]
        curvature_interp = curvature_grid

    # Create figure
    fig = go.Figure()

    # Check if threshold is a list or single value
    if not isinstance(threshold, list):
        threshold = [threshold]

    # Create a different colored isosurface for each threshold
    colors = ['red', 'blue', 'green', 'purple', 'orange']

    for i, thresh in enumerate(threshold):
        color = colors[i % len(colors)]

        # Add isosurface
        fig.add_trace(go.Isosurface(
            x=X.flatten(),
            y=Y.flatten(),
            z=Z.flatten(),
            value=curvature_interp.flatten(),
            isomin=thresh - 0.001,
            isomax=thresh + 0.001,
            opacity=0.7,
            surface_count=1,
            colorscale=[[0, color], [1, color]],
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False)
        ))

    # Configure the layout
    fig.update_layout(
        title=f"Curvature Isosurface(s) at threshold(s): {threshold}",
        scene=dict(
            xaxis_title="Wind Speed (m/s)",
            yaxis_title="Treatment Fuel Height (m)",
            zaxis_title="Circle Radius (m)",
            aspectmode='cube'
        ),
        width=900,
        height=800,
        margin=dict(l=65, r=50, b=65, t=90)
    )

    return fig


# Example usage
data_path = Path("/Volumes/T7 Shield/bandy-circles/grid-search-high-wind")
quantity = "heat_flux"
df = pd.read_csv(data_path / f"curvatures_ideal_{quantity}.csv")

# Scale curvature
line_length = df["circle_radius"] / df["circle_radius"].max()
df["curvature"] = df["curvature"] * line_length

# Plot zero curvature isosurface (transition between concave and convex)
fig = plot_curvature_isosurface(df, threshold=0)
fig.show()

# Plot multiple isosurfaces for different curvature thresholds
fig_multi = plot_curvature_isosurface(df, threshold=[-0.2, 0, 0.2])
fig_multi.show()
