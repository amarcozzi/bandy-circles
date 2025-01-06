import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_curvature_heatmap(radius, data_file='/Volumes/T7 Shield/bandy-circles/grid-search/curvatures.csv'):
    # Read and filter data
    df = pd.read_csv(data_file)
    df_filtered = df[df['circle_radius'] == radius].copy()
    
    # Create pivot table for heatmap
    pivot = df_filtered.pivot(
        index='treatment_fuel_height',
        columns='wind_speed',
        values='curvature'
    )
    
    # Create figure with larger size
    plt.figure(figsize=(12, 8))
    
    # Find symmetric limits for colorbar
    abs_max = max(abs(pivot.min().min()), abs(pivot.max().max()))
    
    # Create heatmap
    sns.heatmap(
        pivot.sort_index(ascending=False),  # Flip y-axis
        cmap='RdBu_r',
        center=0,
        vmin=-abs_max,
        vmax=abs_max,
        annot=True,
        fmt='.2f',
        cbar_kws={'label': 'Curvature'},
        xticklabels=2,
        yticklabels=2
    )
    
    # Customize plot
    plt.title(f'Curvature Heatmap (Circle Radius = {radius})m')
    plt.xlabel('Wind Speed (m/s)')
    plt.ylabel('Treatment Fuel Height (m)')
    
    plt.tight_layout()
    return plt.gcf()

# Example usage:
fig = plot_curvature_heatmap(1.8)  # For radius = 0.9
plt.show()