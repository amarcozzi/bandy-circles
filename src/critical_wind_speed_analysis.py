"""
Critical Wind Speed Analysis: Heat Transfer Model vs FDS Simulations

This script compares the critical wind speed predictions from a simple 2D heat transfer model
with CFD simulation results from FDS wildfire simulations. The analysis examines the formation
of junction fires as an interaction with circular fuel treatments.

The heat transfer model considers:
- Radiative heat transfer between fuel beds
- Convective heat transfer due to wind
- Different view factor formulations

The FDS simulations provide:
- Curvature metrics (positive = junction fire, negative = accelerating fire)
- Critical wind speed where curvature transitions (< 0.1)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import pint
from scipy.optimize import fsolve, brentq
from scipy.interpolate import interp1d
import seaborn as sns

# Set up plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# ============================================================================
# HEAT TRANSFER MODEL FUNCTIONS (from first script)
# ============================================================================

# Initialize unit registry
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# Physical Constants
SIGMA = Q_(5.67e-8, 'watt / (meter**2 * kelvin**4)')  # Stefan-Boltzmann constant

# Temperature conditions
T_UPSTREAM = Q_(800 + 273.15, 'kelvin')  # Flame temperature (800°C)
T_COLD = Q_(295, 'kelvin')  # Ambient temperature (22°C)
T_IGNITION = Q_(350 + 273.15, 'kelvin')  # Ignition temperature (350°C)

# Geometric parameters
H = Q_(0.6, 'meter')  # Height of fuel lines (60 cm)
F_CLOSE = 1.0  # View factor for closely spaced bodies

# Air properties (at average temperature ~500°C)
RHO = Q_(0.5, 'kg/m**3')  # Density
K_FLUID = Q_(55.79, 'mW/(m*K)')  # Thermal conductivity
MU = Q_(3.5e-5, 'Pa*s')  # Dynamic viscosity
CP = Q_(1.092, 'kJ/(kg*K)')  # Specific heat at constant pressure
epsilon = 0.96  # Emissivity of the surface

# Calculate Prandtl number
PR = (MU * CP / K_FLUID).to_base_units()


def view_factor_parallel_plates_perpendicular(h, d):
    """View factor for parallel plates with perpendicular midlines."""
    h_val = h.to('meter').magnitude
    d_val = d.to('meter').magnitude

    if h_val == 0:
        return 0

    ratio = d_val / h_val
    F = np.sqrt(1 + ratio ** 2) - ratio
    return F


def view_factor_small_area_approximation(h, d):
    """Small area approximation for view factor."""
    h_val = h.to('meter').magnitude
    d_val = d.to('meter').magnitude

    if d_val == 0:
        return np.inf

    return h_val / (np.pi * d_val)


def view_factor_infinite_strips(h, d):
    """View factor between infinite parallel strips."""
    h_val = h.to('meter').magnitude
    d_val = d.to('meter').magnitude

    if d_val == 0:
        return 1.0

    F = (1 / np.pi) * np.arctan(h_val / d_val)
    return F


def calculate_h(velocity, D):
    """Calculate convective heat transfer coefficient."""
    # Calculate Reynolds number
    Re = (RHO * velocity * D / MU).to_base_units()

    # Calculate Nusselt number based on flow regime
    if Re.magnitude < 5e4:
        # Laminar flow
        Nu = 0.664 * Re.magnitude ** 0.5 * PR.magnitude ** (1 / 3)
    else:
        # Turbulent flow
        Nu = 0.037 * Re.magnitude ** 0.8 * PR.magnitude ** (1 / 3)

    # Calculate heat transfer coefficient
    h = Nu * K_FLUID / D
    return h


def energy_balance_equation(T_d_kelvin, velocity, D, view_factor_func):
    """Energy balance equation to solve for downstream temperature."""
    # Convert scalar to pint quantity
    T_d = Q_(T_d_kelvin, 'kelvin')

    # Calculate view factor
    F_12 = view_factor_func(H, D)

    # Calculate heat transfer coefficient
    h = calculate_h(velocity, D)

    # Convective heat input
    q_conv_in = h * (T_UPSTREAM - T_d)

    # Radiative heat input from upstream
    q_rad_in = epsilon * SIGMA * F_12 * (T_UPSTREAM ** 4 - T_d ** 4)

    # Radiative heat output to cold body
    q_rad_out = epsilon * SIGMA * F_CLOSE * (T_d ** 4 - T_COLD ** 4)

    # Energy balance: input = output
    residual = (q_conv_in + q_rad_in - q_rad_out).to('W/m**2').magnitude

    return residual


def solve_for_td(velocity, D, view_factor_func):
    """Solve for downstream temperature given velocity and separation distance."""
    # Initial guess (in Kelvin)
    T_d_initial = 500

    # Solve using scipy's fsolve
    T_d_solution = fsolve(energy_balance_equation, T_d_initial,
                          args=(velocity, D, view_factor_func))[0]

    return Q_(T_d_solution, 'kelvin')


def find_critical_velocity(D, view_factor_func):
    """Find critical velocity where downstream temperature reaches ignition."""

    def objective(v):
        velocity = Q_(v, 'm/s')
        T_d = solve_for_td(velocity, D, view_factor_func)
        return (T_d - T_IGNITION).magnitude

    try:
        # Use scipy to find the root
        v_critical = brentq(objective, 0.01, 50)  # Search between 0.01 and 50 m/s
        return Q_(v_critical, 'm/s')
    except ValueError:
        # No solution found in the range
        return None


# ============================================================================
# HEAT TRANSFER MODEL ANALYSIS
# ============================================================================

def run_heat_transfer_model(D_range):
    """Run the heat transfer model over a range of distances using parallel plates view factor."""

    # Use only parallel plates view factor
    vf_func = view_factor_parallel_plates_perpendicular

    distances = []
    critical_velocities = []
    view_factors = []

    for d in D_range:
        D = Q_(d, 'meter')
        v_crit = find_critical_velocity(D, vf_func)

        if v_crit is not None:
            distances.append(d)
            critical_velocities.append(v_crit.magnitude)
            view_factors.append(vf_func(H, D))

    results = {
        'distance': np.array(distances),
        'critical_velocity': np.array(critical_velocities),
        'view_factor': np.array(view_factors)
    }

    return results


# ============================================================================
# FDS SIMULATION DATA ANALYSIS
# ============================================================================

def load_fds_data(data_path):
    """Load FDS simulation results from CSV file using the same approach as heatmap script."""
    try:
        df = pd.read_csv(data_path)

        # Apply the same filtering as in the heatmap script
        df = df[df["wind_speed"] >= 1.5]
        df = df[df["wind_speed"] < 3.2]

        # Apply the same scaling as in the heatmap script
        line_length = df["circle_radius"]
        df["curvature"] = df["curvature"] * line_length

        return df
    except FileNotFoundError:
        print(f"Warning: FDS data file not found at {data_path}")
        return None


def find_critical_wind_speed(df, circle_radius, epsilon=-0.2):
    """
    Find critical wind speed where curvature transitions below epsilon.

    Parameters:
    df: DataFrame with simulation results
    circle_radius: Radius of the circular treatment
    epsilon: Curvature threshold (default 0.1)

    Returns:
    Critical wind speed (m/s) or None if not found
    """
    # Filter for specific circle radius and fuel height = 0.6
    subset = df[(df['circle_radius'] == circle_radius) &
                (df['treatment_fuel_height'] == 0.1)].copy()

    if len(subset) == 0:
        return None

    # Sort by wind speed
    subset = subset.sort_values('wind_speed')

    # Scale curvature by circle radius
    subset['curvature'] = subset['curvature'] * circle_radius

    # Find where curvature drops below epsilon
    below_threshold = subset[subset['curvature'] < epsilon]

    if len(below_threshold) == 0:
        return None

    # Return the minimum wind speed where curvature < epsilon
    return below_threshold['wind_speed'].min()


def process_fds_data(df):
    """Process FDS data to extract critical wind speeds for different circle radii."""
    if df is None:
        return None

    # Use the same radii as defined in the heatmap script
    radii = [0.9, 1.35, 1.8, 2.25, 2.7]

    fds_results = {
        'circle_radius': [],
        'diameter': [],
        'critical_wind_speed': []
    }

    for radius in radii:
        critical_wind = find_critical_wind_speed(df, radius)

        if critical_wind is not None:
            fds_results['circle_radius'].append(radius)
            fds_results['diameter'].append(2 * radius)  # Convert to diameter
            fds_results['critical_wind_speed'].append(critical_wind)
            print(f"Radius {radius} m: Critical wind speed = {critical_wind:.2f} m/s")
        else:
            print(f"Radius {radius} m: No critical wind speed found")

    return fds_results


# ============================================================================
# PLOTTING AND ANALYSIS
# ============================================================================

def create_comparison_plot(heat_transfer_results, fds_results):
    """Create comparison plot of heat transfer model vs FDS simulations."""

    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot: Critical wind speed vs distance
    ax.plot(heat_transfer_results['distance'], heat_transfer_results['critical_velocity'],
            '-o', color='#1f77b4', markersize=6, linewidth=2,
            label='Heat Transfer Model')

    # Add FDS simulation data if available
    if fds_results is not None:
        # ax.plot(fds_results['diameter'], fds_results['critical_wind_speed'],
        #         's', color='red', markersize=10, linewidth=3,
        #         label='FDS Simulations', markerfacecolor='none',
        #         markeredgewidth=2)
        ax.plot(fds_results['diameter'], fds_results['critical_wind_speed'],
                'X', color='red', markersize=10,  # linewidth=3,
                label='FDS Simulations')

    ax.set_xlabel('Separation Distance D (m)', fontsize=12)
    ax.set_ylabel('Critical Wind Speed (m/s)', fontsize=12)
    ax.set_title('Critical Wind Speed for Ignition vs Separation Distance', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    # ax.set_xlim(0, 12)
    # ax.set_ylim(0, 30)

    plt.tight_layout()

    # Add overall title
    fig.suptitle('Heat Transfer Model vs FDS Simulations: Junction Fire Analysis',
                 fontsize=16, y=1.02)

    return fig


def print_analysis_summary(heat_transfer_results, fds_results):
    """Print summary of analysis results."""
    print("=" * 80)
    print("CRITICAL WIND SPEED ANALYSIS SUMMARY")
    print("=" * 80)
    print(f"Fuel bed height: {H.magnitude} m")
    print(f"Ignition temperature: {T_IGNITION.to('degC').magnitude:.0f}°C")
    print(f"Upstream temperature: {T_UPSTREAM.to('degC').magnitude:.0f}°C")
    print()

    print("HEAT TRANSFER MODEL RESULTS (Parallel Plates View Factor):")
    print("-" * 60)
    for d, v in zip(heat_transfer_results['distance'], heat_transfer_results['critical_velocity']):
        print(f"  D = {d:.1f} m: u_critical = {v:.2f} m/s")

    if fds_results is not None:
        print(f"\nFDS SIMULATION RESULTS (H = 0.6m, scaled curvature):")
        print("-" * 60)
        for r, d, v in zip(fds_results['circle_radius'],
                           fds_results['diameter'],
                           fds_results['critical_wind_speed']):
            print(f"  Circle radius = {r:.2f} m (D = {d:.1f} m): u_critical = {v:.2f} m/s")

    print("\n" + "=" * 80)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""

    # Define distance range for heat transfer model
    D_range = np.linspace(2, 8, 50)

    print("Running heat transfer model analysis...")
    heat_transfer_results = run_heat_transfer_model(D_range)

    # Load FDS simulation data
    # Note: Update this path to your actual CSV file location
    csv_path = Path("/Volumes/T7 Shield/bandy-circles/grid-search-combined") / "ideal_curvatures_heat_flux.csv"
    fds_df = load_fds_data(csv_path)

    if fds_df is not None:
        print("Processing FDS simulation data...")
        fds_results = process_fds_data(fds_df)

        if fds_results and len(fds_results['diameter']) > 0:
            print(f"Found {len(fds_results['diameter'])} FDS data points")
        else:
            print("No valid FDS data found for treatment height = 0.6m")
            fds_results = None
    else:
        print("No FDS data loaded - showing heat transfer model only")
        fds_results = None

    # Create comparison plot
    fig = create_comparison_plot(heat_transfer_results, fds_results)

    # Print analysis summary
    print_analysis_summary(heat_transfer_results, fds_results)

    # Show plot
    plt.show()

    # Save plot
    fig.savefig('critical_wind_speed_comparison.png', dpi=300, bbox_inches='tight')
    print("\nPlot saved as 'critical_wind_speed_comparison.png'")


if __name__ == "__main__":
    main()
