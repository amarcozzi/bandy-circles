import numpy as np

# Generate 10 wind samples from 1 to 3 m/s
wind_speeds = np.linspace(1, 3, 10)
wind_speeds_str = ", ".join([f"{speed:.2f}" for speed in wind_speeds])
print(f"Wind speeds: [{wind_speeds_str}]")
