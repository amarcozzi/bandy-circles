import numpy as np

wind_speeds = np.linspace(0.75, 2.5, 10)

print([f"{speed:.2f}" for speed in wind_speeds])

circle_radius = 1.8

circle_radius_samples = [0.5 * circle_radius, 0.75 * circle_radius, circle_radius, 1.25 * circle_radius, 1.5 * circle_radius]
print(circle_radius_samples)