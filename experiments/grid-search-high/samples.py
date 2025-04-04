import numpy as np

wind_speeds = np.linspace(2, 3.4, 15)

print(wind_speeds.round(2).tolist())

circle_radius = 1.8

circle_radius_samples = [0.5 * circle_radius, 0.75 * circle_radius, circle_radius, 1.25 * circle_radius,
                         1.5 * circle_radius]
print(circle_radius_samples)
