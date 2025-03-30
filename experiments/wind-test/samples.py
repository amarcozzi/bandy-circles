import numpy as np

wind_speeds = np.linspace(1.5, 5, 12)

print([float(round(speed, 2)) for speed in wind_speeds])

circle_radius = 1.8

circle_radius_samples = [0.5 * circle_radius, 0.75 * circle_radius, circle_radius, 1.25 * circle_radius,
                         1.5 * circle_radius]
print(circle_radius_samples)
