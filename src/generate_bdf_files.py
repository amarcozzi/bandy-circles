import numpy as np
import matplotlib.pyplot as plt

# constants
Radius = 1.8288  # m

ControlMass = 0.5  # kg/m^3
TreatmentMass = 1.27  # kg/m^3

ControlHeight = 0.6  # m
TreatmentHeight = 0.15  # m


# start of function
xmin = -9  # m
xmax = 9  # m
dx = 0.1  # m

ymin = -12  # m
ymax = 12  # m
dy = 0.1  # m

# x_coords and y_coords cell centered
x_coords = np.arange(xmin + dx/2, xmax + dx/2, dx)
y_coords = np.arange(ymax - dy/2, ymin - dy/2, -dy)

xx, yy = np.meshgrid(x_coords, y_coords, indexing='xy')


inside_circle = np.square(xx) + np.square(yy) <= np.square(Radius)

# make a 3D array (NY, NX, 1)
control_array = np.ones(inside_circle.shape) * ControlMass
control_array *= ~inside_circle
