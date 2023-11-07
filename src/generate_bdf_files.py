import numpy as np
import matplotlib.pyplot as plt

# constants
Radius = 1.8288  # m

ControlMass = 0.5  # kg/m^3
TreatmentMass = 1.27  # kg/m^3

ControlHeight = 0.6  # m
TreatmentHeight = 0.15  # m


# start of function

def generate_bdf_files(xmin=-9, xmax=9, dx=0.1, ymin=-12, ymax=12, dy=0.1):
    '''
    Generate BDF files for the control and treatment areas.

    Parameters
    ----------
        xmin : float
            Minimum x coordinate of the domain.
        xmax : float
            Maximum x coordinate of the domain.
        dx : float
            Spacing between x coordinates.
        ymin : float
            Minimum y coordinate of the domain.
        ymax : float
            Maximum y coordinate of the domain.
        dy : float
            Spacing between y coordinates.

    Returns
    -------
    # TODO: add return value
        None
    '''

    # x_coords and y_coords cell centered
    x_coords = np.arange(xmin + dx/2, xmax + dx/2, dx)
    y_coords = np.arange(ymax - dy/2, ymin - dy/2, -dy)

    xx, yy = np.meshgrid(x_coords, y_coords, indexing='xy')

    # TODO: we need this as a 3D array (?)
    inside_circle = np.square(xx) + np.square(yy) <= np.square(Radius)

    NY = inside_circle.shape[0]
    NX = inside_circle.shape[1]
    # make a 3D array (NY, NX, 1)
    control_array = np.ones((NY, NX, 1)) * ControlMass
    control_array[:, :, 0] *= ~inside_circle

    treatment_array = np.ones((NY, NX, 1)) * TreatmentMass
    treatment_array[:, :, 0] *= inside_circle


if 'main' == __name__:
    generate_bdf_files()
