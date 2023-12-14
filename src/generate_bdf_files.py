import numpy as np
import matplotlib.pyplot as plt
import export_array_to_fds as eaf

# constants
RADIUS = 1.8288  # m

CONTROLMASS = 0.5  # kg/m^3
TREATMENTMASS = 1.27  # kg/m^3

CONTROLHEIGHT = 0.6  # m
TREATMENTHEIGHT = 0.15  # m


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
        None
    '''

    # x_coords and y_coords cell centered
    x_coords = np.arange(xmin + dx/2, xmax + dx/2, dx)
    y_coords = np.arange(ymax - dy/2, ymin - dy/2, -dy)

    xx, yy = np.meshgrid(x_coords, y_coords, indexing='xy')

    inside_circle = np.square(xx) + np.square(yy) <= np.square(RADIUS)

    ny = inside_circle.shape[0]
    nx = inside_circle.shape[1]
    # make a 3D array (ny, nx, 1)
    control_array = np.ones((ny, nx, 1)) * CONTROLMASS
    control_array[:, :, 0] *= ~inside_circle

    treatment_array = np.ones((ny, nx, 1)) * TREATMENTMASS
    treatment_array[:, :, 0] *= inside_circle

    eaf.export_array_to_fds(
        ".", "control", control_array, dx, dy, CONTROLHEIGHT)
    eaf.export_array_to_fds(
        ".", "treatment", treatment_array, dx, dy, TREATMENTHEIGHT)


if '__main__' == __name__:
    # test generate_bdf_files
    generate_bdf_files()
