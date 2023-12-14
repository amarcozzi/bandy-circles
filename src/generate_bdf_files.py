import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import FortranFile

# constants
RADIUS = 1.8288  # m

CONTROL_MASS = 0.5  # kg/m^3
TREATMENT_MASS = 1.27  # kg/m^3

CONTROL_HEIGHT = 0.6  # m
TREATMENT_HEIGHT = 0.15  # m

def export_array_to_fds(output_dir: Path | str, name: str, array: np.ndarray, dx: float, dy: float, dz: float, x_coords: np.ndarray = None, y_coords: np.ndarray = None) -> None:
    """
    Write a 3D numpy array to an FDS bulk density input file
    """
    # Convert the output directory to a Path object if it is a string
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)

    # Get the attributes from the zarr file
    ny, nx, nz = array.shape

    # FastFuels stores the array as (y, x, z) so we need to swap the axes
    # to (x, y, z) for FDS.
    bd_data = np.swapaxes(array, 0, 1).copy()

    # meshgrid of x,y,z voxel centers
    if x_coords is None:
        x_coords = np.arange(dx / 2, nx * dx, dx)
    if y_coords is None:
        y_coords = np.arange(dy / 2, ny * dy, dy)
    z_coords = np.arange(dz / 2, nz * dz, dz)

    # Create a FortranFile object for the output file
    output_path = output_dir / f"{name}.bdf"
    with FortranFile(output_path, 'w') as f:
        # write out global bounding voxel faces.
        # (VXMIN,VXMAX,VYMIN,VYMAX,VZMIN,VZMAX)
        vxbounds = [np.round(np.min(x_coords) - dx / 2, 2), np.round(np.max(x_coords) + dx / 2, 2),
                    np.round(np.min(y_coords) - dy / 2, 2), np.round(np.max(y_coords) + dy / 2, 2),
                    np.round(np.min(z_coords) - dz / 2, 2), np.round(np.max(z_coords) + dz / 2, 2)]
        f.write_record(np.array(vxbounds, dtype=np.float64))

        # write out voxel resolution (VDX, VDY, VDZ)
        f.write_record(np.array([dx, dy, dz], dtype=np.float64))

        # write out number of voxels for relevant fuel class (NVOX)
        nvox = np.count_nonzero(bd_data)
        f.write_record(np.array(nvox, dtype=np.int32))

        # write out center and bulk density for each occupied voxel
        # (VCX, VCY, VCZ)
        # (MASS_PER_VOLUME)
        for i, j, k in np.argwhere(bd_data > 0):
            f.write_record(
                np.array([x_coords[i], y_coords[j], z_coords[k]],
                         dtype=np.float64))
            f.write_record(np.array(bd_data[i, j, k], dtype=np.float64))


def generate_bdf_files(out_path: Path, control_mass: float, treatment_mass: float, xmin=-9, xmax=9, dx=0.1, ymin=-12, ymax=12, dy=0.1):
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

    control_array = np.ones((ny, nx, 1)) * control_mass
    control_array[:, :, 0] *= ~inside_circle

    treatment_array = np.ones((ny, nx, 1)) * treatment_mass
    treatment_array[:, :, 0] *= inside_circle

    export_array_to_fds(
        out_path, "control", control_array, dx, dy, CONTROL_HEIGHT, x_coords, y_coords)
    export_array_to_fds(
        out_path, "treatment", treatment_array, dx, dy, TREATMENT_HEIGHT, x_coords, y_coords)


if '__main__' == __name__:
    # test generate_bdf_files
    generate_bdf_files()
