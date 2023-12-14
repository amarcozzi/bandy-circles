# export_array_to_fds("data", "control", control_array, dx, dy, control_height)
from pathlib import Path
from typing import Union
import numpy as np
from scipy.io import FortranFile


def export_array_to_fds(output_dir: Path | str, name: str, array: np.ndarray, dx: float, dy: float, dz: float) -> None:
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
    x_coords = np.arange(dx / 2, nx * dx, dx)
    y_coords = np.arange(dy / 2, ny * dy, dy)
    z_coords = np.arange(dz / 2, nz * dz, dz)

    # Create a FortranFile object for the output file
    output_path = output_dir / f"{name}.bdf"
    with FortranFile(output_path, 'w') as f:
        # write out global bounding voxel faces.
        # (VXMIN,VXMAX,VYMIN,VYMAX,VZMIN,VZMAX)
        vxbounds = [np.min(x_coords) - dx / 2, np.max(x_coords) + dx / 2,
                    np.min(y_coords) - dy / 2, np.max(y_coords) + dy / 2,
                    np.min(z_coords) - dz / 2, np.max(z_coords) + dz / 2]
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
