import fdsreader
from fdsreader import Simulation
import numpy as np


def get_normalized_ros(data):
    """
    Algorithm:

    Step 1) get a raster of boolean values. 1 is on fire, 0 is not on fire.
            NEED: Some threashold value for on fire/not on fire
    Step 2) Compute velocity at the following points:
            a) Through the center of the circle
            b) Left of the circle.
            c) Right of the circle.
    Step 3) Average left + right velocities and divide by center velocity to
            get normalized rate of spread.

    NOTE: we need to decide if we are going to do this once or for every time step.
    """


def get_curvature(data):
    """
    Algorithm:

    Step 1) get a raster of boolean values. 1 is on fire, 0 is not on fire.
            NEED: Some threashold value for on fire/not on fire
    Step 2) Get the coordinates of the fire front.
    Step 3) Prune list to include coordinates inside circle and within some distance from the circle.
            NEED: Some threshold value for distance from circle.
    Step 4) Fit a 2nd degree polynomial to the coordinates of the fire front.
    Step 5) Curvature = a_2 * 2

    Parameters
    ----------
    data : array-like object
        A 3D array-like object of data of a quantity. Dimensions of the array are (time, y, x).

    Returns
    -------
    curvature: array-like object
        A 1D array-like object of the curvature of the fire front through the circle over time.
        Dimensions of the array are (time).
    """


def get_fire_line(data):
    """
    Takes in a boolean array-like object and returns a list of tuples of the coordinates
    of the fire front. The format of coordinates is (x, y).

    The "fire front" is defined as the furtherst cell in the downwind direction that has
    a value of 1 in the boolean array. In our case the wind will blow from the bottom of
    the domain to the top.

    Note: There is some trickery here between cartesian (xy) and array (coordinates) when
    dealing with notions like "top" and "bottom". We'll need to be on the same page here.

    Parameters
    ----------
    data : array-like object
        A 2D array-like object of boolean values. 1 is on fire, 0 is not on fire.

    Returns
    -------
    fire_line : list of tuples
        A list of tuples of the coordinates of the fire front. The format of coordinates
        is (x, y).
    """


def get_active_fire_array(data, threshold):
    """
    This function will take in a array-like object and return an array-like object of
    boolean values. 1 is on fire, 0 is not on fire.

    The threshold is the minimum value for a cell to be considered on fire for the given
    quantity associated with the data.

    Parameters
    ----------
    data : array-like object
        An array-like object of data of a quantity.
    threshold : float
        The minimum value for a cell to be considered on fire.

    Returns
    -------
    active_fire_array : array-like object
        An array-like object of boolean values. 1 is on fire, 0 is not on fire.
    """


def stitch_mesh_data_to_array(list_of_meshes):
    """
    Takes data from an individual mesh and stitches it to a larger array.

    Step 1) Determine the size of the larger array.
    Step 2) Initialize the larger array.
    Step 3) Iterate over each mesh
        a) Determine the offset of the mesh in the larger array.
        b) Add the mesh to the larger array.

    Parameters
    ----------
    list_of_meshes : list[mesh]
        A list of mesh objects.

    Returns
    -------
    stitched_data : array-like object
        A 3D array-like object of data of a quantity. Dimensions of the array are (time, y, x).
    """

    return np.vstack(list_of_meshes)

    data = []
    meshes = np.vstack(list_of_meshes)
    for row in range(len(meshes)):
        if row == 0:  # add the first row to the data
            data.append(list_of_meshes[row])
        # if the current row is not equal to the previous row, add it to the data
        if not np.array_equal(list_of_meshes[row], list_of_meshes[row-1]):
            data.append(list_of_meshes[row])

    stitched_data = np.array(data)

    return stitched_data


def get_bndf_data(sim, qty):
    """
    This function will take in a simulation file and return a 3D array-like object of
    the data for the given quantity. Dimensions of the array are (time, y, x).

    The data comes from the fdsreader bndf object.

    NOTE: It may be necessary to stitch together multiple meshes to get the full data.

    Parameters
    ----------
    sim : fdsreader.Simulation
        A simulation object.
    qty : string
        The name of the quantity to get data for.

    Returns
    -------
    data : array-like object
        A 3D array-like object of data of a quantity. Dimensions of the array are (time, y, x).
    """
    # get global boundary data arrays for each individual mesh
    data = []

    for mesh in sim.meshes:
        data.append(
            mesh.get_boundary_data(quantity=qty))

    # stitch the data together
    data = stitch_mesh_data_to_array(data)

    return data


def get_slice_data(sim, qty):
    """
    This function will take in a simulation file and return a 3D array-like object of
    the data for the given quantity. Dimensions of the array are (time, y, x).

    The data comes from the fdsreader slice object.

    NOTE: It may be necessary to stitch together multiple meshes to get the full data.

    Parameters
    ----------
    sim : fdsreader.Simulation
        A simulation object.
    qty : string
        The name of the quantity to get data for.

    Returns
    -------
    data : array-like object
        A 3D array-like object of data of a quantity. Dimensions of the array are (time, y, x).
    """
    # get slice data arrays for each individual mesh
    data = []
    coords = []
    times = []

    for slice in sim.slices:

        # Creates a global numpy ndarray from all subslices (of the slice)
        # xyz is the returned matching coordinate for each value on the generated grid (data).
        # can return large sparse arrays for some slices
        slice_data, slice_coords = slice.to_global(
            return_coordinates=True)
        # times for each slice
        slice_times = slice.times

        data.append(slice_data)
        coords.append(slice_coords)
        times.append(slice_times)

    return data


def process_simulation(sim_id):
    pass
    # TODO: create a simulation object with out_simulation_id.smv file
    sim = Simulation(sim_id.smv)

    # TODO: get the simulation data
    hrr_array = get_slice_data(sim, "HRRPUV")
    mass_flux_array = get_slice_data(sim, "MASS FLUX")
    data_dict = {"hrr": hrr_array, "mass_flux": mass_flux_array}

    # TODO: Take the data and get outputs
    for qty, data in data_dict.items():
        min_normalized_ros = get_normalized_ros(data)
        curvature = get_curvature(data)

        # TODO: Do something with the outputs like take the minimum

    # return outputs


def main():
    pass
    # TODO: Parse arguments and get a simulation ID list
    sim_id_list = parse_arguments(args)

    # TODO: Iterate over the simulation ID list and postprocess the simulation
    for sim_id in sim_id_list:
        outputs = process_simulation(sim_id)

    # TODO: Save the outputs to the appropriate file (results.csv or simulation_id.out)
    save_outputs(outputs)


if __name__ == "__main__":
    main()
