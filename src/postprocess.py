import fdsreader as fds
from fdsreader import Simulation
import numpy as np


DX = 0.1  # m
DY = 0.1  # m


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
    print()

    data = np.array(([0, 0, 0, 0, 1],
                    [1, 1, 1, 0, 0],
                    [1, 1, 1, 1, 0],
                    [1, 1, 1, 1, 1],
                    [0, 1, 1, 1, 1]))

    u_ind = []
    v_ind = []

    # grab the top-most array locations where the data = 1
    for column in range(data.shape[1]):  # for each column
        for row in range(data.shape[0]):  # for each row
            if data[row, column] == 1:
                u_ind.append(column)
                v_ind.append(row)
                break

    x_coords = []
    y_coords = []
    # convert array locations into (x, y) cartesian coordinates
    # TODO: determine how to calculate the values of tx and ty
    tx = data.shape[0] * DX / 2
    ty = data.shape[1] * DY / 2

    for u in u_ind:
        x = u * DX - tx
        x_coords.append(x)

    for v in v_ind:
        y = v * -DY + ty
        y_coords.append(y)

    # combine the x and y coordinates into a list of tuples
    fire_line = list(zip(x_coords, y_coords))
    return fire_line


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

    # for a given timestep get active fire array

    active_fire_array = np.zeros(data.shape)

    active_fire_array[data > threshold] = 1
    return active_fire_array


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
    # FIXME: The meshes or slice arrays may not be the same shape. How do we handle this? vstack will not work in this case.
    # determine the size of the stitched array

    # treat as 2D array and iterate over time
    # for each timestep, stitch the meshes together
    # iterate over each timestep and stitch these meshes together

    for timestep in list_of_meshes[0, :, :]:
        continue
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
        if str(slice.quantity.name) == qty:

            # Creates a global numpy ndarray from all subslices (of the slice)
            # xyz is the returned matching coordinate for each value on the generated grid (data).
            # can return large sparse arrays for some slices
            slice_data, slice_coords = slice.to_global(
                return_coordinates=True, masked=True)
            # times for each slice
            slice_times = slice.times

            data.append(slice_data)
            coords.append(slice_coords)
            times.append(slice_times)
            break

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

    sim = fds.Simulation(
        r"./tests/testing_data/test_data/out_crop_circles_cat.smv")

    hrr_array = get_slice_data(sim, "HRRPUV")
    mass_flux_array = get_slice_data(sim, "MASS FLUX")

    # stitched = stitch_mesh_data_to_array(hrr_array)
    data_dict = {"hrr": hrr_array, "mass_flux": mass_flux_array}

    # slice of hrr_array data: (time, y, x)
    # need to get (x, y) coordinates of slice (2D array without time units)

    hrr_x = hrr_array[0][1]  # == hrr_array[0][1, :, :] # THIS IS TRUE
    print(f'HRR X:\n{hrr_x}\n\n')
    hrr_y = hrr_array[0][2]  # == hrr_array[0][2, :, :] # THIS IS TRUE
    print(f'HRR Y:\n{hrr_y}\n\n')
    print((hrr_x) == (hrr_y))
    print(hrr_x[0][-1])
    print(hrr_y[0][0])

    # get active fire array for a given slice  at a given timestep
    # slice 0 timestep 0
    fire = get_active_fire_array(hrr_array[0][0], 0.1)

    # TODO: input data needs to be a 2D array of boolean values
    get_fire_line(fire)
    print()
    return
    # TODO: Parse arguments and get a simulation ID list
    sim_id_list = parse_arguments(args)

    # TODO: Iterate over the simulation ID list and postprocess the simulation
    for sim_id in sim_id_list:
        outputs = process_simulation(sim_id)

    # TODO: Save the outputs to the appropriate file (results.csv or simulation_id.out)
    save_outputs(outputs)


if __name__ == "__main__":
    main()
