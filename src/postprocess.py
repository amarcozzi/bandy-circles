import sys
import numpy as np
import pandas as pd
import fdsreader as fds
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
from numpy import ndarray
from pandas import DataFrame
from fdsreader import Simulation


ON_FIRE_THRESHOLD = 50  # kW/m^2 for TOTAL HEAT FLUX


def get_normalized_ros(data, times, coords, circle_radius):
    """
    Algorithm:

    Step 1) Compute velocity at the following points:
            a) Through the center of the circle
            b) Left of the circle.
            c) Right of the circle.
    Step 2) Average left + right velocities and divide by velocity through the
            treatment to get normalized rate of spread.
    """
    # Get the 1D rate of spread one radius length to the left of the edge of the circle
    left_point_start = (-2 * circle_radius, -circle_radius)
    left_point_end = (-2 * circle_radius, 0)
    left_control_ros = get_1D_ros_along_ray(
        data, times, coords, left_point_start, left_point_end
    )

    # Get the 1D rate of spread one radius length to the right of the edge of the circle
    right_point_start = (2 * circle_radius, -circle_radius)
    right_point_end = (2 * circle_radius, 0)
    right_control_ros = get_1D_ros_along_ray(
        data, times, coords, right_point_start, right_point_end
    )

    # Get the 1D rate of spread through the center of the circle
    center_point_start = (0, -circle_radius)
    center_point_end = (0, 0)
    treatment_ros = get_1D_ros_along_ray(
        data, times, coords, center_point_start, center_point_end
    )

    # Average the left and right ros
    control_ros = (left_control_ros + right_control_ros) / 2

    # Normalize the control ros
    try:
        normalized_ros = treatment_ros / control_ros
    except ZeroDivisionError:
        return 0
    return normalized_ros


def get_control_ros(data, times, coords):
    """
    This function calculates the 1D rate of spread through two rays that are on either side of
    the treatment (circle). The 1D rate of spread is then the average of the two rays.
    """
    # Calculate the 1D ros from (-2, -6) to (-2, 0)
    left_ros = get_1D_ros_along_ray(data, times, coords, (-2, -6), (-2, 0))

    # Calculate the 1D ros from (2, -6) to (2, 0)
    right_ros = get_1D_ros_along_ray(data, times, coords, (2, -6), (2, 0))

    # Average the left and right ros
    ros = (left_ros + right_ros) / 2

    return ros


def get_1D_ros_along_ray(
    fire: ndarray, times: ndarray, coords: dict, start_point: float, end_point: float
) -> float:
    """
    This function calculates the 1D rate of spread along a ray that is on either side of
    the treatment (circle). The 1D rate of spread is then the average of the two rays.
    """
    # Find the cell index of the start and end points
    i_start, j_start = get_cell_index(coords, start_point)
    i_end, j_end = get_cell_index(coords, end_point)

    # Find the first time that the fire front reaches the start point
    start_time = get_time_of_arrival(fire, times, i_start, j_start)

    # Find the first time that the fire front reaches the end point
    end_time = get_time_of_arrival(fire, times, i_end, j_end)

    # This case checks if the fire front never reaches the start or end point
    if start_time is None or end_time is None:
        return 0

    # Calculate the 1D rate of spread
    y_distance = end_point[1] - start_point[1]
    ros = y_distance / (end_time - start_time)

    return ros


def get_cell_index(coords, xy_point):
    """
    This function returns the cell index of a given point.
    """
    col = np.argmin(np.abs(coords["x"] - xy_point[0]))
    row = np.argmin(np.abs(coords["y"] - xy_point[1]))
    return row, col


def get_time_of_arrival(fire: np.ndarray, times: np.ndarray, i: int, j: int):
    """
    This function returns the first time that the fire front reaches a given cell.
    fire is a 3D array-like object of boolean values. 1 is on fire, 0 is not on fire.
    The dimensions of the array are (time, y, x). Times is a 1D array-like object of the
    times associated with the fire array. The indices of the fire array correspond to the
    indices of the times array.
    """
    t = fire[:, i, j].argmax()
    if fire[t, i, j] == 1:
        return times[t]
    else:
        return None


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


def get_fire_line(data, coords):
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

    # TODO: Added this parameter. Is this correct?
    coords : dict
        A dictionary object of the (x, y, z) coordinates from the simulation.

    Returns
    -------
    fire_line : list of tuples
        A list of tuples of the coordinates of the fire front for a given timestep. The format of coordinates
        is (x, y).
    """
    # data is:
    # for testing data returning a list of 337 2D arrays of size (123, 301)

    # TODO: need to determine DX and DY from the simulation

    # TODO: store timestep, column, and row in a dictionary

    u_ind = {}
    v_ind = {}

    # grab the top-most array locations where the data = 1
    # one array location for each column
    for timestep in range(len(data)):
        # for each column in a given timestep
        for column in range(data[timestep].shape[1]):
            # for each row in a given timestep
            for row in range(data[timestep].shape[0]):
                if data[timestep][row, column] == 1:
                    u_ind[timestep] = column
                    v_ind[timestep] = row
                    break
    print()
    # getting dx and dy from the coords dictionary
    dx = round(coords["x"][1] - coords["x"][0], 5)
    dy = round(coords["y"][1] - coords["y"][0], 5)

    x_coords = []
    y_coords = []

    x_coords_data = coords["x"]
    y_coords_data = coords["y"]

    # convert array locations into (x, y) cartesian coordinates
    for u in u_ind:
        x = x_coords_data[u] + DX / 2
        x_coords.append(x)

    for v in v_ind:
        # grabbing -v becuase the y coordinates start in bottom left corner; we need to flip the y coordinates
        y = y_coords_data[-v] - DY / 2
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
    active_fire_array = np.where(data >= threshold, True, False)

    return active_fire_array


def stitch_mesh_data_to_array(list_of_meshes, coords=None):
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
    # treat as 2D array and iterate over time
    # for each timestep, stitch the meshes together
    # iterate over each timestep and stitch these meshes together

    data_array = [arr for arr in list_of_meshes]
    # concatenate time, x, and y
    # TODO: reflect this in docstring?
    stitched_data = np.concatenate([arr[:, :, :] for arr in data_array], axis=1)

    if coords:
        # concatenate x coordinates and add y, z coordinates to coords
        x_coords = [arr for arr in coords]
        coords = {"x": x_coords, "y": coords[0]["y"], "z": coords[0]["z"]}

        return stitched_data, coords

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

    # TODO: Added this return. Is this correct?
    coords : dict
        A dictionary object of the (x, y, z) coordinates from the simulation.
    """

    for mesh in sim.meshes:
        boundary = mesh.get_boundary_data(quantity=qty)

        # fdsreader returns data in (time, x, y) (this makes array operations difficult)
        data = boundary.data[boundary.orientations[0]].data

        # Make the data (time, y, x) so that array operations make sense
        data = np.swapaxes(data, 1, 2)

        times = boundary.times
        coords = boundary.data[boundary.orientations[0]].get_coordinates()

    return data, times, coords


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

    # TODO: Added this return. Is this correct?
    coords : dict
        A dictionary object of the (x, y, z) coordinates from the simulation.
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
                return_coordinates=True, masked=True
            )
            # times for each slice
            slice_times = slice.times

            data.append(slice_data)
            coords.append(slice_coords)
            times.append(slice_times)
            break  # this break only grabs the first quantity

    return data, coords


def process_simulation(sim_dir: Path, sim_params: dict, make_plots: bool = False):
    # create a simulation object with fdsreader
    sim = Simulation(str(sim_dir))

    # Get the boundary data, times, and coordinates
    bndf_array, times, coords = get_bndf_data(sim, "TOTAL HEAT FLUX")
    coords["dx"] = round(coords["x"][1] - coords["x"][0], 1)
    coords["dy"] = round(coords["y"][1] - coords["y"][0], 1)

    # Determine what cells are on fire for a given timestep
    active_fire_mask = get_active_fire_array(bndf_array, ON_FIRE_THRESHOLD)

    # Compute outputs
    control_ros = get_control_ros(active_fire_mask, times, coords)
    normalized_ros = get_normalized_ros(
        active_fire_mask, times, coords, sim_params["circle_radius"]
    )

    # Add outputs to a dictionary
    outputs = {"CONTROL_ROS": control_ros, "NORMALIZED_ROS": normalized_ros}

    return outputs


def postprocess(experiment_directory: str | Path, outputs_directory: str | Path = None):
    """
    This function takes in an experiment directory with a populated output

    Parameters
    ----------
    directory : string
        The path to the directory of simulation output files.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the directory does not exist.
    """
    if isinstance(experiment_directory, str):
        experiment_directory = Path(experiment_directory)
    if not experiment_directory.is_dir():
        raise ValueError(f"Experiment directory {experiment_directory} does not exist.")

    if outputs_directory is None:
        outputs_directory = experiment_directory / "output"
    if isinstance(outputs_directory, str):
        outputs_directory = Path(outputs_directory)
    if not outputs_directory.is_dir():
        raise ValueError(f"Outputs directory {outputs_directory} does not exist.")

    # Load the parameters.csv file into a dataframe
    inputs_df = pd.read_csv(experiment_directory / "parameters.csv")

    # Make a copy of the input parameters dataframe
    outputs_dict = {}

    # Iterate over each simulation in the experiment
    for row in tqdm(inputs_df.iterrows(), total=len(inputs_df)):
        sim_params = row[1].to_dict()
        sim_id = int(sim_params["simulation_id"])

        simulation_directory = outputs_directory / f"simulation_{sim_id}"
        sim_outputs = process_simulation(simulation_directory, sim_params)

        # Add the outputs to the outputs dataframe
        outputs_dict[sim_id] = sim_outputs

    # Combine the inputs and outputs into a single dataframe
    outputs_df = pd.DataFrame(outputs_dict).T
    outputs_df = inputs_df.join(outputs_df)

    # Save the outputs to a csv file
    outputs_df.to_csv(experiment_directory / "outputs.csv")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        postprocess(*sys.argv[1:])
    else:
        postprocess(
            "/home/anthony/Work/UM/bandy-circles/experiments/1D-grid-search-coarse",
            "/mnt/Data/bandy-circles/output",
        )
