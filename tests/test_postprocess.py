"""
Tests for the postprocess module to verify simulation data is processed correctly.
"""

import os
import sys  # noqa
sys.path.append('./')  # noqa

import pytest
import numpy as np
from src import postprocess as pp
import fdsreader as fds

# testing simulation data
sim = fds.Simulation(
    r"./tests/testing_data/3_meshes/test_data/out_crop_circles_cat.smv")


class TestProcessSimulation:

    def test_process_simulation(self):

        pass


class TestStitchedData:

    def test_stitch_np(self):

        # Original array (3, 4, 5)
        original_array = np.array([[[1, 2, 3, 4, 5],
                                    [6, 7, 8, 9, 10],
                                    [11, 12, 13, 14, 15],
                                    [16, 17, 18, 19, 20]],

                                   [[21, 22, 23, 24, 25],
                                    [26, 27, 28, 29, 30],
                                    [31, 32, 33, 34, 35],
                                    [36, 37, 38, 39, 40]],

                                   [[41, 42, 43, 44, 45],
                                    [46, 47, 48, 49, 50],
                                    [51, 52, 53, 54, 55],
                                    [56, 57, 58, 59, 60]]])

        # Create two (3, 2, 5) arrays
        array1 = original_array[:, :2, :]
        array2 = original_array[:, 2:4, :]

        assert np.array_equal(
            original_array, pp.stitch_mesh_data_to_array([array1, array2]))


class TestSliceData:

    def test_number_of_slices(self):

        assert len(sim.slices) == 6

    def test_len_of_slice_data(self):
        return

    def test_get_slice_data_values(self):

        # method_data, _ = pp.get_slice_data(sim, "HRRPUV")
        # slice_data = []
        # data_ind = []
        # idx = 0
        # for slice in sim.slices:

        #     data, coords = slice.to_global(
        #         return_coordinates=True, masked=True)

        #     slice_data.append(data)

        #     data_ind.append(np.array_equal(
        #         data, method_data[idx]))
        #     idx += 1

        # assert all(data_ind)
        return


class TestBNDFData:

    def test_number_of_meshes(self):

        assert len(sim.meshes) == 3

    def test_number_of_boundary_data(self):

        # grab boundary objects
        meshes = []
        for mesh in sim.meshes:
            meshes.append(mesh.get_boundary_data(quantity="TOTAL HEAT FLUX"))

        # grab data from boundary object
        bndf_data = []
        for mesh in range(len(meshes)):
            bndf_data.append(meshes[mesh].data[len(meshes)].data)

        assert len(bndf_data) == len(sim.meshes)

    def test_get_bndf_data_heat_flux(self):

        # grab boundary objects
        meshes = []
        for mesh in sim.meshes:
            meshes.append(mesh.get_boundary_data(quantity='TOTAL HEAT FLUX'))

        # grab data from boundary object
        bndf_data = []
        for mesh in range(len(meshes)):
            bndf_data.append(meshes[mesh].data[len(meshes)].data)

        # TODO: look into if this works once using data with 1 mesh
        if len(bndf_data) == 1:
            assert np.array_equal(bndf_data[0],
                                  pp.get_bndf_data(sim, 'TOTAL HEAT FLUX'))
        # stitch data together
        data_array = [arr for arr in bndf_data]
        # concatenate time, x, and y
        stitched_data = np.concatenate([arr[:, :, :]
                                        for arr in data_array], axis=1)

        postproc_data, _ = pp.get_bndf_data(sim, 'TOTAL HEAT FLUX')

        assert np.array_equal(
            stitched_data, postproc_data)


class TestFireLine:
    pass


class TestActiveFireArray:

    # TODO: test the shape of the returned array is the shape of a slice[i] from sim
    pass


if __name__ == "__main__":
    slice = TestSliceData()
    slice.test_get_slice_data()


# display boolean on fire/ not on fire for validation
