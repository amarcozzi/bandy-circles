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

# sample simulation data
sim = fds.Simulation(
    r"./tests/testing_data/Case_C064_fine_out_cat.smv")


class TestProcessSimulation:

    def test_process_simulation(self):

        pass


class TestStitchedData:

    def test_number_of_meshes(self):

        assert len(sim.meshes) == 36

    def test_stitch_np(self):

        # three 3x3 arrays
        a = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        b = np.array([[6, 7, 8], [9, 10, 11], [12, 13, 14]])
        c = np.array([[12, 13, 14], [15, 16, 17], [18, 19, 20]])

        assert np.array_equal(np.vstack((a, b, c)),
                              pp.stitch_mesh_data_to_array([a, b, c]))

    def test_stitch_mesh_data_to_array(self):

        assert np.array_equal(np.vstack(sim.meshes),
                              pp.stitch_mesh_data_to_array(sim.meshes))


class TestSliceData:

    def test_len_of_slice_data(self):

        assert len(pp.get_slice_data(sim, sim)) == len(sim.slices)

    def test_get_slice_data_values(self):

        method_data = pp.get_slice_data(sim, sim)
        slice_data = []
        data_ind = []
        idx = 0
        for slice in sim.slices:

            data, coords = slice.to_global(
                return_coordinates=True, masked=True)

            slice_data.append(data)

            data_ind.append(np.array_equal(
                data, method_data[idx]))
            idx += 1

        assert all(data_ind)


class TestBNDFData:

    def test_get_bndf_data_wall_temp(self):

        bndf_data = []
        for mesh in sim.meshes:
            bndf_data.append(mesh.get_boundary_data(
                quantity="WALL TEMPERATURE"))
        sim.meshes[0].get_boundary_data(quantity="WALL TEMPERATURE")

        bndf_data = np.vstack(bndf_data)
        assert np.array_equal(
            bndf_data, pp.get_bndf_data(sim, "WALL TEMPERATURE"))

    # def test_get_bndf_data_temp(self):

    #     bndf_data = []
    #     for mesh in sim.meshes:
    #         bndf_data.append(mesh.get_boundary_data(quantity="TEMPERATURE"))

    #     bndf_data = np.vstack(bndf_data)
    #     # assert pp.get_bndf_data(sim, "TEMPERATURE") == bndf_data
    #     pass


class TestFireLine:
    pass


class TestActiveFireArray:
    pass


if __name__ == "__main__":
    slice = TestSliceData()
    slice.test_get_slice_data()
