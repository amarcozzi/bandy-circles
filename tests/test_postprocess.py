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

    def test_stitch_np(self):

        # three 3x3 arrays
        a = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        b = np.array([[6, 7, 8], [9, 10, 11], [12, 13, 14]])
        c = np.array([[12, 13, 14], [15, 16, 17], [18, 19, 20]])

        assert np.array_equal(np.vstack((a, b, c)),
                              pp.stitch_mesh_data_to_array([a, b, c]))

    def test_stitch_mesh_data_to_array(self):
        sim.meshes
        pass


class TestSliceData:

    def test_get_slice_data(self):

        slice_data = []
        for slice in sim.slices:
            slice_data.append(slice.to_global(return_coordinates=False))

        assert np.array_equal(
            slice_data, pp.get_slice_data(sim, sim))


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
