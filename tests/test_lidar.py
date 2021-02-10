#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `lidar` package."""


import unittest
from click.testing import CliRunner

# from lidar import lidar
from lidar import cli

import os
import pkg_resources
import richdem as rd
from scipy import ndimage
import numpy as np
import time
import lidar

class TestLidar(unittest.TestCase):
    """Tests for `lidar` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_000_something(self):
        """Test something."""

    def test_mean_filter(self):

        # identify the sample data directory of the package
        package_name = "lidar"
        data_dir = pkg_resources.resource_filename(package_name, "data/")
        print("Sample data directory: {}".format(data_dir))

        # use the sample dem. Change it to your own dem if needed
        in_dem = os.path.join(data_dir, "dem.tif")
        out_dir = os.path.join(os.path.expanduser("~"), "temp")
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        mean_dem = os.path.join(out_dir, "mean.tif")
        mean = lidar.MeanFilter(in_dem, kernel_size=3)  
        rd.SaveGDAL(mean_dem, mean)

        self.assertTrue(os.path.exists(mean_dem))

    def test_median_filter(self):

        # identify the sample data directory of the package
        package_name = "lidar"
        data_dir = pkg_resources.resource_filename(package_name, "data/")
        print("Sample data directory: {}".format(data_dir))

        # use the sample dem. Change it to your own dem if needed
        in_dem = os.path.join(data_dir, "dem.tif")
        out_dir = os.path.join(os.path.expanduser("~"), "temp")
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        median_dem = os.path.join(out_dir, "median.tif")
        median = lidar.MedianFilter(in_dem, kernel_size=3)  
        rd.SaveGDAL(median_dem, median)

        self.assertTrue(os.path.exists(median_dem))

    def test_gaussian_filter(self):

        # identify the sample data directory of the package
        package_name = "lidar"
        data_dir = pkg_resources.resource_filename(package_name, "data/")
        print("Sample data directory: {}".format(data_dir))

        # use the sample dem. Change it to your own dem if needed
        in_dem = os.path.join(data_dir, "dem.tif")
        out_dir = os.path.join(os.path.expanduser("~"), "temp")
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        gaussian_dem = os.path.join(out_dir, "gaussian.tif")
        gaussian = lidar.GaussianFilter(in_dem,  sigma=1)  
        rd.SaveGDAL(gaussian_dem, gaussian)

        self.assertTrue(os.path.exists(gaussian_dem))

    def test_sink_filling(self):

        # identify the sample data directory of the package
        package_name = "lidar"
        data_dir = pkg_resources.resource_filename(package_name, "data/")
        print("Sample data directory: {}".format(data_dir))

        # use the sample dem. Change it to your own dem if needed
        in_dem = os.path.join(data_dir, "dem.tif")
        # parameters for depression filling
        min_size = 1000  # minimum number of pixels as a depression
        # min_depth = 0.3  # minimum depression depth
        # set output directory
        out_dir = os.path.join(
            os.path.expanduser("~"), "temp"
        )  # create a temp folder under user home directory
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        sink_path = lidar.ExtractSinks(in_dem, min_size=min_size, out_dir=out_dir)
        self.assertTrue(os.path.exists(sink_path))

    def test_slicing(self):

        # set input files
        # identify the sample data directory of the package
        package_name = "lidar"
        data_dir = pkg_resources.resource_filename(package_name, "data/")        
        # in_dem = os.path.join(data_dir, "dem.tif")
        in_sink = os.path.join(data_dir, "sink.tif")
        # parameters for level set method
        min_size = 1000         # minimum number of pixels as a depression
        min_depth = 0.3         # minimum depression depth
        interval = 0.3          # slicing interval, top-down approach
        bool_level_shp = True  # whether or not to extract polygons for each individual level
        # set output directory
        out_dir = os.path.join(os.path.expanduser("~"), "temp")  # create a temp folder under user home directory
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        dep_id_path, dep_level_path = lidar.DelineateDepressions(in_sink, min_size, min_depth, interval, out_dir, bool_level_shp)
        print("Results are saved in: {}".format(out_dir))

        self.assertTrue(os.path.exists(dep_id_path))
        self.assertTrue(os.path.exists(dep_level_path))

    def test_mounts(self):

        # identify the sample data directory of the package
        package_name = "lidar"
        data_dir = pkg_resources.resource_filename(package_name, "data/")

        # use the sample dem. Change it to your own dem if needed
        in_dem = os.path.join(data_dir, "dsm.tif")
        # set output directory. By default, use the temp directory under user's home directory
        out_dir = os.path.join(os.path.expanduser("~"), "temp")

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # parameters for identifying sinks and delineating nested depressions
        min_size = 1000  # minimum number of pixels as a depression
        min_height = 0.3  # minimum depth as a depression
        interval = 0.3  # slicing interval for the level-set method
        bool_shp = False  # output shapefiles for each individual level

        mount_id_path, mount_level_path = lidar.DelineateMounts(
            in_dem, min_size, min_height, interval, out_dir, bool_shp
        )
        self.assertTrue(os.path.exists(mount_id_path))
        self.assertTrue(os.path.exists(mount_level_path))

