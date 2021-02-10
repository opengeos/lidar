# Welcome to the lidar package

[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/lidar-cloud)
[![image](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/giswqs/lidar/master?filepath=examples%2Flidar.ipynb)
[![image](https://img.shields.io/pypi/v/lidar.svg)](https://pypi.python.org/pypi/lidar)
[![image](https://pepy.tech/badge/lidar)](https://pepy.tech/project/lidar)
[![image](https://img.shields.io/conda/vn/conda-forge/lidar.svg)](https://anaconda.org/conda-forge/lidar)
[![image](https://github.com/giswqs/lidar/workflows/build/badge.svg)](https://github.com/giswqs/lidar/actions?query=workflow%3Abuild)
[![image](https://github.com/giswqs/geemap/workflows/docs/badge.svg)](https://lidar.gishub.org)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![image](https://img.shields.io/twitter/follow/giswqs?style=social)](https://twitter.com/giswqs)
[![image](https://img.shields.io/badge/Donate-Buy%20me%20a%20coffee-yellowgreen.svg)](https://www.buymeacoffee.com/giswqs)
[![image](https://joss.theoj.org/papers/005bf3e6f47840e74e71678d8e88facc/status.svg)](https://joss.theoj.org/papers/005bf3e6f47840e74e71678d8e88fac)

**lidar** is Python package for delineating the nested hierarchy of surface depressions in digital elevation models (DEMs). It is
particularly useful for analyzing high-resolution topographic data, such as DEMs derived from Light Detection and Ranging (LiDAR) data.

-   GitHub repo: <https://github.com/giswqs/lidar>
-   Documentation: <https://lidar.gishub.org>
-   PyPI: <https://pypi.org/project/lidar/>
-   Conda-forge: <https://anaconda.org/conda-forge/lidar>
-   Binder: <https://gishub.org/lidar-cloud>
-   Free software: [MIT license](https://opensource.org/licenses/MIT)

## Introduction

**lidar** is a Python package for delineating the nested hierarchy of
surface depressions in digital elevation models (DEMs). In traditional
hydrological modeling, surface depressions in a DEM are commonly treated
as artifacts and thus filled and removed to create a depressionless DEM,
which can then be used to generate continuous stream networks. In
reality, however, surface depressions in DEMs are commonly a combination
of spurious and actual terrain features. Fine-resolution DEMs derived
from Light Detection and Ranging (LiDAR) data can capture and represent
actual surface depressions, especially in glaciated and karst
landscapes. During the past decades, various algorithms have been
developed to identify and delineate surface depressions, such as
depression filling, depression breaching, hybrid breaching-filling, and
contour tree method. More recently, a level-set method based on graph
theory was proposed to delineate the nested hierarchy of surface
depressions. The **lidar** Python package implements the level-set
method and makes it possible for delineating the nested hierarchy of
surface depressions as well as elevated terrain features. It also
provides an interactive Graphical User Interface (GUI) that allows users
to run the program with minimal coding.

## Statement of Need

The **lidar** package is intended for scientists and researchers who
would like to integrate surface depressions into hydrological modeling.
It can also facilitate the identification and delineation of
depressional features, such as sinkholes, detention basins, and prairie
potholes. The detailed topological and geometric properties of surface
depressions can be useful for terrain analysis and hydrological
modeling, including the size, volume, mean depth, maximum depth, lowest
elevation, spill elevation, perimeter, major axis length, minor axis
length, elongatedness.

## Features

-   Smoothing DEMs using mean, median, and Gaussian filters (see [filtering](https://lidar.gishub.org/filtering))
-   Extracting depressions from DEMs (see [filling](https://lidar.gishub.org/filling)).
-   Filtering out small artifact depressions based on user-specified minimum depression size (see [filling](https://lidar.gishub.org/filling/)).
-   Generating refined DEMs with small depressions filled but larger depressions kept intact (see [filling](https://lidar.gishub.org/filling/)).
-   Delineating depression nested hierarchy using the level-set method (see [slicing](https://lidar.gishub.org/slicing/)).
-   Delineating mount nested hierarchy using the level-set method (see [mounts](https://lidar.gishub.org/mounts/)).
-   Computing topological and geometric properties of depressions, including size, volume, mean depth, maximum depth, lowest elevation,
    spill elevation, perimeter, major axis length, minor axis length, elongatedness, eccentricity, orientation, and area-bbox-ratio (see [slicing](https://lidar.gishub.org/slicing/)).
-   Exporting depression properties as a csv file (see [slicing](https://lidar.gishub.org/slicing/)).