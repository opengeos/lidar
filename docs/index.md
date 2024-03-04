# Welcome to the lidar package

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/giswqs/lidar/blob/master/examples/lidar_colab.ipynb)
[![image](https://img.shields.io/pypi/v/lidar.svg)](https://pypi.python.org/pypi/lidar)
[![image](https://pepy.tech/badge/lidar)](https://pepy.tech/project/lidar)
[![image](https://img.shields.io/conda/vn/conda-forge/lidar.svg)](https://anaconda.org/conda-forge/lidar)
[![image](https://github.com/opengeos/lidar/workflows/build/badge.svg)](https://github.com/opengeos/lidar/actions?query=workflow%3Abuild)
[![image](https://github.com/opengeos/lidar/workflows/docs/badge.svg)](https://lidar.gishub.org)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![image](https://img.shields.io/twitter/follow/giswqs?style=social)](https://twitter.com/giswqs)
[![image](https://img.shields.io/badge/Donate-Buy%20me%20a%20coffee-yellowgreen.svg)](https://www.buymeacoffee.com/giswqs)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02965/status.svg)](https://doi.org/10.21105/joss.02965)

**lidar** is Python package for delineating the nested hierarchy of surface depressions in digital elevation models (DEMs). It is
particularly useful for analyzing high-resolution topographic data, such as DEMs derived from Light Detection and Ranging (LiDAR) data.

-   GitHub repo: <https://github.com/opengeos/lidar>
-   Documentation: <https://lidar.gishub.org>
-   PyPI: <https://pypi.org/project/lidar>
-   Conda-forge: <https://anaconda.org/conda-forge/lidar>
-   Open in Colab: <https://gishub.org/lidar-colab>
-   Free software: [MIT license](https://opensource.org/licenses/MIT)

**Citations**

-   **Wu, Q.**, (2021). lidar: A Python package for delineating nested surface depressions from digital elevation data. _Journal of Open Source Software_, 6(59), 2965, <https://doi.org/10.21105/joss.02965>

-   **Wu, Q.**, Lane, C.R., Wang, L., Vanderhoof, M.K., Christensen,
    J.R., & Liu, H. (2019). Efficient Delineation of Nested Depression
    Hierarchy in Digital Elevation Models for Hydrological Analysis
    Using Level-Set Method. _Journal of the American Water Resources
    Association_. <https://doi.org/10.1111/1752-1688.12689> ([PDF](https://spatial.utk.edu/pubs/2019_JAWRA.pdf))

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

## State of the Field

Currently, there are a few open-source Python packages that can perform depression filling on digital elevation data, such as [RichDEM](https://richdem.readthedocs.io/) and [whitebox](https://github.com/giswqs/whitebox-python), the Python frontend for [WhiteboxTools](https://github.com/jblindsay/whitebox-tools). However, there are no Python packages offering tools for delineating the nested hierarchy of surface depressions and catchments as well as simulating inundation dynamics. The **lidar** Python package is intended for filling this gap.

## Key Features

-   Smoothing DEMs using mean, median, and Gaussian filters.
-   Extracting depressions from DEMs.
-   Filtering out small artifact depressions based on user-specified minimum depression size.
-   Generating refined DEMs with small depressions filled but larger depressions kept intact.
-   Delineating depression nested hierarchy using the level-set method.
-   Delineating mount nested hierarchy using the level-set method.
-   Computing topological and geometric properties of depressions, including size, volume, mean depth, maximum depth, lowest elevation,
    spill elevation, perimeter, major axis length, minor axis length, elongatedness, eccentricity, orientation, and area-bbox-ratio.
-   Exporting depression properties as a csv file.
