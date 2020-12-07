---
title: "lidar: A Python package for terrain and hydrological analysis using digital elevation models"
tags:
    - Python
    - Terrain analysis
    - Hydrological analysis
    - Surface depressions
    - Digital elevation models
    - DEM
    - LiDAR
    - Jupyter notebook
authors:
    - name: Qiusheng Wu
      orcid: 0000-0001-5437-4073
      affiliation: "1"
affiliations:
    - name: Department of Geography, University of Tennessee, Knoxville, TN 37996, United States
      index: 1
date: 7 December 2020
bibliography: paper.bib
---

# Summary

**lidar** is a Python package for terrain and hydrological analysis using digital elevation models (DEMs).
In traditional hydrological modeling, surface depressions in a DEM are commonly treated as artifacts and thus filled and removed to create a depressionless DEM, which can then be used to generate continuous stream networks. In reality, however, surface depressions in DEMs are commonly a combination of spurious and actual terrain features [@Lindsay2006]. Fine-resolution DEMs derived from Light Detection and Ranging (LiDAR) data can capture and represent actual surface depressions, especially in glaciated and karst landscapes [@Wu2016, @Wu2016-ub]. During the past decades, various algorithms have been developed to identify and delineate surface depressions, such as depression filling [@Wang2006], depression breaching [@Lindsay2015], hybrid breaching-filling [@Lindsay2016], and contour tree method [@Wu2015]. More recently, a level-set method based on graph theory was proposed to delineate the nested hierarchy of surface depressions [Wu2019]. The **lidar** Python package implements the level-set method and makes it possible for delineating the nested hierarchy of surface depressions as well as elevated terrain features. It also provides an interactive Graphical User Interface (GUI) that allows users to run the program with minimal coding.

# lidar Audience

The **lidar** package is intended for scientists and researchers who would like to integrate surface depressions into hydrological modeling. It can also facilitate the identification and delineation of depressional features, such as sinkholes, detention basins, and prairie potholes. The detailed topological and geometric properties of surface depressions can be useful for terrain analysis and hydrological modeling, including the size, volume, mean depth, maximum depth, lowest elevation, spill elevation, perimeter, major axis length, minor axis length, elongatedness.

# lidar Functionality

The key functionality of the **lidar** package is organized into several modules:

-   [filtering](https://github.com/giswqs/lidar/blob/master/lidar/filtering.py): Smoothing DEMs using mean, median, and Gaussian filters.
-   [filling](https://github.com/giswqs/lidar/blob/master/lidar/filling.py): Delineating surface depressions from DEMs using the traditional depression filling method.
-   [slicing](https://github.com/giswqs/lidar/blob/master/lidar/slicing.py): Delineating the nested hierarchy of surface depressions using the level-set method; computing topological and geometric properties of depressions; and exporting depression properties as a CSV file.
-   [mounts](https://github.com/giswqs/lidar/blob/master/lidar/mounts.py): Delineating the nested hierarchy of elevated features (i.e., mounts) using the level-set method; computing topological and geometric properties of mounts; and exporting mount properties as a CSV file.

# lidar Tutorials

The **lidar** Python package has a C library dependency called [GDAL](https://gdal.org/index.html), which can be challenging for some users to install on their local computer. Alternatively, users can try out the **lidar** package using a browser with Binder or Colab without having to install anything on their computer.

-   Try it out with Binder: <https://gishub.org/lidar-binder>
-   Try it out with Google Colab: <https://gishub.org/lidar-colab>
-   Help documentation: <https://lidar.readthedocs.io>

# Acknowledgments

The author would like to thank the open-source community, especially the developers of the numpy and scikit-image packages, which empower the **lidar** Python package.

# References
