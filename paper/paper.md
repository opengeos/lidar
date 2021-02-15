---
title: "lidar: A Python package for delineating nested surface depressions from digital elevation data"
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
date: 9 February 2021
bibliography: paper.bib
---

# Summary

**lidar** is a Python package for delineating the nested hierarchy of surface depressions in digital elevation models (DEMs).
In traditional hydrological modeling, surface depressions in a DEM are commonly treated as artifacts and thus filled and removed to create a depressionless DEM, which can then be used to generate continuous stream networks. In reality, however, surface depressions in DEMs are commonly a combination of spurious and actual terrain features [@Lindsay2006]. Fine-resolution DEMs derived from Light Detection and Ranging (LiDAR) data can capture and represent actual surface depressions, especially in glaciated [@Wu2016-ub] and karst landscapes [@Wu2016]. During the past decades, various algorithms have been developed to identify and delineate surface depressions, such as depression filling [@Wang2006], depression breaching [@Lindsay2015], hybrid breaching-filling [@Lindsay2016], and contour tree method [@Wu2015]. More recently, a level-set method based on graph theory was proposed to delineate the nested hierarchy of surface depressions [@Wu2019]. The **lidar** Python package implements the level-set method and makes it possible for delineating the nested hierarchy of surface depressions as well as elevated terrain features. It also provides an interactive Graphical User Interface (GUI) that allows users to run the program with minimal coding.

# Statement of Need

The **lidar** package is intended for scientists and researchers who would like to integrate surface depressions into hydrological modeling. It can also facilitate the identification and delineation of depressional features, such as sinkholes, detention basins, and prairie potholes. The detailed topological and geometric properties of surface depressions can be useful for terrain analysis and hydrological modeling, including the size, volume, mean depth, maximum depth, lowest elevation, spill elevation, perimeter, major axis length, minor axis length, elongatedness.

# lidar Functionality

The key functionality of the **lidar** package is organized into several modules:

-   [filtering](https://github.com/giswqs/lidar/blob/master/lidar/filtering.py): Smoothing DEMs using mean, median, and Gaussian filters.
-   [filling](https://github.com/giswqs/lidar/blob/master/lidar/filling.py): Delineating surface depressions from DEMs using the traditional depression filling method.
-   [slicing](https://github.com/giswqs/lidar/blob/master/lidar/slicing.py): Delineating the nested hierarchy of surface depressions using the level-set method; computing topological and geometric properties of depressions; and exporting depression properties as a CSV file.
-   [mounts](https://github.com/giswqs/lidar/blob/master/lidar/mounts.py): Delineating the nested hierarchy of elevated features (i.e., mounts) using the level-set method; computing topological and geometric properties of mounts; and exporting mount properties as a CSV file.
-   [toolbox](https://github.com/giswqs/lidar/blob/master/lidar/toolbox): An [ArcGIS](https://www.esri.com/en-us/arcgis/about-arcgis/overview) toolbox for delineating the nested hierarchy of surface depressions and simulating inundation dynamics.

# lidar Tutorials

The **lidar** Python package has a C library dependency called [GDAL](https://gdal.org/index.html), which can be challenging for some users to install on their computer. Alternatively, users can try out the **lidar** package using just a browser without having to install anything on their computer.

-   Try it out with Binder: <https://gishub.org/lidar-binder>
-   Try it out with Google Colab: <https://gishub.org/lidar-colab>
-   Help documentation: <https://lidar.gishub.org>

# Acknowledgments

The author would like to thank the open-source community, especially the developers of numpy [@Harris2020], scipy [@Virtanen2020], scikit-image [@Van_der_Walt2014], matplotlib [@Hunter2007], and richDEM [@Barnes2018]. These open-source packages empower the **lidar** Python package.

# References
