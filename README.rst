=====
lidar
=====


.. image:: https://img.shields.io/pypi/v/lidar.svg
        :target: https://pypi.python.org/pypi/lidar

.. image:: https://img.shields.io/travis/giswqs/lidar.svg
        :target: https://travis-ci.org/giswqs/lidar

.. image:: https://readthedocs.org/projects/lidar/badge/?version=latest
        :target: https://lidar.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status
.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
        :target: https://opensource.org/licenses/MIT


Author: Qiusheng Wu (https://wetlands.io | wqs@binghamton.edu)


**lidar** is a toolset for terrain and hydrological analysis using digital elevation models (DEMs). It is particularly useful for analyzing high-resolution topographic data, such as DEMs derived from Light Detection and Ranging (LiDAR) data.


* GitHub repo: https://github.com/giswqs/lidar
* Documentation: https://lidar.readthedocs.io.
* PyPI: https://pypi.org/project/lidar/
* Free software: `MIT license`_

.. _`MIT license`: https://en.wikipedia.org/wiki/MIT_License


Features
--------

* Smoothing DEMs using mean, median, and Gaussian filters (see filtering.py_)
* Extracting depressions from DEMs (see filling.py_).
* Filtering out small artifact depressions based on user-specified minimum depression size (see filling.py_).
* Generating refined DEMs with small depressions filled but larger depressions kept intact (see filling.py_).
* Delineating depression nested hierarchy using the level-set method (see slicing.py_).
* Delineating mount nested hierarchy using the level-set method (see mounts.py_).
* Computing topological and geometric properties of depressions, including size, volume, mean depth, maximum depth, lowest elevation, spill elevation, perimeter, major axis length, minor axis length, elongatedness, eccentricity, orientation, and area-bbox-ratio (see slicing.py_).
* Exporting depression properties as a csv file (see slicing.py_).


Using It
--------
Install the Python package using the following command:

.. code:: python

  pip install lidar


And use:

.. code:: python

  import os
  import pkg_resources
  import lidar

  # identify the sample data directory of the package
  package_name = 'lidar'
  data_dir = pkg_resources.resource_filename(package_name, 'data/')

  # use the sample dem. Change it to your own dem if needed
  in_dem = os.path.join(data_dir, 'dem.tif')
  # set output directory. By default, use the temp directory under user's home directory
  out_dir = os.path.join(os.path.expanduser("~"), "temp")

  # parameters for identifying sinks and delineating nested depressions
  min_size = 1000      # minimum number of pixels as a depression
  min_depth = 0.5      # minimum depth as a depression
  interval = 0.3       # slicing interval for the level-set method
  bool_shp = False     # output shapefiles for each individual level

  # extracting sinks based on user-defined minimum depression size
  out_dem = os.path.join(out_dir, "median.tif")
  in_dem = MedianFilter(in_dem, kernel_size=3, out_file=out_dem)
  sink_path = lidar.ExtractSinks(in_dem, min_size, out_dir)
  dep_id_path, dep_level_path = lidar.DelineateDepressions(sink_path, min_size, min_depth, interval, out_dir, bool_shp)

Check the example.py_ for more details.


Examples
--------

The images below show working examples of the level set method for delineating nested depressions in the Cottonwood Lake Study Area (CLSA), North Dakota. More test datasets (e.g., the Pipestem watershed in the Prairie Pothole Region of North Dakota) can be downloaded from http://gishub.org/2018-JAWRA-Data

The following example was conducted on a 64-bit Linux machine with a quad-core Intel i7-7700 CPU and 16 GB RAM. The average running time of the algorithm for this DEM was 0.75 seconds.

.. image:: http://spatial.binghamton.edu/pubs/2018-JAWRA/images/CLSA_DEM.jpg
.. image:: http://spatial.binghamton.edu/pubs/2018-JAWRA/images/CLSA_Result.jpg
.. image:: http://spatial.binghamton.edu/pubs/2018-JAWRA/images/CLSA_Table.jpg


Publications
------------
The level-set algorithm used in the **lidar** package has been published in the following articles.

* **Wu, Q.**, Lane, C.R., Wang, L., Vanderhoof, M.K., Christensen, J.R., & Liu, H. (2018). Efficient Delineation of Nested Depression Hierarchy in Digital Elevation Models for Hydrological Analysis Using Level-Set Methods. *Journal of the American Water Resources Association*. (forthcoming)

Credits
-------
* This algorithms are built on richdem_, numpy_, scipy_, scikit-image_, and pygdal_.

* This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _filtering.py: https://github.com/giswqs/lidar/blob/master/lidar/filtering.py
.. _filling.py: https://github.com/giswqs/lidar/blob/master/lidar/filling.py
.. _slicing.py: https://github.com/giswqs/lidar/blob/master/lidar/slicing.py
.. _mounts.py: https://github.com/giswqs/lidar/blob/master/lidar/mounts.py
.. _example.py: https://github.com/giswqs/lidar/blob/master/lidar/example.py
.. _richdem: https://github.com/r-barnes/richdem
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _scikit-image: http://scikit-image.org/
.. _pygdal: https://github.com/nextgis/pygdal
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
