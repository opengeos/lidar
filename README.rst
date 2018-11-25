=====
lidar
=====

.. image:: https://mybinder.org/badge_logo.svg 
        :target: https://gishub.org/lidar-cloud

.. image:: https://binder.pangeo.io/badge.svg 
        :target: https://binder.pangeo.io/v2/gh/giswqs/lidar/master?filepath=examples%2Flidar.ipynb

.. image:: https://img.shields.io/pypi/v/lidar.svg
        :target: https://pypi.python.org/pypi/lidar

.. image:: https://pepy.tech/badge/lidar
        :target: https://pepy.tech/project/lidar

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
* Binder: https://gishub.org/lidar-cloud
* Free software: `MIT license`_

.. _`MIT license`: https://en.wikipedia.org/wiki/MIT_License


**Contents**

- `Features`_
- `Installation`_
- `Tutorials`_
- `lidar GUI`_
- `Dependencies`_
- `Examples`_
- `References`_
- `Reporting Bugs`_
- `Credits`_


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


Installation
------------
**lidar** supports a variety of platforms, including Microsoft Windows, macOS, and Linux operating systems. Note that you will need to have **Python 3.x** installed. Python 2.x is not supported. The **lidar** Python package can be installed using the following command. If you encounter any errors, please check the Dependencies_ section below.

.. code:: python

  pip install lidar


If you have installed **lidar** before and want to upgrade to the latest version, you can use the following command:

.. code:: python

  pip install lidar -U



Tutorials
---------

Launch the interactive notebook tutorial for the **lidar** Python package with **mybinder.org** or **binder.pangeo.io** now:

.. image:: https://mybinder.org/badge_logo.svg 
        :target: https://gishub.org/lidar-cloud

.. image:: https://binder.pangeo.io/badge.svg 
        :target: https://binder.pangeo.io/v2/gh/giswqs/lidar/master?filepath=examples%2Flidar.ipynb


A Quick Example
===============

.. code:: python

  import os
  import pkg_resources
  from lidar import *

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
  sink_path = ExtractSinks(in_dem, min_size, out_dir)
  dep_id_path, dep_level_path = DelineateDepressions(sink_path, min_size, min_depth, interval, out_dir, bool_shp)

  print('Results are saved in: {}'.format(out_dir))


Check the example.py_ for more details.


An Interactive Jupyter Notebook Tutorial
========================================

This tutorial can be accessed in three ways:

- HTML version: https://gishub.org/whitebox-html
- Viewable Notebook: https://gishub.org/lidar-notebook
- Interactive Notebook: https://gishub.org/lidar-cloud

Launch this tutorial as an interactive Jupyter Notebook on the cloud - https://gishub.org/lidar-cloud.

.. image:: https://i.imgur.com/aIttPVG.gif


lidar GUI
---------

**lidar** also provides a Graphical User Interface (GUI), which can be invoked using the following Python script:

.. code:: python

  import lidar
  lidar.gui()


.. image:: https://i.imgur.com/eSjcSs9.png


Dependencies
------------

lidar's Python dependencies are listed in its requirements.txt file. In addition, lidar has a C library dependency: GDAL >=1.11.2. How to install GDAL in different operating systems will be explained below. More informaton about GDAL can be found here_.

It is highly recommended that you use a Python virtual environment (e.g., conda) to test the lidar package. Please follow the `conda user guide`_ to install conda if necessary. Once you have conda installed, you can use Terminal or an Anaconda Prompt to create a Python virtual environment. Check `managing Python environment`_ for more information.

.. _here: https://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries
.. _`conda user guide`: https://conda.io/docs/user-guide/install/index.html
.. _`managing Python environment`: https://conda.io/docs/user-guide/tasks/manage-environments.html

Once GDAL has been installed, you can then proceed to install the **lidar** Python package using the following command:

.. code:: python

  pip install lidar


Linux
=====

Debian-based Linux
^^^^^^^^^^^^^^^^^^

The following commands can be used to install GDAL for Debian-based Linux distributions (e.g., Ubuntu, Linux Mint).

.. code:: python

  sudo add-apt-repository ppa:ubuntugis/ppa
  sudo apt-get update
  sudo apt-get install gdal-bin libgdal-dev
  pip install lidar


If you encounter any compiling errors, try the following commands. 

.. code:: python

  sudo apt-get install --reinstall build-essential
  sudo apt-get install python3-dev
  pip install wheel


Pacman-based Linux
^^^^^^^^^^^^^^^^^^

The following commands can be used to install GDAL for Pacman-based Linux distributions (e.g., Arch Linux, Manjaro). You might need to use **sudo** if you encounter permission errors.

.. code:: python

  sudo pacman -S yaourt --noconfirm
  yaourt -S gdal --noconfirm
  yaourt -S python-gdal --noconfirm
  pip install lidar


MacOS X
=======
For a Homebrew based Python environment, do the following.

.. code:: python

  brew update
  brew install gdal

Alternatively, you can install GDAL binaries from kyngchaos_. You will then need to add the installed location ``/Library/Frameworks/GDAL.framework/Programs`` to your system path.

.. _kyngchaos: http://www.kyngchaos.com/software/frameworks#gdal_complete


Windows
=======

The instruction below assumes that you have installed Anaconda_. Open **Anaconda Prompt** and enter the following commands to create a conda environment and install required packages

.. code:: python

  conda create -n py36 python=3.6
  activate py36
  conda install -c conda-forge gdal 
  pip install richdem
  pip install lidar

When installing the **richdem** package, if you encounter an error saying 'Microsoft Visual C++ 14.0 is required', please follow the steps below to fix the error and reinstall **richdem**. More infomration can be found at this link `Fix Python 3 on Windows error - Microsoft Visual C++ 14.0 is required`_.  

* Download `Microsoft Build Tools for Visual Studio 2017`_
* Double click to install the downloaded installer - **Microsoft Build Tools for Visual Studio 2017**.
* Open **Microsoft Build Tools for Visual Studio 2017**
* Select **Workloads --> Visual C++ build tools** and click the install button

.. _Anaconda: https://www.anaconda.com/download
.. _`Fix Python 3 on Windows error - Microsoft Visual C++ 14.0 is required`: https://www.scivision.co/python-windows-visual-c++-14-required/
.. _`Microsoft Build Tools for Visual Studio 2017`: https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15


Examples
--------

The images below show working examples of the level set method for delineating nested depressions in the Cottonwood Lake Study Area (CLSA), North Dakota. More test datasets (e.g., the Pipestem watershed in the Prairie Pothole Region of North Dakota) can be downloaded from http://gishub.org/2018-JAWRA-Data

The following example was conducted on a 64-bit Linux machine with a quad-core Intel i7-7700 CPU and 16 GB RAM. The average running time of the algorithm for this DEM was 0.75 seconds.

.. image:: https://wetlands.io/file/images/CLSA_DEM.jpg
.. image:: https://wetlands.io/file/images/CLSA_Result.jpg
.. image:: https://wetlands.io/file/images/CLSA_Table.jpg


References
----------
The level-set algorithm in the **lidar** package has been published in the following article:

* **Wu, Q.**, Lane, C.R., Wang, L., Vanderhoof, M.K., Christensen, J.R., & Liu, H. (2018). Efficient Delineation of Nested Depression Hierarchy in Digital Elevation Models for Hydrological Analysis Using Level-Set Method. *Journal of the American Water Resources Association*. DOI: `10.1111/1752-1688.12689`_ (in press) preprint_

Applications of the level-set and contour-tree methods for feature extraction from LiDAR data:

* **Wu, Q.**, & Lane, C.R. (2017). Delineating wetland catchments and modeling hydrologic connectivity using LiDAR data and aerial imagery. *Hydrology and Earth System Sciences*. 21: 3579-3595. DOI: `10.5194/hess-21-3579-2017`_
* **Wu, Q.**, Deng, C., & Chen, Z. (2016). Automated delineation of karst sinkholes from LiDAR-derived digital elevation models. *Geomorphology*. 266: 1-10. DOI: `10.1016/j.geomorph.2016.05.006`_
* **Wu, Q.**, Su, H., Sherman, D.J., Liu, H., Wozencraft, J.M., Yu, B., & Chen, Z. (2016). A graph-based approach for assessing storm-induced coastal changes. *International Journal of Remote Sensing*. 37:4854-4873. DOI: `10.1080/01431161.2016.1225180`_
* **Wu, Q.**, & Lane, C.R. (2016). Delineation and quantification of wetland depressions in the Prairie Pothole Region of North Dakota. *Wetlands*. 36(2):215–227. DOI: `10.1007/s13157-015-0731-6`_
* **Wu, Q.**, Liu, H., Wang, S., Yu, B., Beck, R., & Hinkel, K. (2015). A localized contour tree method for deriving geometric and topological properties of complex surface depressions based on high-resolution topographic data. *International Journal of Geographical Information Science*. 29(12): 2041-2060. DOI: `10.1080/13658816.2015.1038719`_
* **Wu, Q.**, Lane, C.R., & Liu, H. (2014). An effective method for detecting potential woodland vernal pools using high-resolution LiDAR data and aerial imagery. *Remote Sensing*. 6(11):11444-11467.  DOI: `10.3390/rs61111444`_


Reporting Bugs
--------------
Report bugs at https://github.com/giswqs/lidar/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.


Credits
-------
* The algorithms are built on richdem_, numpy_, scipy_, scikit-image_, and pygdal_.

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
.. _`10.5194/hess-21-3579-2017`: https://doi.org/10.5194/hess-21-3579-2017
.. _`10.1016/j.geomorph.2016.05.006`: http://dx.doi.org/10.1016/j.geomorph.2016.05.006
.. _`10.1007/s13157-015-0731-6`: http://dx.doi.org/10.1007/s13157-015-0731-6
.. _`10.1080/13658816.2015.1038719`: http://dx.doi.org/10.1080/13658816.2015.1038719
.. _`10.1080/01431161.2016.1225180`: http://dx.doi.org/10.1080/01431161.2016.1225180
.. _`10.3390/rs61111444`: http://dx.doi.org/10.3390/rs61111444
.. _preprint: https://www.preprints.org/manuscript/201808.0358/v1
.. _`10.1111/1752-1688.12689`: https://doi.org/10.1111/1752-1688.12689
