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


Author: Qiusheng Wu (https://wetlands.io | wqs@binghamton.edu)


**lidar** is a toolset for terrain and hydrological analysis using digital elevation models (DEMs). It is particularly useful for analyzing high-resolution topographic data, such as DEMs derived from Light Detection and Ranging (LiDAR) data.


* GitHub repo: https://github.com/giswqs/lidar
* Documentation: https://lidar.readthedocs.io.
* PyPI: https://pypi.org/project/lidar/
* Free software: MIT license



Features
--------

* Extracting depressions from DEMs (see filling.py_).
* Filtering out small artifact depressions based on user-specified minimum depression size.
* Delineating depression nested hierarchy using level-set method (see slicing.py_).
* Computing topological and geometric properties of depressions, including size, volume, depth, spill elevation etc.


Using It
--------
Install the Python package with: ``pip install lidar``


And use:

     ``import lidar``

     ``lidar.ExtractSinks(in_dem, min_size, out_dir)``

     ``lidar.DelineateDepressions(sink_path, min_size, min_depth, interval, out_dir, bool_shp)``

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

.. _filling.py: https://github.com/giswqs/lidar/blob/master/lidar/filling.py
.. _slicing.py: https://github.com/giswqs/lidar/blob/master/lidar/slicing.py
.. _example.py: https://github.com/giswqs/lidar/blob/master/lidar/example.py
.. _richdem: https://github.com/r-barnes/richdem
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _scikit-image: http://scikit-image.org/
.. _pygdal: https://github.com/nextgis/pygdal
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
