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

* Extracting depressions from DEMs.
* Filtering out small artifact depressions based on user-specified minimum depression size.
* Delineating depression nested hierarchy using level-set method.
* Computing topological and geometric properties of depressions, including size, volume, depth, spill elevation etc.


Using It
--------
Install the Python package with: ``pip install lidar``


And use:

     ``import lidar``

     ``lidar.ExtractSinks(in_dem, min_size, out_dir)``
     
     ``lidar.DelineateDepressions(sink_path, min_size, min_depth, interval, out_dir, bool_shp)``

Check the example.py_ for more details.

Credits
-------
* This algorithms are built on richdem_, numpy_, scipy_, scikit-image_, and pygdal_.

* This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.


.. _example.py: https://github.com/giswqs/lidar/blob/master/lidar/example.py
.. _richdem: https://github.com/r-barnes/richdem
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _scikit-image: http://scikit-image.org/
.. _pygdal: https://github.com/nextgis/pygdal
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
