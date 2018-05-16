=====
Usage
=====

To use lidar in a project:

.. code:: python

  import os
  import pkg_resources
  import lidar
  import richdem as rd

  # identify the sample data directory of the package
  package_name = 'lidar'
  data_dir = pkg_resources.resource_filename(package_name, 'data/')

  # use the sample dem. Change it to your own dem if needed
  in_dem = os.path.join(data_dir, 'dem.tif')
  # set output directory. By default, use the temp directory under user's home directory
  out_dir = os.path.join(os.path.expanduser("~"), "temp")

  # parameters for identifying sinks and delineating nested depressions
  min_size = 1000             # minimum number of pixels as a depression
  min_depth = 0.3             # minimum depth as a depression
  interval = 0.3      # slicing interval for the level-set method
  bool_shp = False      # output shapefiles for each individual level

  # extracting sinks based on user-defined minimum depression size
  sink_path = lidar.ExtractSinks(in_dem, min_size, out_dir)
  dep_id_path, dep_level_path = lidar.DelineateDepressions(sink_path, min_size, min_depth, interval, out_dir, bool_shp)

  # loading data and results
  dem = rd.LoadGDAL(in_dem)
  sink = rd.LoadGDAL(sink_path)
  dep_id = rd.LoadGDAL(dep_id_path)
  dep_level = rd.LoadGDAL(dep_level_path)

  # plotting results
  dem_fig = rd.rdShow(dem, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
  sink_fig = rd.rdShow(sink, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
  dep_id_fig = rd.rdShow(dep_id, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
  dep_level_path = rd.rdShow(dep_level, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))

Check the example.py_ for more details.
