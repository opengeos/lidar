Launch the interactive notebook tutorial for the **lidar** Python
package with **mybinder.org** or **binder.pangeo.io** now:

[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/lidar-cloud)

[![image](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/giswqs/lidar/master?filepath=examples%2Flidar.ipynb)

## A Quick Example

```python
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
bool_shp = True     # output shapefiles for each individual level

# extracting sinks based on user-defined minimum depression size
out_dem = os.path.join(out_dir, "median.tif")
in_dem = MedianFilter(in_dem, kernel_size=3, out_file=out_dem)
sink_path = ExtractSinks(in_dem, min_size, out_dir)
dep_id_path, dep_level_path = DelineateDepressions(sink_path, min_size, min_depth, interval, out_dir, bool_shp)

print('Results are saved in: {}'.format(out_dir))
```

## An Interactive Notebook Tutorial

This tutorial can be accessed in three ways:

-   HTML version: <https://gishub.org/lidar-html>
-   Viewable Notebook: <https://gishub.org/lidar-notebook>
-   Interactive Notebook: <https://gishub.org/lidar-cloud>

Launch this tutorial as an interactive Jupyter Notebook on the cloud -
<https://gishub.org/lidar-cloud>.

![image](https://i.imgur.com/aIttPVG.gif)

## lidar GUI

**lidar** also provides a Graphical User Interface (GUI), which can be
invoked using the following Python script:

```python
import lidar
lidar.gui()
```

![image](https://i.imgur.com/eSjcSs9.png)

