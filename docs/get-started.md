Launch the interactive notebook tutorial for the **lidar** Python
package with **Google Colab** now:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/giswqs/lidar/blob/master/examples/lidar_colab.ipynb)

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
# set the output directory
out_dir = os.getcwd()

# parameters for identifying sinks and delineating nested depressions
min_size = 1000      # minimum number of pixels as a depression
min_depth = 0.5      # minimum depth as a depression
interval = 0.3       # slicing interval for the level-set method
bool_shp = True      # output shapefiles for each individual level

# extracting sinks based on user-defined minimum depression size
out_dem = os.path.join(out_dir, "median.tif")
in_dem = MedianFilter(in_dem, kernel_size=3, out_file=out_dem)
sink_path = ExtractSinks(in_dem, min_size, out_dir)
dep_id_path, dep_level_path = DelineateDepressions(sink_path,
                                                   min_size,
                                                   min_depth,
                                                   interval,
                                                   out_dir,
                                                   bool_shp)
print('Results are saved in: {}'.format(out_dir))
```
## lidar GUI

**lidar** also provides a Graphical User Interface (GUI), which can be
invoked using the following Python script:

```python
import lidar
lidar.gui()
```

![image](https://i.imgur.com/6hLGeV5.png)


## ArcGIS Toolbox

### Toolbox interface

![toolbox](https://raw.githubusercontent.com/giswqs/lidar/master/images/toolbox_0.png)

![toolbox_ui](https://raw.githubusercontent.com/giswqs/lidar/master/images/toolbox_ui.png)

### Video tutorials

[**Delineating nested surface depressions and catchments using ArcGIS Pro**](https://youtu.be/PpF8sfvCATE)

[![demo](http://img.youtube.com/vi/W9PFHNV3cT0/0.jpg)](http://www.youtube.com/watch?v=W9PFHNV3cT0)

[**Delineating nested surface depressions and catchments using ArcMap**](https://youtu.be/PpF8sfvCATE)

[![demo](http://img.youtube.com/vi/PpF8sfvCATE/0.jpg)](http://www.youtube.com/watch?v=PpF8sfvCATE)


