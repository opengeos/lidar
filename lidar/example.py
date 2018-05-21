import os
import pkg_resources
import richdem as rd
from filtering import MedianFilter
from filling import ExtractSinks
from slicing import DelineateDepressions

# identify the sample data directory of the package
package_name = 'lidar'
data_dir = pkg_resources.resource_filename(package_name, 'data/')

# use the sample dem. Change it to your own dem if needed
in_dem = os.path.join(data_dir, 'dem.tif')
# set output directory. By default, use the temp directory under user's home directory
out_dir = os.path.join(os.path.expanduser("~"), "temp")

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# parameters for identifying sinks and delineating nested depressions
min_size = 1000         # minimum number of pixels as a depression
min_depth = 0.3         # minimum depth as a depression
interval = 0.3          # slicing interval for the level-set method
bool_shp = False        # output shapefiles for each individual level

# extracting sinks based on user-defined minimum depression size
out_dem = os.path.join(out_dir, "median.tif")
in_dem = MedianFilter(in_dem, kernel_size=3, out_file=out_dem)
sink_path = ExtractSinks(in_dem, min_size, out_dir)
dep_id_path, dep_level_path = DelineateDepressions(sink_path, min_size, min_depth, interval, out_dir, bool_shp)

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

