"""Module for filling surface depressions.

"""

import os
import time
import numpy as np
import richdem as rd
from scipy import ndimage
from skimage import measure
from osgeo import gdal, ogr, osr


class Depression:
    """The class for storing depression info."""

    def __init__(
        self,
        id,
        count,
        size,
        volume,
        meanDepth,
        maxDepth,
        minElev,
        bndElev,
        perimeter,
        major_axis,
        minor_axis,
        elongatedness,
        eccentricity,
        orientation,
        area_bbox_ratio,
    ):
        self.id = id
        self.count = count
        self.size = size
        self.volume = volume
        self.meanDepth = meanDepth
        self.maxDepth = maxDepth
        self.minElev = minElev
        self.bndElev = bndElev
        self.perimeter = perimeter
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.elongatedness = elongatedness
        self.eccentricity = eccentricity
        self.orientation = orientation
        self.area_bbox_ratio = area_bbox_ratio


def regionGroup(img_array, min_size, no_data):
    """IdentifIies regions based on region growing method

    Args:
        img_array (np.array): The numpy array containing the image.
        min_size (int): The minimum number of pixels to be considered as a depression.
        no_data (float): The no_data value of the image.

    Returns:
        tuple: The labelled objects and total number of labels.
    """
    img_array[img_array == no_data] = 0
    label_objects, nb_labels = ndimage.label(img_array)
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = sizes > min_size
    mask_sizes[0] = 0
    image_cleaned = mask_sizes[label_objects]
    label_objects, nb_labels = ndimage.label(image_cleaned)
    # nb_labels is the total number of objects. 0 represents background object.
    return label_objects, nb_labels


def np2rdarray(in_array, no_data, projection, geotransform):
    """Converts an numpy array to rdarray.

    Args:
        in_array (np.array): The input numpy array.
        no_data (float): The no_data value of the array.
        projection (str): The projection of the image.
        geotransform (str): The geotransform of the image.

    Returns:
        object: The richDEM array.
    """
    out_array = rd.rdarray(in_array, no_data=no_data)
    out_array.projection = projection
    out_array.geotransform = geotransform
    return out_array


def get_dep_props(objects, resolution):
    """Computes depression attributes.

    Args:
        objects (object): The labeled objects.
        resolution (float): The spatial reoslution of the image.

    Returns:
        list: A list of depression objects with attributes.
    """
    dep_list = []

    for object in objects:
        unique_id = object.label
        count = object.area
        size = count * pow(resolution, 2)  # depression size
        min_elev = np.float(object.min_intensity)  # depression min elevation
        max_elev = np.float(object.max_intensity)  # depression max elevation
        max_depth = max_elev - min_elev  # depression max depth
        mean_depth = np.float(
            (max_elev * count - np.sum(object.intensity_image)) / count
        )  # depression mean depth
        volume = mean_depth * count * pow(resolution, 2)  # depression volume
        perimeter = object.perimeter * resolution
        major_axis = object.major_axis_length * resolution
        minor_axis = object.minor_axis_length * resolution
        if minor_axis == 0:
            minor_axis = resolution
        elongatedness = major_axis * 1.0 / minor_axis
        eccentricity = object.eccentricity
        orientation = object.orientation / 3.1415 * 180
        area_bbox_ratio = object.extent

        dep_list.append(
            Depression(
                unique_id,
                count,
                size,
                volume,
                mean_depth,
                max_depth,
                min_elev,
                max_elev,
                perimeter,
                major_axis,
                minor_axis,
                elongatedness,
                eccentricity,
                orientation,
                area_bbox_ratio,
            )
        )

    return dep_list


def write_dep_csv(dep_list, csv_file):
    """Saves the depression list info to a CSV file.

    Args:
        dep_list (list): A list of depression objects with attributes.
        csv_file (str): File path to the output CSV file.
    """
    csv = open(csv_file, "w")
    header = (
        "region-id"
        + ","
        + "count"
        + ","
        + "area"
        + ","
        + "volume"
        + ","
        + "avg-depth"
        + ","
        + "max-depth"
        + ","
        + "min-elev"
        + ","
        + "max-elev"
        + ","
        + "perimeter"
        + ","
        + "major-axis"
        + ","
        + "minor-axis"
        + ","
        + "elongatedness"
        + ","
        + "eccentricity"
        + ","
        + "orientation"
        + ","
        + "area-bbox-ratio"
    )

    csv.write(header + "\n")
    for dep in dep_list:
        line = "{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}, {:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}".format(
            dep.id,
            dep.count,
            dep.size,
            dep.volume,
            dep.meanDepth,
            dep.maxDepth,
            dep.minElev,
            dep.bndElev,
            dep.perimeter,
            dep.major_axis,
            dep.minor_axis,
            dep.elongatedness,
            dep.eccentricity,
            dep.orientation,
            dep.area_bbox_ratio,
        )
        csv.write(line + "\n")
    csv.close()


def polygonize(img, shp_path):
    """Converts a raster image to vector.

    Args:
        img (str): File path to the input image.
        shp_path (str): File path to the output shapefile.
    """
    # mapping between gdal type and ogr field type
    type_mapping = {
        gdal.GDT_Byte: ogr.OFTInteger,
        gdal.GDT_UInt16: ogr.OFTInteger,
        gdal.GDT_Int16: ogr.OFTInteger,
        gdal.GDT_UInt32: ogr.OFTInteger,
        gdal.GDT_Int32: ogr.OFTInteger,
        gdal.GDT_Float32: ogr.OFTReal,
        gdal.GDT_Float64: ogr.OFTReal,
        gdal.GDT_CInt16: ogr.OFTInteger,
        gdal.GDT_CInt32: ogr.OFTInteger,
        gdal.GDT_CFloat32: ogr.OFTReal,
        gdal.GDT_CFloat64: ogr.OFTReal,
    }

    ds = gdal.Open(img)
    prj = ds.GetProjection()
    srcband = ds.GetRasterBand(1)

    dst_layername = "Shape"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(shp_path)
    srs = osr.SpatialReference(wkt=prj)

    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
    # raster_field = ogr.FieldDefn('id', type_mapping[srcband.DataType])
    raster_field = ogr.FieldDefn("id", type_mapping[gdal.GDT_Int32])
    dst_layer.CreateField(raster_field)
    gdal.Polygonize(srcband, srcband, dst_layer, 0, [], callback=None)
    del img, ds, srcband, dst_ds, dst_layer


def ExtractSinks(in_dem, min_size, out_dir):
    """Extract sinks (e.g., maximum depression extent) from a DEM.

    Args:
        in_dem (str): File path to the input DEM.
        min_size (int): The minimum number of pixels to be considered as a sink.
        out_dir (str): File path to the output directory.

    Returns:
        object: The richDEM array containing sinks.
    """
    start_time = time.time()

    out_dem = os.path.join(out_dir, "dem.tif")
    out_dem_filled = os.path.join(out_dir, "dem_filled.tif")
    out_dem_diff = os.path.join(out_dir, "dem_diff.tif")
    out_sink = os.path.join(out_dir, "sink.tif")
    out_region = os.path.join(out_dir, "region.tif")
    out_depth = os.path.join(out_dir, "depth.tif")
    out_csv_file = os.path.join(out_dir, "regions_info.csv")
    out_vec_file = os.path.join(out_dir, "regions.shp")

    # create output folder if nonexistent
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # load the dem and get dem info
    print("Loading data ...")
    dem = rd.LoadGDAL(in_dem)
    no_data = dem.no_data
    projection = dem.projection
    geotransform = dem.geotransform
    cell_size = np.round(geotransform[1], decimals=2)

    # get min and max elevation of the dem
    max_elev = np.float(np.max(dem))
    min_elev = np.float(np.min(dem[dem > 0]))
    print(
        "min = {:.2f}, max = {:.2f}, no_data = {}, cell_size = {}".format(
            min_elev, max_elev, no_data, cell_size
        )
    )

    # depression filling
    print("Depression filling ...")
    dem_filled = rd.FillDepressions(dem, in_place=False)
    dem_diff = dem_filled - dem
    dem_diff.no_data = 0

    print("Saving filled dem ...")
    rd.SaveGDAL(out_dem_filled, dem_filled)
    rd.SaveGDAL(out_dem_diff, dem_diff)

    # nb_labels is the total number of objects. 0 represents background object.
    print("Region grouping ...")
    label_objects, nb_labels = regionGroup(dem_diff, min_size, no_data)
    dem_diff[label_objects == 0] = 0
    depth = np2rdarray(
        dem_diff, no_data=0, projection=projection, geotransform=geotransform
    )
    rd.SaveGDAL(out_depth, depth)
    del dem_diff, depth

    print("Computing properties ...")
    # objects = measure.regionprops(label_objects, dem, coordinates='xy')
    objects = measure.regionprops(label_objects, dem)
    dep_list = get_dep_props(objects, cell_size)
    write_dep_csv(dep_list, out_csv_file)
    del objects, dep_list

    # convert numpy to richdem data format
    region = np2rdarray(
        label_objects, no_data=0, projection=projection, geotransform=geotransform
    )
    del label_objects

    print("Saving sink dem ...")
    sink = np.copy(dem)
    sink[region == 0] = 0
    sink = np2rdarray(sink, no_data=0, projection=projection, geotransform=geotransform)
    rd.SaveGDAL(out_sink, sink)
    # del sink

    print("Saving refined dem ...")
    dem_refined = dem_filled
    dem_refined[region > 0] = dem[region > 0]
    dem_refined = np2rdarray(
        dem_refined, no_data=no_data, projection=projection, geotransform=geotransform
    )
    rd.SaveGDAL(out_dem, dem_refined)
    rd.SaveGDAL(out_region, region)
    del dem_refined, region, dem

    print("Converting raster to vector ...")
    polygonize(out_region, out_vec_file)

    end_time = time.time()
    print("Total run time:\t\t\t {:.4f} s\n".format(end_time - start_time))

    return out_sink
