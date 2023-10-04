"""Module for filling surface depressions.

"""

import os
import time
import numpy as np
import richdem as rd
from scipy import ndimage
from skimage import measure
from osgeo import gdal, ogr, osr
from .common import *
from .filtering import *


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
        img_array (array): The numpy array containing the image.
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
        in_array (array): The input numpy array.
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
        min_elev = float(object.min_intensity)  # depression min elevation
        max_elev = float(object.max_intensity)  # depression max elevation
        max_depth = max_elev - min_elev  # depression max depth
        mean_depth = float(
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


def ExtractSinks(in_dem, min_size, out_dir, filled_dem=None):
    """Extract sinks (e.g., maximum depression extent) from a DEM.

    Args:
        in_dem (str): File path to the input DEM.
        min_size (int): The minimum number of pixels to be considered as a sink.
        out_dir (str): File path to the output directory.
        fill_dem (str, optional): The filled DEM.

    Returns:
        object: The richDEM array containing sinks.
    """
    start_time = time.time()

    out_dem = os.path.join(out_dir, "dem.tif")
    out_dem_diff = os.path.join(out_dir, "dem_diff.tif")
    out_sink = os.path.join(out_dir, "sink.tif")
    out_region = os.path.join(out_dir, "region.tif")
    out_depth = os.path.join(out_dir, "depth.tif")
    out_csv_file = os.path.join(out_dir, "regions_info.csv")
    out_vec_file = os.path.join(out_dir, "regions.shp")
    if filled_dem is None:
        out_dem_filled = os.path.join(out_dir, "dem_filled.tif")
    else:
        out_dem_filled = filled_dem
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
    max_elev = float(np.max(dem[dem != no_data]))
    min_elev = float(np.min(dem[dem > 0]))
    print(
        "min = {:.2f}, max = {:.2f}, no_data = {}, cell_size = {}".format(
            min_elev, max_elev, no_data, cell_size
        )
    )

    # depression filling
    if filled_dem is None:
        print("Depression filling ...")
        dem_filled = rd.FillDepressions(dem, in_place=False)
    else:
        dem_filled = rd.LoadGDAL(filled_dem)
    dem_diff = dem_filled - dem
    dem_diff.no_data = 0

    if filled_dem is None:
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
    sink = np2rdarray(sink, no_data=0, projection=projection,
                      geotransform=geotransform)
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


def extract_sinks_by_bbox(
    bbox,
    filename,
    min_size=10,
    tmp_dir=None,
    mask=None,
    crs="EPSG:5070",
    kernel_size=3,
    resolution=10,
    to_cog=False,
    keep_files=True,
    ignore_warnings=True,
):
    """Extract sinks from a DEM by HUC8.

    Args:
        bbox (list): The bounding box in the form of [minx, miny, maxx, maxy].
        filename (str, optional): The output depression file name.
        min_size (int, optional): The minimum number of pixels to be considered as a sink. Defaults to 10.
        tmp_dir (str, optional): The temporary directory. Defaults to None, e.g., using the current directory.
        mask (str, optional): The mask file path. Defaults to None.
        crs (str, optional): The coordinate reference system. Defaults to "EPSG:5070".
        kernel_size (int, optional): The kernel size for smoothing the DEM. Defaults to 3.
        resolution (int, optional): The resolution of the DEM. Defaults to 10.
        to_cog (bool, optional): Whether to convert the output to COG. Defaults to False.
        keep_files (bool, optional): Whether to keep the intermediate files. Defaults to True.
        ignore_warnings (bool, optional): Whether to ignore warnings. Defaults to True.
    """
    import shutil
    import warnings

    if ignore_warnings:
        warnings.filterwarnings("ignore")

    start_time = time.time()

    if not filename.endswith(".shp"):
        filename = filename + ".shp"

    filename = os.path.abspath(filename)

    if tmp_dir is None:
        tmp_dir = os.path.join(os.getcwd(), "tmp")

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    merge = os.path.join(tmp_dir, "mosaic.tif")
    clip = os.path.join(tmp_dir, "clip.tif")
    reproj = os.path.join(tmp_dir, "reproj.tif")
    image = os.path.join(tmp_dir, "image.tif")
    median = os.path.join(tmp_dir, "median.tif")
    regions = os.path.join(tmp_dir, "regions.shp")
    regions_info = os.path.join(tmp_dir, "regions_info.csv")

    try:
        download_ned_by_bbox(bbox, out_dir=tmp_dir)

        if not os.path.exists(merge):
            print("Merging NED tiles ...")
            mosaic(tmp_dir, merge)

        if mask is not None:
            clip_image(merge, mask, clip)
        else:
            clip = merge

        reproject_image(clip, reproj, crs)
        resample(reproj, image, resolution)
        MedianFilter(image, kernel_size, median)
        if to_cog:
            image_to_cog(median, median)
        ExtractSinks(median, min_size, tmp_dir)
        join_tables(regions, regions_info, filename)

        for file in [merge, clip, reproj, image]:
            if os.path.exists(file):
                os.remove(file)

        if not keep_files:
            shutil.rmtree(tmp_dir)
    except Exception as e:
        print(e)
        return None

    end_time = time.time()
    print("Total run time:\t\t\t {:.4f} s\n".format(end_time - start_time))


def extract_sinks_by_huc8(
    huc8,
    min_size=10,
    filename=None,
    tmp_dir=None,
    wbd=None,
    crs="EPSG:5070",
    kernel_size=3,
    resolution=10,
    keep_files=True,
    error_file=None,
    ignore_warnings=True,
):
    """Extract sinks from a DEM by HUC8.

    Args:
        huc8 (str): The HUC8 code, e.g., 01010002
        min_size (int, optional): The minimum number of pixels to be considered as a sink. Defaults to 10.
        filename (str, optional): The output depression file name. Defaults to None, e,g., using the HUC8 code.
        tmp_dir (str, optional): The temporary directory. Defaults to None, e.g., using the current directory.
        wbd (str | GeoDataFrame, optional): The watershed boundary file. Defaults to None.
        crs (str, optional): The coordinate reference system. Defaults to "EPSG:5070".
        kernel_size (int, optional): The kernel size for smoothing the DEM. Defaults to 3.
        resolution (int, optional): The resolution of the DEM. Defaults to 10.
        keep_files (bool, optional): Whether to keep the intermediate files. Defaults to True.
        error_file (_type_, optional): The file to save the error IDs. Defaults to None.
        ignore_warnings (bool, optional): Whether to ignore warnings. Defaults to True.
    """
    import shutil
    import warnings
    import geopandas as gpd

    if ignore_warnings:
        warnings.filterwarnings("ignore")

    start_time = time.time()

    if filename is None:
        filename = huc8

    if not filename.endswith(".shp"):
        filename = filename + ".shp"

    filename = os.path.abspath(filename)

    if tmp_dir is None:
        tmp_dir = os.path.join(os.getcwd(), huc8)

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    merge = os.path.join(tmp_dir, "mosaic.tif")
    mask = os.path.join(tmp_dir, "mask.geojson")
    clip = os.path.join(tmp_dir, "clip.tif")
    reproj = os.path.join(tmp_dir, "reproj.tif")
    image = os.path.join(tmp_dir, "image.tif")
    median = os.path.join(tmp_dir, "median.tif")
    regions = os.path.join(tmp_dir, "regions.shp")
    regions_info = os.path.join(tmp_dir, "regions_info.csv")

    try:
        download_ned_by_huc(huc8, out_dir=tmp_dir)

        if wbd is None:
            print("Downloading WBD ...")
            hu8_url = "https://drive.google.com/file/d/1AVBPVVAzsLs8dnF_bCvFvGMCAEgaPthh/view?usp=sharing"
            output = os.path.join(tmp_dir, "WBDHU8_CONUS.zip")
            wbd = download_file(hu8_url, output=output, unzip=False)

        if isinstance(wbd, str):
            print("Reading WBD ...")
            gdf = gpd.read_file(wbd)
        elif isinstance(wbd, gpd.GeoDataFrame):
            gdf = wbd
        else:
            raise ValueError("shp_path must be a filepath or a GeoDataFrame.")

        selected = gdf[gdf["huc8"] == huc8].copy()
        selected.to_crs(epsg=4326, inplace=True)
        selected.to_file(mask)

        if not os.path.exists(merge):
            print("Merging NED tiles ...")
            mosaic(tmp_dir, merge)
        clip_image(merge, mask, clip)
        reproject_image(clip, reproj, crs)
        resample(reproj, image, resolution)
        MedianFilter(image, kernel_size, median)
        ExtractSinks(median, min_size, tmp_dir)
        join_tables(regions, regions_info, filename)

        for file in [merge, mask, clip, reproj, image]:
            if os.path.exists(file):
                os.remove(file)

        if not keep_files:
            shutil.rmtree(tmp_dir)
    except Exception as e:
        if error_file is not None:
            with open(error_file, "a") as f:
                f.write(huc8 + "\n")

        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        print(e)
        return None

    end_time = time.time()
    print("Total run time:\t\t\t {:.4f} s\n".format(end_time - start_time))


def extract_sinks_by_huc8_batch(
    huc_ids=None,
    min_size=10,
    out_dir=None,
    tmp_dir=None,
    wbd=None,
    crs="EPSG:5070",
    kernel_size=3,
    resolution=10,
    keep_files=False,
    reverse=False,
    error_file=None,
    ignore_warnings=True,
    ignored_ids=[],
    overwrite=False,
):
    """Extract sinks from NED by HUC8.

    Args:
        huc8 (str): The HUC8 code, e.g., 01010002
        min_size (int, optional): The minimum number of pixels to be considered as a sink. Defaults to 10.
        filename (str, optional): The output depression file name. Defaults to None, e,g., using the HUC8 code.
        tmp_dir (str, optional): The temporary directory. Defaults to None, e.g., using the current directory.
        wbd (str | GeoDataFrame, optional): The watershed boundary file. Defaults to None.
        crs (str, optional): The coordinate reference system. Defaults to "EPSG:5070".
        kernel_size (int, optional): The kernel size for smoothing the DEM. Defaults to 3.
        resolution (int, optional): The resolution of the DEM. Defaults to 10.
        keep_files (bool, optional): Whether to keep the intermediate files. Defaults to True.
        reverse (bool, optional): Whether to reverse the HUC8 list. Defaults to False.
        error_file (_type_, optional): The file to save the error IDs. Defaults to None.
        ignore_warnings (bool, optional): Whether to ignore warnings. Defaults to True.
        overwrite (bool, optional): Whether to overwrite the existing files. Defaults to False.
    """
    import pandas as pd

    if huc_ids is None:
        url = "https://raw.githubusercontent.com/giswqs/lidar/master/examples/data/huc8.csv"
        df = pd.read_csv(url, dtype=str)
        huc_ids = df["huc8_id"].tolist()

    if not isinstance(huc_ids, list):
        huc_ids = [huc_ids]

    if reverse:
        huc_ids = huc_ids[::-1]

    if out_dir is None:
        out_dir = os.getcwd()

    for index, huc8 in enumerate(huc_ids):
        print(f"Processing {index+1}:{len(huc_ids)}: {huc8} ...")
        if huc8 in ignored_ids:
            continue
        filename = os.path.join(out_dir, str(huc8) + ".shp")
        if not os.path.exists(filename) or (os.path.exists(filename) and overwrite):
            extract_sinks_by_huc8(
                huc8,
                min_size,
                filename,
                tmp_dir,
                wbd,
                crs,
                kernel_size,
                resolution,
                keep_files,
                error_file,
                ignore_warnings,
            )
        else:
            print(f"File already exists: {filename}")


def image_to_cog(source, dst_path=None, profile="deflate", **kwargs):
    """Converts an image to a COG file.

    Args:
        source (str): A dataset path, URL or rasterio.io.DatasetReader object.
        dst_path (str, optional): An output dataset path or or PathLike object. Defaults to None.
        profile (str, optional): COG profile. More at https://cogeotiff.github.io/rio-cogeo/profile. Defaults to "deflate".

    Raises:
        ImportError: If rio-cogeo is not installed.
        FileNotFoundError: If the source file could not be found.
    """
    try:
        from rio_cogeo.cogeo import cog_translate
        from rio_cogeo.profiles import cog_profiles

    except ImportError:
        raise ImportError(
            "The rio-cogeo package is not installed. Please install it with `pip install rio-cogeo` or `conda install rio-cogeo -c conda-forge`."
        )

    if not source.startswith("http"):
        source = check_file_path(source)

        if not os.path.exists(source):
            raise FileNotFoundError(
                "The provided input file could not be found.")

    if dst_path is None:
        if not source.startswith("http"):
            dst_path = os.path.splitext(source)[0] + "_cog.tif"
        else:
            dst_path = temp_file_path(extension=".tif")

    dst_path = check_file_path(dst_path)

    dst_profile = cog_profiles.get(profile)
    cog_translate(source, dst_path, dst_profile, **kwargs)
