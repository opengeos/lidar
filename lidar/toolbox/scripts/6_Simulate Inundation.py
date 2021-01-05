from skimage.external.tifffile import TiffFile
import matplotlib.pyplot as plt
from scipy import ndimage
from skimage import measure
import numpy as np
import math
import time
import os
import shutil
import copy
from osgeo import gdal, ogr, osr
import arcpy

# from joblib import Parallel, delayed
# from skimage.external.tifffile import TiffWriter
# from skimage import segmentation
# import tempfile

## class for true depression
class Depression:
    def __init__(
        self,
        id,
        level,
        size,
        volume,
        meanDepth,
        maxDepth,
        minElev,
        bndElev,
        inNbrId,
        nbrId=0,
    ):
        self.id = id
        self.size = size
        self.level = level
        self.volume = volume
        self.meanDepth = meanDepth
        self.maxDepth = maxDepth
        self.minElev = minElev
        self.bndElev = bndElev
        self.inNbrId = inNbrId
        self.outNbrId = nbrId


# get min and max elevation of a dem
def get_min_max_nodata(img_array):
    max_elev = np.max(img_array)
    nodata = (
        pow(10, math.floor(math.log10(np.max(image))) + 2) - 1
    )  # based on the max value of the image, assign no data value
    image[image <= 0] = nodata  # change no data value
    min_elev = np.min(img_array)
    return min_elev, max_elev, nodata


# identify regions based on region growing method
def regionGroup(img_array, min_size, no_data):
    img_array[img_array == no_data] = 0
    label_objects, nb_labels = ndimage.label(img_array)
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = sizes > min_size
    mask_sizes[0] = 0
    image_cleaned = mask_sizes[label_objects]
    label_objects, nb_labels = ndimage.label(image_cleaned)
    # nb_labels is the total number of objects. 0 represents background object.
    return label_objects, nb_labels


# extract a subset of the raster
def extractByRegion(img_array, obj_array, obj_id, no_data):
    # no_data = 9999
    slice_x, slice_y = ndimage.find_objects(obj_array == obj_id)[0]
    roi = obj_array[slice_x, slice_y]
    img_roi = np.copy(img_array[slice_x, slice_y])
    img_roi[roi != obj_id] = no_data
    return img_roi, roi, slice_x, slice_y


def extractRegionById(obj_array, obj_id):
    slice_x, slice_y = ndimage.find_objects(obj_array == obj_id)[0]
    roi = obj_array[slice_x, slice_y]
    return roi


# write output depression raster
def writeObject(img_array, obj_array, bbox):
    min_row, min_col, max_row, max_col = bbox
    roi = img_array[min_row:max_row, min_col:max_col]
    roi[obj_array > 0] = obj_array[obj_array > 0]
    return img_array


def writeRaster(arr, out_path, template):
    no_data = 0
    # First of all, gather some information from the original file
    data = gdal.Open(template)
    [cols, rows] = arr.shape
    trans = data.GetGeoTransform()
    proj = data.GetProjection()
    # nodatav = 0 #data.GetNoDataValue()
    # Create the file, using the information from the original file
    outdriver = gdal.GetDriverByName("GTiff")
    # http://www.gdal.org/gdal_8h.html
    # GDT_Byte = 1, GDT_UInt16 = 2, GDT_UInt32 = 4, GDT_Int32 = 5, GDT_Float32 = 6,
    outdata = outdriver.Create(str(out_path), rows, cols, 1, gdal.GDT_UInt32)
    # Write the array to the file, which is the original array in this example
    outdata.GetRasterBand(1).WriteArray(arr)
    # Set a no data value if required
    outdata.GetRasterBand(1).SetNoDataValue(no_data)
    # Georeference the image
    outdata.SetGeoTransform(trans)
    # Write projection information
    outdata.SetProjection(proj)
    return arr


# raster to vector
def polygonize(img, shp_path):
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

    tif = os.path.split(img)[1]
    # print("reading {}...".format(tif))
    ds = gdal.Open(img)
    prj = ds.GetProjection()
    srcband = ds.GetRasterBand(1)
    # create shapefile datasource from geotiff file
    shp = os.path.split(shp_path)[1]
    print("creating {}...".format(shp))
    dst_layername = "Shape"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(shp_path)
    srs = osr.SpatialReference(wkt=prj)

    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
    raster_field = ogr.FieldDefn("level", type_mapping[srcband.DataType])
    dst_layer.CreateField(raster_field)
    gdal.Polygonize(srcband, srcband, dst_layer, 0, [], callback=None)


def img_to_shp(in_img_dir, out_shp_dir):
    img_files = os.listdir(in_img_dir)
    for img_file in img_files:
        if img_file.endswith(".tif"):
            img_filename = os.path.join(in_img_dir, img_file)
            shp_filename = os.path.join(out_shp_dir, img_file.replace("tif", "shp"))
            polygonize(img_filename, shp_filename)


def task(region, out_image, no_data, min_size, min_depth, interval, resolution):
    label_id = region.label
    img = region.intensity_image
    bbox = region.bbox
    out_obj = levelSet(
        img, label_id, no_data, min_size, min_depth, interval, resolution
    )
    writeObject(out_image, out_obj, bbox)


def levelSet(
    img,
    obj_id,
    no_data,
    min_size,
    min_depth,
    interval,
    resolution,
    catchment_img,
    rain_intensity,
    rain_time,
):
    level_img = np.copy(img)
    level_img[level_img != 0] = 0
    obj_img = np.copy(level_img)
    flood_img = np.copy(level_img)

    max_elev = np.max(img)
    img[img == 0] = no_data
    min_elev = np.min(img)

    unique_id = 0
    parent_ids = {}  #  store current upper level depressions
    nbr_ids = {}  # store the inner-neighbor ids of current upper level depressions
    dep_list = []  # list for storing depressions

    for elev in np.arange(max_elev, min_elev, interval):
        img[img > elev] = 0
        label_objects, nb_labels = regionGroup(img, min_size, no_data)
        if nb_labels == 0:  # if slicing results in no objects, quit
            break

        objects = measure.regionprops(label_objects, img)
        for i, object in enumerate(objects):
            (row, col) = object.coords[0]
            bbox = object.bbox
            if len(parent_ids) == 0:  # This is the first depression, maximum depression
                cells = object.area
                size = cells * pow(resolution, 2)
                max_depth = object.max_intensity - object.min_intensity
                mean_depth = (
                    object.max_intensity * cells - np.sum(object.intensity_image)
                ) / cells
                volume = mean_depth * cells * pow(resolution, 2)
                spill_elev = object.max_intensity  # to be implemented
                min_elev = object.min_intensity
                max_elev = object.max_intensity
                unique_id += 1
                level = 1
                dep_list.append(
                    Depression(
                        unique_id,
                        level,
                        size,
                        volume,
                        mean_depth,
                        max_depth,
                        min_elev,
                        max_elev,
                        [],
                        0,
                    )
                )
                parent_ids[unique_id] = 0  # number of inner neighbors
                nbr_ids[unique_id] = []  # ids of inner neighbors
                tmp_img = np.zeros(object.image.shape)
                tmp_img[object.image] = unique_id
                writeObject(
                    level_img, tmp_img, bbox
                )  # write the object to the final image

                # determine inundation area
                if catchment_img is not None:
                    catchment_area = catchment_img[row, col]
                    flood_volume = catchment_area * rain_intensity * rain_time
                    if flood_volume >= volume and flood_img[row, col] == 0:
                        tmp_img = np.zeros(object.image.shape)
                        tmp_img[object.image] = 1
                        writeObject(flood_img, tmp_img, bbox)
            else:  # identify inner neighbors of upper level depressions
                parent_id = level_img[row, col]
                parent_ids[parent_id] += 1
                nbr_ids[parent_id].append(i)

                # determine inundation area
                if catchment_img is not None:
                    cells = object.area
                    size = cells * pow(resolution, 2)
                    # max_depth = object.max_intensity - object.min_intensity
                    mean_depth = (
                        object.max_intensity * cells - np.sum(object.intensity_image)
                    ) / cells
                    volume = mean_depth * cells * pow(resolution, 2)
                    # spill_elev = object.max_intensity   # to be implemented
                    # min_elev = object.min_intensity
                    # max_elev = object.max_intensity
                    flood_volume = catchment_area * rain_intensity * rain_time
                    if flood_volume >= volume and flood_img[row, col] == 0:
                        tmp_img = np.zeros(object.image.shape)
                        tmp_img[object.image] = 1
                        writeObject(flood_img, tmp_img, bbox)

        for (
            key
        ) in (
            parent_ids.copy()
        ):  # check how many inner neighbors each upper level depression has
            if parent_ids[key] > 1:  # if the parent has two or more children
                new_parent_keys = nbr_ids[key]
                for new_key in new_parent_keys:
                    object = objects[new_key]
                    cells = object.area
                    size = cells * pow(resolution, 2)
                    max_depth = object.max_intensity - object.min_intensity
                    mean_depth = (
                        object.max_intensity * cells - np.sum(object.intensity_image)
                    ) / cells
                    volume = mean_depth * cells * pow(resolution, 2)
                    # spill_elev = object.max_intensity
                    min_elev = object.min_intensity
                    max_elev = object.max_intensity
                    unique_id += 1
                    level = 1
                    dep_list.append(
                        Depression(
                            unique_id,
                            level,
                            size,
                            volume,
                            mean_depth,
                            max_depth,
                            min_elev,
                            max_elev,
                            [],
                            0,
                        )
                    )
                    dep_list[key - 1].inNbrId.append(unique_id)
                    parent_ids[unique_id] = 0
                    nbr_ids[unique_id] = []
                    bbox = object.bbox
                    tmp_img = np.zeros(object.image.shape)
                    tmp_img[object.image] = unique_id
                    writeObject(level_img, tmp_img, bbox)

                if key in parent_ids.keys():
                    parent_ids.pop(key)
            else:
                parent_ids[key] = 0
                nbr_ids[key] = []
    dep_list = updateLevel(dep_list)
    return level_img, dep_list, flood_img


def updateLevel(dep_list):
    for dep in reversed(dep_list):
        if len(dep.inNbrId) == 0:
            dep.level = 1
        else:
            max_children_level = 0
            for id in dep.inNbrId:
                if dep_list[id - 1].level > max_children_level:
                    max_children_level = dep_list[id - 1].level
            dep.level = max_children_level + 1
    return dep_list


def obj_to_level(obj_img, dep_list):
    level_img = np.copy(obj_img)
    max_id = int(np.max(level_img))
    for i in range(1, max_id + 1):
        level_img[level_img == i] = dep_list[i - 1].level + max_id
    level_img = level_img - max_id
    return level_img


# extracting individual level image
def extract_levels(level_img, min_size, no_data, out_dir, template, bool_comb):
    max_level = int(np.max(level_img))
    combined_images = []
    single_images = []
    img = np.copy(level_img)
    digits = (
        int(math.log10(max_level)) + 1
    )  # determine the level number of output file name
    for i in range(1, max_level + 1):
        img[(img > 0) & (img <= i)] = i
        tmp_img = np.copy(img)
        tmp_img[tmp_img > i] = 0
        if bool_comb == True:
            combined_images.append(np.copy(tmp_img))
            filename_combined = "Combined_level_" + str(i).zfill(digits) + ".tif"
            out_file = os.path.join(out_dir, filename_combined)
            writeRaster(tmp_img, out_file, template)
        lbl_objects, n_labels = regionGroup(tmp_img, min_size, no_data)
        regs = measure.regionprops(lbl_objects, level_img)
        sin_img = np.zeros(img.shape)
        for reg in regs:
            if reg.max_intensity >= i:
                bbox = reg.bbox
                tmp_img = np.zeros(reg.image.shape)
                tmp_img[reg.image] = i
                writeObject(sin_img, tmp_img, bbox)
        single_images.append(np.copy(sin_img))
        filename_single = "Single_level_" + str(i).zfill(digits) + ".tif"
        out_file = os.path.join(out_dir, filename_single)
        writeRaster(sin_img, out_file, template)
    return True


def maxPlotWindow():
    figManager = plt.get_current_fig_manager()  # maximize plot window
    figManager.window.showMaximized()


######################################  main script
if __name__ == "__main__":

    # input files
    # in_sink = "../../data/sink_full.tif"   # input dem image
    # in_catchment = "../../data/catchment.tif"  # input catchment image
    #
    # # output files
    # out_flood_dir = r"C:\temp\ppr\flood"
    # out_flood_file = r"C:\temp\ppr\flood\flood60.tif"

    # # input parameters
    # min_size = 1000         # minimum depression size
    # min_depth = 0.2         # minimum depression depth
    # interval = 0.2 * (-1)     # slicing interval, top-down approach, negative values
    # resolution = 1          # image resolution
    # rain_intensity = 0.05   # rainfall intensity, e.g., 5.0 cm/h
    # rain_time = 100          # rainfall duration, e.g., 2 hours
    # step = 2

    in_sink = arcpy.GetParameterAsText(0)
    in_catchment = arcpy.GetParameterAsText(1)
    min_size = float(arcpy.GetParameterAsText(2))
    min_depth = float(arcpy.GetParameterAsText(3))
    interval = float(arcpy.GetParameterAsText(4)) * -1
    rain_intensity = float(arcpy.GetParameterAsText(5)) / 100
    rain_time = float(arcpy.GetParameterAsText(6))
    step = float(arcpy.GetParameterAsText(7))
    out_flood_dir = arcpy.GetParameterAsText(8)

    init_time = time.time()

    desc = arcpy.Describe(in_sink)
    in_sink = desc.catalogPath

    desc = arcpy.Describe(in_catchment)
    in_catchment = desc.catalogPath

    ras = arcpy.Raster(in_sink)
    resolution = ras.meanCellWidth

    # if os.path.exists(out_flood_dir):
    #     shutil.rmtree(out_flood_dir)
    # os.mkdir(out_flood_dir)

    with TiffFile(in_sink) as tif:  # read dem file as numpy array
        raw_image = tif.asarray()

    image = np.copy(raw_image)  # store original DEM
    min_elev, max_elev, no_data = get_min_max_nodata(
        image
    )  # set nodata value to a large value, e.g., 9999

    with TiffFile(in_catchment) as catchment:  # read catchment file as numpy array
        catchment_img = catchment.asarray()

    # nb_labels is the total number of objects. 0 represents background object.
    label_objects, nb_labels = regionGroup(image, min_size, no_data)
    regions = measure.regionprops(label_objects, image)
    prep_time = time.time()
    arcpy.AddMessage("Data preparation time: {}".format(prep_time - init_time))
    arcpy.AddMessage("Total number of regions: {}".format(nb_labels))

    digits = (
        int(math.log10(rain_time / step)) + 1
    )  # determine the level number of output file name
    for i, rain in enumerate(np.arange(step, rain_time + 0.0001, step)):
        iter_start = time.time()
        iter_number = i + 1
        rain_duration = step * (iter_number)
        arcpy.AddMessage("Iteration: {} ( H = {} )".format(iter_number, rain_duration))
        # print("Time step: H = {} ".format(rain_duration))
        flood_image = np.copy(
            image
        )  # output depression image with unique id for each nested depression
        flood_image[flood_image > 0] = 0  # set default value to 0

        out_flood_name = "flood_" + str(iter_number).zfill(digits) + ".tif"
        out_flood_file = os.path.join(out_flood_dir, out_flood_name)

        for region in copy.deepcopy(regions):  # iterate through each depression region
            label_id = region.label
            img = region.intensity_image  # dem
            bbox = region.bbox
            reg_catchment = catchment_img[
                bbox[0] : bbox[2], bbox[1] : bbox[3]
            ]  # extract the corresponding catchment
            out_obj, dep_list, flood_obj = levelSet(
                img,
                label_id,
                no_data,
                min_size,
                min_depth,
                interval,
                resolution,
                reg_catchment,
                rain_intensity,
                rain_duration,
            )
            flood_image = writeObject(flood_image, flood_obj, bbox)
        arcpy.AddMessage("Iteration time: {:.4f}".format(time.time() - iter_start))
        writeRaster(flood_image, out_flood_file, in_sink)

    end_time = time.time()
    arcpy.AddMessage("Total run time: {:.4f}".format(end_time - init_time))
