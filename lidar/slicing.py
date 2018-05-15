import os
import math
import time
import shutil
import numpy as np
from scipy import ndimage
from skimage import measure
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from osgeo import gdal, ogr, osr
from skimage.external.tifffile import TiffFile


# class for true depression
class Depression:
    def __init__(self, id, level, size, volume, meanDepth, maxDepth, minElev, bndElev, inNbrId, regionId):
        self.id = id
        self.size = size
        self.level = level
        self.volume = volume
        self.meanDepth = meanDepth
        self.maxDepth = maxDepth
        self.minElev = minElev
        self.bndElev = bndElev
        self.inNbrId = inNbrId
        self.regionId = regionId


# get min and max elevation of a dem
def get_min_max_nodata(image):
    max_elev = np.max(image)
    nodata = pow(10, math.floor(math.log10(np.max(image))) + 2) - 1  # assign no data value
    image[image <= 0] = nodata  # change no data value
    min_elev = np.min(image)
    return min_elev, max_elev, nodata


# get cell size of tif file
def get_cell_size(tif):
    cell_size = 0
    try:
        meta = tif.info()  # return text, the cell size looks like "* model_pixel_scale (3d) (1.0, 1.0, 0.0)"
        str_begin_index = meta.find("model_pixel_scale (3d)")
        str_end_index = str_begin_index + len("model_pixel_scale (3d)")
        res_begin_index = meta.find("(", str_end_index) + 1
        res_end_index = meta.find(",", res_begin_index)
        cell_size = float(meta[res_begin_index:res_end_index])
    except:
        print("error getting cell size")
    return cell_size


# set input image parameters for level set method
def set_image_paras(no_data, min_size, min_depth, interval, resolution):
    image_paras = {}
    image_paras["no_data"] = no_data
    image_paras["min_size"] = min_size
    image_paras["min_depth"] = min_depth
    image_paras["interval"] = interval
    image_paras["resolution"] = resolution
    return image_paras


# set rainfall simulation parameters
def set_rain_paras(reg_catchment, rain_intensity, rain_time):
    rain_paras = {}
    rain_paras["reg_catchment"] = reg_catchment
    rain_paras["rain_intensity"] = rain_intensity
    rain_paras["rain_time"] = rain_time
    return rain_paras


# get image parameters
def get_image_paras(image_paras):
    no_data = image_paras["no_data"]
    min_size = image_paras["min_size"]
    min_depth = image_paras["min_depth"]
    interval = image_paras["interval"]
    resolution = image_paras["resolution"]
    return no_data, min_size, min_depth, interval, resolution


# get rainfall simulation parameters
def get_rain_paras(rain_paras):
    reg_catchment = rain_paras["reg_catchment"]
    rain_intensity = rain_paras["rain_intensity"]
    rain_time = rain_paras["rain_time"]
    return reg_catchment, rain_intensity, rain_time


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
    # First of all, gather some information from the template file
    data = gdal.Open(template)
    [cols, rows] = arr.shape
    trans = data.GetGeoTransform()
    proj = data.GetProjection()
    # nodatav = 0 #data.GetNoDataValue()
    # Create the file, using the information from the template file
    outdriver = gdal.GetDriverByName("GTiff")
    # http://www.gdal.org/gdal_8h.html
    # GDT_Byte = 1, GDT_UInt16 = 2, GDT_UInt32 = 4, GDT_Int32 = 5, GDT_Float32 = 6,
    outdata   = outdriver.Create(str(out_path), rows, cols, 1, gdal.GDT_UInt32)
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
def polygonize(img,shp_path):
    # mapping between gdal type and ogr field type
    type_mapping = {gdal.GDT_Byte: ogr.OFTInteger,
                    gdal.GDT_UInt16: ogr.OFTInteger,
                    gdal.GDT_Int16: ogr.OFTInteger,
                    gdal.GDT_UInt32: ogr.OFTInteger,
                    gdal.GDT_Int32: ogr.OFTInteger,
                    gdal.GDT_Float32: ogr.OFTReal,
                    gdal.GDT_Float64: ogr.OFTReal,
                    gdal.GDT_CInt16: ogr.OFTInteger,
                    gdal.GDT_CInt32: ogr.OFTInteger,
                    gdal.GDT_CFloat32: ogr.OFTReal,
                    gdal.GDT_CFloat64: ogr.OFTReal}

    ds = gdal.Open(img)
    prj = ds.GetProjection()
    srcband = ds.GetRasterBand(1)
    dst_layername = "Shape"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(shp_path)
    srs = osr.SpatialReference(wkt=prj)

    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
    raster_field = ogr.FieldDefn('level', type_mapping[srcband.DataType])
    dst_layer.CreateField(raster_field)
    gdal.Polygonize(srcband, srcband, dst_layer, 0, [], callback=None)
    del img, ds, srcband, dst_ds, dst_layer


# convert images in a selected folder to shapefiles
def img_to_shp(in_img_dir, out_shp_dir):
    img_files = os.listdir(in_img_dir)
    for img_file in img_files:
        if img_file.endswith(".tif"):
            img_filename = os.path.join(in_img_dir, img_file)
            shp_filename = os.path.join(out_shp_dir, img_file.replace("tif", "shp"))
            polygonize(img_filename, shp_filename)


# parallel processing
def task(region, out_image, no_data, min_size, min_depth, interval, resolution):
    label_id = region.label
    img = region.intensity_image
    # img[img == 0] = no_data
    bbox = region.bbox
    # out_obj = identifyDepression(img,label_id,no_data,min_size,min_depth)
    # writeObject(out_image,out_obj,bbox)
    out_obj = levelSet(img, label_id, no_data, min_size, min_depth, interval, resolution)
    writeObject(out_image, out_obj, bbox)


# identify nested depressions using level set method
def levelSet(img, region_id, obj_uid, image_paras, rain_paras):

    # unzip input parameters from dict
    no_data, min_size, min_depth, interval, resolution = get_image_paras(image_paras)
    catchment_img, rain_intensity, rain_time = get_rain_paras(rain_paras)

    level_img = np.zeros(img.shape)     # init output level image
    flood_img = np.zeros(img.shape)     # init output flood time image

    max_elev = np.max(img)
    img[img == 0] = no_data
    min_elev = np.min(img)

    print("============================================================================= Region: {}".format(region_id))
    unique_id = obj_uid
    # unique_id = 0
    parent_ids = {}  # store current parent depressions
    nbr_ids = {}  # store the inner-neighbor ids of current parent depressions
    dep_list = []  # list for storing depressions

    for elev in np.arange(max_elev, min_elev, interval):  # slicing operation using top-down approach
        img[img > elev] = 0  # set elevation higher than xy-plane to zero
        label_objects, nb_labels = regionGroup(img, min_size, no_data)
        print('slicing elev = {:.2f}, number of objects = {}'.format(elev, nb_labels))
        if nb_labels == 0:   # if slicing results in no objects, quit
            break

        objects = measure.regionprops(label_objects, img)
        for i, object in enumerate(objects):
            (row, col) = object.coords[0]  # get a boundary cell
            bbox = object.bbox

            if len(parent_ids) == 0:  # This is the first depression, maximum depression
                print("This is the maximum depression extent.")
                cells = object.area
                size = cells * pow(resolution, 2)  # depression size
                max_depth = object.max_intensity - object.min_intensity  # depression max depth
                mean_depth = (object.max_intensity * cells - np.sum(object.intensity_image)) / cells  # depression mean depth
                volume = mean_depth * cells * pow(resolution, 2)  # depression volume
                spill_elev = object.max_intensity   # to be implemented
                min_elev = object.min_intensity   # depression min elevation
                max_elev = object.max_intensity     # depression max elevation
                print("size = {}, max depth = {:.2f}, mean depth = {:.2f}, volume = {:.2f}, spill elev = {:.2f}".format(
                    size, max_depth, mean_depth, volume, spill_elev))
                # plt.imshow(object.intensity_image)
                unique_id += 1
                level = 1
                dep_list.append(Depression(unique_id,level,size,volume,mean_depth,max_depth,min_elev,max_elev,[],region_id))
                parent_ids[unique_id] = 0  # number of inner neighbors
                nbr_ids[unique_id] = []   # ids of inner neighbors
                tmp_img = np.zeros(object.image.shape)
                tmp_img[object.image] = unique_id
                writeObject(level_img, tmp_img, bbox)  # write the object to the final image

                # determine inundation area
                if catchment_img is not None:
                    catchment_area = catchment_img[row,col]   # get the catchment size of the depression
                    flood_volume = catchment_area * rain_intensity * rain_time
                    if flood_volume >= volume and flood_img[row, col] == 0:
                        tmp_img = np.zeros(object.image.shape)
                        tmp_img[object.image] = 1
                        writeObject(flood_img, tmp_img, bbox)

            else:  # identify inner neighbors of parent depressions
                # print("current id: {}".format(parent_ids.keys()))
                # (row, col) = object.coords[0]
                parent_id = level_img[row,col]
                parent_ids[parent_id] += 1
                nbr_ids[parent_id].append(i)

                # determine inundation area
                if catchment_img is not None:
                    cells = object.area
                    size = cells * pow(resolution, 2)
                    # max_depth = object.max_intensity - object.min_intensity
                    mean_depth = (object.max_intensity * cells - np.sum(object.intensity_image)) / cells
                    volume = mean_depth * cells * pow(resolution, 2)
                    # spill_elev = object.max_intensity   # to be implemented
                    # min_elev = object.min_intensity
                    # max_elev = object.max_intensity
                    flood_volume = catchment_area * rain_intensity * rain_time
                    if flood_volume >= volume and flood_img[row, col] == 0:
                        tmp_img = np.zeros(object.image.shape)
                        tmp_img[object.image] = 1
                        writeObject(flood_img, tmp_img, bbox)

        for key in parent_ids.copy():  # check how many inner neighbors each upper level depression has
            if parent_ids[key] > 1:  # if the parent has two or more children
                print("Object id: {} has split into {} objects".format(key, parent_ids[key]))
                new_parent_keys = nbr_ids[key]
                for new_key in new_parent_keys:
                    object = objects[new_key]
                    cells = object.area
                    size = cells * pow(resolution, 2)
                    max_depth = object.max_intensity - object.min_intensity
                    mean_depth = (object.max_intensity * cells - np.sum(object.intensity_image)) / cells
                    volume = mean_depth * cells * pow(resolution, 2)
                    spill_elev = object.max_intensity
                    min_elev = object.min_intensity
                    max_elev = object.max_intensity
                    print(
                        "  --  size = {}, max depth = {:.2f}, mean depth = {:.2f}, volume = {:.2f}, spill elev = {:.2f}".format(
                            size, max_depth, mean_depth, volume, spill_elev))
                    unique_id += 1
                    level = 1
                    dep_list.append(
                        Depression(unique_id, level, size, volume, mean_depth, max_depth, min_elev, max_elev, [], region_id))
                    dep_list[key-1-obj_uid].inNbrId.append(unique_id)
                    parent_ids[unique_id] = 0
                    nbr_ids[unique_id] = []
                    bbox = object.bbox
                    tmp_img = np.zeros(object.image.shape)
                    tmp_img[object.image] = unique_id
                    writeObject(level_img, tmp_img, bbox)

                if key in parent_ids.keys():    # remove parent id that has split
                    parent_ids.pop(key)
            else:
                parent_ids[key] = 0     # if a parent depression has not split, keep it
                nbr_ids[key] = []

    for dep in dep_list:
        print("id: {} has children {}".format(dep.id, dep.inNbrId))
    dep_list = updateLevel(dep_list, obj_uid)   # update the inner neighbors of each depression
    for dep in dep_list:
        print("id: {} is level {}".format(dep.id, dep.level))

    del img
    objects = None
    label_objects = None
    # if objects is not None:
    #     del objects
    # if label_objects is not None:
    #     del label_objects

    return level_img, dep_list, flood_img


# update the inner neighbors of each depression
def updateLevel(dep_list, obj_uid):
    for dep in reversed(dep_list):
        if len(dep.inNbrId) == 0:
            dep.level = 1
        else:
            max_children_level = 0
            for id in dep.inNbrId:
                if dep_list[id-1-obj_uid].level > max_children_level:
                    max_children_level = dep_list[id-1-obj_uid].level
            dep.level = max_children_level + 1
    return dep_list


# derive depression level image based on the depression id image and depression list
def obj_to_level(obj_img, dep_list):
    level_img = np.copy(obj_img)

    max_id = int(np.max(level_img))
    # print("max id = " + str(max_id))
    if max_id > 0:
        min_id = int(np.min(level_img[np.nonzero(level_img)]))
        # print("min_id = " + str(min_id))
        for i in range(min_id, max_id+1):
            level_img[level_img == i] = dep_list[i-1].level + max_id
    level_img = level_img - max_id

    return level_img


# derive depression level image based on the depression id image and depression list
def obj_to_level2(obj_img, dep_list):
    level_img = np.copy(obj_img)
    max_id = int(np.max(level_img))
    for i in range(1, max_id+1):
        level_img[level_img == i] = dep_list[i-1].level + max_id
    level_img = level_img - max_id

    return level_img


# save the depression list info to csv
def write_dep_csv(dep_list, csv_file):
    csv = open(csv_file, "w")
    header = "Depression ID" +","+"Level"+","+"Area"+","+"Volume"+","+"Mean depth"+","+"Maximum depth"+","+\
             "Lowest elevation"+","+"Spill elevation"+","+"Children IDs"+","+"Region ID"
    csv.write(header + "\n")
    for dep in dep_list:
        # id, level, size, volume, meanDepth, maxDepth, minElev, bndElev, inNbrId, nbrId = 0
        line = "{},{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{},{}".format(dep.id, dep.level, dep.size, dep.volume,
                dep.meanDepth, dep.maxDepth, dep.minElev, dep.bndElev, str(dep.inNbrId).replace(",",":"), dep.regionId)
        csv.write(line + "\n")
    csv.close()


# extracting individual level image
def extract_levels_bk(level_img, min_size, no_data, out_dir, template, bool_comb=False):
    max_level = int(np.max(level_img))
    combined_images = []
    single_images = []
    img = np.copy(level_img)

    digits = int(math.log10(max_level)) + 1  # determine the level number of output file name
    for i in range(1, max_level + 1):
        img[(img > 0) & (img <= i) ] = i
        tmp_img = np.copy(img)
        tmp_img[tmp_img > i] = 0
        if bool_comb == True:  # whether to extract combined level image
            combined_images.append(np.copy(tmp_img))
            filename_combined = "Combined_level_" + str(i).zfill(digits) + ".tif"
            out_file = os.path.join(out_dir, filename_combined)
            writeRaster(tmp_img,out_file,template)

        lbl_objects, n_labels = regionGroup(tmp_img, min_size, no_data)
        regs = measure.regionprops(lbl_objects, level_img)
        sin_img = np.zeros(img.shape)

        for reg in regs:
            if reg.max_intensity >= i:
                bbox = reg.bbox
                tmp_img = np.zeros(reg.image.shape)
                tmp_img[reg.image] = i
                writeObject(sin_img, tmp_img, bbox)
        del tmp_img
        # single_images.append(np.copy(sin_img))
        filename_single = "Single_level_" + str(i).zfill(digits) + ".tif"
        out_file = os.path.join(out_dir, filename_single)
        writeRaster(sin_img,out_file,template)
        del sin_img

    del img
    return True


def extract_levels(level_img, min_size, no_data, out_img_dir, out_shp_dir, template, bool_comb=False):
    max_level = int(np.max(level_img))
    combined_images = []
    single_images = []
    img = np.copy(level_img)

    digits = int(math.log10(max_level)) + 1  # determine the level number of output file name
    for i in range(1, max_level + 1):
        img[(img > 0) & (img <= i) ] = i
        tmp_img = np.copy(img)
        tmp_img[tmp_img > i] = 0
        if bool_comb == True:  # whether to extract combined level image
            combined_images.append(np.copy(tmp_img))
            filename_combined = "Combined_level_" + str(i).zfill(digits) + ".tif"
            out_file = os.path.join(out_shp_dir, filename_combined)
            writeRaster(tmp_img,out_file,template)

        lbl_objects, n_labels = regionGroup(tmp_img, min_size, no_data)
        regs = measure.regionprops(lbl_objects, level_img)
        sin_img = np.zeros(img.shape)

        for reg in regs:
            if reg.max_intensity >= i:
                bbox = reg.bbox
                tmp_img = np.zeros(reg.image.shape)
                tmp_img[reg.image] = i
                writeObject(sin_img, tmp_img, bbox)
        del tmp_img
        # single_images.append(np.copy(sin_img))
        filename_single = "Single_level_" + str(i).zfill(digits) + ".shp"
        out_shp_file = os.path.join(out_shp_dir, filename_single)

        out_img_file = os.path.join(out_img_dir, "tmp.tif")
        writeRaster(sin_img, out_img_file, in_sink)
        polygonize(out_img_file, out_shp_file)
        # writeRaster(sin_img,out_file,template)
        del sin_img

    del img
    return True


# maximize plot window
def maxPlotWindow():
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()


# display image and legend
def display_image(img, title, legend ="", max_plot = False):
    if max_plot == True:
        maxPlotWindow()
    im = plt.imshow(img)
    plt.suptitle(title, fontsize = 24, fontweight='bold', color = "black")
    if legend != "":
        values = np.unique(img.ravel())[1:]
        colors = [im.cmap(im.norm(value)) for value in values]
        # create a patch (proxy artist) for every color
        patches = [mpatches.Patch(color=colors[i], label=legend + " {l}".format(l=int(values[i]))) for i in range(len(values))]
        # put those patched as legend-handles into the legend
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()
    return True


# delineate nested depressions
def DelineateDepressions(in_sink, min_size, min_depth, interval, out_dir):

    # The following parameters can be used by default
    in_catchment = None
    resolution = 1              # default image resolution if not specified
    rain_intensity = 0.05       # rainfall intensity, e.g., 5.0 cm/h
    rain_time = 5               # rainfall duration, e.g., 2 hours
    interval = interval * (-1)  # convert slicing interval to negative value

    out_img_dir = os.path.join(out_dir, "img-level")
    out_shp_dir = os.path.join(out_dir, "shp-level")
    out_obj_file = os.path.join(out_dir, "object_id.tif")
    out_level_file = os.path.join(out_dir, "object_level.tif")
    out_flood_file = os.path.join(out_dir, "flood.tif")
    out_vec_file = os.path.join(out_dir, "object_vec.shp")
    out_csv_file = os.path.join(out_dir, "object_info.csv")

    init_time = time.time()

    # delete contents in output folder if existing
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if os.path.exists(out_img_dir):
        shutil.rmtree(out_img_dir)
    os.mkdir(out_img_dir)
    if os.path.exists(out_shp_dir):
        shutil.rmtree(out_shp_dir)
    os.mkdir(out_shp_dir)

    print("Reading data ...")
    read_time = time.time()

    with TiffFile(in_sink) as tif:   # read the dem file as numpy array
        image = tif.asarray()
        rows_cols = image.shape
        print("rows, cols: " + str(rows_cols))
        cell_size = get_cell_size(tif)  # get image cell size
        if cell_size > 0:
            resolution = cell_size
            print("Pixel resolution: " + str(resolution))
    print("Read data time: {}".format(time.time() - read_time))
    # image = np.copy(raw_image)  # store original DEM
    min_elev, max_elev, no_data = get_min_max_nodata(image)  # set nodata value to a large value, e.g., 9999
    # initialize output image
    obj_image = np.zeros(image.shape)  # output depression image with unique id for each nested depression
    level_image = np.zeros(image.shape)  # output depression level image
    flood_image = None
    catchment_img = None

    if in_catchment is not None:  # whether or not to derive inundation area
        with TiffFile(in_catchment) as catchment:  # read catchment file as numpy array
            catchment_img = catchment.asarray()
        flood_image = np.copy(obj_image)  # output inundation image

    # nb_labels is the total number of objects. 0 represents background object.
    label_objects, nb_labels = regionGroup(image, min_size, no_data)
    regions = measure.regionprops(label_objects, image)
    del image  # delete the original image to save memory
    prep_time = time.time()
    print("Data preparation time: {}".format(prep_time - init_time))
    print("Total number of regions: {}".format(nb_labels))
    # plt.imshow(label_objects, cmap=plt.cm.spectral)
    # plt.show()

    identify_time = time.time()

    obj_uid = 0
    global_dep_list = []

    # loop through regions and identify nested depressions in each region using level-set method
    for region in regions:  # iterate through each depression region
        region_id = region.label
        img = region.intensity_image  # dem subset for each region
        bbox = region.bbox

        if in_catchment is not None:
            reg_catchment = catchment_img[bbox[0]:bbox[2], bbox[1]:bbox[3]]   # extract the corresponding catchment
        else:
            reg_catchment = None
            rain_intensity = None
            rain_time = None

        # save all input parameters needed for level set methods as a dict
        image_paras = set_image_paras(no_data, min_size, min_depth, interval, resolution)
        rain_paras = set_rain_paras(reg_catchment, rain_intensity, rain_time)

        # execute level set methods
        out_obj, dep_list, flood_obj = levelSet(img, region_id, obj_uid, image_paras, rain_paras)

        for dep in dep_list:
            global_dep_list.append(dep)

        obj_uid += len(dep_list)

        level_obj = obj_to_level(out_obj, global_dep_list)
        obj_image = writeObject(obj_image, out_obj, bbox)       # write region to whole image
        level_image = writeObject(level_image, level_obj, bbox)

        if in_catchment is not None:
            flood_image = writeObject(flood_image, flood_obj, bbox)

        del out_obj, level_obj, flood_obj, region

    del regions, label_objects

    print("=========== Run time statistics =========== ")
    print("(rows, cols):\t\t\t {0}".format(str(rows_cols)))
    print("Pixel resolution:\t\t {0} m".format(str(resolution)))
    print("Number of regions:\t\t {0}".format(str(nb_labels)))
    print("Data preparation time:\t {:.4f} s".format(prep_time - init_time))
    print("Identify level time:\t {:.4f} s".format(time.time() - identify_time))

    write_time = time.time()
    writeRaster(obj_image, out_obj_file, in_sink)
    writeRaster(level_image, out_level_file, in_sink)
    print("Write image time:\t\t {:.4f} s".format(time.time() - write_time))

    # # del obj_image
    # # extracting single and combined level images
    # level_time = time.time()
    # extract_levels(level_image, min_size, no_data, out_img_dir, in_sink, False)
    # print("Extract level time:\t\t {:.4f} s".format(time.time() - level_time))

    # vector_time = time.time()
    # img_to_shp(out_img_dir,out_shp_dir)
    # print("Vectorizing run time:\t {:.4f} s".format(time.time() - vector_time))

    # del obj_image
    # extracting single and combined level images
    level_time = time.time()
    polygonize(out_obj_file, out_vec_file)
    write_dep_csv(global_dep_list, out_csv_file)
    # extract_levels(level_image, min_size, no_data, out_img_dir, out_shp_dir, in_sink, False)
    print("Extract level time:\t\t {:.4f} s".format(time.time() - level_time))

    # writeRaster(obj_image, out_obj_file, in_sink)
    # writeRaster(level_image, out_level_file, in_sink)
    if in_catchment is not None:
        write_time = time.time()
        writeRaster(flood_image, out_flood_file, in_sink)
        print("Write data time:\t\t {:.4f} s".format(time.time() - write_time))

    end_time = time.time()
    print("Total run time:\t\t\t {:.4f} s".format(end_time - init_time))
    return True


# #####################################  main script
if __name__ == '__main__':

    # ************************ change the following parameters if needed ******************************** #
    # set input files
    in_dem = "../data/dem.tif"
    in_sink = "../data/sink.tif"
    # parameters for level set method
    min_size = 1000         # minimum number of pixels as a depression
    min_depth = 0.3         # minimum depression depth
    interval = 0.2          # slicing interval, top-down approach
    # set output directory
    out_dir = os.path.join(os.path.expanduser("~"), "temp")  # create a temp folder under user home directory
    # **************************************************************************************************#

    results = DelineateDepressions(in_sink, min_size, min_depth, interval, out_dir)

    print("Results are saved in: {}".format(out_dir))
