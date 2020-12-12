import arcpy
import os
import time
import shutil

def DelineateCatchment(DEMRasterPath,flow_direction, SinkPolyPath,OutputPath):
    # arcpy.CheckOutExtension("Spatial")
    workspace = os.path.split(OutputPath)[0]
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True

    if arcpy.Exists(DEMRasterPath) == False:
        print("The input raster does not exist")
        quit()

    if (os.path.splitext(OutputPath)[1].lower() == ".shp"):
        FieldOID = "FID"
        FlowDirection = os.path.join(workspace,"FlowDirection.tif")
        SinkRaster = os.path.join(workspace,"SinkRaster.tif")
        Watershed = os.path.join(workspace,"Watershed.tif")
        Catchment = os.path.join(workspace,"Catchment.shp")
    else:
        FieldOID = "OBJECTID"
        FlowDirection = os.path.join(workspace,"FlowDirection")
        SinkRaster = os.path.join(workspace,"SinkRaster")
        Watershed = os.path.join(workspace,"Watershed")
        Catchment = os.path.join(workspace,"Catchment")

    input_dem = arcpy.Raster(DEMRasterPath)
    # flow_direction = arcpy.sa.FlowDirection(input_dem)
    # flow_direction.save(FlowDirection)

    cell_size = input_dem.meanCellWidth
    arcpy.env.extent = input_dem.extent
    arcpy.PolygonToRaster_conversion(SinkPolyPath,FieldOID,SinkRaster,"CELL_CENTER","NONE",cell_size)

    watershed = arcpy.sa.Watershed(flow_direction,SinkRaster,"Value")
    # watershed.save(Watershed)
    arcpy.RasterToPolygon_conversion(watershed,OutputPath,"NO_SIMPLIFY","Value")
    return OutputPath


def calculateArea(in_shp, fieldname):
    arcpy.AddField_management(in_shp, fieldname,field_type="DOUBLE")
    arcpy.CalculateField_management(in_shp, fieldname, "!shape.area@squaremeters!", "PYTHON_9.3")
    arcpy.DeleteField_management(in_shp,"GRIDCODE")


def spatialJoin(target_shp, join_shp, out_shp):
    arcpy.SpatialJoin_analysis(target_shp, join_shp,out_shp, match_option="INTERSECT")


def mergeShapefiles(in_dir, out_shp):
    arcpy.env.workspace = in_dir
    shapefiles = arcpy.ListFeatureClasses()
    arcpy.Merge_management(shapefiles, out_shp)


def zonalStatistics(in_shp_dir, in_dem):
    arcpy.env.workspace = in_shp_dir
    arcpy.env.overwriteOutput = True
    arcpy.env.snapRaster = in_dem
    dbf_dir = os.path.join(in_shp_dir, "dbf")
    os.mkdir(dbf_dir)
    tif_dir = os.path.join(in_shp_dir, "tif")
    os.mkdir(tif_dir)
    shapefiles = os.listdir(in_shp_dir)
    dem = arcpy.Raster(in_dem)
    cell_size = dem.meanCellHeight
    for shp in shapefiles:
        if shp.endswith(".shp"):
            shp_path = os.path.join(in_shp_dir, shp)
            dbf_path = os.path.join(dbf_dir, "zonal_" + shp.replace("shp", "dbf"))
            tif_path = os.path.join(tif_dir,shp.replace("shp", "tif"))
            arcpy.PolygonToRaster_conversion(shp_path, value_field='FID',out_rasterdataset=tif_path,cell_assignment="CELL_CENTER", priority_field="NONE", cellsize=cell_size)
            arcpy.sa.ZonalStatisticsAsTable(tif_path,"Value",in_dem,dbf_path,"DATA","ALL")
            arcpy.JoinField_management(shp_path,in_field="FID",join_table=dbf_path, join_field="Value",fields="COUNT;AREA;MIN;MAX;RANGE;MEAN;STD;SUM")
            arcpy.AddField_management(shp_path,field_name="dep2catR",field_type="FLOAT")
            arcpy.AddField_management(shp_path,field_name="volume",field_type="FLOAT")
            arcpy.AddField_management(shp_path,field_name="mean_depth",field_type="FLOAT")
            arcpy.CalculateField_management(shp_path,field="dep2catR",expression="!AREA! / !cat_area!", expression_type="PYTHON_9.3")
            arcpy.CalculateField_management(shp_path, field="volume", expression="( !COUNT! * !MAX! - !SUM!) * ( !AREA! / !COUNT! )",expression_type="PYTHON_9.3")
            arcpy.CalculateField_management(shp_path, field="mean_depth", expression="!volume! / !AREA!", expression_type="PYTHON_9.3")
    if os.path.exists(dbf_dir):
        shutil.rmtree(dbf_dir)
    if os.path.exists(tif_dir):
        shutil.rmtree(tif_dir)
    return True


if __name__ == '__main__':

    arcpy.CheckOutExtension("Spatial")
    init_time = time.time()
    in_dem = arcpy.GetParameterAsText(0)
    in_shp_dir = arcpy.GetParameterAsText(1)
    out_img = arcpy.GetParameterAsText(2)

    desc = arcpy.Describe(in_dem)
    in_dem = desc.catalogPath     # get file path

    out_catchment_dir = os.path.split(out_img)[0]
    workspace = out_catchment_dir
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True
    out_shp_dir = os.path.join(out_catchment_dir,"shapefiles")
    out_level_dir = os.path.join(out_shp_dir,"level")

    if os.path.exists(out_shp_dir) == False:
        os.mkdir(out_shp_dir)
    if os.path.exists(out_level_dir) == False:
        os.mkdir(out_level_dir)

    dem = arcpy.Raster(in_dem)
    cell_size = dem.meanCellWidth
    flow_direction = arcpy.sa.FlowDirection(dem)
    flow_direction.save(os.path.join(out_shp_dir, "flowDir.tif"))

    for in_shp in os.listdir(in_shp_dir):
        if in_shp.endswith(".shp")and "Single" in in_shp:
            # print(in_shp)
            # out_shp_name = "Catchment_level_" + in_shp[len(in_shp)-5:]
            out_shp_name = in_shp.replace("Single", "Catchment")
            out_shp_path = os.path.join(out_shp_dir,out_shp_name)
            in_shp_path = os.path.join(in_shp_dir, in_shp)
            arcpy.AddMessage("Generating {} ...".format(out_shp_name))
            DelineateCatchment(in_dem,flow_direction, in_shp_path, out_shp_path)
            calculateArea(out_shp_path, "cat_area")
            out_level_path = os.path.join(out_level_dir,in_shp)
            spatialJoin(in_shp_path,out_shp_path, out_level_path)

    zonalStatistics(out_level_dir,in_dem)
    arcpy.Delete_management(os.path.join(out_shp_dir, "flowDir.tif"))
    arcpy.Delete_management(os.path.join(out_shp_dir, "SinkRaster.tif"))

    out_merge_levels = os.path.join(out_shp_dir, "Merge_levels.shp")
    mergeShapefiles(out_level_dir, out_merge_levels)

    tmp_catchment_img = os.path.join(out_catchment_dir, "catchment_tmp.tif")
    arcpy.PolygonToRaster_conversion(out_merge_levels, value_field="cat_area", out_rasterdataset=tmp_catchment_img, cellsize=cell_size)
    # out_catchment_img = os.path.join(out_catchment_dir, "catchment.tif")
    out_catchment_img = out_img
    arcpy.env.compression = "NONE"
    arcpy.CopyRaster_management(tmp_catchment_img,out_catchment_img,nodata_value="0",pixel_type="32_BIT_UNSIGNED",format="TIFF")
    arcpy.Delete_management(tmp_catchment_img)

    tmp_info_dir = os.path.join(out_shp_dir, "info")
    if os.path.exists(tmp_info_dir):
        shutil.rmtree(tmp_info_dir)
    log_file = os.path.join(out_shp_dir,"log")
    if os.path.exists(log_file):
        os.remove(log_file)

    end_time = time.time()
    arcpy.AddMessage("Total run time: {}".format(end_time - init_time))
