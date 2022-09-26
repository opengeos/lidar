import arcpy
import os
import time
import shutil
from arcpy import env
from arcpy.mp import *

# from arcpy.mapping import *


def extract_sink(in_dem, min_size, min_depth, buffer_dist, out_sink):
    arcpy.CheckOutExtension("Spatial")
    start_time = time.time()
    workspace = os.path.split(out_sink)[0]
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True

    img_ext = ".tif"
    vec_ext = ".shp"
    if os.path.splitext(workspace)[1].lower() == ".gdb":
        img_ext = ""
        vec_ext = ""
    if arcpy.Exists(in_dem) == False:
        arcpy.AddMessage("The input raster does not exist")
        quit()

    ### Mean Focal Statistics
    arcpy.AddMessage("DEM filtering ...")
    ras_mf = arcpy.sa.FocalStatistics(in_dem, "Rectangle 3 3 CELL", "MEAN", "DATA")

    ### Fill depression
    arcpy.AddMessage("Filling sinks ...")
    ras_fill = arcpy.sa.Fill(ras_mf)
    dem_filled_name = "dem_fully_filled" + img_ext
    ras_fill = arcpy.sa.ApplyEnvironment(ras_fill)
    ras_fill.save(os.path.join(workspace, dem_filled_name))
    ### Get sink
    ras_sink = ras_fill - ras_mf

    ### Convert sink to binary image
    arcpy.AddMessage("Creating binary sink image ...")
    ras_sink_bin = arcpy.sa.Con(ras_sink > 0, 1)
    ras_sink_bin_name = "ras_sink_bin" + img_ext
    ras_sink_bin.save(os.path.join(workspace, ras_sink_bin_name))
    ### Region group
    arcpy.AddMessage("Grouping regions ...")
    ras_region_bk = arcpy.sa.RegionGroup(ras_sink_bin, "FOUR", "WITHIN", "ADD_LINK")
    ras_region_bk.save(os.path.join(workspace, "ras_region" + img_ext))

    ras_region_zonal = arcpy.sa.ZonalStatistics(ras_region_bk, "Value", in_dem, "RANGE")
    ras_region_zonal.save("ras_region_zonal" + img_ext)
    ras_region = arcpy.sa.Con(ras_region_zonal > min_depth, 1)
    ras_region = arcpy.sa.ApplyEnvironment(ras_region)
    ras_region.save("ras_region_sub" + img_ext)

    ### Convert raster to polygon
    arcpy.AddMessage("Converting raster to polygon ...")
    region_poly_name = os.path.join(workspace, "region_poly" + vec_ext)
    arcpy.RasterToPolygon_conversion(ras_region, region_poly_name, "NO_SIMPLIFY")

    ### Select polygon based on minimum size
    arcpy.AddMessage("Selecting polygons ...")
    area_field = "Area"
    arcpy.AddField_management(region_poly_name, area_field, "DOUBLE")
    arcpy.CalculateField_management(
        region_poly_name, "Area", "!shape.area@squaremeters!", "PYTHON_9.3", "#"
    )
    sqlExp = area_field + ">=" + str(min_size)

    region_poly_select_name = out_sink
    arcpy.Select_analysis(region_poly_name, region_poly_select_name, sqlExp)
    arcpy.CalculateField_management(region_poly_select_name, "gridcode", "1", "PYTHON")
    region_poly_ras = os.path.join(workspace, "region_poly_ras" + img_ext)
    arcpy.PolygonToRaster_conversion(
        region_poly_select_name, "gridcode", region_poly_ras, "CELL_CENTER", "NONE", "1"
    )

    ### Convert foreground sink to 0
    arcpy.AddMessage("Converting foreground sink ...")
    ras_sink_bg = ras_mf - ras_mf
    ras_sink_bg.save(os.path.join(workspace, "ras_sink_bg" + img_ext))

    # ras_sink_final = "ras_sink_final"
    arcpy.AddMessage("Calculating cell statistics ...")
    in_ras_list = [region_poly_ras, ras_sink_bg]
    ras_sink_final_name = arcpy.sa.CellStatistics(in_ras_list, "SUM", "DATA")
    arcpy.env.extent = ras_mf.extent
    ras_sink_final_name = arcpy.sa.ApplyEnvironment(ras_sink_final_name)

    ### Convert foreground sink
    arcpy.AddMessage("Creating partially filled DEM ...")
    dem_name = arcpy.sa.Con(ras_sink_final_name == 1, ras_mf, ras_fill)
    dem_name = arcpy.sa.ApplyEnvironment(dem_name)
    dem_name.save(os.path.join(workspace, "dem_partially_filled" + img_ext))

    arcpy.AddMessage("Creating sink DEM ...")
    dem_sink = arcpy.sa.Con(ras_sink_final_name == 1, ras_mf, ras_fill)
    dem_sink = arcpy.sa.ApplyEnvironment(dem_sink)

    arcpy.AddMessage("Calculating sink depth ...")
    dem_sink_depth = ras_fill - dem_name
    dem_sink_depth_name = arcpy.sa.Con(dem_sink_depth > 0, dem_sink)
    dem_sink_depth_name = arcpy.sa.ApplyEnvironment(dem_sink_depth_name)
    dem_sink_depth_name.save(os.path.join(workspace, "sink" + img_ext))

    sink_depth = arcpy.sa.Con(dem_sink_depth > 0, dem_sink_depth)
    sink_depth = arcpy.sa.ApplyEnvironment(sink_depth)
    sink_depth.save(os.path.join(workspace, "sink_depth" + img_ext))

    arcpy.AddMessage("Zonal statistics ...")
    zonalStatistics(out_sink, in_dem)

    if buffer_dist > 0:
        arcpy.AddMessage("Creating buffered sink DEM ...")
        sink_buffer_name = os.path.join(workspace, "sink_buffer_poly" + vec_ext)
        sqlExp = str(buffer_dist) + " Meters"
        arcpy.Buffer_analysis(out_sink, sink_buffer_name, sqlExp, "", "", "", "")
        dem_sink_buffer = arcpy.sa.ExtractByMask(dem_name, sink_buffer_name)
        sink_buffer_img = os.path.join(workspace, "sink_buffer" + img_ext)
        arcpy.CopyRaster_management(dem_sink_buffer, sink_buffer_img)

    # add output data to map
    # arcpy.AddMessage("Adding data to map ...")
    # mxd = MapDocument("CURRENT")
    # df = ListDataFrames(mxd, "*")[0]
    # lyr_fully_filled_dem = Layer(os.path.join(workspace, dem_filled_name))
    # AddLayer(df, lyr_fully_filled_dem)
    # lyr_partially_filled_dem = Layer(
    #     os.path.join(workspace, "dem_partially_filled" + img_ext)
    # )
    # AddLayer(df, lyr_partially_filled_dem)
    # lyr_sink_dem = Layer(os.path.join(workspace, "sink" + img_ext))
    # AddLayer(df, lyr_sink_dem)

    # add output data to map
    arcpy.AddMessage("Adding data to map ...")
    p = arcpy.mp.ArcGISProject("CURRENT")
    m = p.listMaps("*")[0]
    lyr_fully_filled_dem = os.path.join(workspace, dem_filled_name)
    m.addDataFromPath(lyr_fully_filled_dem)
    lyr_partially_filled_dem = os.path.join(workspace, "dem_partially_filled" + img_ext)
    m.addDataFromPath(lyr_partially_filled_dem)
    lyr_sink_dem = os.path.join(workspace, "sink" + img_ext)
    m.addDataFromPath(lyr_sink_dem)

    arcpy.AddMessage("Deleting temporary data ...")
    arcpy.Delete_management(region_poly_name)
    arcpy.Delete_management(region_poly_ras)
    arcpy.Delete_management(ras_sink_bin_name)
    arcpy.Delete_management(os.path.join(workspace, "ras_sink_bg" + img_ext))
    arcpy.Delete_management(os.path.join(workspace, "ras_region" + img_ext))
    arcpy.Delete_management(os.path.join(workspace, "ras_region_sub" + img_ext))
    arcpy.Delete_management(os.path.join(workspace, "ras_region_zonal" + img_ext))
    if buffer_dist > 0:
        arcpy.Delete_management(sink_buffer_name)
    arcpy.AddMessage("Extract sink done!")

    end_time = time.time()
    arcpy.AddMessage("Total run time: {:.4f}".format(end_time - start_time))

    return out_sink


def zonalStatistics(in_shp_path, in_dem):
    in_shp_dir = os.path.split(in_shp_path)[0]
    in_shp_name = os.path.split(in_shp_path)[1]
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
        if shp.endswith(".shp") and shp == in_shp_name:
            shp_path = os.path.join(in_shp_dir, shp)
            dbf_path = os.path.join(dbf_dir, "zonal_" + shp.replace("shp", "dbf"))
            tif_path = os.path.join(tif_dir, shp.replace("shp", "tif"))
            arcpy.PolygonToRaster_conversion(
                shp_path,
                value_field="FID",
                out_rasterdataset=tif_path,
                cell_assignment="CELL_CENTER",
                priority_field="NONE",
                cellsize=cell_size,
            )
            arcpy.sa.ZonalStatisticsAsTable(
                tif_path, "Value", in_dem, dbf_path, "DATA", "ALL"
            )
            arcpy.JoinField_management(
                shp_path,
                in_field="FID",
                join_table=dbf_path,
                join_field="Value",
                fields="COUNT;MIN;MAX;RANGE;MEAN;STD;SUM",
            )
            # arcpy.AddField_management(shp_path,field_name="dep2catR",field_type="FLOAT")
            arcpy.AddField_management(shp_path, field_name="volume", field_type="FLOAT")
            arcpy.AddField_management(
                shp_path, field_name="mean_depth", field_type="FLOAT"
            )
            # arcpy.CalculateField_management(shp_path,field="dep2catR",expression="!AREA! / !cat_area!", expression_type="PYTHON_9.3")
            arcpy.CalculateField_management(
                shp_path,
                field="volume",
                expression="( !COUNT! * !MAX! - !SUM!) * ( !AREA! / !COUNT! )",
                expression_type="PYTHON_9.3",
            )
            arcpy.CalculateField_management(
                shp_path,
                field="mean_depth",
                expression="!volume! / !AREA!",
                expression_type="PYTHON_9.3",
            )
            arcpy.CalculateField_management(
                shp_path,
                field="ID",
                expression="!FID! + 1",
                expression_type="PYTHON_9.3",
            )
            arcpy.DeleteField_management(shp_path, drop_field="GRIDCODE")

    if os.path.exists(dbf_dir):
        shutil.rmtree(dbf_dir)
    if os.path.exists(tif_dir):
        shutil.rmtree(tif_dir)
    return True


# main script
if __name__ == "__main__":

    in_dem = arcpy.GetParameterAsText(0)
    min_size = float(arcpy.GetParameterAsText(1))
    min_depth = float(arcpy.GetParameterAsText(2))
    buffer_dist = float(arcpy.GetParameterAsText(3))
    out_sink = arcpy.GetParameterAsText(4)

    extract_sink(in_dem, min_size, min_depth, buffer_dist, out_sink)
