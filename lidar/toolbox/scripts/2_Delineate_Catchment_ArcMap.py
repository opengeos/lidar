import arcpy
import os
import time
from arcpy.mapping import *


def delineate_catchment(in_dem, in_sink, out_catchment):
    arcpy.CheckOutExtension("Spatial")
    workspace = os.path.split(out_catchment)[0]
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True

    if arcpy.Exists(in_dem) == False:
        arcpy.AddMessage("The input raster does not exist")
        quit()

    if os.path.splitext(out_catchment)[1].lower() == ".shp":
        # FieldOID = "FID"
        FieldOID = "ID"
        FlowDirection = os.path.join(workspace, "FlowDirection.tif")
        SinkRaster = os.path.join(workspace, "SinkRaster.tif")
        Watershed = os.path.join(workspace, "Watershed.tif")
        Catchment_tmp = os.path.join(workspace, "Catchment_tmp.shp")
        Catchment_select = os.path.join(workspace, "Catchment_select.shp")
        # Catchment_dissolve = os.path.join(workspace, "Catchment.shp")
    else:
        FieldOID = "OBJECTID"
        FlowDirection = os.path.join(workspace, "FlowDirection")
        SinkRaster = os.path.join(workspace, "SinkRaster")
        Watershed = os.path.join(workspace, "Watershed")
        Catchment_tmp = os.path.join(workspace, "Catchment")

    input_dem = arcpy.Raster(in_dem)
    flow_direction = arcpy.sa.FlowDirection(input_dem)
    flow_direction.save(FlowDirection)

    cell_size = input_dem.meanCellWidth
    arcpy.env.extent = input_dem.extent
    arcpy.PolygonToRaster_conversion(
        in_sink, FieldOID, SinkRaster, "CELL_CENTER", "NONE", cell_size
    )

    watershed = arcpy.sa.Watershed(flow_direction, SinkRaster, "Value")
    watershed.save(Watershed)

    arcpy.RasterToPolygon_conversion(watershed, Catchment_tmp, "NO_SIMPLIFY", "Value")
    field = "GRIDCODE"
    sqlExp = field + ">" + str(0)
    arcpy.Select_analysis(Catchment_tmp, Catchment_select, sqlExp)
    arcpy.Dissolve_management(
        Catchment_select,
        out_catchment,
        dissolve_field="GRIDCODE",
        statistics_fields="",
        multi_part="MULTI_PART",
        unsplit_lines="DISSOLVE_LINES",
    )

    area_field = "cat_area"
    arcpy.AddField_management(out_catchment, area_field, "DOUBLE")
    arcpy.CalculateField_management(
        out_catchment, area_field, "!shape.area@squaremeters!", "PYTHON_9.3", "#"
    )

    arcpy.JoinField_management(
        in_sink, in_field="ID", join_table=Watershed, join_field="Value", fields="Count"
    )
    arcpy.AddField_management(in_sink, field_name="cat_area", field_type="FLOAT")
    arcpy.CalculateField_management(
        in_sink,
        field="cat_area",
        expression="!Count_1! * math.pow(" + str(cell_size) + ",2)",
        expression_type="PYTHON_9.3",
    )
    arcpy.DeleteField_management(in_sink, drop_field="Count_1")
    arcpy.AddField_management(in_sink, field_name="dep2catR", field_type="FLOAT")
    arcpy.CalculateField_management(
        in_sink,
        field="dep2catR",
        expression="!AREA! / !cat_area!",
        expression_type="PYTHON_9.3",
    )
    arcpy.Delete_management(Catchment_tmp)
    arcpy.Delete_management(Catchment_select)

    # #add output data to map
    # mxd = MapDocument("CURRENT")
    # df = ListDataFrames(mxd, "*")[0]
    # lyr_watershed = Layer(Watershed)
    # AddLayer(df, lyr_watershed)

    return out_catchment


# main script
if __name__ == "__main__":

    in_dem = arcpy.GetParameterAsText(0)
    in_sink = arcpy.GetParameterAsText(1)
    out_catchment = arcpy.GetParameterAsText(2)

    start_time = time.time()
    delineate_catchment(in_dem, in_sink, out_catchment)
    end_time = time.time()
    arcpy.AddMessage("Total run time: {:.4f}".format(end_time - start_time))
