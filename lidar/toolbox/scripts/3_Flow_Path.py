import arcpy
import os
import sys
import time
import string
from arcpy import env
import collections


def FlowPath(in_dem, in_sink, rain_intensity, out_flowpath):
    arcpy.CheckOutExtension("Spatial")
    workspace = os.path.split(out_flowpath)[0]
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True
    dem = arcpy.Raster(in_dem)
    cell_size = dem.meanCellWidth

    if arcpy.Exists(in_dem) == False:
        arcpy.AddMessage("The input raster does not exist")
        quit()

    if (os.path.splitext(out_flowpath)[1].lower() == ".shp"):
        FieldOID = "FID"
        FlowDir = os.path.join(workspace,"FlowDir.tif")
        SinkCentroid = os.path.join(workspace,"SinkCentroid.shp")
        CostPath = os.path.join(workspace,"CostPath.tif")
        PathThin = os.path.join(workspace,"PathThin.tif")
        PathLine = os.path.join(workspace,"PathLine.shp")
        PathLineErase = os.path.join(workspace,"PathLineErase.shp")
        Path = os.path.join(workspace,"FlowPath_Raw.shp")
        # Path = out_flowpath
        LineFlip = os.path.join(workspace,"LineFlip.shp")
        LineNoFlip = os.path.join(workspace,"LineNoFlip.shp")
        FlowFrom = os.path.join(workspace,"FlowFrom.shp")
        FlowTo = os.path.join(workspace,"FlowTo.shp")
        PathLineEraseSingle = os.path.join(workspace,"PathLineEraseSingle.shp")
        LineStart = os.path.join(workspace,"LineStart.shp")
        LineEnd = os.path.join(workspace, "LineEnd.shp")
        LineStartElev = os.path.join(workspace,"LineStartElev.shp")
        LineEndElev = os.path.join(workspace, "LineEndElev.shp")
        PathBuffer = os.path.join(workspace,"PathBuffer.shp")
        PathBufferSingle = os.path.join(workspace,"PathBufferSingle.shp")
        FlowFromJoin = os.path.join(workspace, "FlowFromJoin.shp")
        FlowToJoin = os.path.join(workspace, "FlowToJoin.shp")
        FlowFromJoinBuffer = os.path.join(workspace, "FlowFromJoinBuffer.shp")
        FlowToJoinBuffer = os.path.join(workspace, "FlowToJoinBuffer.shp")

    else:
        FieldOID = "OBJECTID"
        FlowDir = os.path.join(workspace,"FlowDir")
        SinkCentroid = os.path.join(workspace,"SinkCentroid")
        CostPath = os.path.join(workspace,"CostPath")
        PathThin = os.path.join(workspace,"PathThin")
        PathLine = os.path.join(workspace,"PathLine")
        PathLineErase = os.path.join(workspace,"PathLineErase")
        Path = os.path.join(workspace,"FlowPath")
        LineFlip = os.path.join(workspace,"LineFlip")
        LineNoFlip = os.path.join(workspace,"LineNoFlip")
        FlowFrom = os.path.join(workspace,"FlowFrom")
        FlowTo = os.path.join(workspace,"FlowTo")
        LineStart = os.path.join(workspace, "LineStart.shp")
        LineEnd = os.path.join(workspace, "LineEnd.shp")
    ### Delineate flow direction
    flow_dir = arcpy.sa.FlowDirection(in_dem)
    flow_dir.save(FlowDir)

    ### Extract the depression polygon centroids
    arcpy.FeatureToPoint_management(in_sink, SinkCentroid, "INSIDE")

    ### Delineate cost path
    cost_path = arcpy.sa.CostPath(SinkCentroid, in_dem, FlowDir, "EACH_CELL", FieldOID)
    cost_path.save(CostPath)

    ### Thin the raster cost path to single-cell width
    path_thin = arcpy.sa.Thin(cost_path,"#","#","#",1)
    path_thin.save(PathThin)

    ### Convert the raster path to vector
    arcpy.RasterToPolyline_conversion(path_thin,PathLine,simplify="NO_SIMPLIFY")

    ### Erase the flow path within depression polygons
    arcpy.Erase_analysis(PathLine, in_sink, PathLineErase)
    arcpy.MultipartToSinglepart_management(PathLineErase, PathLineEraseSingle)
    arcpy.FeatureVerticesToPoints_management(PathLineEraseSingle,LineStart,"START")
    arcpy.FeatureVerticesToPoints_management(PathLineEraseSingle,LineEnd,"END")
    arcpy.sa.ExtractValuesToPoints(LineStart,in_dem,LineStartElev)
    arcpy.sa.ExtractValuesToPoints(LineEnd, in_dem, LineEndElev)
    arcpy.AddField_management(LineStartElev,field_name="FromElev", field_type="FLOAT")
    arcpy.AddField_management(LineEndElev, field_name="ToElev", field_type="FLOAT")
    arcpy.CalculateField_management(in_table=LineStartElev, field="FromElev", expression="!RASTERVALU!",
                                    expression_type="PYTHON", code_block="")
    arcpy.CalculateField_management(in_table=LineEndElev, field="ToElev", expression="!RASTERVALU!",
                                    expression_type="PYTHON", code_block="")
    arcpy.JoinField_management(in_data=PathLineEraseSingle, in_field="FID", join_table=LineStartElev, join_field="FID",
                               fields="FromElev")
    arcpy.JoinField_management(in_data=PathLineEraseSingle, in_field="FID", join_table=LineEndElev, join_field="FID",
                               fields="ToElev")
    arcpy.CopyFeatures_management(PathLineEraseSingle, Path)
    # ExtractElevation(PathLineErase, in_dem, Path)

    arcpy.AddField_management(Path,"Flip","SHORT")

    FromElev = arcpy.AddFieldDelimiters(workspace,"FromElev")
    ToElev = arcpy.AddFieldDelimiters(workspace,"ToElev")
    sql = FromElev + "<" + ToElev
    sql2 = FromElev + ">=" + ToElev

    arcpy.Select_analysis(Path,LineFlip,sql)
    arcpy.CalculateField_management(LineFlip,"Flip","1","PYTHON")
    arcpy.FlipLine_edit(LineFlip)
    arcpy.Select_analysis(Path,LineNoFlip,sql2)

    arcpy.Delete_management(Path)
    MergeList = []
    MergeList.append(LineFlip)
    MergeList.append(LineNoFlip)
    arcpy.Merge_management(MergeList,Path)
    arcpy.AddField_management(Path,field_name="StartElev", field_type="FLOAT")
    arcpy.AddField_management(Path, field_name="EndElev", field_type="FLOAT")
    arcpy.AddField_management(Path,field_name="DiffElev", field_type="FLOAT")
    arcpy.AddField_management(Path, field_name="Length", field_type="FLOAT")
    arcpy.CalculateField_management(in_table=Path, field="StartElev", expression="max( !FromElev! , !ToElev! )",
                                    expression_type="PYTHON", code_block="")
    arcpy.CalculateField_management(in_table=Path, field="EndElev", expression="min( !FromElev! , !ToElev! )",
                                    expression_type="PYTHON", code_block="")
    arcpy.CalculateField_management(in_table=Path, field="DiffElev", expression="!StartElev! - !EndElev!",
                                    expression_type="PYTHON", code_block="")
    arcpy.CalculateField_management(Path, "Length", "!shape.length@meters!", "PYTHON_9.3", "#")
    arcpy.DeleteField_management(in_table=Path,
                                 drop_field="ARCID;GRID_CODE;FROM_NODE;TO_NODE;ORIG_FID;FromElev;ToElev;Flip")
    sql3 = "Length >" + str(2 * cell_size)  # if flow path is shorter than 2 pixels, delete
    arcpy.Select_analysis(Path, out_flowpath, sql3)

    arcpy.FeatureVerticesToPoints_management(out_flowpath,FlowFrom,"START")
    arcpy.FeatureVerticesToPoints_management(out_flowpath,FlowTo,"END")
    arcpy.AddField_management(FlowFrom,field_name="FlowFromID", field_type="Long")
    arcpy.AddField_management(FlowTo, field_name="FlowToID", field_type="Long")
    arcpy.CalculateField_management(in_table=FlowFrom, field="FlowFromID", expression="!FID! + 1",
                                    expression_type="PYTHON", code_block="")
    arcpy.CalculateField_management(in_table=FlowTo, field="FlowToID", expression="!FID! + 1",
                                    expression_type="PYTHON", code_block="")
    # derive sink connectivity

    arcpy.Buffer_analysis(in_features=Path, out_feature_class=PathBuffer,
                          buffer_distance_or_field="0.1 Meters", line_side="FULL", line_end_type="FLAT",
                          dissolve_option="ALL", dissolve_field="", method="PLANAR")
    arcpy.MultipartToSinglepart_management(in_features=PathBuffer,out_feature_class=PathBufferSingle)
    arcpy.AddField_management(PathBufferSingle, field_name="BufferID", field_type="Long")
    arcpy.CalculateField_management(in_table=PathBufferSingle, field="BufferID", expression="!FID! + 1",
                                    expression_type="PYTHON", code_block="")


    search_radius = str(2.1 * cell_size) + " Meters"
    arcpy.SpatialJoin_analysis(target_features=FlowFrom, join_features=in_sink,
                               out_feature_class=FlowFromJoin,
                               join_operation="JOIN_ONE_TO_ONE", join_type="KEEP_COMMON",
                               # field_mapping="""ID "ID" true true false 10 Long 0 10 ,First,#,poly,ID,-1,-1""",
                               match_option="INTERSECT", search_radius=search_radius, distance_field_name="")
    arcpy.SpatialJoin_analysis(target_features=FlowTo, join_features=in_sink,
                               out_feature_class=FlowToJoin,
                               join_operation="JOIN_ONE_TO_ONE", join_type="KEEP_COMMON",
                               # field_mapping="""ID "ID" true true false 10 Long 0 10 ,First,#,poly,ID,-1,-1""",
                               match_option="INTERSECT", search_radius=search_radius, distance_field_name="")
    arcpy.SpatialJoin_analysis(target_features=FlowFromJoin, join_features=PathBufferSingle,
                               out_feature_class=FlowFromJoinBuffer,
                               join_operation="JOIN_ONE_TO_ONE", join_type="KEEP_COMMON",
                               # field_mapping="""ID "ID" true true false 10 Long 0 10 ,First,#,poly,ID,-1,-1""",
                               match_option="INTERSECT", search_radius=search_radius, distance_field_name="")
    arcpy.SpatialJoin_analysis(target_features=FlowToJoin, join_features=PathBufferSingle,
                               out_feature_class=FlowToJoinBuffer,
                               join_operation="JOIN_ONE_TO_ONE", join_type="KEEP_COMMON",
                               # field_mapping="""ID "ID" true true false 10 Long 0 10 ,First,#,poly,ID,-1,-1""",
                               match_option="INTERSECT", search_radius=search_radius, distance_field_name="")
    arcpy.JoinField_management(in_data=FlowFromJoinBuffer, in_field="BufferID",
                               join_table=FlowToJoinBuffer, join_field="BufferID", fields="ID")
    arcpy.JoinField_management(in_data=in_sink, in_field="ID", join_table=FlowFromJoinBuffer, join_field="ID",
                               fields="ID_12")
    arcpy.AddField_management(in_sink, field_name="Downstream", field_type="LONG")
    arcpy.CalculateField_management(in_table=in_sink, field="Downstream", expression="!ID_12!",
                                    expression_type="PYTHON", code_block="")
    arcpy.DeleteField_management(in_table=in_sink, drop_field="ID_12")

    arcpy.AddField_management(in_sink, field_name="simu_depth", field_type="FLOAT")
    arcpy.AddField_management(in_sink, field_name="rain_inten", field_type="FLOAT")
    arcpy.AddField_management(in_sink, field_name="time_inund", field_type="FLOAT")
    arcpy.CalculateField_management(in_table=in_sink, field="simu_depth", expression="!volume! / !cat_area!",
                                    expression_type="PYTHON", code_block="")
    arcpy.CalculateField_management(in_table=in_sink, field="rain_inten", expression=rain_intensity,
                                    expression_type="PYTHON", code_block="")
    arcpy.CalculateField_management(in_table=in_sink, field="time_inund", expression="!simu_depth! / !rain_inten!",
                                    expression_type="PYTHON", code_block="")

    arcpy.JoinField_management(in_data=out_flowpath, in_field="FID", join_table=FlowFromJoin, join_field="ORIG_FID",
                               fields="ID")
    arcpy.AddField_management(in_table=out_flowpath, field_name="start_sink", field_type="LONG")
    arcpy.CalculateField_management(in_table=out_flowpath, field="start_sink", expression="!ID!",
                                    expression_type="PYTHON", code_block="")
    arcpy.DeleteField_management(in_table=out_flowpath, drop_field="ID")
    arcpy.JoinField_management(in_data=out_flowpath, in_field="FID", join_table=FlowToJoin, join_field="ORIG_FID",
                               fields="ID")
    arcpy.AddField_management(in_table=out_flowpath, field_name="end_sink", field_type="LONG")
    arcpy.CalculateField_management(in_table=out_flowpath, field="end_sink", expression="!ID!",
                                    expression_type="PYTHON", code_block="")
    arcpy.DeleteField_management(in_table=out_flowpath, drop_field="ID")
    arcpy.JoinField_management(in_data=out_flowpath, in_field="start_sink", join_table=in_sink, join_field="ID",
                               fields="volume;cat_area;simu_depth;rain_inten;time_inund")



    arcpy.Delete_management(LineFlip)
    arcpy.Delete_management(LineNoFlip)
    arcpy.Delete_management(CostPath)
    arcpy.Delete_management(FlowDir)
    arcpy.Delete_management(PathLineErase)
    arcpy.Delete_management(PathThin)
    arcpy.Delete_management(LineStart)
    arcpy.Delete_management(LineStartElev)
    arcpy.Delete_management(LineEnd)
    arcpy.Delete_management(LineEndElev)
    arcpy.Delete_management(PathLineEraseSingle)
    arcpy.Delete_management(SinkCentroid)
    arcpy.Delete_management(PathLine)
    arcpy.Delete_management(PathBuffer)
    arcpy.Delete_management(PathBufferSingle)
    arcpy.Delete_management(FlowFromJoin)
    arcpy.Delete_management(FlowToJoin)
    arcpy.Delete_management(FlowFromJoinBuffer)
    arcpy.Delete_management(FlowToJoinBuffer)


    arcpy.AddMessage("Flow path delineation done!")
    return out_flowpath


def delete_row(in_shp, field):
    fields = []
    fields.append(field)
    with arcpy.da.UpdateCursor(in_shp,fields) as cursor:
        for row in cursor:
            if row[0] == 0:
                cursor.deleteRow()

def add_rank(in_shp, sort_field, rank_field):
    try:
        arcpy.AddField_management(in_table=in_shp, field_name=rank_field, field_type="LONG")
        inDict = dict()
        with arcpy.da.SearchCursor(in_shp, sort_field) as cursor:
            for row in cursor:
                # arcpy.AddMessage(row[0])
                if row[0] not in inDict.keys():
                    inDict[row[0]] = 0
        od = collections.OrderedDict(sorted(inDict.items()))
        i = 1
        for key, value in od.items():
            od[key] = i
            i += 1
        fields = []
        fields.append(sort_field)
        fields.append(rank_field)
        with arcpy.da.UpdateCursor(in_shp, fields) as cursor:
            for row in cursor:
                row[1] = od[row[0]]
                cursor.updateRow(row)
    except:
        arcpy.GetMessages()



if __name__ == '__main__':

    in_dem = arcpy.GetParameterAsText(0)
    in_sink = arcpy.GetParameterAsText(1)
    rain_intensity = float(arcpy.GetParameterAsText(2)) / 100  # convert cm/h to m/h
    out_flowpath = arcpy.GetParameterAsText(3)

    start_time = time.time()
    FlowPath(in_dem, in_sink, rain_intensity, out_flowpath)
    delete_row(out_flowpath,"volume")
    add_rank(out_flowpath, sort_field="time_inund", rank_field="rank")
    add_rank(in_sink, sort_field="time_inund", rank_field="rank")
    end_time = time.time()
    arcpy.AddMessage("Total run time: {:.4f}".format(end_time - start_time))

