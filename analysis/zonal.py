# -*- coding: utf-8 -*-
"""County statistics using the original ArcMap 10.8 Spatial Analyst tools."""
from __future__ import print_function
import csv
import os
import shutil
import tempfile

import arcpy


def text_path(path):
    return path.decode("mbcs") if isinstance(path, str) else path


def measure(zones, raster, name, statistic):
    table = os.path.join(arcpy.env.workspace, name)
    arcpy.sa.ZonalStatisticsAsTable(zones, "FID", raster, table, "DATA", statistic)
    fields = ["FID", "COUNT"] + (["MEAN"] if statistic == "MEAN" else [])
    return {int(r[0]): r[1:] for r in arcpy.da.SearchCursor(table, fields)}


def main():
    root = os.path.dirname(os.path.dirname(text_path(os.path.abspath(__file__))))
    work = os.path.join(root, "run")
    topic = os.path.join(work, "source", "topic")
    zones = os.path.join(topic, "Gansu-County.shp")
    lucc = os.path.join(topic, "Gansu-lucc-2015.tif")
    rasters = os.path.join(work, "rasters")
    rain = os.path.join(rasters, "precip_OK_log_spherical_T0_N16_prediction.tif")
    evap = os.path.join(rasters, "evap_OK_log_gaussian_T0_N16_prediction.tif")
    if int(arcpy.GetCount_management(zones).getOutput(0)) != 89:
        raise ValueError("Expected 89 source county polygons")
    if arcpy.CheckOutExtension("Spatial") != "CheckedOut":
        raise RuntimeError("Spatial Analyst license required")
    scratch = tempfile.mkdtemp(prefix="zonal_", dir=work)
    workspace = os.path.join(scratch, "work.gdb")
    try:
        arcpy.CreateFileGDB_management(scratch, "work.gdb")
        arcpy.ResetEnvironments()
        arcpy.env.workspace = workspace
        land = arcpy.Describe(lucc)
        counties = arcpy.Describe(zones)
        arcpy.env.outputCoordinateSystem = land.spatialReference
        arcpy.env.cellSize = 1000
        arcpy.env.snapRaster = lucc
        arcpy.env.extent = lucc
        mask = os.path.join(workspace, "vegetation")
        arcpy.sa.ExtractByAttributes(lucc, "VALUE IN (21,22,23,24,31,32,33)").save(mask)

        a, b = counties.extent.projectAs(land.spatialReference), land.extent
        arcpy.ClearEnvironment("snapRaster")
        arcpy.env.extent = arcpy.Extent(max(a.XMin,b.XMin), max(a.YMin,b.YMin),
                                       min(a.XMax,b.XMax), min(a.YMax,b.YMax))
        reference = os.path.join(scratch, "grid.tif")
        arcpy.FeatureToRaster_conversion(zones, "FID", reference, 1000)
        arcpy.env.snapRaster = reference
        arcpy.env.extent = reference
        vegetation = measure(zones, mask, "vegetation_count", "ALL")
        total = measure(zones, lucc, "land_count", "ALL")

        arcpy.env.outputCoordinateSystem = arcpy.Describe(rain).spatialReference
        arcpy.env.cellSize = 5000
        arcpy.env.extent = rain
        arcpy.env.snapRaster = rain
        p = measure(zones, rain, "rainfall", "MEAN")
        e = measure(zones, evap, "evaporation", "MEAN")
        target = os.path.join(work, "tables", "county_zones.csv")
        if not os.path.isdir(os.path.dirname(target)):
            os.makedirs(os.path.dirname(target))
        with open(target + ".tmp", "wb") as stream:
            writer = csv.writer(stream)
            writer.writerow(["ZoneFID","ADCODE99","included","VEG_COUNT","TOTAL_N","percent",
                             "jiangshui","zhengfa","P_COUNT","E_COUNT"])
            for fid, code in arcpy.da.SearchCursor(zones, ["FID","ADCODE99"]):
                count = vegetation.get(fid, (0,))[0]
                n = total[fid][0]
                included = fid in vegetation and fid in p and fid in e
                writer.writerow([fid,code,int(included),count,n,float(count)/n if included else "",
                                 p[fid][1],e[fid][1],p[fid][0],e[fid][0]])
        # MoveFileEx provides replace-existing semantics in ArcMap's Python 2.7.
        import ctypes
        if not ctypes.windll.kernel32.MoveFileExW(target + u".tmp", target, 1):
            raise ctypes.WinError()
        print("County zones:", len(total), "used:", len(vegetation))
    finally:
        arcpy.ResetEnvironments()
        arcpy.ClearWorkspaceCache_management()
        arcpy.CheckInExtension("Spatial")
        if os.path.dirname(os.path.realpath(scratch)) == os.path.realpath(work):
            shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    main()
