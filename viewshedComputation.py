# -*- coding: utf-8 -*-

import arcpy
from arcpy.sa import *
import os

def batch_viewshed(dem_path, observer_points, output_folder,
                   radius=None, observer_height=1.7, target_height=0):
    """
    Batch calculate viewsheds for multiple observation points

    Parameters:
    dem_path: Digital Elevation Model path
    observer_points: Path to observer points feature class
    output_folder: Output folder path
    radius: Analysis radius in meters (None for unlimited)
    observer_height: Observer height in meters
    target_height: Target height in meters
    """

    # Check Spatial Analyst extension
    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
    else:
        raise Exception("Spatial Analyst extension required")

    # Ensure output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Set environment settings
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = output_folder

    # Get DEM as raster object
    dem = Raster(dem_path)

    # Set the analysis environment settings for observer and target heights
    arcpy.env.observerHeight = observer_height
    arcpy.env.targetHeight = target_height

    # Iterate through each observer point
    with arcpy.da.SearchCursor(observer_points, ["SHAPE@", "OID@"]) as cursor:
        for row in cursor:
            point = row[0]
            point_id = row[1]

            try:
                # Calculate viewshed
                if radius is None:
                    # Calculate unlimited viewshed
                    viewshed = Viewshed(dem, point)
                else:
                    # Calculate viewshed with specified radius
                    viewshed = Viewshed(dem, point, outer_radius=radius)

                # Save results
                output_name = "viewshed_{0}.tif".format(point_id)
                output_path = os.path.join(output_folder, output_name)
                viewshed.save(output_path)

                print("Completed viewshed analysis for point {0}".format(point_id))

            except Exception as e:
                print("Error processing point {0}: {1}".format(point_id, str(e)))
                continue

    # Release extension license
    arcpy.CheckInExtension("Spatial")

    print("All viewshed analyses completed")


# Usage example
if __name__ == "__main__":
    # Set input parameters
    dem_path = r"F:\rawData\dem.tif"
    observer_points = r"F:\rawData\GreatWall100.shp"
    output_folder = r"F:\firstStepOut"

    # Run batch viewshed analysis without radius (unlimited viewshed)
    batch_viewshed(dem_path, observer_points, output_folder)

    # Or run with specified radius
    # batch_viewshed(dem_path, observer_points, output_folder, radius=5000)
