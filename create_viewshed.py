# -*- coding: utf-8 -*-

import arcpy
import numpy as np
from arcpy import env
from arcpy.sa import *
import os


def array_to_dem(array_data, output_path, lower_left_corner, cell_size):
    """
    Convert 2D array to DEM raster file with CGCS2000 3-degree GK Zone 39 projection

    Parameters:
    array_data: numpy 2D array containing elevation values
    output_path: full path for output DEM file (.tif format)
    lower_left_corner: coordinates of lower left corner, tuple format (x, y)
    cell_size: raster cell size

    Returns:
    path to output file
    """
    try:
        # Check and create output directory if it doesn't exist
        output_dir = os.path.dirname(output_path)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Check input array
        if not isinstance(array_data, np.ndarray):
            array_data = np.array(array_data)

        if len(array_data.shape) != 2:
            raise ValueError("Input must be a 2D array")

        # Get array dimensions
        rows, cols = array_data.shape

        # Create new raster object
        ll_x, ll_y = lower_left_corner
        refPoint = arcpy.Point(ll_x, ll_y)

        # Create spatial reference for CGCS2000_3_Degree_GK_CM_117E
        # Method 1: Using factory code
        # spatial_ref = arcpy.SpatialReference(4547)  # CGCS2000 code

        # Method 2: Using projection file or WKT
        wkt = '''PROJCS["CGCS2000_3_Degree_GK_CM_117E",
                GEOGCS["GCS_China_Geodetic_Coordinate_System_2000",
                DATUM["D_China_2000",
                SPHEROID["CGCS2000",6378137.0,298.257222101]],
                PRIMEM["Greenwich",0.0],
                UNIT["Degree",0.0174532925199433]],
                PROJECTION["Gauss_Kruger"],
                PARAMETER["False_Easting",500000.0],
                PARAMETER["False_Northing",0.0],
                PARAMETER["Central_Meridian",117.0],
                PARAMETER["Scale_Factor",1.0],
                PARAMETER["Latitude_Of_Origin",0.0],
                UNIT["Meter",1.0]]'''
        spatial_ref = arcpy.SpatialReference()
        spatial_ref.loadFromString(wkt)

        # Create raster object
        raster = arcpy.NumPyArrayToRaster(
            array_data,
            refPoint,
            cell_size,
            cell_size,
            value_to_nodata=None
        )

        # Define projection
        arcpy.DefineProjection_management(raster, spatial_ref)

        # Save to file - use workspace environment
        env.workspace = output_dir
        output_name = os.path.basename(output_path)
        raster.save(output_name)

        return output_path

    except Exception as e:
        print("Error occurred: {0}".format(str(e)))
        raise


# Usage example
if __name__ == "__main__":
    try:
        # Create sample data
        # test_array = np.random.randint(low=0, high=2, size=(1775, 2059))
        test_array = np.loadtxt('F:\compute_verification\location/wall18.txt', dtype=int)
        # Set parameters - using current directory for output
        current_dir = os.path.dirname(os.path.abspath(__file__))
        output_dem = os.path.join(current_dir, "F:\compute_verification\location\wall18Location.tif")

        # Coordinates in CGCS2000_3_Degree_GK_CM_117E system
        lower_left = (513249.020902, 4500154.64174)  # Example coordinates
        cell_size = 12.5  # 100 meters resolution

        # Convert to DEM
        result = array_to_dem(test_array, output_dem, lower_left, cell_size)
        print("DEM file saved to: {0}".format(result))

    except Exception as e:
        print("Error in main: {0}".format(str(e)))