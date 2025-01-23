# -*- coding: utf-8 -*-

import arcpy
from arcpy.sa import *
import os

def batch_mask_extraction(tif_folder, mask_shp, output_folder):
    """
    Batch mask extraction for multiple TIFF files using a single SHP file as mask.

    Parameters:
    tif_folder: Path to folder containing TIFF files
    mask_shp: Path to the shapefile used as mask
    output_folder: Output folder path for masked results
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
    arcpy.env.workspace = tif_folder

    # Get list of TIFF files in the folder
    tif_files = [f for f in os.listdir(tif_folder) if f.endswith('.tif')]

    # Iterate through each TIFF file and apply mask
    for tif_file in tif_files:
        try:
            # Construct full paths
            tif_path = os.path.join(tif_folder, tif_file)
            output_name = os.path.splitext(tif_file)[0] + "_masked.tif"
            output_path = os.path.join(output_folder, output_name)

            # Perform mask extraction
            outExtractByMask = ExtractByMask(tif_path, mask_shp)
            outExtractByMask.save(output_path)

            print("Completed mask extraction for {}".format(tif_file))

        except Exception as e:
            print("Error processing {}: {}".format(tif_file, str(e)))
            continue

    # Release extension license
    arcpy.CheckInExtension("Spatial")

    print("All mask extractions completed")


# Usage example
if __name__ == "__main__":
    # Set input parameters
    tif_folder = r"F:\firstStepOut"  # 替换为包含TIFF文件的文件夹路径
    mask_shp = r"F:\rawData\A100.shp"      # 替换为用作掩膜的SHP文件路径
    output_folder = r"F:\secondStepOut\viewshed100"   # 替换为输出文件夹路径

    # Run batch mask extraction
    batch_mask_extraction(tif_folder, mask_shp, output_folder)
