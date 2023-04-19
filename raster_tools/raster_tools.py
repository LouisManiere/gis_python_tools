# coding: utf-8

"""
DOCME

***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

# packages
import rasterio
from rasterio.merge import merge
import glob
import os

def merge_raster_in_folder(
        input_raster_folder_path: str,
        output_merge_raster_folder_path: str,
        extension: str,
        compress: str = 'deflate',
        zlevel: int = 9):
    """
    Merge all the raster file contain in a folder in Gtiff
    documentation : https://automating-gis-processes.github.io/CSC18/lessons/L6/raster-mosaic.html
    """
    # Make a search criteria to select the tiff raster files
    search_criteria = '*'+extension

    q = os.path.join(input_raster_folder_path, search_criteria)

    # list all raster file
    rasterInFolder = glob.glob(q)

    # create an empty list to gather the list of the raster name
    rasterListToMerge = []

    # read all the raster with rasterio and add it to the rasterListToMerge list
    for raster in rasterInFolder:
        item = rasterio.open(raster)
        rasterListToMerge.append(item)

    # Merge function
    mergeRaster, mergeRaster_trans = merge(datasets = rasterListToMerge, method='max')

    # Copy the metadata
    out_meta = rasterListToMerge[0].meta.copy()

    # Update the metadata
    out_meta.update({"driver": 'GTiff',
                    "height": mergeRaster.shape[1],
                    "width": mergeRaster.shape[2],
                    "transform": mergeRaster_trans,
                    "compress": compress,
                    "zlevel": zlevel
                    }
                    )
    
    # Write the mosaic raster to disk
    with rasterio.open(output_merge_raster_folder_path + 'merge_raster.tif', "w", **out_meta) as dest:
        dest.write(mergeRaster)