import os
from raster_tools.raster_tools import CreateTilesetFromRasters, ExtractRasterTilesFromTileset
from vector_tools.CreateTilesetFromExtent import CreateTilesetFromExtent
from config.config import paths_config, parameters_config
import subprocess
import geopandas as gpd

# parameters
paths = paths_config()
params = parameters_config()

# create dem and landuse tileset

CreateTilesetFromRasters(
    input_dir_path = paths['inputs_dir_landuse_tiles'],
    extension = params['landuse_extension'],
    tileset_path = paths['tileset_landuse'],
    crs = params['crs']
)

CreateTilesetFromRasters(
    input_dir_path = paths['inputs_dir_dem_tiles'],
    extension = params['dem_extension'],
    tileset_path = paths['tileset_dem'],
    crs = params['crs']
)

# create new tileset in the mask extent

CreateTilesetFromExtent(
    tile_size = params['tile_size'],
    study_area_path = paths['mask'],
    tileset_path = paths['tileset_mask_landuse'],
    crs = params['crs']
)

CreateTilesetFromExtent(
    tile_size = params['tile_size'],
    study_area_path = paths['mask'],
    tileset_path = paths['tileset_mask_dem'],
    crs = params['crs']
)

# copy raster tiles

ExtractRasterTilesFromTileset(
    tileset_path = paths['tileset_mask_landuse'],
    raster_dir = paths['inputs_dir_landuse_tiles'],
    dest_dir = paths['outputs_dir_landuse_tiles']
)

ExtractRasterTilesFromTileset(
    tileset_path = paths['tileset_mask_dem'],
    raster_dir = paths['inputs_dir_dem_tiles'],
    dest_dir = paths['outputs_dir_dem_tiles']
)

# create virtual raster

bash_landuse = 'gdalbuildvrt "{}" "{}"*"{}"'.format(paths['landuse_vrt'], paths['outputs_dir_landuse_tiles'], params['landuse_extension'])
vrt_landuse = subprocess.run(bash_landuse, shell=True, capture_output=True, text=True)

bash_dem = 'gdalbuildvrt "{}" "{}"*"{}"'.format(paths['dem_vrt'], paths['outputs_dir_dem_tiles'], params['dem_extension'])
vrt_dem = subprocess.run(bash_dem, shell=True, capture_output=True, text=True)

# Clip reference hydrologic network

# hydro_network = gpd.read_file(paths['hydro_network'])

# hydro_network_mask = hydro_network.clip(mask = paths['mask'], keep_geom_type=True)

# hydro_network_mask.to_file(paths['hydro_network_mask'], driver="GPKG")
