import os
from raster_tools.raster_tools import extract_raster_extent
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()


CreateTilesetFromRasters(
    input_dir_path = paths['inputs_dir'],
    extension = params['raster_extension'],
    tileset_path = paths['vector_in_1'],
    crs = params['crs']
)