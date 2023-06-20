import os
from raster_tools.raster_tools import ExtractRasterTilesFromTileset
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()


ExtractRasterTilesFromTileset(
    tileset_path = paths['vector_in_1'],
    raster_dir = paths['inputs_dir'],
    dest_dir = paths['outputs_dir']
)