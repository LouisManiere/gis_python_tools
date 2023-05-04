import os
from raster_tools.raster_tools import extract_raster_extent
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()


extract_raster_extent(
    input_dir_path = paths['inputs_dir'],
    extension = params['raster_extension'],
    tileset_path = os.path.join(paths['outputs_dir'], paths['tileset']),
    crs = params['crs']
)