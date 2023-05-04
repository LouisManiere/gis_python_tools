
from raster_tools.raster_tools import extract_raster_tiles_from_tileset
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()


extract_raster_tiles_from_tileset(
    tileset_path = './inputs/tileset.gpkg',
    raster_dir = paths['inputs_dir'],
    dest_dir = paths['outputs_dir']
)