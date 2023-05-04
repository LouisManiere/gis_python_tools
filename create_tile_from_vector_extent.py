
from vector_tools.create_tile_from_extent import CreateTileset
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()

CreateTileset(
    tile_size = 10000,
    study_area_path = './inputs/extent_test_fct.gpkg',
    tileset_path = './outputs/tileset.gpkg',
    crs = '2154'
)

