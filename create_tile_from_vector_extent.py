
from vector_tools.create_tile_from_extent import CreateTileset
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()

CreateTileset(
    tile_size = 10000,
    study_area_path = paths['vector_in_1'],
    tileset_path = paths['vector_in_2'],
    crs = '2154'
)

