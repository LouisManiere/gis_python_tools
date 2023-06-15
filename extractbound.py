from vector_tools.vector_tools import ExtractByBoundMask
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()

ExtractByBoundMask(paths['tileset_landuse'], paths['mask'], paths['tileset_mask_landuse'])