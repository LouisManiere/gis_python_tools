
from raster_tools.raster_tools import merge_raster_in_folder
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()


merge_raster_in_folder(
    input_raster_folder_path = paths['inputs_dir'],
    output_merge_raster_folder_path = paths['outputs_dir'],
    extension = params['raster_extension']
)