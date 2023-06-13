from vector_tools.vector_tools import ExtractBylocation
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()

ExtractBylocation(paths['test'], paths['mask_test'], paths['test_extracted'], method = 'contains')