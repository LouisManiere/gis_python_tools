from vector_tools.vector_tools import identify_network_nodes
from config.config import paths_config, parameters_config

# parameters
paths = paths_config()
params = parameters_config()

identify_network_nodes(paths['test'], paths['output_test'], paths['output_test_nodes'])