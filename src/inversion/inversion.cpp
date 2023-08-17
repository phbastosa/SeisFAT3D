# include "inversion.hpp"

void Inversion::set_general_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", file));

    obs_data_folder = catch_parameter("obs_data_folder", file);
    obs_data_prefix = catch_parameter("obs_data_prefix", file);

    convergence_map_folder = catch_parameter("convergence_map_folder", file);
    estimated_model_folder = catch_parameter("estimated_model_folder", file);

    model_per_iteration = str2bool(catch_parameter("model_per_iteration", file));
}