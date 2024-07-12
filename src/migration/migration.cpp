# include "migration.hpp"

void Migration::set_parameters()
{
    set_modeling_type();

    modeling->file = file;

    modeling->set_parameters();

    nr = modeling->total_nodes;
    ns = modeling->total_shots;

    nt = std::stoi(catch_parameter("time_samples", file));
    dt = std::stof(catch_parameter("time_spacing", file));

    input_data_folder = catch_parameter("input_data_folder", file);
    input_data_prefix = catch_parameter("input_data_prefix", file);

    data = new float[nt*nr]();

    Tr = new float[modeling->nPoints]();
    Ts = new float[modeling->nPoints]();
    Im = new float[modeling->nPoints]();

    dTx = new float[modeling->nPoints]();
    dTy = new float[modeling->nPoints]();
    dTz = new float[modeling->nPoints]();

    image = new float[modeling->nPoints]();
}

void Migration::read_input_data()
{
    import_binary_float(input_data_folder + input_data_prefix + std::to_string(modeling->shot_index+1) + ".bin", data, nt*nr);
}

