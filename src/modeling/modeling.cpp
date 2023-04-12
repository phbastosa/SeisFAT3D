# include "modeling.hpp"

void Modeling::set_parameters(std::string file)
{
    std::vector<std::string> splitted;

    splitted = split(catch_parameter("model_samples", file), ',');
    
    nz = std::stoi(splitted[0]);
    nx = std::stoi(splitted[1]);
    ny = std::stoi(splitted[2]);

    splitted = split(catch_parameter("model_spacing", file), ',');

    dz = std::stof(splitted[0]);
    dx = std::stof(splitted[1]);
    dy = std::stof(splitted[2]);

    std::vector<std::string>().swap(splitted);

    export_receiver_output = str2bool(catch_parameter("export_seismogram", file));
    export_wavefield_output = str2bool(catch_parameter("export_snapshots", file));

    receiver_output_folder = catch_parameter("seismogram_folder", file);
    wavefield_output_folder = catch_parameter("snapshots_folder", file);

    Geometry * geometry[] = {new Regular(), new Circular(), new Streamer()};

    int type = std::stoi(catch_parameter("geometry_type", file));

    geometry[type]->set_geometry(file);

    total_shots = geometry[type]->shots.total;
    total_nodes = geometry[type]->nodes.total;
}

void Modeling::set_runtime()
{
    ti = std::chrono::system_clock::now();
}

void Modeling::get_runtime()
{
    tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout<<"\nRun time: "<<elapsed_seconds.count()<<" s."<<std::endl;
}

void Modeling::export_outputs()
{
    if (export_receiver_output) export_binary_float(receiver_output_file, receiver_output, receiver_output_samples);
    if (export_wavefield_output) export_binary_float(wavefield_output_file, wavefield_output, wavefield_output_samples);
}