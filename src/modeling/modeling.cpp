# include "modeling.hpp"

void Modeling::set_parameters(std::string file)
{    
    nx = std::stoi(catch_parameter("x_samples", file));
    ny = std::stoi(catch_parameter("y_samples", file));
    nz = std::stoi(catch_parameter("z_samples", file));

    nPoints = nx*ny*nz;

    dh = std::stof(catch_parameter("model_spacing", file));

    export_receiver_output = str2bool(catch_parameter("export_receiver_output", file));
    export_wavefield_output = str2bool(catch_parameter("export_wavefield_output", file));

    receiver_output_folder = catch_parameter("receiver_output_folder", file); 
    wavefield_output_folder = catch_parameter("wavefield_output_folder", file);

    Geometry * gtypes[] = {new Regular(), new Circular(), new Streamer()};

    int type = std::stoi(catch_parameter("geometry_type", file));

    geometry = gtypes[type];

    total_shots = geometry->shots.total;
    total_nodes = geometry->nodes.total;
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