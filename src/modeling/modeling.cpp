# include "modeling.hpp"

void Modeling::set_parameters(std::string file)
{    
    get_vRAM_init();

    nx = std::stoi(catch_parameter("x_samples", file));
    ny = std::stoi(catch_parameter("y_samples", file));
    nz = std::stoi(catch_parameter("z_samples", file));

    nPoints = nx*ny*nz;

    dh = std::stof(catch_parameter("model_spacing", file));

    export_receiver_output = str2bool(catch_parameter("export_receiver_output", file));
    export_wavefield_output = str2bool(catch_parameter("export_wavefield_output", file));

    receiver_output_folder = catch_parameter("receiver_output_folder", file); 
    wavefield_output_folder = catch_parameter("wavefield_output_folder", file);

    Geometry * gtypes[] = {new Regular(), new Circular()};

    int type = std::stoi(catch_parameter("geometry_type", file));

    geometry = gtypes[type];

    geometry->set_geometry(file);

    // Check if geometry is inside the model dimensions

    total_shots = geometry->shots.total;
    total_nodes = geometry->nodes.total;
}

void Modeling::info_message()
{
    auto clear = system("clear");

    get_RAM_usage();
    get_GPU_usage();
    
    std::cout<<title<<"\n";
    
    std::cout<<"Total x model length = "<<(nx-1)*dh<<" m\n";
    std::cout<<"Total Y model length = "<<(ny-1)*dh<<" m\n";
    std::cout<<"Total Z model length = "<<(nz-1)*dh<<" m\n\n";

    std::cout<<"Shot " << shot_id+1 << " of " << geometry->shots.total;

    std::cout << " at position (x,y,z) = (" << geometry->shots.x[shot_id] << ", " 
                                            << geometry->shots.y[shot_id] << ", " 
                                            << geometry->shots.z[shot_id] << ") m\n\n";

    std::cout<<"Memory usage: \n";
    std::cout<<"RAM = "<<RAM<<" Mb\n";
    std::cout<<"GPU = "<<vRAM<<" Mb\n\n";
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

void Modeling::get_RAM_usage()
{
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    RAM = (int) (usage.ru_maxrss / 1024);
}

void Modeling::get_vRAM_init()
{
	size_t freeMem, totalMem;
	cudaMemGetInfo(&freeMem, &totalMem);
    ivRAM = (int) ((totalMem - freeMem) / (1024 * 1024));
}

void Modeling::get_GPU_usage()
{
	size_t freeMem, totalMem;
	cudaMemGetInfo(&freeMem, &totalMem);
    vRAM = (int) ((totalMem - freeMem) / (1024 * 1024));
    vRAM -= ivRAM;
}

void Modeling::export_outputs()
{
    if (export_receiver_output) export_binary_float(receiver_output_file, receiver_output, receiver_output_samples);
    if (export_wavefield_output) export_binary_float(wavefield_output_file, wavefield_output, wavefield_output_samples);
}