# include "modeling.hpp"

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

void Modeling::expand_boundary(float * input, float * output)
{
    for (int y = nbyl; y < nyy - nbyr; y++)
    {
        for (int x = nbxl; x < nxx - nbyr; x++)
        {
            for (int z = nbzu; z < nzz - nbzd; z++)
            {
                output[z + x*nzz + y*nxx*nzz] = input[(z - nbzu) + (x - nbxl)*nz + (y - nbyl)*nx*nz];       
            }
        }
    }

    for (int y = nbyl; y < nyy - nbyr; y++)
    {
        for (int x = nbxl; x < nxx - nbyr; x++)
        {
            for (int z = 0; z < nbzu; z++)
            {
                output[z + x*nzz + y*nxx*nzz] = input[0 + (x - nbxl)*nz + (y - nbyl)*nx*nz];
            }
        }
    }

    for (int y = nbyl; y < nyy - nbyr; y++)
    {
        for (int x = nbxl; x < nxx - nbyr; x++)
        {
            for (int z = nzz - nbzd; z < nzz; z++)
            {
                output[z + x*nzz + y*nxx*nzz] = input[(nz - 1) + (x - nbxl)*nz + (y - nbyl)*nx*nz];
            }
        }
    }

    for (int y = nbyl; y < nyy - nbyr; y++)
    {
        for (int x = nbxl; x < nxx - nbyr; x++)
        {
            for (int z = nzz - nbzd; z < nzz; z++)
            {
                output[z + x*nzz + y*nxx*nzz] = input[(nz - 1) + (x - nbxl)*nz + (y - nbyl)*nx*nz];
            }
        }
    }

    for (int y = nbyl; y < nyy - nbyr; y++)
    {
        for (int x = 0; x < nbxl; x++)
        {
            for (int z = 0; z < nzz; z++)
            {
                output[z + x*nzz + y*nxx*nzz] = output[z + nbxl*nzz + y*nxx*nzz];
            }
        }
    }

    for (int y = nbyl; y < nyy - nbyr; y++)
    {
        for (int x = nxx-nbxr; x < nxx; x++)
        {
            for (int z = 0; z < nzz; z++)
            {
                output[z + x*nzz + y*nxx*nzz] = output[z + (nxx - nbxr - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    for (int y = 0; y < nbyl; y++)
    {
        for (int x = 0; x < nxx; x++)
        {
            for (int z = 0; z < nzz; z++)
            {
                output[z + x*nzz + y*nxx*nzz] = output[z + x*nzz + nbyl*nxx*nzz];
            }
        }
    }

    for (int y = nyy - nbyr; y < nyy; y++)
    {
        for (int x = 0; x < nxx; x++)
        {
            for (int z = 0; z < nzz; z++)
            {
                output[z + x*nzz + y*nxx*nzz] = output[z + x*nzz + (nyy - nbyr - 1)*nxx*nzz];
            }
        }
    }
}

void Modeling::reduce_boundary(float * input, float * output)
{
    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         
        int x = (int) (index - y*nx*nz) / nz;    
        int z = (int) (index - x*nz - y*nx*nz);  

        output[z + x*nz + y*nx*nz] = input[(z + nbzu) + (x + nbxl)*nzz + (y + nbyl)*nxx*nzz];
    }
}

void Modeling::check_geometry_overflow()
{
    for (int shot = 0; shot < total_shots; shot++)
    {
        if ((geometry->shots.x[shot] < 0) || (geometry->shots.x[shot] > (nx-1)*dx) || 
            (geometry->shots.y[shot] < 0) || (geometry->shots.y[shot] > (ny-1)*dy) ||
            (geometry->shots.z[shot] < 0) || (geometry->shots.z[shot] > (nz-1)*dz))       
        throw std::invalid_argument("\033[31mError: shots geometry overflow!\033[0;0m");
    }

    for (int node = 0; node < total_nodes; node++)
    {
        if ((geometry->nodes.x[node] < 0) || (geometry->nodes.x[node] > (nx-1)*dx) || 
            (geometry->nodes.y[node] < 0) || (geometry->nodes.y[node] > (ny-1)*dy) ||
            (geometry->nodes.z[node] < 0) || (geometry->nodes.z[node] > (nz-1)*dz))       
        throw std::invalid_argument("\033[31mError: nodes geometry overflow!\033[0;0m");
    }
}

void Modeling::set_parameters()
{
    set_generals();

    set_geometry();

    set_specifics();
    
    set_boundary();     
    
    set_vp_model();

    set_volumes();

    set_outputs();
}

void Modeling::set_generals()
{
    threadsPerBlock = 256;

    nx = std::stoi(catch_parameter("x_samples", file));
    ny = std::stoi(catch_parameter("y_samples", file));
    nz = std::stoi(catch_parameter("z_samples", file));
    
    dx = std::stof(catch_parameter("x_spacing", file));
    dy = std::stof(catch_parameter("y_spacing", file));
    dz = std::stof(catch_parameter("z_spacing", file));

    nPoints = nx*ny*nz;

    export_receiver_output = str2bool(catch_parameter("export_receiver_output", file));
    receiver_output_folder = catch_parameter("receiver_output_folder", file); 
}

void Modeling::set_geometry()
{
    geometry = new Geometry(); 

    geometry->file = file;

    geometry->set_geometry();

    total_shots = geometry->shots.total;
    total_nodes = geometry->nodes.total;

    check_geometry_overflow();
}

void Modeling::set_boundary()
{
    nxx = nx + nbxl + nbxr;
    nyy = ny + nbyl + nbyr;
    nzz = nz + nbzu + nbzd;

    volsize = nxx*nyy*nzz;

    blocksPerGrid = (int)(volsize / threadsPerBlock);    
}

void Modeling::set_vp_model()
{
    std::string vp_file = catch_parameter("vp_model_file", file);

    model = new float[nPoints](); 

    import_binary_float(vp_file, model, nPoints);

    V = new float[volsize]();

    expand_boundary(model, V);
}

void Modeling::print_information()
{
    auto clear = system("clear");

    std::cout << "\033[1mSeisFAT3D\033[m - Modeling program\n\n";

    std::cout << "Model dimensions: (z = " << (nz-1)*dz << 
                                  ", x = " << (nx-1)*dx << 
                                  ", y = " << (ny-1)*dy << ") m\n\n";

    std::cout << "Modeling type: \033[1m" << type_message << "\033[m\n\n";

    std::cout << "Running shot " << shot_index+1 << " of " << total_shots;

    std::cout << " at position (z = " << geometry->shots.z[shot_index] << 
                             ", x = " << geometry->shots.x[shot_index] << 
                             ", y = " << geometry->shots.y[shot_index] << ") m\n\n";
}

void Modeling::set_initial_conditions()
{
    sidx = (int)(geometry->shots.x[shot_index] / dx) + nbxl;
    sidy = (int)(geometry->shots.y[shot_index] / dy) + nbyl;
    sidz = (int)(geometry->shots.z[shot_index] / dz) + nbzu;

    source_index = sidz + sidx*nzz + sidy*nxx*nzz;

    initialization();
}

void Modeling::export_outputs()
{
    if (export_receiver_output) 
        export_binary_float(receiver_output_file, receiver_output, receiver_output_samples);
}
