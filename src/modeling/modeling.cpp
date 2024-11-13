# include "modeling.hpp"

void Modeling::set_parameters()
{
    nx = std::stoi(catch_parameter("x_samples", parameters));
    ny = std::stoi(catch_parameter("y_samples", parameters));
    nz = std::stoi(catch_parameter("z_samples", parameters));

    dx = std::stof(catch_parameter("x_spacing", parameters));
    dy = std::stof(catch_parameter("y_spacing", parameters));
    dz = std::stof(catch_parameter("z_spacing", parameters));

    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));

    data_folder = catch_parameter("modeling_output_folder", parameters);

    nPoints = nx*ny*nz;

    geometry = new Geometry();

    geometry->parameters = parameters;
    geometry->set_parameters();

    max_spread = 0;
    for (int index = 0; index < geometry->nrel; index++)
    {   
        if (max_spread < geometry->spread[index])
            max_spread = geometry->spread[index]; 
    }

    set_boundaries();

    nxx = nx + 2*nb;
    nyy = ny + 2*nb;
    nzz = nz + 2*nb;

    volsize = nxx*nyy*nzz;

    set_specifications();
}

void Modeling::expand_boundary(float * input, float * output)
{
    for (int index = 0; index < nPoints; index++)
    {
        int k = (int) (index / (nx*nz));         
        int j = (int) (index - k*nx*nz) / nz;    
        int i = (int) (index - j*nz - k*nx*nz);  

        output[(i + nb) + (j + nb)*nzz + (k + nb)*nxx*nzz] = input[i + j*nz + k*nx*nz];       
    }

    for (int k = nb; k < nyy - nb; k++)
    {   
        for (int j = nb; j < nxx - nb; j++)
        {
            for (int i = 0; i < nb; i++)            
            {
                output[i + j*nzz + k*nxx*nzz] = input[0 + (j - nb)*nz + (k - nb)*nx*nz];
                output[(nzz - i - 1) + j*nzz + k*nxx*nzz] = input[(nz - 1) + (j - nb)*nz + (k - nb)*nx*nz];
            }
        }
    }

    for (int k = 0; k < nyy; k++)
    {   
        for (int j = 0; j < nb; j++)
        {
            for (int i = 0; i < nzz; i++)
            {
                output[i + j*nzz + k*nxx*nzz] = output[i + nb*nzz + k*nxx*nzz];
                output[i + (nxx - j - 1)*nzz + k*nxx*nzz] = output[i + (nxx - nb - 1)*nzz + k*nxx*nzz];
            }
        }
    }

    for (int k = 0; k < nb; k++)
    {   
        for (int j = 0; j < nxx; j++)
        {
            for (int i = 0; i < nzz; i++)
            {
                output[i + j*nzz + k*nxx*nzz] = output[i + j*nzz + nb*nxx*nzz];
                output[i + j*nzz + (nyy - k - 1)*nxx*nzz] = output[i + j*nzz + (nyy - nb - 1)*nxx*nzz];
            }
        }
    }
}

void Modeling::reduce_boundary(float * input, float * output)
{
    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
    {
        int k = (int) (index / (nx*nz));         
        int j = (int) (index - k*nx*nz) / nz;    
        int i = (int) (index - j*nz - k*nx*nz);  

        output[i + j*nz + k*nx*nz] = input[(i + nb) + (j + nb)*nzz + (k + nb)*nxx*nzz];
    }
}

void Modeling::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------------------\n";
    std::cout << "                                 \033[34mSeisFAT3D\033[0;0m\n";
    std::cout << "-------------------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (nz - 1)*dz << ", x = " << (nx - 1) * dx <<", y = " << (ny - 1) * dy << ") m\n\n";

    std::cout << "Running shot " << srcId + 1 << " of " << geometry->nrel << " in total\n\n";

    std::cout << "Current shot position: (z = " << geometry->zsrc[geometry->sInd[srcId]] << 
                                       ", x = " << geometry->xsrc[geometry->sInd[srcId]] << 
                                       ", y = " << geometry->ysrc[geometry->sInd[srcId]] << ") m\n\n";

    std::cout << modeling_name << "\n";
}
