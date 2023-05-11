# include "eikonal.cuh"

void Eikonal::set_parameters(std::string file)
{
    Modeling::set_parameters(file);

    eikonal_type = std::stoi(catch_parameter("eikonal_type", file));

    switch (eikonal_type)
    {
    case 0:
        PAL_parameters();    
        break;    
    case 1:
        FIM_parameters();
        break;
    case 2:
        FSM_parameters();
        break;
    default:
        eikonal_type = 1;
        FIM_parameters();
        break;
    }

    V = new float[nPoints]();

    import_binary_float(catch_parameter("vp_model_file", file), V, nPoints);
}

void Eikonal::set_components()
{
    volsize = nxx*nyy*nzz;

    T = new float[volsize]();
    S = new float[volsize]();
    
    eikonal_type == true ? pad_expansion() : fdm_expansion();

    switch (eikonal_type)
    {
    case 0:
        PAL_components();    
        break;    
    case 1:
        FIM_components();
        break;
    case 2:
        FSM_components();
        break;
    default:
        eikonal_type = 1;
        FIM_components();
        break;
    }
	
    wavefield_output_samples = nPoints;
    receiver_output_samples = geometry->nodes.total;

    receiver_output = new float[receiver_output_samples]();
    wavefield_output = new float[wavefield_output_samples]();
}

void Eikonal::initial_setup()
{
    switch (eikonal_type)
    {
    case 0:
        PAL_init();    
        break;    
    case 1:
        FIM_init();
        break;
    case 2:
        FSM_init();
        break;
    default:
        eikonal_type = 1;
        FIM_init();
        break;
    }
}

void Eikonal::forward_solver()
{
    switch (eikonal_type)
    {
    case 0:
        PAL_solver();    
        break;    
    case 1:
        FIM_solver();
        break;
    case 2:
        FSM_solver();
        break;
    default:
        eikonal_type = 1;
        FIM_solver();
        break;
    }
}

void Eikonal::build_outputs()
{
    get_travelTimes();
    get_firstArrivals();
}

void Eikonal::get_travelTimes()
{
    eikonal_type == true ? pad_reduction() : fdm_reduction();    

    wavefield_output_file = wavefield_output_folder + "travel_times_" + std::to_string(nz) + "x" + std::to_string(nx) + "x" + std::to_string(ny) + "_shot_" + std::to_string(shot_id+1) + ".bin";
}

void Eikonal::get_firstArrivals()
{
    for (int r = 0; r < total_nodes; r++)
    {
        float x = geometry->nodes.x[r];
        float y = geometry->nodes.y[r];
        float z = geometry->nodes.z[r];

        float x0 = floorf(x / dh) * dh;
        float y0 = floorf(y / dh) * dh;
        float z0 = floorf(z / dh) * dh;

        float x1 = floorf(x / dh) * dh + dh;
        float y1 = floorf(y / dh) * dh + dh;
        float z1 = floorf(z / dh) * dh + dh;

        int id = ((int)(z / dh)) + ((int)(x / dh))*nz + ((int)(y / dh))*nx*nz;

        float c000 = wavefield_output[id];
        float c001 = wavefield_output[id + 1];
        float c100 = wavefield_output[id + nz]; 
        float c101 = wavefield_output[id + 1 + nz]; 
        float c010 = wavefield_output[id + nx*nz]; 
        float c011 = wavefield_output[id + 1 + nx*nz]; 
        float c110 = wavefield_output[id + nz + nx*nz]; 
        float c111 = wavefield_output[id + 1 + nz + nx*nz];

        float xd = (x - x0) / (x1 - x0);
        float yd = (y - y0) / (y1 - y0);
        float zd = (z - z0) / (z1 - z0);

        float c00 = c000*(1 - xd) + c100*xd;    
        float c01 = c001*(1 - xd) + c101*xd;    
        float c10 = c010*(1 - xd) + c110*xd;    
        float c11 = c011*(1 - xd) + c111*xd;    

        float c0 = c00*(1 - yd) + c10*yd;
        float c1 = c01*(1 - yd) + c11*yd;

        receiver_output[r] = c0*(1 - zd) + c1*zd;
    }

    receiver_output_file = receiver_output_folder + "first_arrivals_" + std::to_string(geometry->nodes.total) + "_shot_" + std::to_string(shot_id+1) + ".bin";
}

void Eikonal::pad_expansion()
{
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = 1.0f / V[z + x*nz + y*nx*nz];
            }
        }
    }

    for (int z = 0; z < pdz; z++)
    {
        for (int y = 0; y < nyy - pdy; y++)
        {
            for (int x = 0; x < nxx - pdx; x++)
            {
                S[(nzz - z - 1) + x*nzz + y*nxx*nzz] = 1.0f / V[(nz - 1) + x*nz + y*nx*nz];
            }
        }
    }

    for (int x = 0; x < pdx; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = 0; y < nyy - pdy; y++)
            {
                S[z + (nxx - x - 1)*nzz + y*nxx*nzz] = S[z + (nxx - pdx - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    for (int y = 0; y < pdy; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                S[z + x*nzz + (nyy - y - 1)*nxx*nzz] = S[z + x*nzz + (nyy - pdy - 1)*nxx*nzz];
            }
        }
    }
}

void Eikonal::fdm_expansion()
{
    // Centering
    for (int z = padb; z < nzz - padb; z++)
    {
        for (int y = padb; y < nyy - padb; y++)
        {
            for (int x = padb; x < nxx - padb; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = 1.0f / V[(z - padb) + (x - padb)*nz + (y - padb)*nx*nz];
            }
        }
    }

    // Z direction
    for (int z = 0; z < padb; z++)
    {
        for (int y = padb; y < nyy - padb; y++)
        {
            for (int x = padb; x < nxx - padb; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = 1.0f / V[0 + (x - padb)*nz + (y - padb)*nx*nz];
                S[(nzz - z - 1) + x*nzz + y*nxx*nzz] = 1.0f / V[(nz - 1) + (x - padb)*nz + (y - padb)*nx*nz];
            }
        }
    }

    // X direction
    for (int x = 0; x < padb; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = padb; y < nyy - padb; y++)
            {
                S[z + x*nzz + y*nxx*nzz] = S[z + padb*nzz + y*nxx*nzz];
                S[z + (nxx - x - 1)*nzz + y*nxx*nzz] = S[z + (nxx - padb - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    // Y direction
    for (int y = 0; y < padb; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = S[z + x*nzz + padb*nxx*nzz];
                S[z + x*nzz + (nyy - y - 1)*nxx*nzz] = S[z + x*nzz + (nyy - padb - 1)*nxx*nzz];
            }
        }
    }
}

void Eikonal::pad_reduction()
{
    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         
        int x = (int) (index - y*nx*nz) / nz;    
        int z = (int) (index - x*nz - y*nx*nz);  

        wavefield_output[z + x*nz + y*nx*nz] = T[z + x*nzz + y*nxx*nzz];
    }
}

void Eikonal::fdm_reduction()
{
    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         
        int x = (int) (index - y*nx*nz) / nz;    
        int z = (int) (index - x*nz - y*nx*nz);  

        wavefield_output[z + x*nz + y*nx*nz] = T[(z + padb) + (x + padb)*nzz + (y + padb)*nxx*nzz];
    }
}

void Eikonal::free_space()
{
    delete[] T;
    delete[] S;

    switch (eikonal_type)
    {
    case 0:
        PAL_free_space();    
        break;    
    case 1:
        FIM_free_space();
        break;
    case 2:
        FSM_free_space();
        break;
    default:
        eikonal_type = 1;
        FIM_free_space();
        break;
    }
}


