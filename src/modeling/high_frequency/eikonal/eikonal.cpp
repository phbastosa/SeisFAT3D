# include "eikonal.hpp"

void Eikonal::set_specifications()
{
    V = new float[nPoints]();

    // import_binary_float(catch_parameter("vp_model_file", file), V, nPoints);

    for (int index = 0; index < nPoints; index++) V[index] = 1500.0f; 

    set_model_boundaries();

    nxx = nx + nbxl + nbxr;
    nyy = ny + nbyl + nbyr;
    nzz = nz + nbzu + nbzd;

    volsize = nxx*nyy*nzz;

    S = new float[volsize]();
    T = new float[volsize]();

    expand_boundary(V, S);

    for (int index = 0; index < volsize; index++) S[index] = 1.0f / S[index];

    set_preconditioners();

    wavefield_output_samples = nPoints;
    receiver_output_samples = total_nodes;

    receiver_output = new float[receiver_output_samples]();
    wavefield_output = new float[wavefield_output_samples]();
}

void Eikonal::build_outputs()
{
    get_travel_times();
    get_first_arrivals();
}

void Eikonal::get_travel_times()
{
    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         
        int x = (int) (index - y*nx*nz) / nz;    
        int z = (int) (index - x*nz - y*nx*nz);  

        wavefield_output[z + x*nz + y*nx*nz] = T[(z + nbzu) + (x + nbxl)*nzz + (y + nbyl)*nxx*nzz];
    }

    wavefield_output_file = wavefield_output_folder + modeling_method +"_time_volume_" + std::to_string(nz) + "x" + std::to_string(nx) + "x" + std::to_string(ny) + "_shot_" + std::to_string(shot_id+1) + ".bin";
}

void Eikonal::get_first_arrivals()
{
    for (int r = 0; r < total_nodes; r++)
    {
        float x = geometry->nodes.x[r];
        float y = geometry->nodes.y[r];
        float z = geometry->nodes.z[r];

        float x0 = floorf(x / dx) * dx;
        float y0 = floorf(y / dy) * dy;
        float z0 = floorf(z / dz) * dz;

        float x1 = floorf(x / dx) * dx + dx;
        float y1 = floorf(y / dy) * dy + dy;
        float z1 = floorf(z / dz) * dz + dz;

        int id = ((int)(z / dz)) + ((int)(x / dx))*nz + ((int)(y / dy))*nx*nz;

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

    receiver_output_file = receiver_output_folder + modeling_method + "_data_" + std::to_string(geometry->nodes.total) + "_shot_" + std::to_string(shot_id+1) + ".bin";
}


























