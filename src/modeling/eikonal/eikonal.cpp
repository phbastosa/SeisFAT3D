# include "eikonal.hpp"

void Eikonal::set_parameters(std::string file)
{
    Modeling::set_parameters(file);

    parameters();

    V = new float[nPoints]();

    import_binary_float(catch_parameter("vp_model_file", file), V, nPoints);
}

void Eikonal::set_components()
{
    volsize = nxx*nyy*nzz;

    T = new float[volsize]();
    S = new float[volsize]();

    expansion();
    components();

    wavefield_output_samples = nPoints;
    receiver_output_samples = geometry->nodes.total;

    receiver_output = new float[receiver_output_samples]();
    wavefield_output = new float[wavefield_output_samples]();
}

void Eikonal::build_outputs()
{
    get_travelTimes();
    get_firstArrivals();
}

void Eikonal::get_travelTimes()
{
    reduction();

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