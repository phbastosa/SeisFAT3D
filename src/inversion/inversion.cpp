# include "inversion.hpp"

void Inversion::set_general_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", file));

    smooth = str2bool(catch_parameter("smooth_per_iteration", file));
    samples = std::stoi(catch_parameter("gaussian_filter_samples", file));
    stdv = std::stoi(catch_parameter("gaussian_filter_stdv", file));

    obs_data_folder = catch_parameter("obs_data_folder", file);
    obs_data_prefix = catch_parameter("obs_data_prefix", file);

    convergence_map_folder = catch_parameter("convergence_map_folder", file);
    estimated_model_folder = catch_parameter("estimated_model_folder", file);

    gradient_folder = catch_parameter("gradient_folder", file);

    write_gradient_per_iteration = str2bool(catch_parameter("export_gradient", file));
    write_model_per_iteration = str2bool(catch_parameter("export_model_per_iteration", file));
}

void Inversion::smoothing(float * input, float * output, int nx, int ny, int nz)
{
    int init = samples / 2;
    int nPoints = nx * ny * nz;
    int nKernel = samples * samples * samples;

    float pi = 4.0f * atanf(1.0f); 

    float * kernel = new float[nKernel]();

    for (int i = 0; i < nPoints; i++) 
        output[i] = input[i];

    int mid = (int)(samples / 2); 

    kernel[mid + mid*samples + mid*samples*samples] = 1.0f;

    if (stdv != 0.0f)
    {
        float sum = 0.0f;

        for (int y = -init; y <= init; y++)
        {
            for (int x = -init; x <= init; x++)
            {
                for (int z = -init; z <= init; z++)
                {          
                    int index = (z+init) + (x+init)*samples + (y+init)*samples*samples; 

                    kernel[index] = 1.0f / (2.0f*pi*stdv*stdv) * expf(-((x*x + y*y + z*z) / (2.0f*stdv*stdv)));
        
                    sum += kernel[index]; 
                }
            }
        }

        for (int i = 0; i < nKernel; i++) 
            kernel[i] /= sum;
    }
        
    for (int k = init; k < ny - init; k++)
    {   
        for (int j = init; j < nx - init; j++)
        {
            for (int i = init; i < nz - init; i++)
            {       
                float accum = 0.0f;
                
                for (int yk = 0; yk < samples; yk++)
                {      
                    for (int xk = 0; xk < samples; xk++)
                    {      
                        for (int zk = 0; zk < samples; zk++)
                        {   
                            int index = zk + xk*samples + yk*samples*samples;   
                            int partial = (i-init+zk) + (j-init+xk)*nz + (k-init+yk)*nx*nz; 

                            accum += input[partial] * kernel[index];
                        }        
                    }
                }
                
                output[i + j*nz + k*nx*nz] = accum;
            }
        }   
    }

    delete[] kernel;
}