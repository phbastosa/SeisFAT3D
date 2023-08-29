# include "inversion.hpp"

void Inversion::set_general_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", file));

    smooth = str2bool(catch_parameter("smooth_per_iteration", file));
    smoother_samples = std::stoi(catch_parameter("gaussian_filter_samples", file));
    smoother_stdv = std::stoi(catch_parameter("gaussian_filter_stdv", file));

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
    int init = smoother_samples / 2;
    int nPoints = nx * ny * nz;
    int nKernel = smoother_samples * smoother_samples * smoother_samples;

    float pi = 4.0f * atanf(1.0f); 

    float * kernel = new float[nKernel]();

    # pragma omp parallel for
    for (int i = 0; i < nPoints; i++) 
        output[i] = input[i];

    int mid = (int)(smoother_samples / 2); 

    kernel[mid + mid*smoother_samples + mid*smoother_samples*smoother_samples] = 1.0f;

    if (smoother_stdv != 0.0f)
    {
        float sum = 0.0f;

        for (int y = -init; y <= init; y++)
        {
            for (int x = -init; x <= init; x++)
            {
                for (int z = -init; z <= init; z++)
                {          
                    int index = (z+init) + (x+init)*smoother_samples + (y+init)*smoother_samples*smoother_samples; 
                    
                    float r = sqrtf(x*x + y*y + z*z);

                    kernel[index] = 1.0f / (pi*smoother_stdv) * expf(-((r*r)/(2.0f*smoother_stdv*smoother_stdv)));
        
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
                
                for (int yk = 0; yk < smoother_samples; yk++)
                {      
                    for (int xk = 0; xk < smoother_samples; xk++)
                    {      
                        for (int zk = 0; zk < smoother_samples; zk++)
                        {   
                            int index = zk + xk*smoother_samples + yk*smoother_samples*smoother_samples;   
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