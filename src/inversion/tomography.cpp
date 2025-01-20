# include "tomography.hpp"

void Tomography::set_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", parameters));

    obs_data_folder = catch_parameter("obs_data_folder", parameters);
    obs_data_prefix = catch_parameter("obs_data_prefix", parameters);

    smooth_model_per_iteration = str2bool(catch_parameter("smooth_per_iteration", parameters));
    smoother_samples = std::stoi(catch_parameter("gaussian_filter_samples", parameters));
    smoother_stdv = std::stoi(catch_parameter("gaussian_filter_stdv", parameters));
    
    convergence_map_folder = catch_parameter("convergence_folder", parameters);
    estimated_model_folder = catch_parameter("inversion_output_folder", parameters);

    write_model_per_iteration = str2bool(catch_parameter("export_model_per_iteration", parameters));

    set_forward_modeling();
    set_inversion_elements();
    set_specifications();
} 

void Tomography::set_forward_modeling()
{
    modeling = new Eikonal_ISO();
    modeling->parameters = parameters;
    modeling->set_parameters();
}

void Tomography::set_inversion_elements()
{    
    n_data = 0;
    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
        n_data += modeling->geometry->spread[shot];

    dcal = new float[n_data]();
    dobs = new float[n_data]();

    perturbation = new float[modeling->nPoints]();
}

void Tomography::import_obsData()
{
    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        float * data = new float[modeling->geometry->spread[modeling->srcId]]();

        std::string path = obs_data_folder + obs_data_prefix + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin";

        import_binary_float(path, data, modeling->geometry->spread[modeling->srcId]);

        int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];    
        
        for (int i = 0; i < modeling->geometry->spread[modeling->srcId]; i++) 
            dobs[i + skipped] = data[i];

        delete[] data;
    }
}

void Tomography::forward_modeling()
{
    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        show_information();

        modeling->forward_solver();

        concatenate_data();
        
        if (iteration != max_iteration)
            apply_inversion_technique();
    }
}

void Tomography::show_information()
{
    modeling->show_information();    
    
    std::cout << "\nInversion type: " << inversion_method << "\n\n";

    if (iteration == max_iteration) 
        std::cout << "-------- Checking final residuo --------\n\n";
    else
    {    
        std::cout << "-------- Computing iteration " << iteration + 1 << " of " << max_iteration << " --------\n\n";

        if (iteration > 0) std::cout << "Previous residuo: " << residuo.back() << "\n\n";   
    }
}

void Tomography::concatenate_data()
{
    int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];

    for (int i = 0; i < modeling->geometry->spread[modeling->srcId]; i++) 
        dcal[i + skipped] = modeling->synthetic_data[i];    
}

void Tomography::check_convergence()
{
    float square_difference = 0.0f;

    for (int i = 0; i < n_data; i++)
        square_difference += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(square_difference));

    if ((iteration >= max_iteration))
    {
        std::cout << "Final residuo: "<< residuo.back() <<"\n";
        converged = true;
    }
    else
    {
        iteration += 1;
        converged = false;
    }
}

void Tomography::model_update()
{
    if (smooth_model_per_iteration)
    {
        int aux_nx = modeling->nx + 2*smoother_samples;
        int aux_ny = modeling->ny + 2*smoother_samples;
        int aux_nz = modeling->nz + 2*smoother_samples;

        int aux_nPoints = aux_nx*aux_ny*aux_nz;

        float * dm_aux = new float[aux_nPoints]();
        float * dm_smooth = new float[aux_nPoints]();

        # pragma omp parallel for
        for (int index = 0; index < modeling->nPoints; index++)
        {
            int k = (int) (index / (modeling->nx*modeling->nz));         
            int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   
            int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); 

            int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz + (k + smoother_samples)*aux_nx*aux_nz;

            dm_aux[ind_filt] = perturbation[i + j*modeling->nz + k*modeling->nx*modeling->nz];
        }

        smooth_volume(dm_aux, dm_smooth, aux_nx, aux_ny, aux_nz);

        # pragma omp parallel for    
        for (int index = 0; index < modeling->nPoints; index++)
        {
            int k = (int) (index / (modeling->nx*modeling->nz));         
            int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   
            int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); 

            int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz + (k + smoother_samples)*aux_nx*aux_nz;

            perturbation[i + j*modeling->nz + k*modeling->nx*modeling->nz] = dm_smooth[ind_filt];
        }
    
        delete[] dm_aux;
        delete[] dm_smooth;
    }   
    
    for (int index = 0; index < modeling->nPoints; index++)
    {
        int k = (int) (index / (modeling->nx*modeling->nz));         
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); 

        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz + (k + modeling->nb)*modeling->nxx*modeling->nzz;

        modeling->S[indb] += perturbation[index];
        modeling->Vp[index] = 1.0f / modeling->S[indb];
    }

    if (write_model_per_iteration)
    {
        std::string model_iteration_path = estimated_model_folder + inversion_name + "model_iteration_" + std::to_string(iteration) + "_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "x" + std::to_string(modeling->ny) + ".bin";

        export_binary_float(model_iteration_path, modeling->Vp, modeling->nPoints);
    }
}

void Tomography::smooth_volume(float * input, float * output, int nx, int ny, int nz)
{
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

        for (int y = -mid; y <= mid; y++)
        {
            for (int x = -mid; x <= mid; x++)
            {
                for (int z = -mid; z <= mid; z++)
                {          
                    int index = (z + mid) + (x + mid)*smoother_samples + (y + mid)*smoother_samples*smoother_samples; 
                    
                    float r = sqrtf(x*x + y*y + z*z);

                    kernel[index] = 1.0f / (pi*smoother_stdv) * expf(-((r*r)/(2.0f*smoother_stdv*smoother_stdv)));
        
                    sum += kernel[index]; 
                }
            }
        }

        for (int i = 0; i < nKernel; i++) 
            kernel[i] /= sum;
    }
        
    for (int k = mid; k < ny - mid; k++)
    {   
        for (int j = mid; j < nx - mid; j++)
        {
            for (int i = mid; i < nz - mid; i++)
            {       
                float accum = 0.0f;
                
                for (int yk = 0; yk < smoother_samples; yk++)
                {      
                    for (int xk = 0; xk < smoother_samples; xk++)
                    {      
                        for (int zk = 0; zk < smoother_samples; zk++)
                        {   
                            int index = zk + xk*smoother_samples + yk*smoother_samples*smoother_samples;   
                            int partial = (i - mid + zk) + (j - mid + xk)*nz + (k - mid + yk)*nx*nz; 

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

void Tomography::export_results()
{    
    std::string estimated_model_path = estimated_model_folder + inversion_name + "final_model_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "x" + std::to_string(modeling->ny) + ".bin";
    std::string convergence_map_path = convergence_map_folder + inversion_name + "convergence_" + std::to_string(iteration) + "_iterations.txt"; 

    export_binary_float(estimated_model_path, modeling->Vp, modeling->nPoints);

    std::ofstream resFile(convergence_map_path, std::ios::out);
    
    for (int r = 0; r < residuo.size(); r++) 
        resFile << residuo[r] << "\n";

    resFile.close();

    std::cout << "Text file \033[34m" << convergence_map_path << "\033[0;0m was successfully written." << std::endl;
}
