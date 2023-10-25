# include "tomography.hpp"

void Tomography::set_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", file));

    obs_data_folder = catch_parameter("obs_data_folder", file);
    obs_data_prefix = catch_parameter("obs_data_prefix", file);

    smooth = str2bool(catch_parameter("smooth_per_iteration", file));
    smoother_samples = std::stoi(catch_parameter("gaussian_filter_samples", file));
    smoother_stdv = std::stoi(catch_parameter("gaussian_filter_stdv", file));
    
    convergence_map_folder = catch_parameter("convergence_folder", file);
    estimated_model_folder = catch_parameter("estimated_model_folder", file);

    gradient_folder = catch_parameter("gradient_folder", file);

    write_model_per_iteration = str2bool(catch_parameter("export_model_per_iteration", file));
    write_gradient_per_iteration = str2bool(catch_parameter("export_gradient_per_iteration", file));

    set_forward_modeling();

    set_inversion_volumes();

    set_specific_parameters();
} 

void Tomography::set_forward_modeling()
{
    std::vector<Eikonal *> possibilities = 
    {
        new Podvin_and_Lecomte(),
        new Block_FIM(),
        new Accurate_FSM()
    };
    
    auto type = std::stoi(catch_parameter("modeling_type", file));

    modeling = possibilities[type];

    modeling->file = file;

    modeling->set_parameters();

    modeling->set_runtime();
}

void Tomography::set_inversion_volumes()
{
    n_data = modeling->total_shots * modeling->total_nodes;
    
    dcal = new float[n_data]();
    dobs = new float[n_data]();

    dm = new float[modeling->nPoints]();
    model = new float[modeling->nPoints]();

    gradient = new float[modeling->nPoints]();

    modeling->reduce_boundary(modeling->S, model);
}

void Tomography::import_obs_data()
{
    int ptr = 0; 
    
    float * data = new float[modeling->total_nodes]();

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        import_binary_float(obs_data_folder + obs_data_prefix + std::to_string(shot+1) + ".bin", data, modeling->total_nodes);

        for (int d = ptr; d < ptr + modeling->total_nodes; d++) 
            dobs[d] = data[d - ptr];

        ptr += modeling->total_nodes;        
    }

    delete[] data;    
}

void Tomography::forward_modeling()
{
    modeling->expand_boundary(model, modeling->S);

    for (int index = 0; index < modeling->nPoints; index++)    
        gradient[index] = 0.0f;

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_id = shot;
    
        modeling->info_message();

        set_tomography_message();

        modeling->initial_setup();
        modeling->forward_solver();
        modeling->build_outputs();

        extract_calculated_data();
        
        if (iteration != max_iteration)
            apply_inversion_technique();
    }

    gradient_preconditioning();

    export_gradient();
}

void Tomography::set_tomography_message()
{
    std::cout<<"Tomography:"<<"\n";
    std::cout<<inversion_method<<"\n\n";

    if (iteration == max_iteration)
    { 
        std::cout<<"------- Checking final residuo ------------\n\n";
    }
    else
        std::cout<<"------- Computing iteration "<<iteration+1<<" of "<<max_iteration<<" ------------\n\n";

    if (iteration > 0) std::cout<<"Previous residuo: "<<residuo.back()<<"\n\n";    
}

void Tomography::extract_calculated_data()
{
    int skipped = modeling->shot_id * modeling->total_nodes;

    for (int i = 0; i < modeling->total_nodes; i++) 
        dcal[i + skipped] = modeling->receiver_output[i];
}

void Tomography::check_convergence()
{
    float square_difference = 0.0f;

    for (int i = 0; i < n_data; i++)
        square_difference += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(square_difference));

    if ((iteration >= max_iteration))
    {
        std::cout << "\nFinal residuo: "<< residuo.back() <<"\n\n";
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
    float * velocity = new float[modeling->nPoints]();

    if (smooth)
    {
        int aux_nx = modeling->nx + 2*smoother_samples;
        int aux_ny = modeling->ny + 2*smoother_samples;
        int aux_nz = modeling->nz + 2*smoother_samples;

        int aux_nPoints = aux_nx*aux_ny*aux_nz;

        float * dm_aux = new float[aux_nPoints]();
        float * dm_smooth = new float[aux_nPoints]();

        for (int index = 0; index < modeling->nPoints; index++)
        {
            int k = (int) (index / (modeling->nx*modeling->nz));        
            int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
            int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);          

            int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz + (k + smoother_samples)*aux_nx*aux_nz;

            dm_aux[ind_filt] = dm[i + j*modeling->nz + k*modeling->nx*modeling->nz];
        }

        smooth_volume(dm_aux, dm_smooth, aux_nx, aux_ny, aux_nz);

        for (int index = 0; index < modeling->nPoints; index++)
        {
            int k = (int) (index / (modeling->nx*modeling->nz));        
            int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
            int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);          

            int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz + (k + smoother_samples)*aux_nx*aux_nz;

            dm[i + j*modeling->nz + k*modeling->nx*modeling->nz] = dm_smooth[ind_filt];
        }
    
        delete[] dm_aux;
        delete[] dm_smooth;
    }

    for (int index = 0; index < modeling->nPoints; index++)
    {
        model[index] += dm[index];

        velocity[index] = 1.0f / model[index]; 
    }

    if (write_model_per_iteration)
    {
        std::string model_iteration_path = estimated_model_folder + "model_iteration_" + std::to_string(iteration) + "_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "x" + std::to_string(modeling->ny) + ".bin";

        export_binary_float(model_iteration_path, velocity, modeling->nPoints);
    }

    delete[] velocity;
}

void Tomography::export_results()
{
    float * final_model = new float[modeling->nPoints]();

    for (int index = 0; index < modeling->nPoints; index++)
    {
        final_model[index] = 1.0f / model[index];
    }
    
    std::string estimated_model_path = estimated_model_folder + "final_model_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "x" + std::to_string(modeling->ny) + ".bin";
    std::string convergence_map_path = convergence_map_folder + "convergence_" + std::to_string(iteration) + "_iterations.txt"; 

    export_binary_float(estimated_model_path, final_model, modeling->nPoints);

    std::ofstream resFile(convergence_map_path, std::ios::out);
    
    for (int r = 0; r < residuo.size(); r++) 
        resFile << residuo[r] << "\n";

    resFile.close();

    std::cout<<"Text file "<<convergence_map_path<<" was successfully written."<<std::endl;

    modeling->get_runtime();

    delete[] final_model;
}

void Tomography::export_gradient()
{
    if (write_gradient_per_iteration)
    {
        std::string gradient_path = gradient_folder + "gradient_iteration_" + std::to_string(iteration+1) + "_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "x" + std::to_string(modeling->ny) + ".bin";

        export_binary_float(gradient_path, gradient, modeling->nPoints);
    }
}

void Tomography::smooth_volume(float * input, float * output, int nx, int ny, int nz)
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









