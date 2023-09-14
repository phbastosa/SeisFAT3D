# include "tomography.hpp"

void Tomography::set_forward_modeling()
{
    std::vector<Eikonal *> possibilities = 
    {
        new Podvin_and_Lecomte(),
        new Fast_Iterative_Method(),
        new Fast_Sweeping_Method()
    };
    
    auto type = std::stoi(catch_parameter("modeling_type", file));

    modeling = possibilities[type];

    modeling->file = file;

    modeling->set_parameters();

    modeling->set_runtime();
}

void Tomography::set_main_components()
{
    n_data = modeling->total_shots * modeling->total_nodes;
    
    dcal = new float[n_data]();
    dobs = new float[n_data]();

    dm = new float[modeling->nPoints]();
    model = new float[modeling->nPoints]();

    gradient = new float[modeling->nPoints]();

    for (int index = 0; index < modeling->nPoints; index++)
    {
        int k = (int) (index / (modeling->nx*modeling->nz));
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;  
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);

        int indB = (i + modeling->nbzu) + (j + modeling->nbxl)*modeling->nzz + (k + modeling->nbyl)*modeling->nxx*modeling->nzz;

        model[index] = modeling->S[indB];
    }
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

void Tomography::init_modeling()
{
    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
    {    
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);          

        int indB = (i+modeling->nbzu) + (j+modeling->nbxl)*modeling->nzz + (k+modeling->nbyl)*modeling->nxx*modeling->nzz;

        modeling->S[indB] = model[index];

        gradient[index] = 0.0f;
    }
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

void Tomography::tomography_message()
{
    std::cout<<"Inversion:"<<"\n";
    std::cout<<inversion_method<<"\n\n";

    if (iteration == max_iteration)
    { 
        std::cout<<"------- Checking final residuo ------------\n\n";
    }
    else
        std::cout<<"------- Computing iteration "<<iteration+1<<" of "<<max_iteration<<" ------------\n\n";

    if (iteration > 0) std::cout<<"Previous iteration residuo: "<<residuo.back()<<"\n\n";    
}

void Tomography::extract_calculated_data()
{
    int skipped = modeling->shot_id * modeling->total_nodes;

    for (int i = 0; i < modeling->total_nodes; i++) 
        dcal[i + skipped] = modeling->receiver_output[i];
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
        std::string gradient_path = gradient_folder + "gradient_iteration_" + std::to_string(iteration) + "_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "x" + std::to_string(modeling->ny) + ".bin";

        export_binary_float(gradient_path, gradient, modeling->nPoints);
    }
}

