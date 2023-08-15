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
}

void Tomography::set_main_components()
{
    n_data = modeling->total_shots * modeling->total_nodes;
    
    dcal = new float[n_data]();
    dobs = new float[n_data]();

    dm = new float[modeling->nPoints]();
    model = new float[modeling->nPoints]();

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

void Tomography::check_convergence()
{
    float r = 0.0f;

    for (int i = 0; i < n_data; i++)
        r += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(r));

    if ((iteration >= max_iteration) || (residuo.back() <= tolerance))
    {
        std::cout << "\nFinal residuo: "<< residuo.back() << std::endl;
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
    modeling->info_message();

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
    for (int index = 0; index < modeling->nPoints; index++)
    {
        model[index] += dm[index];

        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); 

        int indB = (i + modeling->nbzu) + (j + modeling->nbxl)*modeling->nzz + (k + modeling->nbyl)*modeling->nxx*modeling->nzz;

        modeling->S[indB] = model[index];
    }
}

void Tomography::export_results()
{
    float * final_model = new float[modeling->nPoints]();

    for (int index = 0; index < modeling->nPoints; index++)
    {
        final_model[index] = 1.0f / model[index];
    }
    
    std::string estimated_model_path = estimated_model_folder + "final_model_iteration_" + std::to_string(iteration) + ".bin";
    std::string convergence_map_path = convergence_map_folder + "convergency.txt"; 

    export_binary_float(estimated_model_path, final_model, modeling->nPoints);

    std::ofstream resFile(convergence_map_path, std::ios::out);
    
    for (int r = 0; r < residuo.size(); r++) 
        resFile << residuo[r] << "\n";

    resFile.close();

    delete[] final_model;
}



