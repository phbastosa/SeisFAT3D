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
    n_model = modeling->nPoints;

    dcal = new float[n_data]();
    dobs = new float[n_data]();

    dm = new float[n_model]();
    model = new float[n_model]();
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

void Tomography::model_update()
{


}

void Tomography::export_results()
{

    
}



