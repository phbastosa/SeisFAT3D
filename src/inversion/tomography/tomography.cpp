# include "tomography.hpp"

void Tomography::set_parameters()
{
    modeling = 
    {
        new Podvin_and_Lecomte(),
        new Fast_Iterative_Method(),
        new Fast_Sweeping_Method()
    };
    
    if (!isInteger(catch_parameter("modeling_type", file)))
        throw std::invalid_argument("\033[31mError: Wrong modeling type! \033[0;0m");

    auto type = std::stoi(catch_parameter("modeling_type", file));

    if ((type < 0) || (type >= modeling.size()))
        throw std::invalid_argument("\033[31mError: Modeling type must to be an eikonal equation solver! \033[0;0m");

    modeling[type]->file = file;

    modeling[type]->set_parameters();

    set_general_parameters();    

    n_data = modeling[type]->total_shots * modeling[type]->total_nodes;
    n_model = modeling[type]->nx * modeling[type]->ny * modeling[type]->nz;

    std::cout<<n_data<<std::endl;
    std::cout<<n_model<<std::endl;
} 

void Tomography::import_obs_data()
{


}

void Tomography::forward_modeling()
{


}

void Tomography::check_convergence()
{


}

void Tomography::optimization()
{


}

void Tomography::model_update()
{


}

void Tomography::export_results()
{

    
}



