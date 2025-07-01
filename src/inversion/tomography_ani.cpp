# include "tomography_ani.hpp"

void Tomography_ANI::set_modeling_type()
{
    modeling = new Eikonal_ANI();
    modeling->parameters = parameters;
    modeling->set_parameters();

    inversion_name = "tomography_ani";
    inversion_method = "Anisotropic First-Arrival Tomography";
}

void Tomography_ANI::set_objective_function()
{


}

void Tomography_ANI::set_sensitivity_matrix()
{


    
}

void Tomography_ANI::set_regularization()
{


    
}