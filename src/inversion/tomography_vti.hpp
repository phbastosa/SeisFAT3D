# ifndef TOMOGRAPHY_VTI_HPP
# define TOMOGRAPHY_VTI_HPP

# include "inversion.hpp"

class Tomography_VTI : public Inversion
{
private:

    float * E = nullptr;
    float * D = nullptr;

    float * dE = nullptr;
    float * dD = nullptr;

    Eikonal_ANI * eikonal = nullptr;

    void set_modeling_type();
    void set_sensitivity_matrix();
    void get_parameter_variation();
    void export_estimated_models();    

public:

    void model_update();
};

# endif