# ifndef TOMOGRAPHY_ISO_HPP
# define TOMOGRAPHY_ISO_HPP

# include "inversion.hpp"

class Tomography_ISO : public Inversion
{
private:

    void set_modeling_type();
    void set_sensitivity_matrix();
    void get_parameter_variation();
    void export_estimated_models();    

public:

    void model_update();
};

# endif