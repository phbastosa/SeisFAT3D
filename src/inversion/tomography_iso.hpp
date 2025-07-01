# ifndef TOMOGRAPHY_ISO_HPP
# define TOMOGRAPHY_ISO_HPP

# include "inversion.hpp"

class Tomography_ISO : public Inversion
{
    void set_modeling_type();
    void set_objective_function();
    void set_sensitivity_matrix();
    void set_regularization();
};

# endif