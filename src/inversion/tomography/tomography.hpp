# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../inversion.hpp"

# include "../../modeling/eikonal_equation/podvin_and_lecomte/podvin_and_lecomte.cuh"
# include "../../modeling/eikonal_equation/fast_iterative_method/fast_iterative_method.cuh"
# include "../../modeling/eikonal_equation/fast_sweeping_method/fast_sweeping_method.cuh"

class Tomography : public Inversion
{
protected:

    Eikonal * modeling;

    void set_forward_modeling();
    void set_main_components();
    void tomography_message();
    void export_gradient();
    void init_modeling();
    
    void extract_calculated_data();

public:

    void import_obs_data();
    void check_convergence();

    void model_update();
    void export_results();

    virtual void optimization() = 0;
    virtual void set_parameters() = 0;
    virtual void forward_modeling() = 0;
};

# endif