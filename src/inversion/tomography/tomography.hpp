# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../inversion.hpp"

# include "../../modeling/eikonal_equation/podvin_and_lecomte/podvin_and_lecomte.cuh"
# include "../../modeling/eikonal_equation/fast_iterative_method/fast_iterative_method.cuh"
# include "../../modeling/eikonal_equation/fast_sweeping_method/fast_sweeping_method.cuh"

class Tomography : public Inversion
{
private:


protected:

    std::vector<Eikonal *> modeling;

public:

    void optimization();
    void model_update();
    void export_results();
    void set_parameters();
    void import_obs_data();
    void forward_modeling();
    void check_convergence();
};

# endif