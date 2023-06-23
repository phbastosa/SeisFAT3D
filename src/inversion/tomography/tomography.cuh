# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../inversion.hpp"
# include "../../modeling/eikonal/eikonal.cuh"

class Tomography : public Inversion
{
private:

    float * source = nullptr;       
    float * adjoint = nullptr;

    void adjoint_state_solver();
    void fill_calculated_data();

public:

    void set_parameters();
    void import_obs_data();
    void forward_modeling();
    void compute_residuals();
    void check_convergence();
    void optimization();
    void model_update();
};

__global__ void adjoint_state_kernel();

# endif