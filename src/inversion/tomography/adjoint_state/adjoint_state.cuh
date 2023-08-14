# ifndef ADJOINT_STATE_HPP
# define ADJOINT_STATE_HPP

# include "../tomography.hpp"

class Adjoint_State : public Tomography
{
private:


protected:


public:


};

# endif

//     int n_data;
//     int n_model;


//     float * gradient = nullptr;
//     float * source = nullptr;       
//     float * adjoint = nullptr;
    
//     float * d_source = nullptr;
//     float * d_adjoint = nullptr;



//     void adjoint_state_solver();
//     void fill_calculated_data();

//     bool converged;

//     void optimization();
//     void model_update();
//     void export_results();
//     void set_parameters();
//     void import_obs_data();
//     void forward_modeling();
//     void compute_residuals();
//     void check_convergence();

// __global__ void adjoint_state_kernel(float * adjoint, float * source, float * T, int level, int xOffset, int yOffset, 
//                                      int xSweepOffset, int ySweepOffset, int zSweepOffset, int nxx, int nyy, int nzz, 
//                                      float dx, float dy, float dz);

