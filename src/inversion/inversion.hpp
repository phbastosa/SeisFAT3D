# ifndef INVERSION_HPP
# define INVERSION_HPP

# include <vector>
# include <string>

class Inversion
{
private:


protected:


public:

    bool converged;

    int iteration;
    int max_iteration;

    float tolerance;

    std::string file;

    virtual void optimization() = 0;
    virtual void model_update() = 0;
    virtual void export_results() = 0;
    virtual void set_parameters() = 0;
    virtual void import_obs_data() = 0;
    virtual void forward_modeling() = 0;
    virtual void check_convergence() = 0;
};

# endif

//     int n_data;
//     int n_model;

//     float * dobs = nullptr;
//     float * dcal = nullptr;

//     float * gradient = nullptr;
//     float * source = nullptr;       
//     float * adjoint = nullptr;
    
//     float * d_source = nullptr;
//     float * d_adjoint = nullptr;

//     std::string obs_data_folder;
//     std::string obs_data_prefix;

//     std::vector<float> residuo;

//     void adjoint_state_solver();
//     void fill_calculated_data();

//     std::string file;

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

