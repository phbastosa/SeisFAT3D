# ifndef INVERSION_HPP
# define INVERSION_HPP

# include "../modeling/modeling.cuh"

class Inversion : public Modeling
{
private:


protected:

    int n_data;
    int n_model;

    float tolerance;

    int iteration;
    int max_iteration;

    float * dobs = nullptr;
    float * dcal = nullptr;

    float * gradient = nullptr;
    float * source = nullptr;       
    float * adjoint = nullptr;
    
    std::string obs_data_folder;
    std::string obs_data_prefix;

    std::vector<float> residuo;

    void adjoint_state_solver();
    void fill_calculated_data();

public:

    std::string file;

    bool converged;

    void optimization();
    void model_update();
    void export_results();
    void set_parameters();
    void import_obs_data();
    void forward_modeling();
    void compute_residuals();
    void check_convergence();
};

__global__ void adjoint_state_kernel();

# endif