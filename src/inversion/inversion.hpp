# ifndef INVERSION_HPP
# define INVERSION_HPP

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

class Inversion
{
private:

    int iteration;
    int max_iteration;
    int smoother_samples;

    float smoother_stdv;

    std::string obs_data_folder;
    std::string obs_data_prefix;
    std::string convergence_map_folder;
    std::string estimated_model_folder;

    bool write_model_per_iteration;
    bool smooth_model_per_iteration;

    void show_information();
    void concatenate_data();
    void gradient_ray_tracing();

    void solve_linear_system_lscg();

    void smooth_volume(float * input, float * output, int nx, int ny, int nz);

protected:

    int n_data;
    int n_model;

    int tk_order;
    float tk_param;

    float * dcal = nullptr;
    float * dobs = nullptr;

    int M, N, NNZ;

    int * iA = nullptr;
    int * jA = nullptr;
    float * vA = nullptr;
    float * B = nullptr;
    float * x = nullptr; 

    std::vector< int > iG;
    std::vector< int > jG;
    std::vector<float> vG;

    Modeling * modeling = nullptr;

    std::vector<float> residuo;

    std::string inversion_name;
    std::string inversion_method;

    virtual void set_modeling_type() = 0;
    virtual void set_objective_function() = 0;
    virtual void set_sensitivity_matrix() = 0;
    virtual void set_regularization() = 0;

public:
    
    bool converged;

    std::string parameters;

    void set_parameters();
    void import_obsData();

    void forward_modeling();
    void check_convergence();

    void optimization();
    void model_update();

    void export_results();
};

# endif