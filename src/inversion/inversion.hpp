# ifndef INVERSION_HPP
# define INVERSION_HPP

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

class Inversion
{
private:

    int smoother_samples;

    float smoother_stdv;

    std::string obs_data_folder;
    std::string obs_data_prefix;
    std::string convergence_map_folder;

    bool write_model_per_iteration;
    bool smooth_model_per_iteration;

    void show_information();
    void concatenate_data();
    void gradient_ray_tracing();
    void solve_linear_system_lscg();
    void set_regularization_matrix();

protected:

    int n_data;
    int n_model;

    int iteration;
    int max_iteration;

    int tk_order;
    float tk_param;

    float * dcal = nullptr;
    float * dobs = nullptr;

    int M, N, NNZ;

    int * iA = nullptr;
    int * jA = nullptr;
    float * vA = nullptr;

    int * iR = nullptr;
    int * jR = nullptr;
    float * vR = nullptr;

    float * B = nullptr;
    float * x = nullptr; 

    float * W = nullptr;
    float * R = nullptr;    

    float * dS = nullptr;

    std::vector< int > iG;
    std::vector< int > jG;
    std::vector<float> vG;

    Modeling * modeling = nullptr;

    std::vector<float> residuo;

    std::string inversion_name;
    std::string inversion_method;
    std::string estimated_model_folder;

    virtual void set_modeling_type() = 0;
    virtual void set_sensitivity_matrix() = 0;
    virtual void get_parameter_variation() = 0;
    virtual void export_estimated_models() = 0;

    void model_smoothing(float * model);
    void smooth_volume(float * input, float * output, int nx, int ny, int nz);

public:
    
    bool converged;

    std::string parameters;

    void set_parameters();
    void import_obsData();

    void forward_modeling();
    void check_convergence();

    void optimization();
    
    virtual void model_update() = 0;

    void export_results();
};

# endif