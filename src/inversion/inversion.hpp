# ifndef INVERSION_HPP
# define INVERSION_HPP

# include "../modeling/modeling.cuh"

class Inversion
{
protected:

    int n_data, n_model, max_iteration;

    bool write_model_per_iteration;
    bool smooth_model_per_iteration;
    int iteration, smoother_samples;
    float smoother_stdv;

    float * dcal = nullptr;
    float * dobs = nullptr;

    float * perturbation = nullptr;

    Modeling * modeling;

    std::vector<float> residuo;

    std::string inversion_name;
    std::string inversion_method;

    std::string obs_data_folder;
    std::string obs_data_prefix;
    std::string convergence_map_folder;
    std::string estimated_model_folder;

    void set_forward_modeling();
    void set_inversion_elements();

    void show_information();
    void concatenate_data();

    virtual void set_specifications() = 0;
    virtual void apply_inversion_technique() = 0;

    void smooth_volume(float * input, float * output, int nx, int ny, int nz);

public:
    
    bool converged;

    std::string parameters;

    void set_parameters();
    void import_obsData();

    void forward_modeling();
    void check_convergence();

    virtual void optimization() = 0;
    
    void model_update();

    void export_results();
};

# endif