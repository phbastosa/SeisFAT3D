# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../modeling/podvin_and_lecomte/podvin_and_lecomte.cuh"
# include "../modeling/fast_iterative_method/block_FIM.cuh"
# include "../modeling/fast_sweeping_method/accurate_FSM.cuh"

class Tomography
{
protected:

    int iteration;
    int max_iteration;
    int n_data, n_model;

    float tolerance;

    bool smooth;
    int smoother_samples;
    float smoother_stdv;

    float * dobs = nullptr;
    float * dcal = nullptr;

    float * dm = nullptr;
    float * model = nullptr;
    float * gradient = nullptr;

    bool write_model_per_iteration;
    bool write_gradient_per_iteration;

    std::vector<float> residuo;

    std::string obs_data_folder;
    std::string obs_data_prefix;

    std::string gradient_folder;
    std::string convergence_map_folder;
    std::string estimated_model_folder;

    std::string inversion_method;

    Eikonal * modeling;

    virtual void set_specific_parameters() = 0;
    virtual void apply_inversion_technique() = 0;
    virtual void gradient_preconditioning() = 0;

    void extract_calculated_data();
    void set_tomography_message();
    void set_inversion_volumes();
    void set_forward_modeling();
    void export_gradient();    

    void smooth_volume(float * input, float * output, int nx, int ny, int nz);

public:

    bool converged;

    std::string file;

    void set_parameters();
    void import_obs_data();

    void forward_modeling();
    void check_convergence();

    virtual void optimization() = 0;
    
    void model_update();
    void export_results();
};

# endif
