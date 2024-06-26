# ifndef INVERSION_CPP
# define INVERSION_CPP

# include "../modeling/modeling.hpp"

class Inversion
{
private:

    void print_information();

protected:

    int iteration;
    int max_iteration;
    int n_data, n_model;

    bool update_smooth;
    int smoother_samples;
    float smoother_stdv;

    float * dobs = nullptr;
    float * dcal = nullptr;

    float * gradient = nullptr;
    float * slowness = nullptr;
    float * variation = nullptr;

    float max_variation;

    bool export_model_per_iteration;

    std::vector<float> residuo;

    std::string type_name;
    std::string type_message;

    std::string obs_data_folder;
    std::string obs_data_prefix;

    std::string convergence_map_folder;
    std::string estimated_model_folder;

    Modeling * modeling;

    virtual void get_objective_function() = 0;

    virtual void set_forward_modeling() = 0;
    virtual void set_inversion_volumes() = 0;

    virtual void extract_calculated_data() = 0;
    virtual void adjoint_propagation() = 0;

    virtual void update_specifications() = 0;

    void update_smoothing();

    void gaussian_smoothing(float * input, float * output, int nx, int ny, int nz);

public:

    std::string file;

    bool converged;

    void set_parameters();
    
    virtual void import_obs_data() = 0;

    void forward_modeling();
    void check_convergence();

    void optimization();
    void model_update();

    void export_results();
};



# endif